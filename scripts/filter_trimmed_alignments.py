import sys
import gzip
import subprocess
import multiprocessing

from intervaltree import Interval, IntervalTree
import pysam
import pandas as pd
import numpy as np



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
COMP_NTS = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
MAX_MISMATCHES = 5
MAX_MISMATCH_PERCENT = 0.1
MAX_GAP_LENGTH = 3
ALIGNMENTS_INITIAL_COLUMNS = ('read_id', 'trim_seq', 'chrom', 'align_is_reverse', 'align_start', 'align_end')
FAILED_ALIGNMENTS_COLS = ['read_id',
						  'read_is_reverse',
						  'trim_seq',
						  'chrom',
						  'align_is_reverse',
						  'align_start',
						  'align_end', 
						  'fivep_pos',
						  'read_bp_pos',
						  'read_bp_nt',
						  'bp_pos',
						  'threep_pos',
						  'bp_dist_to_threep',
						  'genomic_bp_context',
						  'genomic_bp_nt',
						  'filter_failed']

TEMPLATE_SWITCHING_COLS = ['read_id',
						'fivep_sites',
						'temp_switch_sites',
						'read_seq', 
						'fivep_seq',
						'trim_seq']

FINAL_INFO_TABLE_COLS = ['read_id', 
						'read_is_reverse', 
						'read_seq', 
						'chrom', 
						'strand', 
						'align_start',
						'align_end', 
						'align_is_reverse',
						'fivep_pos', 
						'bp_pos', 
						'read_bp_pos',
						'read_bp_nt', 
						'genomic_bp_context', 
						'genomic_bp_nt', 
						'threep_pos', 
						'bp_dist_to_threep',
						'gene_id', 
						'gene_name', 
						'gene_type', 
						]

filtered_out_lock = multiprocessing.Lock()
failed_out_lock = multiprocessing.Lock()
temp_switch_lock = multiprocessing.Lock()



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def comma_join(x): 
	return ','.join(set(x))


def reverse_complement(seq):
	return ''.join([COMP_NTS[seq[i]] for i in range(len(seq)-1,-1,-1)])


def parse_intron_info(ref_introns:str) -> tuple:
	'''
	Returns a dict formatted as follows:
	{Chromosome: {Strand(+ or -): Intervaltree(StartPosition(int), EndPosition(int), {"gene_id": GeneID})}}
	'''
	if ref_introns.split('.')[-1] == 'gz':
		intron_file = gzip.open(ref_introns, 'rt')
	else:
		intron_file = open(ref_introns)

	introns, introns_done = {}, set()
	fivep_genes = {}
	for line in intron_file:
		chrom, start, end, intron_info, _, strand = line.strip().split('\t')
		if chrom not in introns:
			introns[chrom] = {s: IntervalTree() for s in ['+', '-']}
		intron_id = '{}_{}_{}_{}'.format(chrom, strand, start, end)
		if intron_id not in introns_done:
			start, end = int(start), int(end)
			if end-start > 20:
				intron_genes = set([i.split('-')[0] for i in intron_info.split(';')[-1].split('|')])
				introns[chrom][strand].add(Interval(start, end, {'gene_id':intron_genes}))
				introns_done.add(intron_id)
				fivep_pos = start if strand=='+' else end-1
				fivep_site = f'{chrom};{fivep_pos};{strand}'
				if fivep_site not in fivep_genes:
					fivep_genes[fivep_site] = set()
				fivep_genes[fivep_site].update(intron_genes)

	intron_file.close()
	return introns, fivep_genes


def parse_attributes(attribute_string:str, anno_type:str) -> dict:
	if anno_type == 'gtf':
		attributes = attribute_string.rstrip('";').split('; ')
		attributes = [attr.split(' ') for attr in attributes]
		tags = [attr_val.strip('"') for attr_name, attr_val in attributes if attr_name=='tag']
		attributes = {attr_name: attr_val.strip('"') for attr_name, attr_val in attributes if attr_name!='tag'}
		attributes['tags'] = tags
	else:
		attributes = [attr.split('=') for attr in attribute_string.split(';')]
		attributes = [(attr[0].lstrip(), attr[1]) for attr in attributes]
		attributes = dict(attributes)
		if 'tag' in attributes:
			attributes['tags'] = attributes['tag'].split(',')

	return attributes


def parse_gene_info(ref_anno:str) -> dict:
	prev_ext, last_ext = ref_anno.split('.')[-2:]
	if last_ext == 'gz':
		in_file, anno_type = gzip.open(ref_anno, 'rt'), prev_ext
	else:
		in_file, anno_type = open(ref_anno), last_ext

	genes = {}
	for line in in_file:
		if line[0] != '#':
			_, _, feature, _, _, _, _, _, attributes = line.strip().split('\t')
			if feature == 'transcript':
				attributes = parse_attributes(attributes, anno_type)
				genes[attributes['gene_id']] = [attributes[e] for e in ['gene_name', 'gene_type']]
	in_file.close()
	
	return genes


def parse_fivep_alignments(output_base:str, fivep_genes:dict):
	fivep_alignments = pd.read_csv(f'{output_base}fivep_info_table.tsv', sep='\t')

	# Unpack fivep_sites
	fivep_alignments.fivep_sites = fivep_alignments.fivep_sites.str.split(',')
	fivep_alignments = fivep_alignments.explode('fivep_sites')
	fivep_alignments['gene_id'] = fivep_alignments.apply(lambda row: tuple(fivep_genes[row['fivep_sites']]), axis=1)
	fivep_alignments[['fivep_chrom', 'fivep_pos', 'strand']] = fivep_alignments.fivep_sites.str.split(';', expand=True)
	fivep_alignments.fivep_pos = fivep_alignments.fivep_pos.astype(int)

	return fivep_alignments


def parse_trimmed_alignments(sam_file:str) -> pd.DataFrame:
	alignments = []
	for align in pysam.AlignmentFile(sam_file, 'r'):		
		mismatches = align.get_tag('XM')
		mismatch_percent = mismatches/align.query_length
		gaps = [length for tag_num, length in align.cigartuples if tag_num in (1,2)]
		
		# Check alignment quality and fail alignment if it doesn't meet all thresholds
		if mismatches > MAX_MISMATCHES:
			continue
		if mismatch_percent > MAX_MISMATCH_PERCENT:
			continue
		if len(gaps) > 1 or (len(gaps) == 1 and gaps[0] > MAX_GAP_LENGTH):
			continue
		
		# If the alignment passed the quality thresholds, add it to the list
		align_info = [
					align.query_name,
					align.query_sequence,
					align.reference_name, 
					align.is_reverse, 
					align.reference_start, 
					align.reference_end,
					]
		alignments.append(align_info)
	
	# Flatten the passed alignments into a list of lists, then covert to a DataFrame and return
	alignments = pd.DataFrame(alignments, columns=ALIGNMENTS_INITIAL_COLUMNS)
	return alignments
	

def get_bp_seqs(alignments:pd.DataFrame, genome_fasta:str):
	alignments = alignments.copy()
	alignments['window_start'] = alignments.apply(lambda row: row['bp_pos']-4 if row['strand']=='+' else row['bp_pos']-5, axis=1).astype(str)
	alignments['window_end'] = alignments.apply(lambda row: row['bp_pos']+6 if row['strand']=='+' else row['bp_pos']+5, axis=1).astype(str)
	alignments['align_num'] = alignments.index.to_series().astype(str)
	alignments['zero'] = '0'
	# alignments['bedtools_line'] = alignments.apply(lambda row: '\t'.join([str(e) for e in row[['chrom', 'window_start', 'window_end', 'align_num', 'zero', 'strand']]]), axis=1)
	alignments['bedtools_line'] = alignments[['chrom', 'window_start', 'window_end', 'align_num', 'zero', 'strand']].agg('\t'.join, axis=1)
	
	bedtools_input = '\n'.join(alignments.bedtools_line) + '\n'
	bedtools_call = f'bedtools getfasta -nameOnly -s -tab -fi {genome_fasta} -bed -'
	bedtools_lines = subprocess.run(bedtools_call.split(' '), input=bedtools_input, check=True, capture_output=True, text=True).stdout.strip().split('\n')
	
	bp_seqs = pd.DataFrame([l.split('\t') for l in bedtools_lines], columns=['align_num', 'genomic_bp_context'])
	bp_seqs.align_num = bp_seqs.align_num.str.slice(0,-3)
	bp_seqs.genomic_bp_context = bp_seqs.genomic_bp_context.transform(lambda genomic_bp_context: genomic_bp_context.upper())
	
	alignments = pd.merge(alignments, bp_seqs, on='align_num')
	return alignments


def add_read_bp_nt(row:pd.Series) -> str:
	if not row['read_is_reverse'] and not row['align_is_reverse']:
		return row['trim_seq'][-1]
	elif not row['read_is_reverse'] and row['align_is_reverse']:
		return reverse_complement(row['trim_seq'][0])
	elif row['read_is_reverse'] and not row['align_is_reverse']:
		return reverse_complement(row['trim_seq'][0])
	elif row['read_is_reverse'] and row['align_is_reverse']:
		return row['trim_seq'][-1]


def add_nearest_threep(row:pd.Series):
	if row['strand'] == '+':
		threep_pos = min(row['overlap_introns'], key=lambda s: s.end-row['bp_pos']).end - 1
	else:
		threep_pos = min(row['overlap_introns'], key=lambda s: row['bp_pos']-s.begin).begin

	return threep_pos


def final_filters(row:pd.Series):
	# Check if bad alignment orientation combination, which do NOT leave the branchpoint adjacent to the 5'ss in the read as is expected of lariats
	if row['read_is_reverse'] is False:
		if row['align_is_reverse'] is False and row['strand']=='-':
			return 'wrong_orientation'
		elif row['align_is_reverse'] is True and row['strand']=='+':
			return 'wrong_orientation'
	if row['read_is_reverse'] is True:
		if row['align_is_reverse'] is False and row['strand']=='+':
			return 'wrong_orientation'
		elif row['align_is_reverse'] is True and row['strand']=='-':
			return 'wrong_orientation'

	# Check if the 5'ss is at or downstream of the BP
	if row['strand'] == '+' and row['fivep_pos'] >= row['bp_pos']:
		return '5p_bp_order'
	if row['strand'] == '-' and row['fivep_pos'] <= row['bp_pos']:
		return '5p_bp_order'

	return pd.NA


def drop_failed_alignments(trim_alignments:pd.DataFrame, output_base:str) -> pd.DataFrame:
		trim_alignments = trim_alignments.copy()

		# Get alignments that failed one of the filters
		failed_alignments = trim_alignments.loc[trim_alignments.filter_failed.notna()].copy()

		# Add in any missing cols as needed
		for col in FAILED_ALIGNMENTS_COLS:
			if col not in failed_alignments.columns:
				failed_alignments[col] = ''
		
		# Arrange cols and a write to file
		failed_alignments = failed_alignments[FAILED_ALIGNMENTS_COLS].drop_duplicates()
		with failed_out_lock:
			failed_alignments.to_csv(f'{output_base}failed_trimmed_alignments.tsv', mode='a', sep='\t', header=False, index=False)

		# Remove failed alignments from DataFrame
		trim_alignments = trim_alignments.loc[trim_alignments.filter_failed.isna()]
		return trim_alignments


def filter_trimmed_alignments(trim_alignments:pd.DataFrame, genes:dict, introns:dict, genome_fasta:str, fivep_alignments:pd.DataFrame, output_base:str) -> pd.DataFrame:
	# Merge
	trim_alignments = pd.merge(trim_alignments, fivep_alignments, 'left', on=['read_id'])
	if trim_alignments.fivep_pos.isna().any():
		raise RuntimeError(f"{trim_alignments.fivep_pos.isna().sum()} alignments didn't match any fivep_info_table read IDs, this shouldn't be possible")
	
	# Infer info
	trim_alignments.fivep_pos = trim_alignments.fivep_pos.astype(int)
	trim_alignments['bp_pos'] = trim_alignments.apply(lambda row: row['align_end']-1 if row['strand']=='+' else row['align_start'], axis=1)
	trim_alignments = get_bp_seqs(trim_alignments, genome_fasta)
	trim_alignments['genomic_bp_nt'] = trim_alignments.genomic_bp_context.str.get(4)
	trim_alignments['read_bp_nt'] = trim_alignments.apply(add_read_bp_nt, axis=1)
	
	# Identify template-switching reads
	trim_alignments['template_switching'] = trim_alignments.apply(lambda row: row['genomic_bp_context'][5:]==row['fivep_seq'][:5].upper(), axis=1)

	# Output template-switching reads
	temp_switch_alignments = trim_alignments.loc[trim_alignments.template_switching].copy()
	temp_switch_alignments = temp_switch_alignments.astype(str)
	temp_switch_alignments['read_id'] = temp_switch_alignments.read_id.transform(lambda rid: rid.split('/')[0])
	temp_switch_alignments['fivep_sites'] = temp_switch_alignments[['fivep_chrom', 'strand', 'fivep_pos']].agg(';'.join, axis=1)
	temp_switch_alignments['temp_switch_sites'] = temp_switch_alignments[['chrom', 'bp_pos']].agg(';'.join, axis=1)
	temp_switch_alignments = temp_switch_alignments[['read_id', 'read_seq', 'fivep_seq', 'trim_seq', 'fivep_sites', 'temp_switch_sites']]
	temp_switch_alignments = temp_switch_alignments.groupby(['read_id', 'read_seq']).agg({'fivep_sites':comma_join, 'temp_switch_sites':comma_join, 'fivep_seq':comma_join, 'trim_seq':comma_join}).reset_index()
	temp_switch_alignments = temp_switch_alignments[TEMPLATE_SWITCHING_COLS]
	with temp_switch_lock:
		temp_switch_alignments.to_csv(f'{output_base}template_switching_reads.tsv', mode='a', sep='\t', header=False, index=False)

	# Filter out template-switching reads
	trim_alignments = trim_alignments.loc[~trim_alignments.template_switching].drop(columns='template_switching')
	if trim_alignments.empty:
		return trim_alignments
	
	# Identify introns that overlap the alignment
	trim_alignments['overlap_introns'] = trim_alignments.apply(lambda row: introns[row['chrom']][row['strand']].overlap(row['align_start'], row['align_end']), axis=1)
	
	# Filter out alignments that don't overlap an alignment
	trim_alignments.loc[trim_alignments.overlap_introns.transform(len)==0, 'filter_failed'] = 'overlaps_intron'
	trim_alignments = drop_failed_alignments(trim_alignments, output_base)
	if trim_alignments.empty:
		return trim_alignments
	
	# Filter out alignments where 5'ss and BP segments aren't in the same gene
	trim_alignments = trim_alignments.explode('gene_id')
	trim_alignments.overlap_introns = trim_alignments.apply(lambda row: tuple(intron for intron in row['overlap_introns'] if row['gene_id'] in intron.data['gene_id']), axis=1)
	trim_alignments.loc[trim_alignments.overlap_introns.transform(len).eq(0), 'filter_failed'] = 'fivep_intron_match'
	trim_alignments = drop_failed_alignments(trim_alignments, output_base)
	if trim_alignments.empty:
		return trim_alignments
	
	# Filter alignments based on proper read orientation and 5'ss-BP ordering
	trim_alignments['filter_failed'] = trim_alignments.apply(final_filters, axis=1, result_type='reduce')
	trim_alignments = drop_failed_alignments(trim_alignments, output_base)
	if trim_alignments.empty:
		return trim_alignments
	
	# Infer more info
	trim_alignments[['gene_name', 'gene_type']] = trim_alignments.apply(lambda row:genes[row['gene_id']], axis=1, result_type='expand')
	trim_alignments['threep_pos'] = trim_alignments.apply(add_nearest_threep, axis=1)
	trim_alignments['bp_dist_to_threep'] = trim_alignments.apply(lambda row: -abs(row['bp_pos']-row['threep_pos']) if pd.notna(row['threep_pos']) else pd.NA, axis=1)

	# Drop all the uneeded columns
	trim_alignments = trim_alignments[FINAL_INFO_TABLE_COLS]
	
	# Some reads get mapped to coordinates with multiple overlapping gene annotations
	# We resolve this by collapsing the duplicated rows and concatenating the gene_id, gene_name, and gene_type columns
	# Code adapted from https://stackoverflow.com/questions/27298178/concatenate-strings-from-several-rows-using-pandas-groupby
	trim_alignments = (trim_alignments.groupby([col for col in trim_alignments.columns if col not in ('gene_id', 'gene_name', 'gene_type')], as_index=False)
								.agg({'gene_id': comma_join, 'gene_name': comma_join, 'gene_type': comma_join})
					)

	with filtered_out_lock:
		trim_alignments.to_csv(f'{output_base}trimmed_info_table.tsv', mode='a', sep='\t', index=False, header=False)



# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	threads, ref_anno, ref_introns, genome_fasta, output_base = sys.argv[1:]
	threads = int(threads)

	# Load auxiliary data for processing alignments
	introns, fivep_genes = parse_intron_info(ref_introns)
	genes = parse_gene_info(ref_anno)
	fivep_alignments = parse_fivep_alignments(output_base, fivep_genes)
	
	# Load trim_alignments, skipping bad-quality alignments
	trim_alignments = parse_trimmed_alignments(f'{output_base}trimmed_reads_to_genome.sam')

	#TODO: change parse_intron_info to have all possible chromosomes, even those with 0 annotated genes, so this fix isn't needed to prevent a key-not-found error later
	for chrom in trim_alignments.chrom.unique():
		if chrom not in introns.keys():
			introns[chrom] ={s: IntervalTree() for s in ['+', '-']}

	# Record read count
	rids = set(rid.split('/')[0] for rid in fivep_alignments.read_id.values)
	with open(f'{output_base}read_counts.tsv', 'a') as a:
		a.write(f'fivep_filter_passed\t{len(rids)}\n')	

	# Write headers for the failed_trimmed_alignments and template_switching_alignments out-files
	# The rows will get appended in chunks during filter_thread()
	with open(f'{output_base}trimmed_info_table.tsv', 'w') as w:
		w.write('\t'.join(FINAL_INFO_TABLE_COLS) + '\n')
	with open(f'{output_base}template_switching_reads.tsv', 'w') as w:
		w.write('\t'.join(TEMPLATE_SWITCHING_COLS) + '\n')
	with open(f'{output_base}failed_trimmed_alignments.tsv', 'w') as w:
		w.write('\t'.join(FAILED_ALIGNMENTS_COLS) + '\n')

	# Split the alignments DataFrame into a list of (almost) equal-sized chunks
	trim_alignments = np.array_split(trim_alignments, threads)

	# Make a wrapper function that processes in the multiprocessing pool will use
	# We need to do this because imap_unordered only passes 1 positional arg to the specified function
	def filter_chunk(trim_alignments_chunk):
		return filter_trimmed_alignments(trim_alignments_chunk, genes=genes, introns=introns, genome_fasta=genome_fasta, fivep_alignments=fivep_alignments, output_base=output_base)
	
	# Enter the Pool context manager 
	with multiprocessing.Pool(processes=threads) as pool:
		# Iteratively call filter_chunk with each chunk in trim_alignments
		# Each iteration creates a process that filters the alignments in its chunk and appends the results to final_info_table.tsv 
		for _ in pool.imap_unordered(filter_chunk, trim_alignments):
			pass

	# Get rid of duplicated rows in template_switching_alignments
	# They show up because alignments to the same original read can end up in different threads 
	# temp_switch_alignments = pd.read_csv(f'{output_base}template_switching_reads.tsv', sep='\t')
	# temp_switch_alignments.drop_duplicates().to_csv(f'{output_base}template_switching_reads.tsv', sep='\t', index=False)


