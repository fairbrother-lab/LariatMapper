import sys
import gzip
import subprocess
import pysam
import pandas as pd
import numpy as np
import itertools
from intervaltree import Interval, IntervalTree
import threading
import os

# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
COMP_NTS = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
MAX_MISMATCHES = 5
MAX_MISMATCH_PERCENT = 0.1
MAX_GAP_LENGTH = 3
ALIGNMENTS_INITIAL_COLUMNS = ('align_num', 'read_id', 'trim_seq', 'chrom', 'align_is_reverse', 'align_start', 'align_end')
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
						'gene_id', 
						'gene_name', 
						'gene_type', 
						'fivep_pos', 
						'bp_pos', 
						'read_bp_pos',
						'read_bp_nt', 
						'genomic_bp_context', 
						'genomic_bp_nt', 
						'threep_pos', 
						'bp_dist_to_threep']
failed_out_lock = threading.Lock()

# =============================================================================#
#                                   Classes                                    #
# =============================================================================#
class filter_thread (threading.Thread):
	def __init__(self, genes, introns, genome_fasta, fivep_alignments, trim_alignments, output_base):
		threading.Thread.__init__(self)
		self.genes, self.introns, self.genome_fasta, self.fivep_alignments, self.trim_alignments, self.output_base= \
		genes, introns, genome_fasta, fivep_alignments, trim_alignments, output_base
	

	def drop_failed_alignments(trim_alignments:pd.DataFrame) -> pd.DataFrame:
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
	

	def run(self):
		genes, introns, genome_fasta, fivep_alignments, trim_alignments, output_base = \
		self.genes, self.introns, self.genome_fasta, self.fivep_alignments, self.trim_alignments, self.output_base
		
		# Merge
		trim_alignments = pd.merge(trim_alignments, fivep_alignments, 'left', on=['read_id'])
		# failed_alignments = trim_alignments.loc[trim_alignments.fivep_pos.isna()].copy()
		# failed_alignments = failed_alignments.assign(c='', d='', e='', f='', g='', filter_failed = 'fivep_match')
		# failed_alignments = failed_alignments[[col for col in FAILED_ALIGNMENTS_COLS if col in failed_alignments.columns]].drop_duplicates()
		# with failed_out_lock:
		# 	failed_alignments.to_csv(f'{output_base}failed_trimmed_alignments.tsv', mode='a', sep='\t', header=False, index=False)
		# trim_alignments = trim_alignments.loc[trim_alignments.fivep_pos.notna()]
		# trim_alignments.loc[trim_alignments.fivep_pos.isna(), 'filter_failed'] = 'fivep_match'
		# trim_alignments = filter_thread.drop_failed_alignments(trim_alignments)
		# if trim_alignments.empty:
			# self.trim_alignments = trim_alignments
			# return
		if trim_alignments.fivep_pos.isna().any():
			raise RuntimeError(f"{trim_alignments.fivep_pos.isna().sum()} alignments didn't match any fivep_info_table read IDs, this shouldn't be possible")
		
		# Infer info
		trim_alignments.fivep_pos = trim_alignments.fivep_pos.astype(int)
		trim_alignments['bp_pos'] = trim_alignments.apply(lambda row: row['align_end']-1 if row['strand']=='+' else row['align_start'], axis=1)
		trim_alignments = get_bp_seqs(trim_alignments, genome_fasta)
		trim_alignments['genomic_bp_nt'] = trim_alignments.genomic_bp_context.str.get(4)
		trim_alignments['read_bp_nt'] = trim_alignments.apply(add_read_bp_nt, axis=1)
		
		# Filter out template switching alignments
		trim_alignments['template_switching'] = trim_alignments.apply(lambda row: row['genomic_bp_context'][5:]==row['fivep_seq'][:5].upper(), axis=1)
		# failed_alignments = trim_alignments.loc[trim_alignments.template_switching].copy()
		# with failed_out_lock:
		# 	output_header = not os.path.isfile(f'{output_base}template_switching_alignments.tsv')
		# 	failed_alignments['read_id'] = failed_alignments.apply(lambda row: row['read_id'].split('/')[0], axis=1)
		# 	failed_alignments['fivep_sites'] = failed_alignments.apply(lambda row:';'.join([str(e) for e in row[['fivep_chrom', 'strand', 'fivep_pos']]]), axis=1)
		# 	failed_alignments['temp_switch_sites'] = failed_alignments.apply(lambda row:';'.join([str(e) for e in row[['chrom', 'bp_pos']]]), axis=1)
		# 	failed_alignments = failed_alignments[['read_id', 'read_seq', 'fivep_seq', 'trim_seq', 'fivep_sites', 'temp_switch_sites']]
		# 	def comma_join(x): return ','.join(set(x))
		# 	failed_alignments = failed_alignments.groupby(['read_id', 'read_seq']).agg({'fivep_sites':comma_join, 'temp_switch_sites':comma_join, 'fivep_seq':comma_join, 'trim_seq':comma_join}).reset_index()
		# 	failed_alignments.to_csv(f'{output_base}template_switching_alignments.tsv', mode='a', sep='\t', header=output_header, index=False)
		temp_switch_alignments = trim_alignments.loc[trim_alignments.template_switching].copy()
		temp_switch_alignments['read_id'] = temp_switch_alignments.apply(lambda row: row['read_id'].split('/')[0], axis=1)
		temp_switch_alignments['fivep_sites'] = temp_switch_alignments.apply(lambda row:';'.join([str(e) for e in row[['fivep_chrom', 'strand', 'fivep_pos']]]), axis=1)
		temp_switch_alignments['temp_switch_sites'] = temp_switch_alignments.apply(lambda row:';'.join([str(e) for e in row[['chrom', 'bp_pos']]]), axis=1)
		temp_switch_alignments = temp_switch_alignments[['read_id', 'read_seq', 'fivep_seq', 'trim_seq', 'fivep_sites', 'temp_switch_sites']]
		temp_switch_alignments = temp_switch_alignments.groupby(['read_id', 'read_seq']).agg({'fivep_sites':comma_join, 'temp_switch_sites':comma_join, 'fivep_seq':comma_join, 'trim_seq':comma_join}).reset_index()
		temp_switch_alignments = temp_switch_alignments[TEMPLATE_SWITCHING_COLS]
		temp_switch_alignments.to_csv(f'{output_base}template_switching_alignments.tsv', mode='a', sep='\t', header=False, index=False)

		trim_alignments = trim_alignments.loc[~trim_alignments.template_switching].drop(columns='template_switching')
		if trim_alignments.empty:
			self.trim_alignments = trim_alignments
			return
		
		# Filter out alignments that aren't within an intron
		# trim_alignments['overlaps_intron'] = trim_alignments.apply(lambda row: False if row['chrom'] not in introns else introns[row['chrom']][row['strand']].overlaps(row['align_start'], row['align_end']), axis=1)
		# failed_alignments = trim_alignments.loc[trim_alignments.overlaps_intron==False].copy().drop(columns='overlaps_intron')
		# failed_alignments = failed_alignments[[col for col in FAILED_ALIGNMENTS_COLS if col in failed_alignments.columns]].drop_duplicates()
		# with failed_out_lock:
		# 	failed_alignments.to_csv(f'{output_base}failed_trimmed_alignments.tsv', mode='a', sep='\t', header=False, index=False)
		# trim_alignments['filter_failed'] = trim_alignments.apply(intron_overlap_filter, introns=introns, axis=1)
		trim_alignments['overlap_introns'] = trim_alignments.apply(lambda row: introns[row['chrom']][row['strand']].overlaps(row['align_start', 'align_end']), axis=1)
		trim_alignments.loc[trim_alignments.overlap_introns.transform(len)==0, 'filter_failed'] = 'overlaps_intron'
		trim_alignments = filter_thread.drop_failed_alignments(trim_alignments)
		if trim_alignments.empty:
			self.trim_alignments = trim_alignments
			return
		
		# Filter out alignments where 5'ss and BP segments aren't in the same gene
		# trim_alignments['trimmed_gene_id'] = trim_alignments.apply(lambda row: list(itertools.chain(*[i.data['gene_id'] for i in introns[row['chrom']][row['strand']].overlap(row['align_start'], row['align_end']) if i.begin<=row['align_start'] and i.end>row['align_end']])), axis=1)
		# trim_alignments = trim_alignments.explode('trimmed_gene_id').drop_duplicates()
		# trim_alignments = trim_alignments.loc[trim_alignments.trimmed_gene_id==trim_alignments.gene_id].drop(columns='trimmed_gene_id')
		trim_alignments = trim_alignments.explode('gene_id').drop_duplicates()
		trim_alignments.overlap_introns = trim_alignments.apply(lambda row: tuple(intron for intron in row['overlap_introns'] if intron.data[0]['gene_id']==row['gene_id']), axis=1)
		trim_alignments.loc[trim_alignments.overlap_introns.transform(len)==0, 'filter_failed'] = 'fivep_intron_match'
		trim_alignments = filter_thread.drop_failed_alignments(trim_alignments)
		if trim_alignments.empty:
			self.trim_alignments = trim_alignments
			return
		
		# Filter alignments based on proper read orientation and 5'ss-BP ordering
		trim_alignments['filter_failed'] = trim_alignments.apply(final_filters, axis=1)
		trim_alignments = filter_thread.drop_failed_alignments(trim_alignments)
		if trim_alignments.empty:
			self.trim_alignments = trim_alignments
			return
		
		# Infer more info
		trim_alignments[['gene_name', 'gene_type']] = trim_alignments.apply(lambda row:genes[row['gene_id']], axis=1, result_type='expand')
		# trim_alignments['threep_pos'] = trim_alignments.apply(add_nearest_threep, introns=introns, axis=1)
		trim_alignments['threep_pos'] = trim_alignments.apply(add_nearest_threep, axis=1)
		trim_alignments['bp_dist_to_threep'] = trim_alignments.apply(lambda row: -abs(row['bp_pos']-row['threep_pos']) if pd.notna(row['threep_pos']) else pd.NA, axis=1)

		# Code adapted from https://stackoverflow.com/questions/27298178/concatenate-strings-from-several-rows-using-pandas-groupby
		# Some reads get mapped to coordinates with multiple overlapping gene annotations
		# We resolve this by collapsing the duplicated rows and concatenating the gene_id, gene_name, and gene_type columns
		trim_alignments = (trim_alignments.groupby([col for col in trim_alignments.columns if col not in ('gene_id', 'gene_name', 'gene_type')])
									.agg({'gene_id': comma_join, 'gene_name': comma_join, 'gene_type': comma_join})
									.reset_index()
						)
		
		# Return the alignments that passed filtering
		self.trim_alignments = trim_alignments[FINAL_INFO_TABLE_COLS]
		

# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def comma_join(x): 
	return ','.join(set(x))


def reverse_complement(seq):
	return ''.join([COMP_NTS[seq[i]] for i in range(len(seq)-1,-1,-1)])


# def intron_overlap_filter(row, introns):
# 	has_overlaps = introns[row['chrom']][row['strand']].overlaps(row['align_start'], row['align_end'])
# 	if has_overlaps is True:
# 		return np.nan
# 	else:
# 		return 'overlaps_intron'


def parse_intron_info(ref_introns:str) -> tuple:
	'''
	Returns a dict formatted as follows:
	{Chromosome: {Strand(+ or -): Intervaltree(StartPosition(int), EndPosition(int))}}
	for 3' splice sites (+/- 2 bases), 5' splice sites (+/- 2 bases), and the introns they come from (start to end)
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
				intron_genes = [i.split('-')[0] for i in intron_info.split(';')[-1].split('|')]
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
	fivep_alignments['gene_id'] = fivep_alignments.apply(lambda row: fivep_genes[row['fivep_sites']], axis=1)
	fivep_alignments[['fivep_chrom', 'fivep_pos', 'strand']] = fivep_alignments.fivep_sites.str.split(';', expand=True)
	fivep_alignments.fivep_pos = fivep_alignments.fivep_pos.astype(int)

	return fivep_alignments


def parse_trimmed_alignments(bam_file:str):
	# alignments = {}
	alignments = []
	i = 0
	for align in pysam.AlignmentFile(bam_file, 'rb'):
		i += 1
		
		# read_is_reverse = align.query_name.endswith('_rev')
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
		align_info = [i,
					align.query_name,
					# read_is_reverse,
					align.query_sequence,
					align.reference_name, 
					align.is_reverse, 
					align.reference_start, 
					align.reference_end,]
		alignments.append(align_info)
	
	# Flatten the passed alignments into a list of lists, then covert to a DataFrame and return
	alignments = pd.DataFrame(alignments, columns=ALIGNMENTS_INITIAL_COLUMNS)
	return alignments


def get_bp_seqs(alignments:pd.DataFrame, genome_fasta:str):
	alignments['window_start'] = alignments.apply(lambda row: row['bp_pos']-4 if row['strand']=='+' else row['bp_pos']-5, axis=1)
	alignments['window_end'] = alignments.apply(lambda row: row['bp_pos']+6 if row['strand']=='+' else row['bp_pos']+5, axis=1)
	alignments['zero'] = 0
	alignments['bedtools_line'] = alignments.apply(lambda row: '\t'.join([str(e) for e in row[['chrom', 'window_start', 'window_end', 'align_num', 'zero', 'strand']]]), axis=1)
	bedtools_input = '\n'.join(alignments.bedtools_line) + '\n'
	bedtools_call = f'bedtools getfasta -nameOnly -s -tab -fi {genome_fasta} -bed -'
	bedtools_lines = subprocess.run(bedtools_call.split(' '), input=bedtools_input, check=True, capture_output=True, text=True).stdout.strip().split('\n')
	bp_seqs = pd.DataFrame([l.split('\t') for l in bedtools_lines], columns=['align_num', 'genomic_bp_context'])
	bp_seqs.align_num = bp_seqs.align_num.transform(lambda align_num: int(align_num[:-3]))
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
# def add_nearest_threep(row:pd.Series, introns:dict):
	# enveloping_introns = list(introns[row['chrom']][row['strand']].at(row['bp_pos']))
	# if len(enveloping_introns) == 0:
	# 	raise RuntimeError(f"No introns enveloped the BP in the following:\n{row}")
	
	# if row['strand'] == '+':
	# 	threep_pos = min(enveloping_introns, key=lambda s: s.end-row['bp_pos']).end - 1
	# else:
	# 	threep_pos = min(enveloping_introns, key=lambda s: row['bp_pos']-s.begin).begin
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



# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':

	threads, ref_anno, ref_introns, genome_fasta, output_base = sys.argv[1:]

	# Load auxiliary data for processing alignments
	introns, fivep_genes = parse_intron_info(ref_introns)
	genes = parse_gene_info(ref_anno)
	fivep_alignments = parse_fivep_alignments(output_base, fivep_genes)
	
	# Load trim_alignments, skipping bad-quality alignments
	trim_alignments = parse_trimmed_alignments(f'{output_base}trimmed_reads_to_genome.bam')

	# DELETE THIS AFTER TESTING, NOT INTENDED BEHAVIOR
	trim_alignments = trim_alignments.loc[~trim_alignments.chrom.str.contains('\.')].reset_index(drop=True)

	# Write headers for the failed_trimmed_alignments and template_switching_alignments out-files
	# The rows will get appended in chunks during filter_thread()
	with open(f'{output_base}failed_trimmed_alignments.tsv', 'w') as w:
		w.write('\t'.join(FAILED_ALIGNMENTS_COLS) + '\n')
	with open(f'{output_base}template_switching_alignments.tsv', 'w') as w:
		w.write('\t'.join(TEMPLATE_SWITCHING_COLS) + '\n')

	thread_list, threads = [], int(threads)
	for thread_num in range(threads):
		new_thread = filter_thread(genes, introns, genome_fasta, fivep_alignments, trim_alignments[thread_num::threads].copy(), output_base)
		new_thread.start()
		thread_list.append(new_thread)

	for t in thread_list:
		t.join()
	
	# Write filtered alignments to final info file
	trim_alignments = pd.concat([t.trim_alignments for t in thread_list if not t.trim_alignments.empty])
	trim_alignments.to_csv(f'{output_base}final_info_table.tsv', sep='\t', index=False)

	# Get rid of duplicated rows in template_switching_alignments
	# They show up because alignments to the same original read can end up in different threads 
	temp_switch_alignments = pd.read_csv(f'{output_base}template_switching_alignments.tsv', sep='\t')
	temp_switch_alignments.drop_duplicates().to_csv(f'{output_base}template_switching_alignments.tsv', sep='\t', index=False)

	

