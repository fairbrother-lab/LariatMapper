import sys
import os
import multiprocessing as mp
import collections
import tempfile

from intervaltree import Interval, IntervalTree
import pandas as pd

import functions



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
# In files
HEADS_TO_GENOME_FILE = "{}heads_to_genome.sam"
TAILS_FILE = "{}tails.tsv"
# Out files
FAILED_HEADS_FILE = "{}failed_head_alignments.tsv"
TEMP_SWITCH_FILE = "{}template_switching_reads.tsv"
CIRCULARS_FILE = "{}circularized_intron_reads.tsv"
PUTATITVE_LARIATS_FILE = "{}putative_lariats.tsv"
RUN_DATA_FILE = "{}read_counts.tsv"

CIGAR_OPERATORS = ('M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X')
MAX_MISMATCHES = 5
MAX_MISMATCH_PERCENT = 0.1
MAX_GAP_LENGTH = 3
ALIGN_CHUNKSIZE = 100_000
BP_CONTEXT_LENGTH = 8

TEMPLATE_SWITCHING_COLS = ['read_id',
						'fivep_sites',
						'temp_switch_sites',
						'read_seq', 
						'fivep_seq',
						'genomic_bp_context',
						'read_bp_pos',]

FINAL_INFO_TABLE_COLS = ['read_id', 
						'read_is_reverse', 
						'read_seq', 
						'chrom', 
						'strand', 
						'align_start',
						'align_end', 
						'align_is_reverse',
					  	'quality',
						'fivep_pos', 
						'bp_pos', 
						'read_bp_pos',
						'read_bp_nt', 
						'genomic_bp_context', 
						'genomic_bp_nt', 
						'threep_pos', 
						'bp_dist_to_threep',
						'gene_id', 
						]

filtered_out_lock = mp.Lock()
failed_out_lock = mp.Lock()
temp_switch_lock = mp.Lock()
circulars_lock = mp.Lock()



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def parse_intron_info(ref_introns):
	introns_df = pd.read_csv(ref_introns, sep='\t')
	introns_df['gene_id_dict'] = introns_df.gene_id.transform(lambda gene_id: {'gene_id': set(gene_id.split(','))})
	
	introns = {}
	for chrom in introns_df.chrom.unique():
		if chrom not in introns:
			introns[chrom] = {'+': IntervalTree(), '-': IntervalTree()}
		pos_subset = introns_df.loc[(introns_df.chrom==chrom) & (introns_df.strand=='+'), ['start', 'end', 'gene_id_dict']]
		introns[chrom]['+'] = IntervalTree.from_tuples([row for i, row in pos_subset.iterrows()])
		neg_subset = introns_df.loc[(introns_df.chrom==chrom) & (introns_df.strand=='-'), ['start', 'end', 'gene_id_dict']]
		introns[chrom]['-'] = IntervalTree.from_tuples([row for i, row in neg_subset.iterrows()])

	fivep_genes = introns_df.copy()
	fivep_genes['fivep_site'] = fivep_genes.apply(lambda row: f"{row['chrom']};{row['start']};{row['strand']}" if row['strand']=='+' else f"{row['chrom']};{row['end']-1};{row['strand']}", axis=1)
	fivep_genes = fivep_genes.groupby('fivep_site').agg({'gene_id': set}, as_index=False).gene_id.to_dict()

	return introns, fivep_genes	


def parse_tails(output_base:str, fivep_genes:dict):
	tails = pd.read_csv(TAILS_FILE.format(output_base), sep='\t', dtype={'fivep_chrom': 'category', 'strand': 'category'})
	tails = tails.drop(columns=['read_is_reverse'])

	# Unpack fivep_sites
	tails.fivep_sites = tails.fivep_sites.str.split(',')
	tails = tails.explode('fivep_sites')
	tails['gene_id'] = tails.apply(lambda row: tuple(fivep_genes[row['fivep_sites']]), axis=1)
	tails[['fivep_chrom', 'fivep_pos', 'strand']] = tails.fivep_sites.str.split(';', expand=True)
	tails.fivep_pos = tails.fivep_pos.astype(int)
	tails.fivep_chrom = tails.fivep_chrom.astype('category')
	tails.strand = tails.strand.astype('category')

	return tails


def cigar_str_to_tup(cigar:str) -> tuple[tuple[str, int]]:
	out = []
	length = ''
	for char in cigar:
		if char in CIGAR_OPERATORS:
			out.append((char, int(length)))
			length = ''
		else:
			length += char

	return tuple(out)

def infer_mismatches(tags) -> int:
	mismatches_tag = [tag for tag in tags.values if tag.startswith('XM:i:')][0]
	mismatches = int(mismatches_tag.lstrip('XM:i:'))
	return mismatches


def infer_read_bp(row:pd.Series) -> str:
	if not row['read_is_reverse'] and not row['align_is_reverse']:
		return row['head_seq'][-1]
	elif not row['read_is_reverse'] and row['align_is_reverse']:
		return functions.reverse_complement(row['head_seq'][0])
	elif row['read_is_reverse'] and not row['align_is_reverse']:
		return functions.reverse_complement(row['head_seq'][0])
	elif row['read_is_reverse'] and row['align_is_reverse']:
		return row['head_seq'][-1]


def parse_alignments_chunk(alignments_sam:str, chunk_start:int, chunk_end:int, n_aligns:int):
	# We get mismatch number from the XM tag, but it can end up 3rd, 4th, or 5th in the tags depending on whether or not ZS and YS are included
	# So we load all three columns it could be in and grab it
	# We grab the line directly before the starting line and the line directly after the ending line 
	# to check that we're keeping all of each read's alignments in 1 chunk 
	alignments = pd.read_csv(alignments_sam, 
								sep='\t',
								comment='@',
								header=None,
								usecols=[0, 1, 2, 3, 4, 5, 9, 13, 14, 15],
								names=['read_id', 'flag', 'chrom', 'align_start', 'quality', 'cigar', 'head_seq', 'tag1', 'tag2', 'tag3'],
								dtype={'read_id': 'string', 'chrom': 'category', 'align_start': 'UInt64', 'quality': 'UInt16'},
								skiprows=chunk_start-2,
								nrows=chunk_end-chunk_start+2,)
	
	# Same chunk start check as in filter_fivep_aligns.py, just implemented with pandas
	if chunk_start != 1:
		previous_line_rid = alignments.iloc[0]['read_id'][:-6]
		start_line_rid = alignments.iloc[1]['read_id'][:-6]
		if start_line_rid == previous_line_rid:
			alignments = alignments.loc[alignments.read_id.str.slice(0,-6)!=previous_line_rid].reset_index()

	# The chunk end check is also functionally the same
	# Well, it SHOUlD be
	if chunk_end != n_aligns:
		last_rid = alignments.iloc[-2]['read_id'][:-6]
		next_rid = alignments.iloc[-2]['read_id'][:-6]
		while last_rid == next_rid:
			next_line = pd.read_csv(alignments_sam, 
									sep='\t',
									comment='@',
									header=None,
									usecols=[0, 1, 2, 3, 4, 5, 9, 13, 14, 15],
									names=['read_id', 'flag', 'chrom', 'align_start', 'quality', 'cigar', 'head_seq', 'tag1', 'tag2', 'tag3'],
									dtype={'read_id': 'string', 'chrom': 'category', 'align_start': 'UInt64', 'quality': 'UInt16'},
									skiprows=chunk_end-chunk_start+2,
									nrows=1)
			next_rid = next_line.iloc[0]['read_id']
			alignments = pd.concat([alignments, next_line])
			chunk_end += 1
		alignments = alignments.iloc[:-1]
		
	# Convert 1-based inclusive to 0-based inclusive
	alignments.align_start = (alignments.align_start-1).astype('UInt64')
	# Infer some more info
	alignments['read_is_reverse'] = alignments.read_id.transform(lambda rid: {'_rev': True, '_for': False}[rid[-4:]]).astype('bool')
	alignments['length'] = alignments.head_seq.str.len()
	alignments['align_end'] = (alignments.align_start + alignments.length).astype('UInt64')
	alignments['align_is_reverse'] = alignments.flag.transform(functions.align_is_reverse)
	alignments.cigar = alignments.cigar.transform(cigar_str_to_tup)
	alignments['mismatches'] = alignments[['tag1', 'tag2', 'tag3']].agg(infer_mismatches, axis=1)
	alignments['mismatches_p'] = alignments.mismatches / alignments.length
	alignments['read_bp_nt'] = alignments.apply(infer_read_bp, axis=1).astype('category')
	# Discard uneeded columns
	alignments = alignments.drop(columns=['flag', 'tag1', 'tag2', 'tag3', 'length'])

	return alignments


def check_gaps(cigar:tuple) -> tuple:
	gaps = [length for op, length in cigar if op in ('I', 'D')]
	if len(gaps) > 1 or (len(gaps) == 1 and gaps[0] > MAX_GAP_LENGTH):
		return False
	else:
		return True


def get_bp_seqs(alignments:pd.DataFrame, genome_fasta:str):
	alignments = alignments.copy()
	alignments['window_start'] = pd.Series(alignments.bp_pos - BP_CONTEXT_LENGTH, dtype='str')
	alignments['window_end'] = pd.Series(alignments.bp_pos + BP_CONTEXT_LENGTH+1, dtype='str')
	alignments['align_num'] = alignments.index.to_series().astype(str)
	alignments['zero'] = '0'
	alignments['bedtools_line'] = alignments[['chrom', 'window_start', 'window_end', 'align_num', 'zero', 'strand']].agg('\t'.join, axis=1)

	bedtools_input = '\n'.join(alignments.bedtools_line) + '\n'
	# We can't parse the standard output for the sequences because warnings will be included in some lines
	# in a non-deterministic pattern 
	bp_seqs = functions.getfasta(genome_fasta, bedtools_input, log=None)
	bp_seqs.columns = ['align_num', 'genomic_bp_context']

	bp_seqs.genomic_bp_context = bp_seqs.genomic_bp_context.transform(lambda genomic_bp_context: genomic_bp_context.upper())
	
	alignments = pd.merge(alignments, bp_seqs, on='align_num')
	return alignments


def is_template_switch(row):
	bp_adj_seq = row['genomic_bp_context'][BP_CONTEXT_LENGTH+1:]
	base_matches = sum([bp_adj_seq[i]==row['fivep_seq'][i] for i in range(5)])
	return base_matches == 5


def	enveloping_introns(row, introns) -> list:
	if row['chrom'] not in introns:
		return []
	
	overlaps = introns[row['chrom']][row['strand']].overlap(row['align_start'], row['align_end'])
	envelops = [intron for intron in overlaps if intron.begin<=row['align_start'] and intron.end>=row['align_end']]

	return envelops


def filter_introns(row):
	matched_introns = []
	for intron in row['overlap_introns']:
		if len(intron.data['gene_id'].intersection(set(row['gene_id'])))>0:
			matched_introns.append(intron)

	return matched_introns


def add_nearest_threep(row:pd.Series):
	if row['strand'] == '+':
		candidate_introns = [intron for intron in row['overlap_introns'] if intron.end>row['bp_pos']]
		threep_pos = min(candidate_introns, key=lambda i: i.end-row['bp_pos']).end - 1
	else:
		candidate_introns = [intron for intron in row['overlap_introns'] if intron.begin<=row['bp_pos']]
		threep_pos = min(candidate_introns, key=lambda i: row['bp_pos']-i.begin).begin

	return threep_pos


def more_filters(row:pd.Series):
	# Check if bad alignment orientation combination, which do NOT leave the branchpoint adjacent to the 5'ss in the read as is expected of lariats
	if row['read_is_reverse'] is False:
		if row['align_is_reverse'] is False and row['strand']=='-':
			return 'wrong_orient'
		elif row['align_is_reverse'] is True and row['strand']=='+':
			return 'wrong_orient'
	if row['read_is_reverse'] is True:
		if row['align_is_reverse'] is False and row['strand']=='+':
			return 'wrong_orient'
		elif row['align_is_reverse'] is True and row['strand']=='-':
			return 'wrong_orient'

	# Check if the 5'ss is at or downstream of the tail's start
	if row['strand'] == '+' and row['fivep_pos'] > row['align_start']:
		return '5p_bp_order'
	if row['strand'] == '-' and row['fivep_pos'] < row['align_end']-1:
		return '5p_bp_order'

	return pd.NA


def drop_failed_alignments(alignments:pd.DataFrame, output_base:str) -> pd.DataFrame:
		# Get alignments that failed one of the filters
		failed_aligns = alignments.loc[alignments.filter_failed.notna()].copy()

		# Add in any missing cols as needed
		for col in FINAL_INFO_TABLE_COLS:
			if col not in failed_aligns.columns:
				failed_aligns[col] = ''
		
		failed_aligns.gene_id = failed_aligns.gene_id.transform(lambda gids: gids if isinstance(gids, str) else functions.str_join(gids))

		# Arrange cols and a write to file
		failed_aligns = failed_aligns[[*FINAL_INFO_TABLE_COLS, 'filter_failed']].drop_duplicates()
		with failed_out_lock:
			failed_aligns.to_csv(FAILED_HEADS_FILE.format(output_base), mode='a', sep='\t', header=False, index=False)

		# Remove failed alignments from DataFrame
		alignments = alignments.loc[alignments.filter_failed.isna()]
		return alignments


def filter_alignments_chunk(chunk_start, chunk_end, n_aligns, tails, introns, output_base, log_level) -> None:
	# We have to set the log level in each process because the children don't inherit the log level from their parent,
	# even if you pass the log object itself
	log = functions.get_logger(log_level)
	log.debug(f'Process {os.getpid()}: Born and assigned lines {chunk_start:,}-{chunk_end:,}')

	# Load in the assigned chunk of alignments, excluding skipping low-quality alignments
	alignments = parse_alignments_chunk(HEADS_TO_GENOME_FILE.format(output_base), chunk_start, chunk_end, n_aligns)

	# Filter out low-quality alignments
	pass_mismatch_filter = ((alignments.mismatches<=MAX_MISMATCHES) & (alignments.mismatches_p<=MAX_MISMATCH_PERCENT))
	pass_gap_filter = alignments.cigar.transform(check_gaps)
	alignments.loc[~pass_mismatch_filter, 'filter_failed'] = 'mismatches'
	alignments.loc[(alignments.filter_failed.isna()) & (~pass_gap_filter), 'filter_failed'] = 'gaps'
	alignments = drop_failed_alignments(alignments, output_base)
	if alignments.empty:
		log.debug(f'Process {os.getpid()}: Chunk exhausted after alignment quality filter')
		return 
	
	# Merge alignments with tails
	# This expands each alignment row into alignment-5'ss-combination rows
	alignments = pd.merge(alignments, tails, 'left', on=['read_id'])
	if alignments.fivep_pos.isna().any():
		raise RuntimeError(f"{alignments.fivep_pos.isna().sum()} alignments didn't match any tails read IDs, this shouldn't be possible")
	
	# Infer info
	alignments.fivep_pos = alignments.fivep_pos.astype('UInt64')
	alignments['bp_pos'] = alignments.apply(lambda row: row['align_end']-1 if row['strand']=='+' else row['align_start'], axis=1)
	alignments = get_bp_seqs(alignments, genome_fasta)
	alignments['genomic_bp_nt'] = alignments.genomic_bp_context.str.get(8)
	
	# Identify template-switching reads
	alignments['template_switching'] = alignments.apply(is_template_switch, axis=1)

	# Output template-switching reads
	temp_switches = alignments.loc[alignments.template_switching].copy()
	if not temp_switches.empty:
		temp_switches = temp_switches.astype(str)
		temp_switches.read_id = temp_switches.read_id.str.slice(0,-6)
		temp_switches['fivep_sites'] = temp_switches[['fivep_chrom', 'strand', 'fivep_pos']].agg(';'.join, axis=1)
		temp_switches['temp_switch_sites'] = temp_switches[['chrom', 'bp_pos']].agg(';'.join, axis=1)
		temp_switches = temp_switches[TEMPLATE_SWITCHING_COLS]
		temp_switches = temp_switches.groupby('read_id', as_index=False).agg({col: functions.str_join for col in temp_switches.columns if col != 'read_id'})
		with temp_switch_lock:
			temp_switches.to_csv(TEMP_SWITCH_FILE.format(output_base), mode='a', sep='\t', header=False, index=False)

	# Filter out template-switching reads
	alignments = alignments.loc[~alignments.template_switching].drop(columns='template_switching')
	if alignments.empty:
		log.debug(f'Process {os.getpid()}: Chunk exhausted after template-switching filter')
		return 
	
	# Identify introns that envelop the alignment
	# alignments['overlap_introns'] = alignments.apply(lambda row: introns[row['chrom']][row['strand']].overlap(row['align_start'], row['align_end']), axis=1)
	# alignments['overlap_introns'] = alignments.apply(lambda row: [intron for intron in row['overlap_introns'] if intron.begin<=row['align_start'] and intron.end>=row['align_end']], axis=1)
	alignments['overlap_introns'] = alignments.apply(enveloping_introns, introns=introns, axis=1)
	
	# Filter out alignments that don't overlap any introns
	alignments.loc[alignments.overlap_introns.transform(len)==0, 'filter_failed'] = 'overlap_introns'
	alignments = drop_failed_alignments(alignments, output_base)
	if alignments.empty:
		log.debug(f'Process {os.getpid()}: Chunk exhausted after overlap_introns filter')
		return 
	
	# Filter out alignments where 5'ss and BP segments aren't in the same gene
	# alignments = alignments.explode('gene_id')
	# alignments.overlap_introns = alignments.apply(lambda row: tuple(intron for intron in row['overlap_introns'] if row['gene_id'] in intron.data['gene_id']), axis=1)
	alignments.overlap_introns = alignments.apply(filter_introns, axis=1)
	alignments.gene_id = alignments.gene_id.transform(functions.str_join)
	alignments.loc[alignments.overlap_introns.transform(len).eq(0), 'filter_failed'] = 'fivep_intron_match'
	alignments = drop_failed_alignments(alignments, output_base)
	if alignments.empty:
		log.debug(f'Process {os.getpid()}: Chunk exhausted after fivep_intron_match filter')
		return 
	
	# Infer more info
	alignments['threep_pos'] = alignments.apply(add_nearest_threep, axis=1)
	alignments['bp_dist_to_threep'] = alignments.apply(lambda row: -abs(row['bp_pos']-row['threep_pos']) if pd.notna(row['threep_pos']) else pd.NA, axis=1)
	
	# Filter alignments based on proper read orientation and 5'ss-BP ordering
	alignments['filter_failed'] = alignments.apply(more_filters, axis=1, result_type='reduce')
	alignments = drop_failed_alignments(alignments, output_base)
	if alignments.empty:
		log.debug(f'Process {os.getpid()}: Chunk exhausted after final filters')
		return 
	
	# Identify circularized intron reads
	alignments['circular'] = alignments.bp_dist_to_threep.isin((0, -1, -2))

	# Output circularized intron reads
	circulars = alignments.loc[alignments.circular].copy()
	if not circulars.empty:
		circulars = circulars.astype(str)
		circulars.read_id = circulars.read_id.str.slice(0,-6)
		circulars['fivep_sites'] = circulars[['fivep_chrom', 'strand', 'fivep_pos']].agg(functions.str_join, axis=1)
		circulars = circulars[FINAL_INFO_TABLE_COLS]
		circulars = circulars.groupby('read_id', as_index=False).agg({col: functions.str_join for col in circulars.columns if col != 'read_id'})
		with circulars_lock:
			circulars.to_csv(CIRCULARS_FILE.format(output_base), mode='a', sep='\t', header=False, index=False)

	# Filter out circularized intron reads
	alignments = alignments.loc[~alignments.circular]
	if alignments.empty:
		log.debug(f'Process {os.getpid()}: Chunk exhausted after circularized intron filter')
		return 

	# Drop all the uneeded columns
	alignments = alignments[FINAL_INFO_TABLE_COLS]
	
	# Some reads get mapped to coordinates with multiple overlapping gene annotations
	# We resolve this by collapsing the duplicated rows and concatenating the gene_id column
	alignments = (alignments.groupby([col for col in alignments.columns if col != 'gene_id'], as_index=False, observed=True)
									.gene_id
									.agg(functions.str_join)
					)

	with filtered_out_lock:
		alignments.to_csv(PUTATITVE_LARIATS_FILE.format(output_base), mode='a', sep='\t', index=False, header=False)

	log.debug(f'Process {os.getpid()}: Finished')


# If any processes hit an error, raise an error
def error_callback(exception):
	raise exception



# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	# Get args
	threads, ref_introns, genome_fasta, output_base, log_level = sys.argv[1:]

	# Get logger
	log = functions.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}')

	threads = int(threads)
	
	with open(HEADS_TO_GENOME_FILE.format(output_base)) as sam:
		n_aligns = sum(1 for _ in sam)
	log.debug(f'{n_aligns:,} head alignments')
	chunk_ranges = [[chunk_start, chunk_start+ALIGN_CHUNKSIZE] for chunk_start in range(1, n_aligns+1, ALIGN_CHUNKSIZE)]
	chunk_ranges[-1][-1] = n_aligns
	log.debug(f'chunk_ranges: {chunk_ranges}')

	if n_aligns == 0:
		log.info('No reads remaining')
		exit()

	# Load reference data for processing alignments
	introns, fivep_genes = parse_intron_info(ref_introns)
	tails = parse_tails(output_base, fivep_genes)

	# Record # of reads that passed fivep_filtering
	count = tails.read_id.str.rstrip('/1').str.rstrip('/2').nunique()
	with open(RUN_DATA_FILE.format(output_base), 'a') as a:
		a.write(f'fivep_filter_passed\t{count}\n')

	# Write headers for the outfiles
	# The rows will get appended in chunks during filter_alignments_chunk()
	with open(PUTATITVE_LARIATS_FILE.format(output_base), 'w') as w:
		w.write('\t'.join(FINAL_INFO_TABLE_COLS) + '\n')
	with open(FAILED_HEADS_FILE.format(output_base), 'w') as w:
		w.write('\t'.join(FINAL_INFO_TABLE_COLS) + '\tfilter_failed' + '\n')
	with open(TEMP_SWITCH_FILE.format(output_base), 'w') as w:
		w.write('\t'.join(TEMPLATE_SWITCHING_COLS) + '\n')
	with open(CIRCULARS_FILE.format(output_base), 'w') as w:
		w.write('\t'.join(FINAL_INFO_TABLE_COLS) + '\n')

	# multiprocessing won't run correctly with just 1 chunk for some reason
	if len(chunk_ranges) == 1:
		log.debug(f'Parallel processing {len(chunk_ranges):,} chunks...')
		filter_alignments_chunk(chunk_ranges[0][0], chunk_ranges[0][1], n_aligns, tails, introns, output_base, log_level)
	else:
		log.debug(f'Parallel processing {len(chunk_ranges):,} chunks...')
		pool = mp.Pool(processes=threads)
		for chunk_start, chunk_end in chunk_ranges:
			pool.apply_async(filter_alignments_chunk, 
							args=(chunk_start, chunk_end, n_aligns, tails, introns, output_base, log_level,), 
							error_callback=error_callback)
		
		# Don't create any more processes
		pool.close()
		# Wait until all processes are finished
		pool.join()

	log.debug('End of script')
