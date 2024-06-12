import multiprocessing.shared_memory
from dataclasses import dataclass, field
import sys
import os
import gzip
import subprocess
import multiprocessing
import time
import collections
import logging

from intervaltree import Interval, IntervalTree
import pandas as pd
import numpy as np
# import psutil

import functions



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
CIGAR_OPERATORS = ('M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X')
MAX_MISMATCHES = 5
MAX_MISMATCH_PERCENT = 0.1
MAX_GAP_LENGTH = 3
ALIGN_CHUNKSIZE = 1_000_000

ALIGN_INITIAL_COLS = ['read_id', 
					  'read_is_reverse',
					  'chrom', 
					  'align_start', 
					  'align_end', 
					  'align_is_reverse', 
					  'read_bp_nt']

FIVEP_INFO_COLS = ['read_id',
				'read_seq',
				'fivep_seq',
				'fivep_sites',
				'read_fivep_start',
				'read_fivep_end',
				'read_bp_pos',
				'gene_id',
				'fivep_chrom',
				'fivep_pos',
				'strand']

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

filtered_out_lock = multiprocessing.Lock()
failed_out_lock = multiprocessing.Lock()
temp_switch_lock = multiprocessing.Lock()



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def parse_intron_info(ref_introns):
	introns_df = pd.read_csv(ref_introns, sep='\t')
	introns_df['gene_id_dict'] = introns_df.gene_id.transform(lambda gene_id: {'gene_id': gene_id})
	
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


def parse_fivep_info(output_base:str, fivep_genes:dict):
	fivep_info = pd.read_csv(f'{output_base}fivep_info_table.tsv', sep='\t', dtype={'fivep_chrom': 'category', 'strand': 'category'})
	fivep_info = fivep_info.drop(columns=['read_is_reverse'])

	# Unpack fivep_sites
	fivep_info.fivep_sites = fivep_info.fivep_sites.str.split(',')
	fivep_info = fivep_info.explode('fivep_sites')
	fivep_info['gene_id'] = fivep_info.apply(lambda row: tuple(fivep_genes[row['fivep_sites']]), axis=1)
	fivep_info[['fivep_chrom', 'fivep_pos', 'strand']] = fivep_info.fivep_sites.str.split(';', expand=True)
	fivep_info.fivep_pos = fivep_info.fivep_pos.astype(int)
	fivep_info.fivep_chrom = fivep_info.fivep_chrom.astype('category')
	fivep_info.strand = fivep_info.strand.astype('category')

	return fivep_info


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


def pass_gap_filter(cigar:tuple) -> tuple:
	gaps = [length for op, length in cigar if op in ('I', 'D')]
	if len(gaps) > 1 or (len(gaps) == 1 and gaps[0] > MAX_GAP_LENGTH):
		return False
	else:
		return True


def infer_mismatches(tags) -> int:
	mismatches_tag = [tag for tag in tags.values if tag.startswith('XM:i:')][0]
	mismatches = int(mismatches_tag.lstrip('XM:i:'))
	return mismatches


def infer_read_bp(row:pd.Series) -> str:
	if not row['read_is_reverse'] and not row['align_is_reverse']:
		return row['trim_seq'][-1]
	elif not row['read_is_reverse'] and row['align_is_reverse']:
		return functions.reverse_complement(row['trim_seq'][0])
	elif row['read_is_reverse'] and not row['align_is_reverse']:
		return functions.reverse_complement(row['trim_seq'][0])
	elif row['read_is_reverse'] and row['align_is_reverse']:
		return row['trim_seq'][-1]


def parse_alignments_chunk(alignments_sam:str, chunk_start:int):
	# We get mismatch number from the XM tag, but it can end up 3rd, 4th, or 5th in the tags depending on whether or not ZS and YS are included
	# So we load all three columns it could be in and grab it
	alignments = pd.read_csv(alignments_sam, 
								sep='\t',
								comment='@',
								header=None,
								usecols=[0, 1, 2, 3, 5, 9, 13, 14, 15],
								names=['read_id', 'flag', 'chrom', 'align_start', 'cigar', 'trim_seq', 'tag1', 'tag2', 'tag3'],
								dtype={'read_id': 'string', 'chrom': 'category', 'align_start': 'UInt64'},
								skiprows=chunk_start-1,
								nrows=ALIGN_CHUNKSIZE,)

	# Convert 1-based inclusive to 0-based inclusive
	alignments.align_start = (alignments.align_start-1).astype('UInt64')
	# Infer some more info
	alignments['read_is_reverse'] = alignments.read_id.transform(lambda rid: {'_rev': True, '_for': False}[rid[-4:]]).astype('bool')
	alignments['length'] = alignments.trim_seq.str.len()
	alignments['align_end'] = (alignments.align_start + alignments.length).astype('UInt64')
	alignments['align_is_reverse'] = alignments.flag.transform(functions.align_is_reverse)
	alignments.cigar = alignments.cigar.transform(cigar_str_to_tup)
	alignments['mismatches'] = alignments[['tag1', 'tag2', 'tag3']].agg(infer_mismatches, axis=1)
	alignments['mismatches_p'] = alignments.mismatches / alignments.length
	alignments['read_bp_nt'] = alignments.apply(infer_read_bp, axis=1).astype('category')

	# Exclude alignments which are too low-quality 
	alignments['pass_gap_filter'] = alignments.cigar.transform(pass_gap_filter)
	alignments['pass_mismatch_filter'] = ((alignments.mismatches<=MAX_MISMATCHES) & (alignments.mismatches_p<=MAX_MISMATCH_PERCENT))
	alignments = alignments.loc[(alignments.pass_gap_filter) & (alignments.pass_mismatch_filter)]
	# Discard uneeded columns
	alignments = alignments[ALIGN_INITIAL_COLS]

	return alignments


def get_bp_seqs(alignments:pd.DataFrame, genome_fasta:str):
	alignments = alignments.copy()
	alignments['window_start'] = alignments.apply(lambda row: row['bp_pos']-4 if row['strand']=='+' else row['bp_pos']-5, axis=1).astype(str)
	alignments['window_end'] = alignments.apply(lambda row: row['bp_pos']+6 if row['strand']=='+' else row['bp_pos']+5, axis=1).astype(str)
	alignments['align_num'] = alignments.index.to_series().astype(str)
	alignments['zero'] = '0'
	alignments['bedtools_line'] = alignments[['chrom', 'window_start', 'window_end', 'align_num', 'zero', 'strand']].agg('\t'.join, axis=1)
	
	bedtools_input = '\n'.join(alignments.bedtools_line) + '\n'
	bedtools_call = f'bedtools getfasta -nameOnly -s -tab -fi {genome_fasta} -bed -'
	bedtools_lines = subprocess.run(bedtools_call.split(' '), input=bedtools_input, check=True, capture_output=True, text=True).stdout.strip().split('\n')
	
	bp_seqs = pd.DataFrame([l.split('\t') for l in bedtools_lines], columns=['align_num', 'genomic_bp_context'])
	bp_seqs.align_num = bp_seqs.align_num.str.slice(0,-3)
	bp_seqs.genomic_bp_context = bp_seqs.genomic_bp_context.transform(lambda genomic_bp_context: genomic_bp_context.upper())
	
	alignments = pd.merge(alignments, bp_seqs, on='align_num')
	return alignments


def is_template_switch(row):
	bp_adj_seq = row['genomic_bp_context'][5:]
	base_matches = sum([bp_adj_seq[i]==row['fivep_seq'][i] for i in range(5)])
	return base_matches == 5


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


def drop_failed_alignments(alignments:pd.DataFrame, output_base:str) -> pd.DataFrame:
		alignments = alignments.copy()

		# Get alignments that failed one of the filters
		failed_alignments = alignments.loc[alignments.filter_failed.notna()].copy()

		# Add in any missing cols as needed
		for col in FAILED_ALIGNMENTS_COLS:
			if col not in failed_alignments.columns:
				failed_alignments[col] = ''
		
		# Arrange cols and a write to file
		failed_alignments = failed_alignments[FAILED_ALIGNMENTS_COLS].drop_duplicates()
		with failed_out_lock:
			failed_alignments.to_csv(f'{output_base}failed_trimmed_alignments.tsv', mode='a', sep='\t', header=False, index=False)

		# Remove failed alignments from DataFrame
		alignments = alignments.loc[alignments.filter_failed.isna()]
		return alignments


def filter_alignments_chunk(chunk_start, fivep_info, introns_shared, output_base, log) -> None:
	log.debug(f'Process {os.getpid()} created, assigned lines {chunk_start:,}-{chunk_start+ALIGN_CHUNKSIZE-1:,}')

	# Load in the assigned chunk of alignments, excluding skipping low-quality alignments
	alignments = parse_alignments_chunk(f'{output_base}trimmed_reads_to_genome.sam', chunk_start)
	mem_mb = alignments.memory_usage(index=False, deep=True).divide(1_000_000).sum()
	log.debug(f'Process {os.getpid()}: Alignments memory use is {mem_mb:.3} MB')
	# log.debug(f'Current memory use is {psutil.virtual_memory().used/1_000_000:0.3} MB ({psutil.virtual_memory().percent}%)')

	# Merge alignments with fivep_info
	# This expands each alignment row into alignment-5'ss-combination rows
	alignments = pd.merge(alignments, fivep_info, 'left', on=['read_id'])
	if alignments.fivep_pos.isna().any():
		raise RuntimeError(f"{alignments.fivep_pos.isna().sum()} alignments didn't match any fivep_info_table read IDs, this shouldn't be possible")
	
	# Infer info
	alignments.fivep_pos = alignments.fivep_pos.astype('UInt64')
	alignments['bp_pos'] = alignments.apply(lambda row: row['align_end']-1 if row['strand']=='+' else row['align_start'], axis=1)
	alignments = get_bp_seqs(alignments, genome_fasta)
	alignments['genomic_bp_nt'] = alignments.genomic_bp_context.str.get(4)
	
	# Identify template-switching reads
	alignments['template_switching'] = alignments.apply(is_template_switch, axis=1)

	# Output template-switching reads
	temp_switches = alignments.loc[alignments.template_switching].copy()
	if not temp_switches.empty:
		temp_switches = temp_switches.astype(str)
		temp_switches['fivep_sites'] = temp_switches[['fivep_chrom', 'strand', 'fivep_pos']].agg(';'.join, axis=1)
		temp_switches['temp_switch_sites'] = temp_switches[['chrom', 'bp_pos']].agg(';'.join, axis=1)
		temp_switches = temp_switches[TEMPLATE_SWITCHING_COLS]
		with temp_switch_lock:
			temp_switches.to_csv(f'{output_base}template_switching_reads.tsv', mode='a', sep='\t', header=False, index=False)

	# Filter out template-switching reads
	alignments = alignments.loc[~alignments.template_switching].drop(columns='template_switching')
	if alignments.empty:
		return 
	
	# Identify introns that overlap the alignment
	alignments['overlap_introns'] = alignments.apply(lambda row: introns_shared[row['chrom']][row['strand']].overlap(row['align_start'], row['align_end']), axis=1)
	
	# Filter out alignments that don't overlap any introns
	alignments.loc[alignments.overlap_introns.transform(len)==0, 'filter_failed'] = 'overlap_introns'
	alignments = drop_failed_alignments(alignments, output_base)
	if alignments.empty:
		return 
	
	# Filter out alignments where 5'ss and BP segments aren't in the same gene
	alignments = alignments.explode('gene_id')
	alignments.overlap_introns = alignments.apply(lambda row: tuple(intron for intron in row['overlap_introns'] if row['gene_id'] in intron.data['gene_id']), axis=1)
	alignments.loc[alignments.overlap_introns.transform(len).eq(0), 'filter_failed'] = 'fivep_intron_match'
	alignments = drop_failed_alignments(alignments, output_base)
	if alignments.empty:
		return 
	
	# Filter alignments based on proper read orientation and 5'ss-BP ordering
	alignments['filter_failed'] = alignments.apply(final_filters, axis=1, result_type='reduce')
	alignments = drop_failed_alignments(alignments, output_base)
	if alignments.empty:
		return 
	
	# Infer more info
	alignments['threep_pos'] = alignments.apply(add_nearest_threep, axis=1)
	alignments['bp_dist_to_threep'] = alignments.apply(lambda row: -abs(row['bp_pos']-row['threep_pos']) if pd.notna(row['threep_pos']) else pd.NA, axis=1)

	# Drop all the uneeded columns
	alignments = alignments[FINAL_INFO_TABLE_COLS]
	
	# Some reads get mapped to coordinates with multiple overlapping gene annotations
	# We resolve this by collapsing the duplicated rows and concatenating the gene_id, gene_name, and gene_type columns
	# Code adapted from https://stackoverflow.com/questions/27298178/concatenate-strings-from-several-rows-using-pandas-groupby
	alignments = (alignments.groupby([col for col in alignments.columns if col != 'gene_id'], as_index=False, observed=True)
									.gene_id
									.agg(functions.comma_join)
					)

	with filtered_out_lock:
		alignments.to_csv(f'{output_base}trimmed_info_table.tsv', mode='a', sep='\t', index=False, header=False)

	log.debug(f'Process {os.getpid()}: Chunk finished')
	# log.debug(f'Current memory use is {psutil.virtual_memory().used/1_000_000:0.3} MB ({psutil.virtual_memory().percent:%})')



# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	# Get logger
	log = logging.getLogger()
	# log.setLevel('DEBUG')
	handler = logging.StreamHandler(sys.stdout)
	handler.setLevel('DEBUG')
	log.addHandler(handler)

	# Get args
	threads, ref_introns, genome_fasta, output_base = sys.argv[1:]
	log.debug(f'Args recieved: {sys.argv[1:]}')
	threads = int(threads)
	
	with open(f'{output_base}trimmed_reads_to_genome.sam') as sam:
		n_aligns = sum(1 for _ in sam)
	log.debug(f'{n_aligns:,} trimmed alignments')
	chunk_starts = [start for start in range(1, n_aligns+1, ALIGN_CHUNKSIZE)]
	log.debug(f'chunk starts: {chunk_starts}')

	if n_aligns == 0:
		print(time.strftime('%m/%d/%y - %H:%M:%S') + '| No reads remaining')
		exit()

	# Load reference data for processing alignments
	introns, fivep_genes = parse_intron_info(ref_introns)
	fivep_info = parse_fivep_info(output_base, fivep_genes)

	# Record # of reads that passed fivep_filtering
	count = fivep_info.read_id.str.rstrip('/1').str.rstrip('/2').nunique()
	with open(f'{output_base}read_counts.tsv', 'a') as a:
		a.write(f'fivep_filter_passed\t{count}\n')

	# Write headers for the failed_trimmed_alignments and template_switching_alignments out-files
	# The rows will get appended in chunks during filter_alignments_chunk()
	with open(f'{output_base}trimmed_info_table.tsv', 'w') as w:
		w.write('\t'.join(FINAL_INFO_TABLE_COLS) + '\n')
	with open(f'{output_base}template_switching_reads.tsv', 'w') as w:
		w.write('\t'.join(TEMPLATE_SWITCHING_COLS) + '\n')
	with open(f'{output_base}failed_trimmed_alignments.tsv', 'w') as w:
		w.write('\t'.join(FAILED_ALIGNMENTS_COLS) + '\n')

	# 
	log.debug('Processing')
	if len(chunk_starts) ==1:
		filter_alignments_chunk(chunk_starts[0], fivep_info, introns, output_base, log)
	else:
		pool = multiprocessing.Pool(processes=threads)
		for start in chunk_starts:
			pool.apply_async(filter_alignments_chunk, args=(start, fivep_info, introns, output_base, log,))
		
		# Don't create any more processes
		pool.close()
		# Wait until all processes are finished
		pool.join()

	# Collapse the template-switching reads rows so each row is one read
	# We have to do this at the end because alignments to the same read can end up in different alignment chunks,
	# as can alignments to different mates from the same read if paired-end sequencing data
	log.debug('Collapsing temp switch')
	temp_switches = pd.read_csv(f'{output_base}template_switching_reads.tsv', sep='\t', dtype='string')
	temp_switches['read_id'] = temp_switches.read_id.transform(lambda rid: rid[:-4].split('/')[0])
	temp_switches = (temp_switches
								.groupby(['read_id'], as_index=False)
								.agg({col: functions.comma_join for col in TEMPLATE_SWITCHING_COLS if col != 'read_id'})
							)
	temp_switches.to_csv(f'{output_base}template_switching_reads.tsv', sep='\t', index=False)

	log.info('End of script')

