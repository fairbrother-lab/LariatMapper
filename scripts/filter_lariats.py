import sys
import gzip
import random
import os
from collections import Counter, defaultdict
import subprocess

import pandas as pd
import numpy as np
from intervaltree import Interval, IntervalTree

import utils





# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
# In files
OUTPUT_SJ_FILE = "{}SJ.out.tab"
TEMPLATE_SWITCH_FILE = "{}template_switching_reads.tsv"
CIRCULARS_FILE = "{}circularized_intron_reads.tsv"
PUTATITVE_LARIATS_FILE = "{}putative_lariats.tsv"
# Out files
FAILED_LARIATS_FILE = "{}failed_lariat_alignments.tsv"
LARIATS_FILE = "{}lariat_reads.tsv"
RUN_DATA_FILE = "{}read_counts.tsv"


# =============================================================================#
#                                  Functions                                   #
# =============================================================================#

def parse_genes(ref_genes:str, log) -> dict[str, dict[str, IntervalTree]]:
	'''
	Load gene coordinates
	'''
	genes_df, genes = pd.read_csv(ref_genes, sep='\t'), {}
	for _, gene_info in genes_df.iterrows():
		gene_id, chrom, strand, start, end = [gene_info[c] for c in ['gene_id', 'chrom', 'strand', 'start', 'end']]
		if chrom not in genes:
			genes[chrom] = {'+': IntervalTree(), '-': IntervalTree()}
		genes[chrom][strand].add(Interval(start, end, {'gene_id':gene_id}))
	
	return genes

def parse_introns(ref_introns:str, log) -> dict[str, dict[str, IntervalTree]]:
	'''
	Load intron coordinates and genes
	'''
											
	introns_df = pd.read_csv(ref_introns, sep='\t', na_filter=False)
	introns_df['gene_id'] = introns_df.gene_id.transform(lambda gene_id: set(gene_id.split(',')))

	chromosomes = set(introns_df.chrom.unique())
	introns = {chrom:{'+': IntervalTree(), '-': IntervalTree()} for chrom in chromosomes}
	for _, intron_info in introns_df.iterrows():
		chrom, strand, start, end = [intron_info[c] for c in ['chrom', 'strand', 'start', 'end']]
		intron_id = f'{chrom}_{strand}_{start}_{end}'
		introns[chrom][strand].add(Interval(start, end, {'gene_id': intron_info['gene_id'], 'intron_id':intron_id}))
	
	return introns


#TODO: Implement a more elaborate check. We only need to check for UBC-type repeats as possible false-positives
def check_repeat_overlap(lariat_reads: pd.DataFrame, ref_repeatmasker:str, log) -> set:
	''' 
    Check if both the 5'SS and the BP overlap with a repetitive region
    '''
	if not os.path.isfile(ref_repeatmasker):
		log.info('Repeatmasker file not found, skipping repeats check...')
		return set()

	# Write the 5'ss and BP coordinates to bedtools input strings
	bedtools_fivep_input, bedtools_bp_input = '', ''
	for _, row in lariat_reads.iterrows():
		bedtools_fivep_input += f"{row['chrom']}\t{row['fivep_pos']-1}\t{row['fivep_pos']+1}\t{row['read_id']}\n"
		bedtools_bp_input += f"{row['chrom']}\t{row['bp_pos']-1}\t{row['bp_pos']+1}\t{row['read_id']}\n"

	# Identify 5'ss that overlap a repeat region
	bedtools_call = f'bedtools intersect -u -a - -b {ref_repeatmasker}'
	bedtools_fivep_output = utils.run_command(bedtools_call, input=bedtools_fivep_input, log=log).split('\n')
	bedtools_bp_output = utils.run_command(bedtools_call, input=bedtools_bp_input, log=log).split('\n')

	# Add reads where both sites overlapped to repeat_rids for removal
	fivep_repeat_rids, bp_repeat_rids = set(), set()
	for fivep_line in bedtools_fivep_output:
		fivep_repeat_rids.add(fivep_line.split('\t')[-1])
	for bp_line in bedtools_bp_output:
		bp_repeat_rids.add(bp_line.split('\t')[-1])
	repeat_rids = fivep_repeat_rids.intersection(bp_repeat_rids)

	return repeat_rids


def filter_lariats(row:pd.Series, repeat_rids:set, temp_switch_rids:set, circular_rids:set):
	'''
	Filter the candidate lariat reads to EXCLUDE any that meet the following criteria:
			- Read maps to UBB or UBC (likely false positive due to the repetitive nature of the genes)
			- Both the 5'SS and the BP overlap with repetitive regions from RepeatMasker (likely false positive due to sequence repetition)
			- BP is within 2bp of a splice site (likely an intron circle, not a lariat)
	'''

	if row['read_id'] in temp_switch_rids:
		return 'template_switching'
	
	if row['read_id'] in circular_rids:
		return 'circularized_intron'

	if row['read_id'] in repeat_rids:
		return 'in_repeat'
	
	return np.nan


def parse_line(sam_line:str):
	alignment_info = sam_line.split('\t')
	read_id, flag, chrom, start_pos, _, cigar_str = alignment_info[:6]
	tags = [tag_info.split(':') for tag_info in alignment_info[11:]]
	tags = {tag_info[0]:tag_info[2] for tag_info in tags}
	alignment_score, mismatch_count = int(tags['AS']), int(tags['NM'])

	cigar_tuples = utils.get_cigar_tuples(cigar_str)
	start_pos = int(start_pos)-1
	end_pos = utils.get_align_end(start_pos, cigar_tuples)
	gaps = [length for length, op in cigar_tuples if op in 'ID']
	query_length = sum([length for length, op in cigar_tuples if op in 'MIS=X'])
	mismatch_p = mismatch_count/float(query_length)
	read_orient_to_chrom = utils.align_orient(flag)

	return read_id, chrom, start_pos, end_pos, read_orient_to_chrom, cigar_str, mismatch_count, mismatch_p, gaps

def verify_fivep_alignment(lariat_reads:pd.DataFrame, output_base:str, hisat2_index:str, genes:dict, threads:str, log):
	'''
	Perform full alignment of tail segment to confirm high-quality mapping to the 5'ss identified
	when aligning 5'ss sequences to the unmapped reads 
	'''

	tmp_fasta = f"{output_base}tail_full_seq.tmp.fa"
	read_info_dict = {}
	with open(tmp_fasta, 'w') as out_file:
		for _, read_info in lariat_reads.iterrows():
			if pd.isna(read_info['filter_failed']):
				fivep_id = f"{read_info['chrom']};{read_info['fivep_pos']};{read_info['strand']}"
				tmp_read_id = f"{read_info['read_id']}|{fivep_id}"
				tmp_read_seq = read_info['read_seq_forward'][(read_info['read_bp_pos']+1):]
				out_file.write(f'>{tmp_read_id}\n{tmp_read_seq}\n')
				read_info_dict[tmp_read_id] = read_info
		
	hisat2_cmd = f"hisat2 --no-softclip --no-spliced-alignment --no-unal --sensitive -k 5 --threads {threads} -f -x {hisat2_index} -U {tmp_fasta}"
	hisat2_output = subprocess.run(hisat2_cmd, shell=True, capture_output=True, text=True).stdout.strip().split('\n')
	passed_alignments = []
	read_failures = defaultdict(list)
	for sam_line in hisat2_output:
		if sam_line[0] != '@':
			tmp_read_id, chrom, start_pos, end_pos, read_orient_to_chrom, cigar_str, mismatch_count, mismatch_p, gaps = parse_line(sam_line)
			read_id, fivep_id = tmp_read_id.split('|')
			fivep_chrom, fivep_pos, fivep_strand = fivep_id.split(';')

			if mismatch_p > 0.1:
				read_failures[read_id].append('mismatch_rate')
				continue

			if len(gaps) > 1:
				read_failures[read_id].append('multiple_gaps')
				continue

			if sum(gaps) > 3:
				read_failures[read_id].append('long_gap')
				continue

			if fivep_chrom != chrom or (fivep_strand == '+' and read_orient_to_chrom == 'Reverse') or (fivep_strand == '-' and read_orient_to_chrom == 'Forward'):
				read_failures[read_id].append('chrom_strand_mismatch')
				continue
			
			align_genes = set([g.data['gene_id'] for g in genes[chrom][fivep_strand].overlap(start_pos, end_pos)])
			orig_genes = set(read_info_dict[tmp_read_id]['gene_id'].split(','))
			overlap_genes = align_genes.intersection(orig_genes)
			if len(overlap_genes) == 0:
				read_failures[read_id].append('gene_mismatch')
				continue

			full_align_fivep_pos = start_pos if fivep_strand == '+' else end_pos-1
			if full_align_fivep_pos != int(fivep_pos):
				read_failures[read_id].append('fivep_align_pos_mismatch')
				continue

			gaps = '0' if len(gaps)==0 else str(gaps[0])

			passed_alignments.append([read_id, chrom, fivep_strand, int(fivep_pos), start_pos, end_pos, cigar_str, mismatch_count, gaps])

	passed_rids = set([a[0] for a in passed_alignments])
	failure_score = {'mismatch_rate':0, 'multiple_gaps':1, 'long_gap':2,
				     'chrom_strand_mismatch':3, 'gene_mismatch':4, 'fivep_align_pos_mismatch':5}
	failed_alignments = []
	for i, read_info in lariat_reads.iterrows():
		read_id, chrom, strand, fivep_pos, filter_failed = [read_info[k] for k in ['read_id', 'chrom', 'strand', 'fivep_pos', 'filter_failed']]
		if pd.isna(filter_failed) and not read_id in passed_rids:
			filter_failed_new = max(read_failures[read_id], key=lambda f:failure_score[f]) if read_id in read_failures else 'no_alignment'
			lariat_reads.at[i, 'filter_failed'] = f'tail_align_{filter_failed_new}'
			failed_alignments.append([read_id, chrom, strand, fivep_pos, 'NA', 'NA', 'NA', 'NA', 'NA'])
	
	all_alignments = pd.DataFrame(passed_alignments+failed_alignments, columns=['read_id', 'chrom', 'strand', 'fivep_pos', 'tail_align_start_pos', 'tail_align_end_pos',
																			 	'tail_cigar_str', 'tail_mismatches', 'tail_gaps'])
	
	lariat_reads = pd.merge(lariat_reads, all_alignments, on=['read_id', 'chrom', 'strand', 'fivep_pos'], how='right')

	os.remove(tmp_fasta)

	return lariat_reads

def filter_nested_reads(lariat_reads:pd.DataFrame, log):
	'''
	Filter reads where the head and tail alignments are nested within each other
	(potential template-switching event)
	'''
	for i, read_info in lariat_reads.iterrows():
		if pd.isna(read_info['filter_failed']):
			head_align_start_pos, head_align_end_pos, tail_align_start_pos, tail_align_end_pos = [read_info[c] for c in ['head_align_start_pos', 'head_align_end_pos', 'tail_align_start_pos', 'tail_align_end_pos']]

			head_in_tail = head_align_start_pos >= tail_align_start_pos and head_align_end_pos <= tail_align_end_pos
			tail_in_head = tail_align_start_pos >= head_align_start_pos and tail_align_end_pos <= head_align_end_pos
			if head_in_tail or tail_in_head:
				lariat_reads.at[i, 'filter_failed'] = f'nested_read_segments'
			
	return lariat_reads


def filter_align_score(lariat_reads:pd.DataFrame):

	lariat_reads.read_id = lariat_reads.read_id.str.slice(0,-4)
	lariat_reads[['read_id', 'read_num']] = lariat_reads.read_id.str.split('/', expand=True)
	lariat_reads.read_num = lariat_reads.read_num.astype(int)
	lariat_reads.insert(1, 'read_num', lariat_reads.pop('read_num'))

	for i, read_info in lariat_reads.iterrows():
		rid_aligns = lariat_reads[lariat_reads['read_id']==read_info['read_id']]
		rid_score_min_diff = (rid_aligns['head_max_score']-rid_aligns['head_align_score']).abs().min()
		if abs(read_info['head_max_score']-read_info['head_align_score']) > rid_score_min_diff:
			lariat_reads.at[i, 'filter_failed'] = 'align_score'

	return lariat_reads


def choose_read_mapping(lariat_reads:pd.DataFrame):
	'''
	For reads with multiple lariat mappings that have passed all filters, choose just one to assign to the read and fail the others
	'''

	for rid in lariat_reads.read_id.unique():
		valid_maps = lariat_reads[(lariat_reads.read_id==rid) & (lariat_reads.filter_failed.isna())]
		# If either 1 or 0 valid lariat mappings, no need to choose
		if len(valid_maps) < 2:
			continue
		valid_read_one_lariat_mappings = valid_maps[valid_maps.read_num==1]
		valid_read_two_lariat_mappings = valid_maps[valid_maps.read_num==2]
		if len(valid_read_one_lariat_mappings) > 0:
			valid_maps = valid_read_one_lariat_mappings
			lariat_reads.loc[valid_read_two_lariat_mappings.index, 'filter_failed'] = 'not_chosen'
		else:
			valid_maps = valid_read_two_lariat_mappings

		# Prioritize a mapping based on whether or not the BP base is a mismatch and the alignment orientation
		# If multiple valid mappings in the same category, choose one at random
		mismatch_and_forward = valid_maps[(valid_maps.bp_mismatch) & (valid_maps.read_orient_to_gene=="Forward")]
		mismatch_and_reverse = valid_maps[(valid_maps.bp_mismatch) & (valid_maps.read_orient_to_gene=="Reverse")]
		match_and_forward = valid_maps[(~valid_maps.bp_mismatch) & (valid_maps.read_orient_to_gene=="Forward")]
		match_and_reverse = valid_maps[(~valid_maps.bp_mismatch) & (valid_maps.read_orient_to_gene=="Reverse")]

		random.seed(1)		# For consistent output
		for category in (mismatch_and_forward, mismatch_and_reverse, match_and_forward, match_and_reverse):
			if category.empty is True:
				continue
			
			chosen_index = random.sample(category.index.to_list(), 1)[0]
			rejected_indices = [ind for ind in valid_maps.index if ind!=chosen_index]
			lariat_reads.loc[rejected_indices, 'filter_failed'] = 'not_chosen'
			break

	return lariat_reads

def parse_splice_junction_count(output_base:str, log) -> int:
	'''
	Get the total count of splice junctions from the SJ.out.tab file output by STAR
	'''
	total_sj_count = 0
	with open(OUTPUT_SJ_FILE.format(output_base)) as in_file:
		for line in in_file:
			uniq_read_count = line.strip().split('\t')[6]
			total_sj_count += int(uniq_read_count)

	return total_sj_count


# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	# Get args
	threads, output_base, hisat2_index, ref_genes, ref_introns, ref_repeatmasker, log_level = sys.argv[1:]

	# Get logger
	log = utils.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}')

	# Load putatitve lariat alignments
	log.debug('Parsing lariat reads...')
	lariat_reads = pd.read_csv(PUTATITVE_LARIATS_FILE.format(output_base), sep='\t', dtype={'filter_failed': 'object'}, na_filter=False)
	if lariat_reads.empty:
		log.info('No reads remaining')
		sys.exit(4)
	read_id_count = lariat_reads['read_id'].str.slice(stop=-6).unique().size
	log.info(f'Pre-filter read count = {read_id_count}')

	# Load gene and intron coordinates
	genes = parse_genes(ref_genes, log)
	introns = parse_introns(ref_introns, log)

	# Check for reads aligned to annotated repetitive region 
	log.debug('Checking repeat regions')
	repeat_rids = check_repeat_overlap(lariat_reads, ref_repeatmasker, log)

	# Check for template-switching reads
	log.debug('Getting template-switching read ids')
	temp_switch_rids = set(pd.read_csv(TEMPLATE_SWITCH_FILE.format(output_base), sep='\t', usecols=['read_id'], na_filter=False).read_id)
	# Check for circularized intron reads
	log.debug('Getting circularized intron read ids')
	circular_rids = set(pd.read_csv(CIRCULARS_FILE.format(output_base), sep='\t', usecols=['read_id'], na_filter=False).read_id)

	# Filter lariat reads
	log.debug('Filtering lariat reads')
	lariat_reads['filter_failed'] = lariat_reads.apply(filter_lariats, repeat_rids=repeat_rids, temp_switch_rids=temp_switch_rids, circular_rids=circular_rids, axis=1).astype('object')

	# Perform full alignment of tail segment and filter out reads where tail doesn't completely
	# align to the expected position
	lariat_reads = verify_fivep_alignment(lariat_reads, output_base, hisat2_index, genes, threads, log)

	# Filter out reads where head and tail segments are nested 
	lariat_reads = filter_nested_reads(lariat_reads, log)

	# Filter remaining alignment to those with the least difference between
	# the head alignment score and the maximum head score
	lariat_reads = filter_align_score(lariat_reads)

	# Choose 1 lariat mapping per read id 
	lariat_reads = choose_read_mapping(lariat_reads)

	# Get splice junction count
	lariat_reads['sj_count_total'] = parse_splice_junction_count(output_base, log)
	
	# Separate failed mappings from passed mappings
	failed_aligns = lariat_reads[lariat_reads.filter_failed.notna()].copy()
	filtered_lariats = lariat_reads.loc[lariat_reads.filter_failed.isna()]
	log.info(f'Post-filter read count = {filtered_lariats.read_id.nunique()}')

	# Now write it all to file
	log.info('Writing results to output files...')
	failed_aligns.to_csv(FAILED_LARIATS_FILE.format(output_base), sep='\t', index=False)

	filtered_lariats = filtered_lariats.drop(columns=['head_align_score', 'head_max_score', 'filter_failed'])
	filtered_lariats.to_csv(LARIATS_FILE.format(output_base), sep='\t', index=False)

	log.debug('End of script')