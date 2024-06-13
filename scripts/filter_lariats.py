import sys
import subprocess 
import gzip
import random
import logging

import pandas as pd
import numpy as np

import functions



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
FINAL_RESULTS_COLS = ['read_id',
                        'gene_id',
                        'chrom',
                        'strand',
                        'fivep_pos',
                        'bp_pos',
                        'threep_pos',
                        'bp_dist_to_threep',
						'read_alignment',
						'read_bp_pos',
                        'read_seq',
                        'read_bp_nt',
						'genomic_bp_nt',
                        'genomic_bp_context',
                        'total_mapped_reads',
						]



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def load_lariat_table(output_base:str, log) -> pd.DataFrame:
	'''
	For a given lariat-mapping of a fastq file, retrieve all the lariat reads from the XXX_lariat_info_table.tsv and put them in a dict, which
	can then be added to the experiment-wide superset dict
	'''
	lariat_reads = pd.read_csv(f'{output_base}trimmed_info_table.tsv', sep='\t')
	lariat_reads = lariat_reads.rename(columns={'align_start': 'head_start', 'align_end': 'head_end'})

	if lariat_reads.empty:
		log.info('No reads remaining')
		with open(f'{output_base}lariat_reads.tsv', 'w') as w:
			w.write('\t'.join(FINAL_RESULTS_COLS))
		with open(f'{output_base}failed_lariat_alignments.tsv', 'w') as w:
			w.write('\t'.join(FINAL_RESULTS_COLS) + '\tfilter_failed')
		with open(f'{output_base}circularized_intron_reads.tsv', 'w') as w:
			w.write('\t'.join(FINAL_RESULTS_COLS))
		exit()
	
	lariat_reads.read_id = lariat_reads.read_id.str.slice(0,-4)
	lariat_reads[['read_id', 'read_num']] = lariat_reads.read_id.str.split('/', expand=True)
	lariat_reads.read_num = lariat_reads.read_num.astype(int)
	lariat_reads['align_mismatch'] = lariat_reads.read_bp_nt != lariat_reads.genomic_bp_nt
	lariat_reads['read_alignment'] = lariat_reads.read_is_reverse.map({True: 'reverse', False: 'forward'})
	lariat_reads['max_quality'] = lariat_reads.groupby('read_id').quality.agg('max')
	
	return lariat_reads


def add_mapped_reads(output_base:str) -> int:
	'''
	Get the number of reads that mapped to the reference genome from the *read_counts.tsv file 
	'''
	with open(f'{output_base}read_counts.tsv', 'r') as file:
		line = file.readline()
		sample_read_count = int(line.split('\t')[1])

	return sample_read_count


#TODO: Implement a more elaborate check. We only need to check for UBC-type repeats as possible false-positives
def check_repeat_overlap(lariat_reads: pd.DataFrame, ref_repeatmasker:str) -> set:
	''' 
    Check if both the 5'SS and the BP overlap with a repetitive region
    '''
	# Write the 5'ss and BP coordinates to bedtools input strings
	bedtools_fivep_input, bedtools_bp_input = '', ''
	for _, row in lariat_reads.iterrows():
		bedtools_fivep_input += f"{row['chrom']}\t{row['fivep_pos']-1}\t{row['fivep_pos']+1}\t{row['read_id']}\n"
		bedtools_bp_input += f"{row['chrom']}\t{row['bp_pos']-1}\t{row['bp_pos']+1}\t{row['read_id']}\n"

	# Identify 5'ss that overlap a repeat region
	bedtools_call = f'bedtools intersect -u -a - -b {ref_repeatmasker}'
	bedtools_fivep_output = subprocess.run(bedtools_call.split(' '), input=bedtools_fivep_input, capture_output=True, text=True).stdout.strip().split('\n')
	bedtools_bp_output = subprocess.run(bedtools_call.split(' '), input=bedtools_bp_input, capture_output=True, text=True).stdout.strip().split('\n')

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
	if row['quality'] < row['max_quality']:
		return 'align_quality'

	if row['read_id'] in repeat_rids:
		return 'in_repeat'
	
	if row['read_id'] in temp_switch_rids:
		return 'template_switching'
	
	if row['read_id'] in circular_rids:
		return 'circularized_intron'

	return np.nan


def choose_read_mapping(lariat_reads):
	'''
	For reads with multiple lariat mappings that have passed all filters, choose just one to assign to the read and fail the others
	'''
	for rid in lariat_reads.read_id.unique():
		valid_lariat_mappings = lariat_reads[(lariat_reads.read_id==rid) & (lariat_reads.filter_failed.isna())]
		# If either 1 or 0 valid lariat mappings, no need to choose
		if len(valid_lariat_mappings) < 2:
			continue
		valid_read_one_lariat_mappings = valid_lariat_mappings[valid_lariat_mappings.read_num==1]
		valid_read_two_lariat_mappings = valid_lariat_mappings[valid_lariat_mappings.read_num==2]
		if len(valid_read_one_lariat_mappings) > 0:
			valid_lariat_mappings = valid_read_one_lariat_mappings
			lariat_reads.loc[valid_read_two_lariat_mappings.index, 'filter_failed'] = 'not_chosen'
		else:
			valid_lariat_mappings = valid_read_two_lariat_mappings

		# Prioritize a mapping based on whether or not the BP base is a mismatch and the alignment orientation
		# If multiple valid mappings in the same category, choose one at random
		mismatch_and_forward = valid_lariat_mappings[(valid_lariat_mappings.align_mismatch) & (~valid_lariat_mappings.read_is_reverse)]
		mismatch_and_reverse = valid_lariat_mappings[(valid_lariat_mappings.align_mismatch) & (valid_lariat_mappings.read_is_reverse)]
		match_and_forward = valid_lariat_mappings[(~valid_lariat_mappings.align_mismatch) & (~valid_lariat_mappings.read_is_reverse)]
		match_and_reverse = valid_lariat_mappings[(~valid_lariat_mappings.align_mismatch) & (valid_lariat_mappings.read_is_reverse)]

		random.seed(1)		# For consistent output
		for category in (mismatch_and_forward, mismatch_and_reverse, match_and_forward, match_and_reverse):
			if category.empty is True:
				continue
			
			chosen_index = random.sample(category.index.to_list(), 1)[0]
			rejected_indices = [ind for ind in valid_lariat_mappings.index if ind!=chosen_index]
			lariat_reads.loc[rejected_indices, 'filter_failed'] = 'not_chosen'
			break

	return lariat_reads



# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	# Get args
	ref_fasta, ref_repeatmasker, output_base, log_level = sys.argv[1:]

	# Get logger
	log = functions.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}')

	# Load putatitve lariat alignments
	log.debug('Parsing lariat reads...')
	lariat_reads = load_lariat_table(output_base, log)
	log.info(f'Pre-filter read count = {len(lariat_reads.read_id.unique())}')

	# Get linear mapped reads count 
	lariat_reads['total_mapped_reads'] = add_mapped_reads(output_base)

	# Check for reads aligned to annotated repetitive region 
	repeat_rids = check_repeat_overlap(lariat_reads, ref_repeatmasker)

	# Check for reads which were probably created from the reverse-transcriptase switching RNA templates at the branchpoint
	temp_switch_rids = set(pd.read_csv(f'{output_base}template_switching_reads.tsv', sep='\t', usecols=['read_id']).read_id)

	# Check for circularized intron reads
	circular_rids = set(pd.read_csv(f'{output_base}circularized_intron_reads.tsv', sep='\t', usecols=['read_id']).read_id)

	# Filter lariat reads
	lariat_reads['filter_failed'] = lariat_reads.apply(filter_lariats, repeat_rids=repeat_rids, temp_switch_rids=temp_switch_rids, circular_rids=circular_rids, axis=1)

	# Choose 1 lariat mapping per read id 
	lariat_reads = choose_read_mapping(lariat_reads)
	
	# Seperate failed mappings from passed mappings
	failed_aligns = lariat_reads[lariat_reads.filter_failed.notna()].copy()
	filtered_lariats = lariat_reads.loc[lariat_reads.filter_failed.isna(), FINAL_RESULTS_COLS]
	log.info(f'Post-filter read count = {filtered_lariats.read_id.nunique()}')

	# Now write it all to file
	log.info('Writing results to output files...')
	failed_aligns.to_csv(f'{output_base}failed_lariat_alignments.tsv', sep='\t', index=False)
	filtered_lariats.to_csv(f'{output_base}lariat_reads.tsv', sep='\t', index=False)

	# Record read count
	with open(f'{output_base}read_counts.tsv', 'a') as a:
		a.write(f'trimmed_filter_passed\t{lariat_reads.read_id.nunique()}\n')	
		a.write(f'lariat\t{filtered_lariats.read_id.nunique()}\n')
		a.write(f'template-switch\t{len(temp_switch_rids)}\n')
		a.write(f'circularized_intron\t{len(circular_rids)}\n')

	log.debug('End of script')