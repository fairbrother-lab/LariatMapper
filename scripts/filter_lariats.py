import sys
import gzip
import random
import os

import pandas as pd
import numpy as np

import functions



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
# In files
OUTPUT_BAM_FILE = "{}output.bam"
TEMP_SWITCH_FILE = "{}template_switching_reads.tsv"
CIRCULARS_FILE = "{}circularized_intron_reads.tsv"
PUTATITVE_LARIATS_FILE = "{}putative_lariats.tsv"
# Out files
FAILED_LARIATS_FILE = "{}failed_lariat_alignments.tsv"
LARIATS_FILE = "{}lariat_reads.tsv"
RUN_DATA_FILE = "{}read_counts.tsv"

FINAL_RESULTS_COLS = ['read_id',
                        'gene_id',
                        'chrom',
                        'strand',
                        'fivep_pos',
                        'bp_pos',
                        'threep_pos',
                        'bp_dist_to_threep',
						'read_is_reverse',
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
	lariat_reads = pd.read_csv(PUTATITVE_LARIATS_FILE.format(output_base), sep='\t', dtype={'filter_failed': 'object'}, na_filter=False)
	lariat_reads = lariat_reads.rename(columns={'align_start': 'head_start', 'align_end': 'head_end'})

	# If no reads are left, end the run early
	if lariat_reads.empty:
		log.info('No reads remaining')
		with open(LARIATS_FILE.format(output_base), 'w') as w:
			w.write('\t'.join(FINAL_RESULTS_COLS))
		with open(FAILED_LARIATS_FILE.format(output_base), 'w') as w:
			w.write('\t'.join(FINAL_RESULTS_COLS) + '\tfilter_failed')

		sys.exit(4)
	
	lariat_reads.read_id = lariat_reads.read_id.str.slice(0,-4)
	lariat_reads[['read_id', 'read_num']] = lariat_reads.read_id.str.split('/', expand=True)
	lariat_reads.read_num = lariat_reads.read_num.astype(int)
	lariat_reads['align_mismatch'] = lariat_reads.read_bp_nt != lariat_reads.genomic_bp_nt
	lariat_reads['max_quality'] = lariat_reads.groupby('read_id').quality.agg('max')
	
	return lariat_reads


def add_mapped_reads(output_base:str, seq_type:str, log) -> int:
	'''
	Get the number of reads that mapped to the reference genome
	'''
	command = f'samtools view --count --exclude-flags 12 {OUTPUT_BAM_FILE.format(output_base)}'
	count = int(functions.run_command(command, log=log))
	
	if seq_type == 'paired':
		count = count//2

	return count


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
	bedtools_fivep_output = functions.run_command(bedtools_call, input=bedtools_fivep_input, log=log).split('\n')
	bedtools_bp_output = functions.run_command(bedtools_call, input=bedtools_bp_input, log=log).split('\n')

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

	if row['read_id'] in temp_switch_rids:
		return 'template_switching'
	
	if row['read_id'] in circular_rids:
		return 'circularized_intron'

	if row['read_id'] in repeat_rids:
		return 'in_repeat'
	
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
	output_base, log_level, seq_type, ref_fasta, ref_repeatmasker = sys.argv[1:]

	# Get logger
	log = functions.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}')

	# Load putatitve lariat alignments
	log.debug('Parsing lariat reads...')
	lariat_reads = load_lariat_table(output_base, log)
	log.info(f'Pre-filter read count = {len(lariat_reads.read_id.unique())}')

	# Get linear mapped reads count 
	lariat_reads['total_mapped_reads'] = add_mapped_reads(output_base, seq_type, log)

	# Check for reads aligned to annotated repetitive region 
	log.debug('Checking repeat regions')
	repeat_rids = check_repeat_overlap(lariat_reads, ref_repeatmasker, log)

	# Check for reads which were probably created from the reverse-transcriptase switching RNA templates at the branchpoint
	# Check for circularized intron reads
	log.debug('Getting template-switching and circularized intron reads')
	temp_switch_rids = set(pd.read_csv(TEMP_SWITCH_FILE.format(output_base), sep='\t', usecols=['read_id'], na_filter=False).read_id)
	circular_rids = set(pd.read_csv(CIRCULARS_FILE.format(output_base), sep='\t', usecols=['read_id'], na_filter=False).read_id)

	# Filter lariat reads
	log.debug('Filtering lariat reads')
	lariat_reads['filter_failed'] = lariat_reads.apply(filter_lariats, repeat_rids=repeat_rids, temp_switch_rids=temp_switch_rids, circular_rids=circular_rids, axis=1).astype('object')

	# Choose 1 lariat mapping per read id 
	lariat_reads = choose_read_mapping(lariat_reads)
	
	# Seperate failed mappings from passed mappings
	failed_aligns = lariat_reads[lariat_reads.filter_failed.notna()].copy()
	filtered_lariats = lariat_reads.loc[lariat_reads.filter_failed.isna(), FINAL_RESULTS_COLS]
	log.info(f'Post-filter read count = {filtered_lariats.read_id.nunique()}')

	# Now write it all to file
	log.info('Writing results to output files...')
	failed_aligns.to_csv(FAILED_LARIATS_FILE.format(output_base), sep='\t', index=False)
	filtered_lariats.to_csv(LARIATS_FILE.format(output_base), sep='\t', index=False)

	log.debug('End of script')