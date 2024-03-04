import sys
import os
import subprocess 
import gzip
import time
import random

from intervaltree import Interval, IntervalTree
import pandas as pd


# =============================================================================#
#                                  Constants                                   #
# =============================================================================#
FINAL_RESULTS_COLUMNS = ('gene_name',
                        'gene_id',
                        'gene_type',
                        'read_id',
                        'read_seq',
                        'chrom',
                        'strand',
                        'fivep_pos',
                        'threep_pos',
                        'bp_pos',
                        'read_bp_nt',
                        'genomic_bp_nt',
                        'genomic_bp_context',
                        'bp_dist_to_threep',
                        'total_mapped_reads')



# =============================================================================#
#                                  Functions                                  #
# =============================================================================#
def load_splice_site_info(ref_introns) -> tuple:
    '''
    Returns a dict formatted as follows:
    {Chromosome: {Strand(+ or -): Intervaltree(StartPosition(int), EndPosition(int))}}
    for 3' splice sites (+/- 2 bases), 5' splice sites (+/- 2 bases), and the introns they come from (start to end)
    '''
    threep_sites, fivep_sites = {}, {}

    if ref_introns[-2:] == 'gz':
        intron_file = gzip.open(ref_introns, 'rt')
    else:
        intron_file = open(ref_introns)

    introns_done = set()
    for line in intron_file:
        chrom, start, end, _, _, strand = line.strip().split()
        if chrom not in threep_sites:
            threep_sites[chrom] = {s: IntervalTree() for s in ['+', '-']}
            fivep_sites[chrom] = {s: IntervalTree() for s in ['+', '-']}
        intron_id = '{}_{}_{}_{}'.format(chrom, strand, start, end)
        if intron_id not in introns_done:
            start, end = int(start), int(end)
            if strand == '+':
                threep_sites[chrom][strand].add(Interval(end-2, end+2))
                fivep_sites[chrom][strand].add(Interval(start-2, start+2))
            else:
                threep_sites[chrom][strand].add(Interval(start-2, start+2))
                fivep_sites[chrom][strand].add(Interval(end-2, end+2))
            introns_done.add(intron_id)

    intron_file.close()

    return threep_sites, fivep_sites


def load_lariat_table(output_base: str) -> pd.DataFrame:
	'''
	For a given lariat-mapping of a fastq file, retrieve all the lariat reads from the XXX_lariat_info_table.tsv and put them in a dict, which
	can then be added to the experiment-wide superset dict
	'''
	lariat_reads = pd.read_csv(f'{output_base}final_info_table.tsv', sep='\t')
	lariat_reads = lariat_reads.rename(columns={'fivep_site': 'fivep_pos', 'threep_site': 'threep_pos', 'bp_site': 'bp_pos'})

	if len(lariat_reads) == 0:
		print(time.strftime('%m/%d/%y - %H:%M:%S') + '| No reads remaining')
		with open(f'{output_base}lariat_reads.tsv', 'w') as w:
			w.write('\t'.join(FINAL_RESULTS_COLUMNS))
		with open(f'{output_base}failed_lariat_mappings.tsv', 'w') as w:
			w.write('\t'.join(FINAL_RESULTS_COLUMNS) + '\tfilter_failed')
		exit()

	# Code adapted from https://stackoverflow.com/questions/27298178/concatenate-strings-from-several-rows-using-pandas-groupby
	# Some reads get mapped to coordinates with multiple overlapping gene annotations
	# We resolve this by collapsing the duplicated rows and concatenating the gene_id, gene_name, and gene_type columns
	lariat_reads = (lariat_reads.groupby([col for col in lariat_reads.columns if col not in ('gene_id', 'gene_name', 'gene_type')])
				 				.agg({'gene_id': ','.join, 'gene_name': ','.join, 'gene_type': ','.join})
								.reset_index()
					)

	return lariat_reads


def add_mapped_reads(output_base:str) -> int:
	'''
	Get the number of reads that mapped to the reference genome from the *_run_data.tsv file 
	'''
	with open(f'{output_base}run_data.tsv', 'r') as file:
		line = file.readline()
		sample_read_count = int(line.split('\t')[1])

	return sample_read_count


def check_repeat_overlap(lariat_reads: pd.DataFrame, output_base:str, delete_temp:bool=True) -> set:
	''' 
    Check if both the 5'SS and the BP overlap with a repetitive region
    '''
	# Write the 5'ss and BP coordinates to BED files
	fivep_tmp_bed, bp_tmp_bed = f'{output_base}fivep_tmp.bed', f'{output_base}bp_tmp.bed'
	with open(fivep_tmp_bed, 'w') as fivep_out, open(bp_tmp_bed, 'w') as bp_out:
		for i, row in lariat_reads.iterrows():
			fivep_out.write(f"{row['chrom']}\t{row['fivep_pos']-1}\t{row['fivep_pos']+1}\t{row['read_id']}\n")
			bp_out.write(f"{row['chrom']}\t{row['bp_pos']-1}\t{row['bp_pos']+1}\t{row['read_id']}\n")

	# Identify 5'ss's that overlap a repeat region
	fivep_overlap_bed, bp_overlap_bed = f'{output_base}fivep_repeat_overlaps.bed', f'{output_base}bp_repeat_overlaps.bed'
	with open(fivep_overlap_bed, 'w') as out_file:
		subprocess.run(f'bedtools intersect -u -a {fivep_tmp_bed} -b {ref_repeatmasker}'.split(' '),
			stdout=out_file)
	# Identify BPs that overlap a repeat region
	with open(bp_overlap_bed, 'w') as out_file:
		subprocess.run(f'bedtools intersect -u -a {bp_tmp_bed} -b {ref_repeatmasker}'.split(' '),
			stdout=out_file)

	# Add reads where both sites overlapped to repeat_rids for removal
	fivep_repeat_rids, bp_repeat_rids = set(), set()
	with open(fivep_overlap_bed) as in_file:
		for line in in_file:
			_, _, _, rid = line.strip().split()
			fivep_repeat_rids.add(rid)
	with open(bp_overlap_bed) as in_file:
		for line in in_file:
			_, _, _, rid = line.strip().split()
			bp_repeat_rids.add(rid)
	repeat_rids = fivep_repeat_rids.intersection(bp_repeat_rids)

	# Delete the temporary files
	if delete_temp is True:
		for temp_file in (fivep_tmp_bed, bp_tmp_bed, fivep_overlap_bed, bp_overlap_bed):
			os.remove(temp_file)

	return repeat_rids


def filter_lariats(row:pd.Series, fivep_sites:dict, threep_sites:dict, repeat_rids:set):
	'''
	Filter the candidate lariat reads to EXCLUDE any that meet the following criteria:
			- BP is within 2bp of a splice site (likely an intron circle, not a lariat)
			- Read maps to UBB or UBC (likely false positive due to the repetitive nature of the genes)
			- Both the 5'SS and the BP overlap with repetitive regions from RepeatMasker (likely false positive due to sequence repetition)
	'''
	# Check if BP is within 2bp of an annotated splice site
	if fivep_sites[row['chrom']][row['strand']].overlaps(row['bp_pos']) or threep_sites[row['chrom']][row['strand']].overlaps(row['bp_pos']):
		return 'near_ss'

	# Check if read mapped to a ubiquitin gene
	if row['gene_name'] in ('UBC', 'UBB'):
		return 'ubiquitin_gene'
	
	# Check if read mapped to a repetitive region
	if row['read_id'] in repeat_rids:
		return 'in_repeat'

	return pd.NA


def choose_read_mapping(lariat_reads):
	'''
	For reads with multiple lariat mappings that have passed all filters, choose just one to assign to the read and fail the others
	'''
	lariat_reads['align_mismatch'] = lariat_reads.read_bp_nt != lariat_reads.genomic_bp_nt
	lariat_reads.read_id = lariat_reads.read_id.str.slice(0,-4)

	for rid in lariat_reads.read_id.unique():
		valid_lariat_mappings = lariat_reads[(lariat_reads.read_id==rid) & (lariat_reads.filter_failed.isna())]
		# If either 1 or 0 valid lariat mappings, no need to choose
		if len(valid_lariat_mappings) < 2:
			continue

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
			
			chosen_index = random.sample(category.index.to_list(), 1)
			rejected_indices = [ind for ind in valid_lariat_mappings.index if ind!=chosen_index]
			lariat_reads.loc[rejected_indices, 'filter_failed'] = 'not_chosen'
			break

	return lariat_reads



# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	ref_gtf, ref_introns, ref_repeatmasker, output_base = sys.argv[1:]
	# ref_gtf = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.gencode.v44.basic.annotation.gtf'
	# ref_introns = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.gencode.v44.basic.introns.bed'
	# ref_repeatmasker = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.repeat_masker.bed.gz'
	# output_base = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/output/C22_R1_lariat_mapping/'

	print(ref_gtf, ref_introns, ref_repeatmasker, output_base)

	# Load splice site coordinates
	print(time.strftime('%m/%d/%y - %H:%M:%S | Parsing splice site info...'))
	threep_sites, fivep_sites = load_splice_site_info(ref_introns)

	print(time.strftime('%m/%d/%y - %H:%M:%S | Parsing lariat reads...'))
	lariat_reads = load_lariat_table(output_base)

	print(time.strftime('%m/%d/%y - %H:%M:%S | Adding total mapped read count...'))
	lariat_reads['total_mapped_reads'] = add_mapped_reads(output_base)

	print(time.strftime('%m/%d/%y - %H:%M:%S | Checking for overlaps with repeat regions...'))
	repeat_rids = check_repeat_overlap(lariat_reads, output_base, delete_temp=False)

	# Filter lariat reads
	print(time.strftime('%m/%d/%y - %H:%M:%S | Filtering lariat reads...'))
	lariat_reads['filter_failed'] = lariat_reads.apply(filter_lariats, repeat_rids=repeat_rids, fivep_sites=fivep_sites, threep_sites=threep_sites, axis=1)

	# Choose 1 lariat mapping per read id and remove the _for/_rev suffix
	lariat_reads = choose_read_mapping(lariat_reads)
	
	# Seperate failed mappings from passed mappings
	failed_mappings = lariat_reads[lariat_reads.filter_failed.notna()].copy()
	filtered_lariats = lariat_reads.loc[lariat_reads.filter_failed.isna(), FINAL_RESULTS_COLUMNS]

	print(time.strftime('%m/%d/%y - %H:%M:%S') + f' | Pre-filter read count = {len(lariat_reads.read_id.unique())}')
	print(time.strftime('%m/%d/%y - %H:%M:%S') + f' | Post-filter read count = {len(filtered_lariats.read_id.unique())}')

	# Now write it all to file
	print(time.strftime('%m/%d/%y - %H:%M:%S | Writing results to output files...'))
	failed_mappings.to_csv(f'{output_base}failed_lariat_mappings.tsv', sep='\t', index=False)
	filtered_lariats.to_csv(f'{output_base}lariat_reads.tsv', sep='\t', index=False)
	
	# Record final lariat read count
	with open(f'{output_base}run_data.tsv', 'a') as a:
		a.write(f'filtered_lariats\t{len(filtered_lariats)}\n')