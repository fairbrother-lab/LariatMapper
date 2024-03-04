import sys
import os
import gzip
import subprocess
import time

import pysam
import pandas as pd
from intervaltree import Interval, IntervalTree



# =============================================================================#
#                                  Constants                                   #
# =============================================================================#
COMP_NTS = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
MAX_MISMATCHES = 5
MAX_MISMATCH_PERCENT = 0.1
MAX_GAP_LENGTH = 3

FINAL_INFO_TABLE_COLS = ['read_id', 
						'read_is_reverse', 
						'read_seq', 
						'chrom', 
						'strand', 
						'align_start',
						'align_end', 
						'gene_id', 
						'gene_name', 
						'gene_type', 
						'fivep_site', 
						'bp_site', 
						'read_bp_nt', 
						'genomic_bp_context', 
						'genomic_bp_nt', 
						'threep_site', 
						'bp_dist_to_threep']



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def load_intron_info(ref_introns:str) -> tuple:
	'''
	Returns a dict formatted as follows:
	{Chromosome: {Strand(+ or -): Intervaltree(StartPosition(int), EndPosition(int))}}
	for 3' splice sites (+/- 2 bases), 5' splice sites (+/- 2 bases), and the introns they come from (start to end)
	'''
	introns = {}

	if ref_introns[-2:] == 'gz':
		intron_file = gzip.open(ref_introns, 'rt')
	else:
		intron_file = open(ref_introns)

	introns_done = set()
	for line in intron_file:
		chrom, start, end, _, _, strand = line.strip().split('\t')
		if chrom not in introns:
			introns[chrom] = {s: IntervalTree() for s in ['+', '-']}
		intron_id = '{}_{}_{}_{}'.format(chrom, strand, start, end)
		if intron_id not in introns_done:
			start, end = int(start), int(end)
			introns[chrom][strand].add(Interval(start, end))
			introns_done.add(intron_id)

	intron_file.close()
	return introns


def load_fivep_info(output_base:str):
	fivep_info = pd.read_csv(f'{output_base}fivep_info_table.tsv', sep='\t')

	# Unpack fivep_sites
	fivep_info.fivep_sites = fivep_info.fivep_sites.str.split(',')
	fivep_info = fivep_info.explode('fivep_sites')
	fivep_info[['chrom', 'fivep_site', '_', 'strand']] = fivep_info.fivep_sites.str.split(';', expand=True)
	fivep_info.fivep_site = fivep_info.fivep_site.astype(int)

	return fivep_info


def load_gene_ranges(gtf_file):
	'''
	Read through the GTF annotation file, return a dict with all annotated gene coordinates grouped by chromosome, then strand
	Output format = { chromosome: {strand: IntervalTree(gene start position, gene end position, gene name) } }
	'''
	# open file
	if gtf_file[-2:] == 'gz':
		in_file = gzip.open(gtf_file, 'rt')
	else:
		in_file = open(gtf_file)

	gene_ranges = {}		
	for line in in_file:
		if line[0] != '#':
			chrom, _, feat, start, end, _, strand, _, attributes = line.strip().split('\t')
			if feat == 'gene':
				# Convert attributes string to dict of {tag name: value}
				attributes = attributes[:-1].split('; ')
				attributes = {a.split(' ')[0]:a.split(' ')[1].replace('\"', '') for a in attributes}
				if chrom not in gene_ranges:
					gene_ranges[chrom] = {s:IntervalTree() for s in ('-', '+')}
				gene_ranges[chrom][strand].add(Interval(int(start)-1, int(end), {'gene_id':attributes['gene_id']}))

	in_file.close()
	return gene_ranges


def add_gene_ids(row:pd.Series, gene_ranges):
	return [g.data['gene_id'] for g in gene_ranges[row['chrom']][row['strand']].overlap(row['fivep_site'], row['fivep_site']+1)]
	

# ALIGNMENT_COLUMNS = ('read_id', 'read_is_reverse', 'chrom', 'align_is_reverse', 'align_start', 'align_end', 'quality', 'mismatches', 'mismatch_percent', 'gaps', 'introns')
# def get_read_aligns(bam_file):
# 	'''
# 	Generator
# 	'''
# 	alignments = []
# 	current_read_id = next(pysam.AlignmentFile(bam_file)).query_name

# 	for align in pysam.AlignmentFile(bam_file, 'rb'):
# 		# Parse alignment info
# 		read_is_reverse = True if align.query_name.endswith('_rev') else False
# 		mismatches = align.get_tag('XM')
# 		mismatch_percent = mismatches/align.query_length
# 		gaps = [length for tag_num, length in align.cigartuples if tag_num in (1,2)]
# 		if align.has_tag('YB'):
# 			introns = align.get_tag('YB').split(',')[1:]
# 			introns = [intron.split(';') for intron in introns]
# 		else:
# 			introns = []

# 		# Collect the info in a list which will become a DataFrame row
# 		next_align = [align.query_name, 
# 					read_is_reverse, 
# 					align.reference_name, 
# 					align.is_reverse, 
# 					align.reference_start, 
# 					align.reference_end, 
# 					align.mapping_quality, 
# 					mismatches, 
# 					mismatch_percent, 
# 					gaps, 
# 					introns]
		
# 		# If still processing alignments from the same read, add next_align to alignments
# 		if align.query_name == current_read_id:
# 			alignments.append(next_align)

# 		# If reached alignments for a new read, yield the accumulated alignments and reset 
# 		else:
# 			trim_seq = align.get_forward_sequence()
# 			alignments = pd.DataFrame(alignments, columns=ALIGNMENT_COLUMNS)
# 			yield current_read_id, trim_seq, alignments

# 			current_read_id = align.query_name
# 			alignments = [next_align]

# 	# Upon reaching the end of the BAM file, yield the last alignments DataFrame
# 	trim_seq = align.get_forward_sequence()
# 	alignments = pd.DataFrame(alignments, columns=ALIGNMENT_COLUMNS)
# 	yield current_read_id, trim_seq, alignments


def parse_align_info(align:pysam.AlignedSegment) -> list:
	# Parse alignment info
	read_is_reverse = True if align.query_name.endswith('_rev') else False
	mismatches = align.get_tag('XM')
	mismatch_percent = mismatches/align.query_length
	gaps = [length for tag_num, length in align.cigartuples if tag_num in (1,2)]
	if align.has_tag('YB'):
		introns = align.get_tag('YB').split(',')[1:]
		introns = [intron.split(';') for intron in introns]
	else:
		introns = []

	return [align.query_name, 
			read_is_reverse, 
			align.query_sequence,
			align.reference_name, 
			align.is_reverse, 
			align.reference_start, 
			align.reference_end, 
			align.mapping_quality, 
			mismatches, 
			mismatch_percent, 
			gaps, 
			introns]


ALIGNMENT_COLUMNS = ('read_id', 'read_is_reverse', 'trim_seq', 'chrom', 'align_is_reverse', 'align_start', 'align_end', 'quality', 'mismatches', 'mismatch_percent', 'gaps', 'introns')
def get_read_aligns(bam_file:str, batch_size:int):
	'''
	Generator
	'''
	i = 0
	alignment_batch = []
	bam_iterator = pysam.AlignmentFile(bam_file, 'rb')
	for align in bam_iterator:
		i += 1

		# Collect the info in a list which will become a DataFrame row
		align_row = parse_align_info(align)
		alignment_batch.append(align_row)
		
		if i == batch_size:
			# Check that you have all of the alignments for the last read 
			# If not, add them
			next_align = next(bam_iterator)
			next_align_row = parse_align_info(next_align)
			while align.query_name == next_align.query_name:
				alignment_batch.append(next_align_row)
				next_align = next(bam_iterator)
				next_align_row = parse_align_info(next_align)

			# Convert list of lists to DataFrame, add max_quality, and yield
			alignment_batch = pd.DataFrame(alignment_batch, columns=ALIGNMENT_COLUMNS)
			alignment_batch['max_quality'] = alignment_batch.read_id.map(alignment_batch.groupby('read_id').quality.agg('max'))
			yield alignment_batch

			# Reset for next batch
			i = 0
			alignment_batch = [next_align_row]

	# Upon reaching the end of the BAM file, yield the last alignments DataFrame
	alignment_batch = pd.DataFrame(alignment_batch, columns=ALIGNMENT_COLUMNS)
	alignment_batch['max_quality'] = alignment_batch.read_id.map(alignment_batch.groupby('read_id').quality.agg('max'))
	yield alignment_batch


MAX_MISMATCHES = 5
MAX_MISMATCH_PERCENT = 0.1
MAX_GAP_LENGTH = 3
def first_filters(row:pd.Series):
	# Check if mapping quality is not max
	if row['quality'] < row['max_quality']:
		return 'map_quality'

	# Check if too many mismatches
	if row['mismatches'] > MAX_MISMATCHES:
		return 'mismatch_num'
	if row['mismatch_percent'] > MAX_MISMATCH_PERCENT:
		return 'mismatch_percent'
	
	# Check if more than one insertion/deletion OR one insertion/deletion that's too long
	if len(row['gaps']) > 1 or (len(row['gaps']) == 1 and row['gaps'][0] > MAX_GAP_LENGTH):
		return 'indel'

	# # Check if the 5'ss alignment orientation doesn't match the genomic alignment orientation
	# if row['read_is_reverse'] != row['align_is_reverse']:
	# 	return 'same_orientation'
	
	# Check if not within an intron
	if row['introns'] == []:
		return 'within_intron'
	
	return pd.NA

COMP_NTS = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
def reverse_complement(seq):
	return ''.join([COMP_NTS[seq[i]] for i in range(len(seq)-1,-1,-1)])


def add_read_bp_nt(row:pd.Series) -> str:
	if row['strand']=='+' and not row['align_is_reverse']:
		return row['trim_seq'][-1]
	elif row['strand']=='+' and row['align_is_reverse']:
		return reverse_complement(row['trim_seq'][0])
	elif row['strand']=='-' and not row['align_is_reverse']:
		return reverse_complement(row['trim_seq'][0])
	elif row['strand']=='-' and row['align_is_reverse']:
		return row['trim_seq'][-1]


def add_nearest_threep(row:pd.Series, introns:dict):
	overlap_introns = list(introns[row['chrom']][row['strand']].overlap(row['bp_site'], row['bp_site']+1))
	if len(overlap_introns) == 0:
		raise RuntimeError(f"No introns overlapped the following:\n{row}")
	
	if row['strand'] == '+':
		threep_site = min(overlap_introns, key=lambda s: s.end-row['bp_site']).end - 1
	else:
		threep_site = min(overlap_introns, key=lambda s: row['bp_site']-s.begin).begin

	return threep_site


def add_bp_seq(row:pd.Series, output_base:str, genome_fasta:str):
	temp_bp_bed, temp_bp_seq = output_base+'_temp_bp_seqs.bed', output_base+'temp_bp_seqs.txt'
	temp_file = open(temp_bp_bed, 'w')
	if row['strand'] == '+':
		bp_start, bp_end = row['bp_site']-4, row['bp_site']+6
	else:
		bp_start, bp_end = row['bp_site']-5, row['bp_site']+5
	temp_file.write(f"{row['chrom']}\t{bp_start}\t{bp_end}\t{row['chrom']};{row['bp_site']};{row['strand']}\t0\t{row['strand']}\n")
	temp_file.close()
	subprocess.run(f'bedtools getfasta -fi {genome_fasta} -bed {temp_bp_bed} -fo {temp_bp_seq} -nameOnly -s -tab'.split(' '))

	temp_file = open(temp_bp_seq)
	name, genomic_bp_window = temp_file.readline().strip().split()
	genomic_bp_window = genomic_bp_window.upper()
	temp_file.close()


	return genomic_bp_window



def second_filters(row:pd.Series):
	# Check if no intron overlaps BP 
	if pd.isna(row['threep_site']):
		return 'no_intron'
	# Check if no fivep match
	if pd.isna(row['fivep_site']):
		return 'no_fivep_matches'

	# Check if the 5'ss is at or downstream of the BP
	if row['strand'] == '+' and row['fivep_site'] >= row['bp_site']:
		return '5p_bp_order'
	if row['strand'] == '-' and row['fivep_site'] <= row['bp_site']:
		return '5p_bp_order'

	return pd.NA



# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	ref_gtf, ref_introns, genome_fasta, output_base= sys.argv[1:]
	# ref_gtf = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.gencode.v44.basic.annotation.gtf'
	# ref_transcripts = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.gencode.v44.basic.transcripts.bed'
	# ref_introns = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.gencode.v44.basic.introns_mod.bed'
	# genome_fasta	= '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.gencode.v44.primary_assembly.fa'
	# output_base = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/output/C22_R1_lariat_mapping/'
	batch_size = 1_000
	# batch_size = 100_000


	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Beginning trimmed alignment filtering')
	# print(ref_gtf, ref_introns, genome_fasta, output_base)

	# 
	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Load')
	introns = load_intron_info(ref_introns)
	fivep_info = load_fivep_info(output_base)
	gene_ranges = load_gene_ranges(ref_gtf)
	fivep_info['gene_id'] = fivep_info.apply(add_gene_ids, gene_ranges=gene_ranges, axis=1, result_type='reduce')
	fivep_info = fivep_info.explode('gene_id')[['read_id', 'read_seq', 'fivep_site', 'gene_id']]
	fivep_info

	# 
	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Write')
	with open(f'{output_base}final_info_table.tsv', 'w') as w:
		w.write('\t'.join(FINAL_INFO_TABLE_COLS) + '\n')
	with open(f'{output_base}failed_trimmed_alignments.tsv', 'w') as w:
		w.write('\n')

	# 
	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Iterate')
	for alignments in get_read_aligns(f'{output_base}trimmed_reads_to_genome.bam', batch_size):
		alignments['filter_failed'] = alignments.apply(first_filters, axis=1)
		alignments[alignments.filter_failed.notna()].to_csv(f'{output_base}failed_trimmed_alignments.tsv', mode='a', sep='\t', index=False, header=False)
		alignments = alignments[alignments.filter_failed.isna()]
		if alignments.empty is True:
			continue

		# Explode introns
		alignments = alignments.explode('introns', ignore_index=True)
		print(alignments)
		alignments[['_', 'transcript_id', 'intron_num', 'intron_start', 'intron_end', 'strand', 'gene_id', 'gene_name', 'gene_type']] = alignments.introns.to_list()
		alignments.intron_start = alignments.intron_start.astype(int)
		alignments.intron_end = alignments.intron_end.astype(int)

		# Merge
		alignments = pd.merge(alignments, fivep_info, 'left', on=['read_id', 'gene_id'])

		# Add more info
		alignments['read_bp_nt'] = alignments.apply(add_read_bp_nt, axis=1)
		alignments['bp_site'] = alignments.apply(lambda row: row['align_end']-1 if row['strand']=='+' else row['align_start'], axis=1)
		alignments['threep_site'] = alignments.apply(add_nearest_threep, introns=introns, axis=1)
		alignments['bp_dist_to_threep'] = alignments.apply(lambda row: row['bp_site']-row['threep_site'] if pd.notna(row['threep_site']) else pd.NA, axis=1)
		alignments['genomic_bp_context'] = alignments.apply(add_bp_seq, output_base=output_base, genome_fasta=genome_fasta, axis=1)
		alignments['genomic_bp_nt'] = alignments.genomic_bp_context.str.get(4)
		
		# Second filters
		alignments['filter_failed'] = alignments.apply(second_filters, axis=1)
		alignments[alignments.filter_failed.notna()].to_csv(f'{output_base}failed_trimmed_alignments.tsv', mode='a', sep='\t', index=False, header=False)
		alignments = alignments[alignments.filter_failed.isna()]
		if alignments.empty is True:
			continue

		# Some genes have a bunch of transcripts that overlap the same alignment, so duplicates show up
		alignments.fivep_site = alignments.fivep_site.astype(int)	# For some reason one read that makes it this far (NGSNJ-086:229:GW200110425th:1:1101:2166:13448_for in C22-1_R1) has a float value for fivep_site instead of int
		alignments = alignments[FINAL_INFO_TABLE_COLS]
		alignments = alignments.drop_duplicates()

		# Write alignments to file
		alignments.to_csv(f'{output_base}final_info_table.tsv', mode='a', sep='\t', index=False, header=False)

	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Trimmed alignment filtering complete')

