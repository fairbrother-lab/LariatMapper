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
	

ALIGNMENT_COLUMNS = ('read_id', 'read_is_reverse', 'chrom', 'align_is_reverse', 'align_start', 'align_end', 'quality', 'mismatches', 'mismatch_percent', 'gaps', 'introns')
def get_read_aligns(bam_file):
	'''
	Generator
	'''
	alignments = []
	current_read_id = next(pysam.AlignmentFile(bam_file)).query_name

	for align in pysam.AlignmentFile(bam_file, 'rb'):
		# Parse alignment info
		read_is_reverse = True if align.query_name.endswith('_rev') else False
		mismatches = align.get_tag('XM')
		mismatch_percent = mismatches/align.query_length
		gaps = [length for tag_num, length in align.cigartuples if tag_num in (1,2)]
		if align.has_tag('YB'):
			introns = align.get_tag('YB').split(',')
			introns = [intron.split(';') for intron in introns]
		else:
			introns = []

		# Collect the info in a list which will become a DataFrame row
		next_align = [align.query_name, 
					read_is_reverse, 
					align.reference_name, 
					align.is_reverse, 
					align.reference_start, 
					align.reference_end, 
					align.mapping_quality, 
					mismatches, 
					mismatch_percent, 
					gaps, 
					introns]
		
		# If still processing alignments from the same read, add next_align to alignments
		if align.query_name == current_read_id:
			alignments.append(next_align)

		# If reached alignments for a new read, yield the accumulated alignments and reset 
		else:
			trim_seq = align.get_forward_sequence()
			alignments = pd.DataFrame(alignments, columns=ALIGNMENT_COLUMNS)
			yield current_read_id, trim_seq, alignments

			current_read_id = align.query_name
			alignments = [next_align]


MAX_MISMATCHES = 5
MAX_MISMATCH_PERCENT = 0.1
MAX_GAP_LENGTH = 3
def first_filters(row:pd.Series, max_quality:int):
	# Check if mapping quality is not max
	if row['quality'] < max_quality:
		return 'map_quality'

	# Check if too many mismatches
	if row['mismatches'] > MAX_MISMATCHES:
		return 'mismatch_num'
	if row['mismatch_percent'] > MAX_MISMATCH_PERCENT:
		return 'mismatch_percent'
	
	# Check if more than one insertion/deletion OR one insertion/deletion that's too long
	if len(row['gaps']) > 1 or (len(row['gaps']) == 1 and row['gaps'][0] > MAX_GAP_LENGTH):
		return 'indel'

	# Check if the 5'ss alignment orientation doesn't match the genomic alignment orientation
	if row['read_is_reverse'] != row['align_is_reverse']:
		return 'same_orientation'
	
	# Check if not within an intron
	if row['introns'] == []:
		return 'within_intron'
	
	return pd.NA

COMP_NTS = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
def reverse_complement(seq):
	return ''.join([COMP_NTS[seq[i]] for i in range(len(seq)-1,-1,-1)])


def add_read_bp_nt(row:pd.Series, trim_seq:str) -> str:
	if row['strand']=='+' and not row['align_is_reverse']:
		return trim_seq[-1]
	elif row['strand']=='+' and row['align_is_reverse']:
		return reverse_complement(trim_seq[0])
	elif row['strand']=='-' and not row['align_is_reverse']:
		return reverse_complement(trim_seq[0])
	elif row['strand']=='-' and row['align_is_reverse']:
		return trim_seq[-1]
	

def add_nearest_threep(row:pd.Series, introns:dict):
	overlap_introns = list(introns[row['chrom']][row['strand']].overlap(row['bp_site'], row['bp_site']+1))
	if len(overlap_introns) == 0:
		raise RuntimeError(f"No introns overlapped the following:\n{row}")
	
	if row['strand'] == '+':
		threep_site = min(overlap_introns, key=lambda s: s.end-row['bp_site']-1).end
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
	# ref_gtf, ref_introns, genome_fasta, output_base = sys.argv[1:]
	ref_gtf = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.gencode.v44.basic.annotation.gtf'
	ref_transcripts = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.gencode.v44.basic.transcripts.bed'
	ref_introns = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.gencode.v44.basic.introns_mod.bed'
	genome_fasta	= '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.gencode.v44.primary_assembly.fa'
	output_base = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/output/C22_R1_lariat_mapping/'


	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Beginning trimmed alignment filtering')
	print(ref_gtf, ref_introns, genome_fasta, output_base)

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
	i = 0
	for read_id, trim_seq, alignments in get_read_aligns('/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/output/C22_R1_lariat_mapping/trimmed_reads_to_genome.bam'):
		i += 1
		if i % 1_000 == 0:
			print(time.strftime('%m/%d/%y - %H:%M:%S') + f' | {i:,} reads')
			print(alignments.info())

		max_quality = alignments.quality.max()
		
		alignments['filter_failed'] = alignments.apply(first_filters, max_quality=max_quality, axis=1)
		alignments[alignments.filter_failed.notna()].to_csv(f'{output_base}failed_trimmed_alignments.tsv', mode='a', sep='\t', index=False, header=False)
		alignments = alignments[alignments.filter_failed.isna()]
		if alignments.empty is True:
			continue
		
		# Explode introns
		alignments = alignments.explode('introns', ignore_index=True)
		alignments[['transcript_id', 'intron_num', 'intron_start', 'intron_end', 'strand', 'gene_id', 'gene_name', 'gene_type']] = alignments.introns.to_list()
		alignments.intron_start = alignments.intron_start.astype(int)
		alignments.intron_end = alignments.intron_end.astype(int)

		# 
		alignments['read_bp_nt'] = alignments.apply(add_read_bp_nt, trim_seq=trim_seq, axis=1)
		alignments['bp_site'] = alignments.apply(lambda row: row['align_end'] if row['strand']=='+' else row['align_start'], axis=1)
		# alignments['threep_site'] = alignments.apply(lambda row: row['intron_start'] if row['strand']=='+' else row['intron_end'], axis=1)
		alignments['threep_site'] = alignments.apply(add_nearest_threep, introns=introns, axis=1)
		alignments['bp_dist_to_threep'] = alignments.apply(lambda row: row['bp_site']-row['threep_site'], axis=1)
		alignments['genomic_bp_context'] = alignments.apply(add_bp_seq, output_base=output_base, genome_fasta=genome_fasta, axis=1)
		alignments['genomic_bp_nt'] = alignments.genomic_bp_context.str.get(4)

		# Merge
		alignments = pd.merge(alignments, fivep_info, 'left', on=['read_id', 'gene_id'])

		# Second filters
		alignments['filter_failed'] = alignments.apply(second_filters, axis=1)
		alignments[alignments.filter_failed.notna()].to_csv(f'{output_base}failed_trimmed_alignments.tsv', mode='a', sep='\t', index=False, header=False)
		alignments = alignments[alignments.filter_failed.isna()]

		# Some genes have a bunch of transcripts that overlap the same alignment, so duplicates show up
		alignments = alignments[FINAL_INFO_TABLE_COLS]
		alignments = alignments.drop_duplicates()
		if alignments.empty is False:
			print(alignments.info())
			print(alignments)

		# Write alignments to file
		alignments.to_csv(f'{output_base}final_info_table.tsv', mode='a', sep='\t', index=False, header=False)

	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Trimmed alignment filtering complete')





















# import sys
# import os
# import gzip
# import subprocess
# import time

# import pysam
# import pandas as pd
# from intervaltree import Interval, IntervalTree



# # =============================================================================#
# #                                  Constants                                   #
# # =============================================================================#
# COMP_NTS = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
# MAX_MISMATCHES = 5
# MAX_MISMATCH_PERCENT = 0.1
# MAX_GAP_LENGTH = 3
# TRANSCRIPT_OVERLAPS_COLUMNS = ('chrom', 'align_start', 'align_end', 'read_id', '_', 'align_orient', '_', '_', '_', '_', '_', '_', '_', 'transcript_start', 'transcript_end', 'transcript_id', '_', 'transcript_strand', '_', '_', '_', '_', '_', '_')
# ALIGNMENTS_INITIAL_COLUMNS = ('read_id',
# 							'trimmed_seq',
# 							'read_is_reverse',
# 							'chrom',
# 							'align_start',
# 							'align_end',
# 							'align_is_reverse',
# 							'mapping_quality',
# 							'mismatches',
# 							'mismatch_percent',
# 							'gap_lengths')

# FINAL_INFO_TABLE_COLS = ['read_id', 
# 						'read_is_reverse', 
# 						'read_seq', 
# 						'chrom', 
# 						'strand', 
# 						'align_start',
# 						'align_end', 
# 						'gene_id', 
# 						'gene_name', 
# 						'gene_type', 
# 						'fivep_site', 
# 						'bp_site', 
# 						'read_bp_nt', 
# 						'genomic_bp_context', 
# 						'genomic_bp_nt', 
# 						'threep_site', 
# 						'bp_dist_to_threep']



# # =============================================================================#
# #                                  Functions                                   #
# # =============================================================================#
# def reverse_complement(seq):
# 	return ''.join([COMP_NTS[seq[i]] for i in range(len(seq)-1,-1,-1)])


# def load_alignments(output_base:str) -> dict:
# 	rows = []
# 	for align in pysam.AlignmentFile(f'{output_base}trimmed_reads_to_genome.bam', 'rb'):	
# 		if align.query_name.endswith('_for'):
# 			read_is_reverse = False
# 		elif align.query_name.endswith('_rev'):
# 			read_is_reverse = True
# 		else:
# 			raise RuntimeError(f'{align}')

# 		for tag, val in align.get_tags():
# 			if tag == 'XM':
# 				mismatches = val
# 				break
# 		mismatch_percent = mismatches/align.query_length

# 		gap_lengths = []
# 		for tag_num, length in align.cigartuples:
# 			if tag_num in (1, 2):
# 				gap_lengths.append(length)
# 			assert tag_num in (0, 1, 2), align

# 		rows.append([align.query_name, align.query_sequence, read_is_reverse, align.reference_name, align.reference_start, align.reference_end, align.is_reverse, align.mapping_quality, mismatches, mismatch_percent, gap_lengths])

# 	alignments = pd.DataFrame(rows, columns=ALIGNMENTS_INITIAL_COLUMNS)
	
# 	max_quality = alignments.groupby('read_id').mapping_quality.aggregate('max')
# 	max_quality.name = 'max_quality'
# 	alignments = pd.merge(alignments, max_quality, how='left', left_on='read_id', right_index=True)
	
# 	return alignments


# def parse_gtf_line(line:str) -> dict:
# 	'''
# 	Returns {'Chromosome': str, 'Source': str, 'Feature': str, 'Start': int, 'End': int, 'Score': str, 'Strand': str, 'Frame': str,
# 			[attribute 1]: [attribute 1 value], [attribute 2]: [attribute 2 value]...}
# 	'''
# 	info = {}
# 	chrom, source, feature, start, end, score, strand, frame, attributes = line.strip().split('\t')

# 	for key, val in (('Chromosome', chrom),
# 					('Source', source),
# 					('Feature', feature),
# 					('Start', int(start)),
# 					('End', int(end)),
# 					('Score', score),
# 					('Strand', strand),
# 					('Frame', frame)
# 					):
# 		info[key] = val

# 	attributes = attributes.split('; ')
# 	attribute_dir = {attr.split(' ')[0]: attr.split(' ')[1].strip('"') for attr in attributes if not attr.startswith('tag "')}
# 	tags = [attr.lstrip('tag "').rstrip('"') for attr in attributes if attr.startswith('tag "')]
# 	attribute_dir['tags'] = tuple(tags)

# 	info.update(attribute_dir)

# 	return info


# def get_annotations(gtf_path:str, number = 'all', chromosome : str = None, feature : str = None, start : str = None, end : str = None, strand : str = None, none_except : bool = True) -> list:
# 	'''
# 	number = 'all' or an int >= 1
# 	'''
# 	assert isinstance(number, (str, int)), f'{number}, {type(number)}'
# 	assert number == 'all' or number >= 1, number
# 	if [chromosome, feature, start, end, strand] == [None]*5:
# 		raise RuntimeError(f'You must specify at least one argument from chromosome, start, end, strand, or feature')

# 	annotations = []
# 	with open(gtf_path, 'r') as gtf:
# 		for line in gtf:
# 			if line.startswith('##'):
# 				continue

# 			line_info = parse_gtf_line(line)
			
# 			if feature is not None and line_info['Feature'] != feature:
# 				continue
# 			if start is not None and line_info['Start'] != start:
# 				continue
# 			if end is not None and line_info['Snd'] != end:
# 				continue
# 			if strand is not None and line_info['Strand'] != strand:
# 				continue

# 			annotations.append(line_info)
# 			if number != 'all' and len(annotations) == number:
# 				break
	
# 	if annotations == [] and none_except:
# 		raise RuntimeError(f'Annotations not found in {gtf_path} matching args "{chromosome} {feature} {start} {end} {strand}"')

# 	return annotations


# def process_exon_lengths(row: pd.Series):
# 	out = [s for s in row['exon_lengths'].split(',') if s != '']
# 	out = [int(s) for s in out]
# 	return out

# def process_exon_starts(row: pd.Series):
# 	out = [s for s in row['exon_starts'].split(',') if s != '']
# 	out = [int(s)+row['transcript_start'] for s in out]
# 	out = [int(s) for s in out]
# 	return out

# def add_exon_ends(row:pd.Series):
# 	out = [row['exon_starts'][i]+row['exon_lengths'][i] for i in range(int(row['exon_count']))]
# 	return out

# def add_fiveps(row: pd.Series):
# 	if row['transcript_strand'] == '+':
# 		return row['exon_ends'][:-1]
# 	else:
# 		return row['exon_starts'][1:]

# def add_gtf_info(row:pd.Series, gtf_info:dict):
# 	transcript = gtf_info[row['transcript_id']]
# 	return [transcript['gene_id'], transcript['gene_name'],transcript['gene_type']]


# def load_transcript_overlaps(output_base:str, ref_transcripts:str, transcript_gtf_info:dict):
# 	transcript_overlaps = pd.read_csv(f'{output_base}transcript_overlaps.bed', sep='\t', header=None)
# 	transcript_overlaps.columns = TRANSCRIPT_OVERLAPS_COLUMNS
	
# 	# transcript_overlaps.align_end = transcript_overlaps.align_end - 1			# Make both start and end 0-based inclusive 
# 	# transcript_overlaps.transcript_end = transcript_overlaps.transcript_end - 1
# 	transcript_overlaps['align_is_reverse'] = transcript_overlaps.align_orient.map({'+': False, '-': True})
# 	transcript_overlaps = transcript_overlaps.drop(columns=['_', 'align_orient'])

# 	# We have to get the exon counts, lengths, and starts from the reference file because "bedtools intersect" truncates output lines for transcripts with lots of exons 
# 	transcripts_bed = pd.read_csv(ref_transcripts, sep='\t', skiprows=1, header=None)
# 	transcripts_bed.columns = ('_', '_', '_', 'transcript_id', '_', '_', '_', '_', '_', 'exon_count', 'exon_lengths', 'exon_starts')
# 	transcripts_bed = transcripts_bed.drop(columns='_')
# 	transcript_overlaps = pd.merge(transcript_overlaps, transcripts_bed, 'left', on='transcript_id')

# 	transcript_overlaps.exon_lengths = transcript_overlaps.apply(process_exon_lengths, axis=1)
# 	transcript_overlaps.exon_starts = transcript_overlaps.apply(process_exon_starts, axis=1)
# 	transcript_overlaps['exon_ends'] = transcript_overlaps.apply(add_exon_ends, axis=1)
# 	# transcript_overlaps['transcript_fiveps'] = transcript_overlaps.apply(add_fiveps, axis=1)
# 	transcript_overlaps[['gene_id', 'gene_name', 'gene_type']] = transcript_overlaps.apply(add_gtf_info, axis=1, gtf_info=transcript_gtf_info, result_type='expand')

# 	return transcript_overlaps


# def load_gene_ranges(gtf_file):
# 	'''
# 	Read through the GTF annotation file, return a dict with all annotated gene coordinates grouped by chromosome, then strand
# 	Output format = { chromosome: {strand: IntervalTree(gene start position, gene end position, gene name) } }
# 	'''
# 	# open file
# 	if gtf_file[-2:] == 'gz':
# 		in_file = gzip.open(gtf_file, 'rt')
# 	else:
# 		in_file = open(gtf_file)

# 	gene_ranges = {}		
# 	for line in in_file:
# 		if line[0] != '#':
# 			chrom, _, feat, start, end, _, strand, _, attributes = line.strip().split('\t')
# 			if feat == 'gene':
# 				# Convert attributes string to dict of {tag name: value}
# 				attributes = attributes[:-1].split('; ')
# 				attributes = {a.split(' ')[0]:a.split(' ')[1].replace('\"', '') for a in attributes}
# 				if chrom not in gene_ranges:
# 					gene_ranges[chrom] = {s:IntervalTree() for s in ('-', '+')}
# 				gene_ranges[chrom][strand].add(Interval(int(start)-1, int(end), {'gene_id':attributes['gene_id']}))

# 	in_file.close()
# 	return gene_ranges


# def add_gene_ids(row:pd.Series, gene_ranges):
# 	return [g.data['gene_id'] for g in gene_ranges[row['chrom']][row['strand']].overlap(row['fivep_site'], row['fivep_site']+1)]


# def load_fivep_info(output_base:str, ref_gtf:str):
# 	fivep_info = pd.read_csv(f'{output_base}fivep_info_table.tsv', sep='\t')

# 	# Unpack fivep_sites
# 	fivep_info.fivep_sites = fivep_info.fivep_sites.str.split(',')
# 	fivep_info = fivep_info.explode('fivep_sites')
# 	fivep_info[['chrom', 'fivep_site', 'fivep_end', 'strand']] = fivep_info.fivep_sites.str.split(';', expand=True)
# 	fivep_info.fivep_site = fivep_info.fivep_site.astype(int)
# 	fivep_info.fivep_end = fivep_info.fivep_end.astype(int)

# 	# Add gene ids
# 	gene_ranges = load_gene_ranges(ref_gtf)
# 	fivep_info['gene_id'] = fivep_info.apply(add_gene_ids, gene_ranges=gene_ranges, axis=1, result_type='reduce')
# 	fivep_info = fivep_info.explode('gene_id')

# 	return fivep_info

# # add gene ranges, change fivep matching


# def load_intron_info(ref_introns:str) -> tuple:
# 	'''
# 	Returns a dict formatted as follows:
# 	{Chromosome: {Strand(+ or -): Intervaltree(StartPosition(int), EndPosition(int))}}
# 	for 3' splice sites (+/- 2 bases), 5' splice sites (+/- 2 bases), and the introns they come from (start to end)
# 	'''
# 	introns = {}

# 	if ref_introns[-2:] == 'gz':
# 		intron_file = gzip.open(ref_introns, 'rt')
# 	else:
# 		intron_file = open(ref_introns)

# 	introns_done = set()
# 	for line in intron_file:
# 		chrom, start, end, _, _, strand = line.strip().split('\t')
# 		if chrom not in introns:
# 			introns[chrom] = {s: IntervalTree() for s in ['+', '-']}
# 		intron_id = '{}_{}_{}_{}'.format(chrom, strand, start, end)
# 		if intron_id not in introns_done:
# 			start, end = int(start), int(end)
# 			introns[chrom][strand].add(Interval(start, end))
# 			introns_done.add(intron_id)

# 	intron_file.close()
# 	return introns


# # test from https://stackoverflow.com/questions/3269434/whats-the-most-efficient-way-to-test-if-two-ranges-overlap
# def check_exon_overlap(row:pd.Series) -> bool:
# 	for i in range(len(row['exon_starts'])):
# 		exon_start = row['exon_starts'][i]
# 		exon_end = row['exon_ends'][i]
# 		if row['align_start'] <= exon_end and row['align_end'] >= exon_start:
# 			return True
		
# 	return False


# def add_read_bp_nt(row:pd.Series) -> str:
# 	if row['strand']=='+' and row['align_is_reverse']:
# 		return reverse_complement(row['trimmed_seq'][0])
# 	elif row['strand']=='+' and not row['align_is_reverse']:
# 		return row['trimmed_seq'][-1]
# 	elif row['strand']=='-' and row['align_is_reverse']:
# 		return row['trimmed_seq'][-1]
# 	elif row['strand']=='-' and not row['align_is_reverse']:
# 		return reverse_complement(row['trimmed_seq'][0])


# def merge_transcript_overlaps(alignments:pd.DataFrame, transcript_overlaps:pd.DataFrame) -> pd.DataFrame:
# 	alignments = pd.merge(alignments, transcript_overlaps, 'left', ['read_id', 'chrom', 'align_start', 'align_end', 'align_is_reverse'])

# 	failed_alignments = alignments[alignments.transcript_id.isna()].copy()
# 	failed_alignments['filter_failed'] = 'transcript_overlap'
# 	alignments = alignments[alignments.transcript_id.notna()]
	
# 	alignments = alignments.rename(columns={'transcript_strand': 'strand'})
# 	alignments.exon_count = alignments.exon_count.astype(int, errors='ignore')
# 	alignments['exon_overlap'] = alignments.apply(check_exon_overlap, axis=1)
# 	alignments['read_bp_nt'] = alignments.apply(add_read_bp_nt, axis=1)
# 	alignments['bp_site'] = alignments.apply(lambda row: row['align_end'] if row['strand']=='+' else row['align_start'], axis=1)

# 	return alignments, failed_alignments


# # def match_fivep_sites(row:pd.Series, fivep_info:pd.DataFrame) -> pd.DataFrame:
# # 	matches = fivep_info.loc[(fivep_info.read_id==row['read_id']) & (fivep_info.gene_id==row['gene_id'])]
# # 	return matches.fivep_site.to_list()


# def first_filters(row:pd.Series):
# 	# Check if mapping quality is not max
# 	if row['mapping_quality'] < row['max_quality']:
# 		return 'map_quality'

# 	# Check if too many mismatches
# 	if row['mismatches'] > MAX_MISMATCHES:
# 		return 'mismatch_num'
# 	if row['mismatch_percent'] > MAX_MISMATCH_PERCENT:
# 		return 'mismatch_percent'
	
# 	# Check if more than one insertion/deletion OR one insertion/deletion that's too long
# 	if len(row['gap_lengths']) > 1 or (len(row['gap_lengths']) == 1 and row['gap_lengths'][0] > MAX_GAP_LENGTH):
# 		return 'indel'

# 	# Check if the 5'ss alignment orientation doesn't match the genomic alignment orientation
# 	if row['read_is_reverse'] != row['align_is_reverse']:
# 		return 'same_orientation'

# 	# Check if the alignment overlaps any exons
# 	if row['exon_overlap'] is True:
# 		return 'exon_overlap'
	
# 	return pd.NA
	

# def second_filters(row:pd.Series):
# 	# Check if not one fivep match
# 	if pd.isna(row['fivep_site']):
# 		return 'no_matches'

# 	# Check if the 5'ss is at or downstream of the BP
# 	if row['strand'] == '+' and row['fivep_site'] >= row['bp_site']:
# 		return '5p_bp_order'
# 	elif row['strand'] == '-' and row['fivep_site'] <= row['bp_site']:
# 		return '5p_bp_order'

# 	return pd.NA


# def add_nearest_threep(row:pd.Series, introns:dict):
# 	overlap_introns = list(introns[row['chrom']][row['strand']].overlap(row['bp_site'], row['bp_site']+1))
# 	if len(overlap_introns) == 0:
# 		raise RuntimeError(f"No introns overlapped the following:\n{row}")
	
# 	if row['strand'] == '+':
# 		threep_site = min(overlap_introns, key=lambda s: s.end-row['bp_site']).end
# 	else:
# 		threep_site = min(overlap_introns, key=lambda s: row['bp_site']-s.begin).begin

# 	return threep_site


# def add_bp_seq(row:pd.Series, output_base:str, genome_fasta:str):
# 	temp_bp_bed, temp_bp_seq = output_base+'_temp_bp_seqs.bed', output_base+'temp_bp_seqs.txt'
# 	temp_file = open(temp_bp_bed, 'w')
# 	if row['strand'] == '+':
# 		bp_start, bp_end = row['bp_site']-4, row['bp_site']+6
# 	else:
# 		bp_start, bp_end = row['bp_site']-5, row['bp_site']+5
# 	temp_file.write(f"{row['chrom']}\t{bp_start}\t{bp_end}\t{row['chrom']};{row['bp_site']};{row['strand']}\t0\t{row['strand']}\n")
# 	temp_file.close()
# 	subprocess.run(f'bedtools getfasta -fi {genome_fasta} -bed {temp_bp_bed} -fo {temp_bp_seq} -nameOnly -s -tab'.split(' '))

# 	temp_file = open(temp_bp_seq)
# 	name, genomic_bp_window = temp_file.readline().strip().split()
# 	genomic_bp_window = genomic_bp_window.upper()
# 	temp_file.close()

# 	# os.remove(temp_bp_bed)
# 	# os.remove(temp_bp_seq)

# 	return genomic_bp_window



# # =============================================================================#
# #                                    Main                                      #
# # =============================================================================#
# if __name__ == '__main__':
# 	ref_gtf, ref_transcripts, ref_introns, genome_fasta, output_base = sys.argv[1:]
# 	# ref_gtf = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.gencode.v44.basic.annotation.gtf'
# 	# ref_transcripts = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.gencode.v44.basic.transcripts.bed'
# 	# ref_introns = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.gencode.v44.basic.introns.bed'
# 	# genome_fasta	= '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.gencode.v44.primary_assembly.fa'
# 	# output_base = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/output/HEK293T_R2_lariat_mapping/'

# 	# Load data
# 	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' load')
# 	alignments = load_alignments(output_base)
# 	transcript_info = get_annotations(ref_gtf, feature='transcript')
# 	transcript_info = {info['transcript_id']: info for info in transcript_info}
# 	transcript_overlaps = load_transcript_overlaps(output_base, ref_transcripts, transcript_info)
# 	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' more load')
# 	fivep_info = load_fivep_info(output_base, ref_gtf)
# 	introns = load_intron_info(ref_introns)

# 	# Record the number of reads that mapped
# 	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' record')
# 	mapped_rids = alignments.read_id.str.slice(0, -4)
# 	mapped_rids = len(mapped_rids.unique())
# 	with open(f'{output_base}run_data.tsv', 'a') as a:
# 		a.write(f'trimmed_mapped_reads\t{mapped_rids}\n')

# 	# Merge read_seq from fivep_info
# 	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' merge and then merge')
# 	read_seqs = fivep_info[['read_id', 'read_seq']].drop_duplicates()
# 	alignments = pd.merge(alignments, read_seqs, 'left', on='read_id')

# 	# Merge transcript overlaps from a "bedtools intersect" and move alignments that didn't overlap any transcripts to failed_alignments 
# 	alignments, failed_alignments = merge_transcript_overlaps(alignments, transcript_overlaps)

# 	# Run alignments through first set of filters and move alignments that failed a filter into failed_alignments
# 	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' filter')
# 	alignments['filter_failed'] = alignments.apply(first_filters, axis=1)
# 	failed_alignments = pd.concat([failed_alignments, alignments[alignments.filter_failed.notna()]])
# 	alignments = alignments.loc[alignments.filter_failed.isna()]

# 	# Match 5'ss's to the alignments
# 	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' match')
# 	# fivep_info = fivep_info[fivep_info.read_id.isin(alignments.read_id)]
# 	# alignments = match_fivep_sites(alignments, fivep_info)
# 	# alignments['fivep_matches'] = alignments.apply(match_fivep_sites, fivep_info=fivep_info, axis=1, result_type='reduce')
# 	# alignments['fivep_site'] = alignments.apply(lambda row: row['fivep_matches'][0] if len(row['fivep_matches'])==1 else pd.NA, axis=1)
# 	alignments = pd.merge(alignments, fivep_info[['read_id', 'gene_id', 'fivep_site', 'read_fivep_start', 'read_fivep_end']], 'left', on=['read_id', 'gene_id'])

# 	# Run alignments through second set of filters and move alignments that failed a filter into failed_alignments
# 	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' remove')
# 	# alignments.loc[alignments.fivep_site.isna(), 'filter_failed'] = 'no matches'
# 	alignments['filter_failed'] = alignments.apply(second_filters, axis=1)
# 	failed_alignments = pd.concat([failed_alignments, alignments[alignments.filter_failed.notna()]])
# 	alignments = alignments.loc[alignments.filter_failed.isna()]
# 	alignments.fivep_site = alignments.fivep_site.astype(int)

# 	# Move filter_failed column to the back for ease of reading, then write failed alignments to file
# 	failed_alignments = failed_alignments[[col for col in failed_alignments if col != 'filter_failed'] + ['filter_failed']]
# 	failed_alignments.to_csv(f'{output_base}failed_trimmed_alignments.tsv', sep='\t', index=False)

# 	# Assign the nearest downstream(?) 3'ss to the alignment
# 	# We couldn't do this before removing alignments that overlap exons since the whole alignment could be in an exon
# 	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' infer')
# 	if len(alignments) == 0:
# 		print(time.strftime('%m/%d/%y - %H:%M:%S') + '| No alignments remaining')
# 		alignments.to_csv(f'{output_base}final_info_table.tsv', sep='\t', index=False)
# 		exit()
# 	alignments['threep_site'] = alignments.apply(add_nearest_threep, axis=1, introns=introns)
# 	alignments['bp_dist_to_threep'] = alignments.apply(lambda row: abs(row['bp_site']-row['threep_site']), axis=1)

# 	# Add the BP's base identity from the reference genome and the 9nt window covering everything within +/- 4nt of the BP
# 	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' index file time')
# 	alignments['genomic_bp_context'] = alignments.apply(add_bp_seq, axis=1, output_base=output_base, genome_fasta=genome_fasta)
# 	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' index file time is over')
# 	alignments['genomic_bp_nt'] = alignments.genomic_bp_context.str.get(4)

# 	# Some genes have a bunch of transcripts that overlap the same alignment, so duplicates show up
# 	alignments = alignments[FINAL_INFO_TABLE_COLS]
# 	alignments = alignments.drop_duplicates()

# 	# Write alignments to file
# 	alignments.to_csv(f'{output_base}final_info_table.tsv', sep='\t', index=False)

# 	# Record filtered read count
# 	filtered_rids = alignments.read_id.str.slice(0, -4)
# 	filtered_rids = len(filtered_rids.unique())
# 	with open(f'{output_base}run_data.tsv', 'a') as a:
# 		a.write(f'trimmed_filtered_reads\t{filtered_rids}\n')
	
# 	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' done with trimmed sequence alignment filtering')