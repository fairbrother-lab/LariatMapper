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
ALIGNMENTS_INITIAL_COLUMNS = ('align_num', 'read_id', 'read_is_reverse', 'trim_seq', 'chrom', 'align_is_reverse', 'align_start', 'align_end', 'quality', 'introns')
FAILED_ALIGNMENTS_COLS = ['align_num', 
						  'read_id',
						  'read_is_reverse',
						  'trim_seq',
						  'chrom',
						  'align_is_reverse',
						  'align_start',
						  'align_end', 
						  'intron',
						  'fivep_pos',
						  'read_bp_nt',
						  'bp_pos',
						  'threep_pos',
						  'bp_dist_to_threep',
						  'genomic_bp_context',
						  'genomic_bp_nt',
						  'filter_failed']
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
						'fivep_pos', 
						'bp_pos', 
						'read_bp_nt', 
						'genomic_bp_context', 
						'genomic_bp_nt', 
						'threep_pos', 
						'bp_dist_to_threep']



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def parse_intron_info(ref_introns:str) -> tuple:
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


# def parse_gene_ranges(gtf_file):
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
# 	return [g.data['gene_id'] for g in gene_ranges[row['chrom']][row['strand']].overlap(row['fivep_pos'], row['fivep_pos']+1)]


def parse_attributes(attribute_string:str) -> dict:
	attributes = attribute_string.rstrip('";').split('; ')
	attributes = [attr.split(' ') for attr in attributes]
	tags = [attr_val.strip('"') for attr_name, attr_val in attributes if attr_name=='tag']
	attributes = {attr_name: attr_val.strip('"') for attr_name, attr_val in attributes if attr_name!='tag'}
	attributes['tags'] = tags

	return attributes

def parse_transcript_info(ref_gtf:str):
	'''
	
	'''
	# Count header lines in GTF
	header_lines = 0
	with open(ref_gtf) as r:
		for line in r:
			if line.startswith('##'):
				header_lines += 1
			else:
				break

	# Load GTF
	transcripts = pd.read_csv(ref_gtf, names=['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'], skiprows=header_lines, sep='\t')
	transcripts = transcripts.loc[transcripts.feature=='transcript'].reset_index(drop=True)
	
	# Pull out transcript x gene info
	transcripts.attributes = transcripts.attributes.transform(parse_attributes)
	transcripts['transcript_id'] = transcripts.attributes.transform(lambda attributes: attributes['transcript_id'])
	transcripts['gene_id'] = transcripts.attributes.transform(lambda attributes: attributes['gene_id'])
	transcripts['gene_name'] = transcripts.attributes.transform(lambda attributes: attributes['gene_name'])
	transcripts['gene_type'] = transcripts.attributes.transform(lambda attributes: attributes['gene_type'])

	return transcripts[['transcript_id', 'gene_id', 'gene_name', 'gene_type']]


# def parse_fivep_info(output_base:str):
# 	fivep_info = pd.read_csv(f'{output_base}fivep_info_table.tsv', sep='\t')

# 	# Unpack fivep_sites
# 	fivep_info.fivep_sites = fivep_info.fivep_sites.str.split(',')
# 	fivep_info = fivep_info.explode('fivep_sites')
# 	fivep_info[['chrom', 'fivep_start', 'fivep_end', 'strand']] = fivep_info.fivep_sites.str.split(';', expand=True)
# 	fivep_info['fivep_pos'] = fivep_info.apply(lambda row: int(row['fivep_start']) if row['strand']=='+' else int(row['fivep_end'])-1, axis=1)

# 	# Load gene ranges and annotate the 5'ss with gene ids
# 	fivep_info['gene_id'] = fivep_info.apply(add_gene_ids, gene_ranges=gene_ranges, axis=1, result_type='reduce')
# 	fivep_info = fivep_info.explode('gene_id')[['read_id', 'read_seq', 'fivep_pos', 'gene_id']]

# 	return fivep_info
def parse_fivep_info(output_base:str):
	fivep_info = pd.read_csv(f'{output_base}fivep_info_table.tsv', sep='\t')

	# Unpack fivep_sites
	fivep_info.fivep_sites = fivep_info.fivep_sites.str.split(',')
	fivep_info = fivep_info.explode('fivep_sites')
	fivep_info[['chrom', 'fivep_pos', 'strand', 'gene_id']] = fivep_info.fivep_sites.str.split(';', expand=True)
	fivep_info.fivep_pos = fivep_info.fivep_pos.astype(int)

	# Unpack each fivep_site's associated genes
	fivep_info.gene_id = fivep_info.gene_id.str.split('|')
	fivep_info = fivep_info.explode('gene_id')

	return fivep_info[['read_id', 'read_seq', 'fivep_pos', 'gene_id']]


# def parse_trimmed_alignments(bam_file:str, failed_alignments_file:str):
def parse_trimmed_alignments(bam_file:str):
	alignments = {}
	# failed_alignments = []
	i = 0
	for align in pysam.AlignmentFile(bam_file, 'rb'):
		i += 1
		if i % 100_000 == 0:
			print(f'{i:,} alignments parsed')

		if align.query_name not in alignments:
			alignments[align.query_name] = []

		read_is_reverse = True if align.query_name.endswith('_rev') else False
		mismatches = align.get_tag('XM')
		mismatch_percent = mismatches/align.query_length
		gaps = [length for tag_num, length in align.cigartuples if tag_num in (1,2)]
		if align.has_tag('YB'):
			introns = align.get_tag('YB').split(',')
		else:
			introns = []
		
		align_info = [i,
					align.query_name,
					read_is_reverse,
					align.query_sequence,
					align.reference_name, 
					align.is_reverse, 
					align.reference_start, 
					align.reference_end,
					align.mapping_quality,
					introns]
		
		# Check alignment quality and fail alignment if it doesn't meet all thresholds
		if mismatches > MAX_MISMATCHES:
			# failed_alignments.append(align_info + [pd.NA]*17 + ['mismatch_number'])
			continue
		if mismatch_percent > MAX_MISMATCH_PERCENT:
			# failed_alignments.append(align_info + [pd.NA]*17 + ['mismatch_percent'])
			continue
		if len(gaps) > 1 or (len(gaps) == 1 and gaps[0] > MAX_GAP_LENGTH):
			# failed_alignments.append(align_info + [pd.NA]*17 + ['indel'])
			continue
		if introns == []:
			# failed_alignments.append(align_info + [pd.NA]*17 + ['within_intron'])
			continue
		
		# If the alignment passed the quality thresholds, add it to the list
		alignments[align.query_name].append(align_info)
	
	# # Convert the failed alignments to a DataFrame and write to the failed file
	# pd.DataFrame(failed_alignments, columns=FAILED_ALIGNMENTS_COLS).to_csv(failed_alignments_file, sep='\t', index=False)

	# Flatten the passed alignments into a list of lists, then covert to a DataFrame and return
	alignments = [align_info for aligns_list in alignments.values() for align_info in aligns_list]
	alignments = pd.DataFrame(alignments, columns=ALIGNMENTS_INITIAL_COLUMNS)
	return alignments


def reverse_complement(seq):
	return ''.join([COMP_NTS[seq[i]] for i in range(len(seq)-1,-1,-1)])


def add_read_bp_nt(row:pd.Series) -> str:
	# if row['strand']=='+' and not row['align_is_reverse']:
	# 	return row['trim_seq'][-1]
	# elif row['strand']=='+' and row['align_is_reverse']:
	# 	return reverse_complement(row['trim_seq'][0])
	# elif row['strand']=='-' and not row['align_is_reverse']:
	# 	return reverse_complement(row['trim_seq'][0])
	# elif row['strand']=='-' and row['align_is_reverse']:
	# 	return row['trim_seq'][-1]
	if not row['read_is_reverse'] and not row['align_is_reverse']:
		return row['trim_seq'][-1]
	elif not row['read_is_reverse'] and row['align_is_reverse']:
		return reverse_complement(row['trim_seq'][0])
	elif row['read_is_reverse'] and not row['align_is_reverse']:
		return reverse_complement(row['trim_seq'][0])
	elif row['read_is_reverse'] and row['align_is_reverse']:
		return row['trim_seq'][-1]


def add_nearest_threep(row:pd.Series, introns:dict):
	enveloping_introns = list(introns[row['chrom']][row['strand']].at(row['bp_pos']))
	if len(enveloping_introns) == 0:
		raise RuntimeError(f"No introns enveloped the BP in the following:\n{row}")
	
	if row['strand'] == '+':
		threep_pos = min(enveloping_introns, key=lambda s: s.end-row['bp_pos']).end - 1
	else:
		threep_pos = min(enveloping_introns, key=lambda s: row['bp_pos']-s.begin).begin

	return threep_pos


# def add_bp_seq(row:pd.Series, output_base:str, genome_fasta:str, temp_bp_bed:str, temp_bp_seq:str):
# 	temp_file = open(temp_bp_bed, 'w')
# 	if row['strand'] == '+':
# 		bp_start, bp_end = row['bp_pos']-4, row['bp_pos']+6
# 	else:
# 		bp_start, bp_end = row['bp_pos']-5, row['bp_pos']+5
# 	temp_file.write(f"{row['chrom']}\t{bp_start}\t{bp_end}\t{row['chrom']};{row['bp_pos']};{row['strand']}\t0\t{row['strand']}\n")
# 	temp_file.close()
# 	subprocess.run(f'bedtools getfasta -fi {genome_fasta} -bed {temp_bp_bed} -fo {temp_bp_seq} -nameOnly -s -tab'.split(' '))

# 	temp_file = open(temp_bp_seq)
# 	name, genomic_bp_window = temp_file.readline().strip().split()
# 	genomic_bp_window = genomic_bp_window.upper()
# 	temp_file.close()

# 	return genomic_bp_window

		# for i, row in alignments.iterrows():
		# alignments['window_start'] = alignments.apply(lambda row: row['bp_pos']-4 if row['strand']=='+' else row['bp_pos']-5, axis=1)
		# alignments['window_end'] = alignments.apply(lambda row: row['bp_pos']+6 if row['strand']=='+' else row['bp_pos']+5, axis=1)
# 
def get_bp_seqs(alignments:pd.DataFrame, output_base:str, genome_fasta:str, keep_intermediates:bool):
	temp_bp_bed = output_base + 'temp_bp_seqs.bed'
	temp_bp_seq = output_base + 'temp_bp_seqs.tsv'
	with open(temp_bp_bed, 'w') as w:
		for i, row in alignments.iterrows():
			if row['strand'] == '+':
				window_start = row['bp_pos']-4
				window_end = row['bp_pos']+6
			else:
				window_start = row['bp_pos']-5
				window_end = row['bp_pos']+5

			w.write(f"{row['chrom']}\t{window_start}\t{window_end}\t{i}\t0\t{row['strand']}\n")

	command = f'bedtools getfasta -s -tab -fi {genome_fasta} -bed {temp_bp_bed} -fo {temp_bp_seq}'
	subprocess.run(command.split(' '))
	bp_seqs = pd.read_csv(temp_bp_seq, sep='\t', names=['_', 'genomic_bp_context'])

	if keep_intermediates is False:
		os.remove(temp_bp_bed)
		os.remove(temp_bp_seq)

	return bp_seqs['genomic_bp_context']


def filter_alignments(row:pd.Series):
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

	# Check if no fivep match
	if pd.isna(row['fivep_pos']):
		return 'no_fivep_matches'

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
	print(time.strftime('%m/%d/%y - %H:%M:%S') + f' | Arguments recieved: {sys.argv[1:]}')
	ref_gtf, ref_introns, genome_fasta, output_base, keep_intermediates = sys.argv[1:]
	keep_intermediates = {'True':True, 'False':False}[keep_intermediates]

	# 
	# with open(f'{output_base}final_info_table.tsv', 'w') as w:
	# 	w.write('\t'.join(FINAL_INFO_TABLE_COLS) + '\n')
	# with open(f'{output_base}failed_trimmed_alignments.tsv', 'w') as w:
	# 	w.write('\t'.join(FAILED_ALIGNMENTS_COLS) + '\n')

	# Load auxiliary data for processing alignments 
	introns = parse_intron_info(ref_introns)
	# gene_ranges = parse_gene_ranges(ref_gtf)
	transcripts = parse_transcript_info(ref_gtf)
	fivep_info = parse_fivep_info(output_base)
	
	# Load alignments, pre-filtering out bad-quality alignments and alignments not within an intron
	# failed_alignments_file = f'{output_base}failed_trimmed_alignments.tsv'
	# alignments = parse_trimmed_alignments(f'{output_base}trimmed_reads_to_genome.bam', failed_alignments_file)
	alignments = parse_trimmed_alignments(f'{output_base}trimmed_reads_to_genome.bam')

	# Explode introns
	alignments = alignments.explode('introns', ignore_index=True).rename(columns={'introns': 'intron'})
	# alignments[['_', 'transcript_id', 'intron_num', 'intron_start', 'intron_end', 'strand', 'gene_id', 'gene_name', 'gene_type']] = alignments.introns.str.split(';', expand=True)
	alignments[['_', '_', 'intron_start', 'intron_end', 'strand', 'transcript_and_nums']] = alignments.intron.str.split(';', expand=True)
	alignments.intron_start = alignments.intron_start.astype(int)
	alignments.intron_end = alignments.intron_end.astype(int)

	# Unpack "ENST00000469289.1-1|ENST00000473358.1-2" into "ENST00000469289.1" "1"
	# 														"ENST00000473358.1" "2"
	alignments.transcript_and_nums = alignments.transcript_and_nums.transform(lambda t_ns: [t_n.split('-') for t_n in t_ns.split('|')])
	alignments = alignments.explode('transcript_and_nums')
	alignments[['transcript_id', 'intron_num']] = alignments.transcript_and_nums.to_list()

	# Swap out transcript ids for gene ids to use for merging with fivep_info df
	alignments = pd.merge(alignments, transcripts, how='left', on='transcript_id')
	alignments = alignments.drop(columns=['transcript_and_nums', 'transcript_id', 'intron_num']).drop_duplicates()

	# Merge
	alignments = pd.merge(alignments, fivep_info, 'left', on=['read_id', 'gene_id'])

	# Add more info
	alignments['read_bp_nt'] = alignments.apply(add_read_bp_nt, axis=1)
	alignments['bp_pos'] = alignments.apply(lambda row: row['align_end']-1 if row['strand']=='+' else row['align_start'], axis=1)
	alignments['threep_pos'] = alignments.apply(add_nearest_threep, introns=introns, axis=1)
	alignments['bp_dist_to_threep'] = alignments.apply(lambda row: -abs(row['bp_pos']-row['threep_pos']) if pd.notna(row['threep_pos']) else pd.NA, axis=1)
	
	# Add BP seqs
	# alignments.to_csv(f'{output_base}alignments.tsv', sep='\t')
	# temp_bp_bed, temp_bp_seq = output_base+'_temp_bp_seqs.bed', output_base+'temp_bp_seqs.txt'
	# alignments['genomic_bp_context'] = alignments.apply(add_bp_seq, output_base=output_base, genome_fasta=genome_fasta, temp_bp_bed=temp_bp_bed, temp_bp_seq=temp_bp_seq, axis=1)
	# alignments['genomic_bp_context'] = alignments.apply(lambda row: get_sequence(genome_fasta, row['chrom'], row['bp_pos']-4, row['bp_pos']+5, False) if row['strand']=='+' else get_sequence(genome_fasta, row['chrom'], row['bp_pos']-5, row['bp_pos']+4, True), axis=1)
	alignments['genomic_bp_context'] = get_bp_seqs(alignments, output_base, genome_fasta, keep_intermediates)
	alignments['genomic_bp_nt'] = alignments.genomic_bp_context.str.get(4)
	# if keep_intermediates is False:
	# 	os.remove(temp_bp_bed)
	# 	os.remove(temp_bp_seq)
	
	# Filtering
	alignments['filter_failed'] = alignments.apply(filter_alignments, axis=1)
	
	# failed_alignments = alignments.loc[alignments.filter_failed.notna(), ]
	failed_alignments = alignments.loc[alignments.filter_failed.notna()].drop_duplicates()
	# failed_alignments['intron'] = failed_alignments.apply(lambda row: f"{row['transcript_id']};{row['intron_num']};{row['intron_start']};{row['intron_end']};{row['strand']};{row['gene_id']};{row['gene_name']};{row['gene_type']}", axis=1)
	# failed_alignments = failed_alignments.groupby([col for col in FAILED_ALIGNMENTS_COLS if col!='intron'], as_index=False).intron.agg(','.join)
	# failed_alignments[FAILED_ALIGNMENTS_COLS].to_csv(f'{output_base}failed_trimmed_alignments.tsv', mode='a', sep='\t', index=False, header=False)
	failed_alignments[FAILED_ALIGNMENTS_COLS].to_csv(f'{output_base}failed_trimmed_alignments.tsv', sep='\t', index=False)
	
	alignments = alignments.loc[alignments.filter_failed.isna(), FINAL_INFO_TABLE_COLS]
	alignments.fivep_pos = alignments.fivep_pos.astype(int)	

	# Some genes have a bunch of transcripts that overlap the same alignment, so duplicates show up
	alignments = alignments.drop_duplicates()

	# Write alignments to file
	alignments.to_csv(f'{output_base}final_info_table.tsv', sep='\t', index=False)

	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Trimmed alignment filtering complete\n')
