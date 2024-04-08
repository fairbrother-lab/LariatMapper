import sys
import os
import gzip
import subprocess
import time
import site

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
ALIGNMENTS_INITIAL_COLUMNS = ('align_num', 'read_id', 'read_is_reverse', 'trim_seq', 'chrom', 'align_is_reverse', 'align_start', 'align_end', 'introns')
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


def parse_attributes(attribute_string:str) -> dict:
	attributes = attribute_string.rstrip('";').split('; ')
	attributes = [attr.split(' ') for attr in attributes]
	tags = [attr_val.strip('"') for attr_name, attr_val in attributes if attr_name=='tag']
	attributes = {attr_name: attr_val.strip('"') for attr_name, attr_val in attributes if attr_name!='tag'}
	attributes['tags'] = tags

	return attributes

def parse_transcript_info(ref_gtf:str):
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


def parse_trimmed_alignments(bam_file:str):
	# alignments = {}
	alignments = []
	i = 0
	for align in pysam.AlignmentFile(bam_file, 'rb'):
		i += 1
		# if align.query_name not in alignments:
			# alignments[align.query_name] = []

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
					introns]
		
		# Check alignment quality and fail alignment if it doesn't meet all thresholds
		if mismatches > MAX_MISMATCHES:
			continue
		if mismatch_percent > MAX_MISMATCH_PERCENT:
			continue
		if len(gaps) > 1 or (len(gaps) == 1 and gaps[0] > MAX_GAP_LENGTH):
			continue
		
		# If the alignment passed the quality thresholds, add it to the list
		# alignments[align.query_name].append(align_info)
		alignments.append(align_info)
	
	# Flatten the passed alignments into a list of lists, then covert to a DataFrame and return
	# alignments = [align_info for aligns_list in alignments.values() for align_info in aligns_list]
	alignments = pd.DataFrame(alignments, columns=ALIGNMENTS_INITIAL_COLUMNS)
	return alignments


def reverse_complement(seq):
	return ''.join([COMP_NTS[seq[i]] for i in range(len(seq)-1,-1,-1)])


def add_read_bp_nt(row:pd.Series) -> str:
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


def get_bp_seqs(alignments:pd.DataFrame, output_base:str, genome_fasta:str, keep_intermediates:bool):
	temp_bp_bed = output_base + 'temp_bp_seqs.bed'
	temp_bp_seq = output_base + 'temp_bp_seqs.tsv'
	alignments['window_start'] = alignments.apply(lambda row: row['bp_pos']-4 if row['strand']=='+' else row['bp_pos']-5, axis=1)
	alignments['window_end'] = alignments.apply(lambda row: row['bp_pos']+6 if row['strand']=='+' else row['bp_pos']+5, axis=1)
	alignments['zero'] = 0
	alignments[['chrom', 'window_start', 'window_end', 'align_num', 'zero', 'strand']].drop_duplicates().to_csv(temp_bp_bed, sep='\t', index=False, header=False)

	command = f'bedtools getfasta -nameOnly -s -tab -fi {genome_fasta} -bed {temp_bp_bed} -fo {temp_bp_seq}'
	subprocess.run(command.split(' '))
	bp_seqs = pd.read_csv(temp_bp_seq, sep='\t', names=['align_num', 'genomic_bp_context'])
	bp_seqs.align_num = bp_seqs.align_num.transform(lambda align_num: int(align_num[:-3]))
	alignments = pd.merge(alignments, bp_seqs, on='align_num')

	if keep_intermediates is False:
		os.remove(temp_bp_bed)
		os.remove(temp_bp_seq)

	return alignments


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

	# # Check if no fivep match
	# if pd.isna(row['fivep_pos']):
	# 	return 'no_fivep_matches'

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

	# Load auxiliary data for processing alignments 
	introns = parse_intron_info(ref_introns)
	transcripts = parse_transcript_info(ref_gtf)
	fivep_info = parse_fivep_info(output_base)
	
	# Load alignments, skipping bad-quality alignments
	alignments = parse_trimmed_alignments(f'{output_base}trimmed_reads_to_genome.bam')

	# Filter out alignments that aren't within an intron
	with open(f'{output_base}failed_trimmed_alignments.tsv', 'w') as w:
		w.write('\t'.join(FAILED_ALIGNMENTS_COLS) + '\n')
	failed_alignments = alignments.loc[alignments.introns.transform(len).eq(0)].copy()
	failed_alignments = failed_alignments.assign(a='', b='', c='', d='', e='', f='', g='', filter_failed = 'within_intron')
	failed_alignments.to_csv(f'{output_base}failed_trimmed_alignments.tsv', mode='a', sep='\t', header=False, index=False)
	alignments = alignments.loc[alignments.introns.transform(len).gt(0)]

	# Explode introns
	alignments = alignments.explode('introns', ignore_index=True).rename(columns={'introns': 'intron'})
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
	failed_alignments = alignments.loc[alignments.fivep_pos.isna()].copy()
	failed_alignments = failed_alignments.assign(c='', d='', e='', f='', g='', filter_failed = 'fivep_match')
	failed_alignments = failed_alignments[[col for col in FAILED_ALIGNMENTS_COLS if col in failed_alignments.columns]].drop_duplicates()
	failed_alignments.to_csv(f'{output_base}failed_trimmed_alignments.tsv', mode='a', sep='\t', header=False, index=False)
	alignments = alignments.loc[alignments.fivep_pos.notna()]

	# Add more info
	alignments.fivep_pos = alignments.fivep_pos.astype(int)
	alignments['read_bp_nt'] = alignments.apply(add_read_bp_nt, axis=1)
	alignments['bp_pos'] = alignments.apply(lambda row: row['align_end']-1 if row['strand']=='+' else row['align_start'], axis=1)
	alignments['threep_pos'] = alignments.apply(add_nearest_threep, introns=introns, axis=1)
	alignments['bp_dist_to_threep'] = alignments.apply(lambda row: -abs(row['bp_pos']-row['threep_pos']) if pd.notna(row['threep_pos']) else pd.NA, axis=1)
	# alignments['genomic_bp_context'] = get_bp_seqs(alignments, output_base, genome_fasta, keep_intermediates)
	alignments = get_bp_seqs(alignments, output_base, genome_fasta, keep_intermediates)
	alignments['genomic_bp_nt'] = alignments.genomic_bp_context.str.get(4)
	
	# Filter alignments
	alignments['filter_failed'] = alignments.apply(filter_alignments, axis=1)
	failed_alignments = alignments.loc[alignments.filter_failed.notna()]
	failed_alignments[FAILED_ALIGNMENTS_COLS].to_csv(f'{output_base}failed_trimmed_alignments.tsv', mode='a', sep='\t', header=False, index=False)
	alignments = alignments.loc[alignments.filter_failed.isna(), FINAL_INFO_TABLE_COLS]
	# print(alignments.info())	

	# Some genes have a bunch of transcripts that overlap the same alignment, so duplicates show up
	alignments = alignments.drop_duplicates()

	# Write alignments to file
	alignments.to_csv(f'{output_base}final_info_table.tsv', sep='\t', index=False)

	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Trimmed alignment filtering complete\n')

