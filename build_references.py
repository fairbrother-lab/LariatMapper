import sys
import os
import shutil
import subprocess

import pandas as pd



#=============================================================================#
#                                  Functions                                  #
#=============================================================================#
def parse_attributes(attribute_string:str) -> dict:
	attributes = attribute_string.rstrip('";').split('; ')
	attributes = [attr.split(' ') for attr in attributes]
	tags = [attr_val.strip('"') for attr_name, attr_val in attributes if attr_name=='tag']
	attributes = {attr_name: attr_val.strip('"') for attr_name, attr_val in attributes if attr_name!='tag'}
	attributes['tags'] = tags

	return attributes


def build_transcript_info(ref_gtf:str, out_transcript_info_tsv:str):
	'''
	
	'''
	print('Processing transcripts...')
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
	transcripts['gene_id'] = transcripts.attributes.transform(lambda attributes: attributes['gene_id'])
	transcripts['gene_name'] = transcripts.attributes.transform(lambda attributes: attributes['gene_name'])
	transcripts['gene_type'] = transcripts.attributes.transform(lambda attributes: attributes['gene_type'])
	transcripts['transcript_id'] = transcripts.attributes.transform(lambda attributes: attributes['transcript_id'])
	
	# Write to file
	print('Writing transcript information to TSV file...')
	transcripts[['transcript_id', 'gene_id', 'gene_name', 'gene_type']].to_csv(out_transcript_info_tsv, sep='\t', header=False, index=False)

	print('Transcript file done')


def build_exons_introns(ref_gtf:str, out_exons_bed:str, out_introns_bed:str) -> None:
	'''
	
	'''
	print('Processing exons and introns...')
	# Count header lines in GTF
	header_lines = 0
	with open(ref_gtf) as r:
		for line in r:
			if line.startswith('##'):
				header_lines += 1
			else:
				break

	# Load GTF
	exons = pd.read_csv(ref_gtf, names=['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'], skiprows=header_lines, sep='\t')
	exons = exons.loc[exons.feature=='exon'].reset_index(drop=True)
	
	# Parse info
	exons.start = exons.start.astype(int) - 1
	exons.end = exons.end.astype(int)
	exons.attributes = exons.attributes.transform(parse_attributes)
	exons['transcript_id'] = exons.attributes.transform(lambda attributes: attributes['transcript_id'])
	exons['num'] = exons.attributes.transform(lambda attributes: int(attributes['exon_number']))
	exons['transcript_and_num'] = exons.transcript_id + '-' + exons.num.astype(str)
	exons['terminal_exon'] = exons.transcript_id.map(exons.groupby('transcript_id').num.agg('max'))
	exons['terminal_exon'] = exons.num == exons.terminal_exon

	# Copy exons DF at this point to derive introns from it
	pos_introns = exons.loc[exons.strand=='+'].copy()
	neg_introns = exons.loc[exons.strand=='-'].copy()

	# Collapse exons that have the same genomic coordinates but different annotated transcripts
	# We'll treat such exons as 1 exon
	exons = exons.groupby(['chrom', 'start', 'end', 'strand'], as_index=False).transcript_and_num.agg('|'.join)
	exons['bed_name'] = exons.apply(lambda row: f"exon;{row['chrom']};{row['start']};{row['end']};{row['strand']};{row['transcript_and_num']}", axis=1)
	
	# Write exons to file
	print('Writing exons to BED file...')
	exons['zero'] = 0
	exons[['chrom', 'start', 'end', 'bed_name', 'zero', 'strand']].to_csv(out_exons_bed, sep='\t', header=False, index=False)

	# Derive plus-strand introns from exons
	pos_introns = pos_introns.sort_values(['chrom', 'transcript_id', 'num'], ignore_index=True)
	pos_introns = pos_introns.rename(columns={'start': 'exon_start', 'end': 'exon_end'})
	pos_introns['start'] = pos_introns.exon_end 
	pos_introns['end'] = pos_introns.exon_start.shift(-1)
	assert pos_introns.iloc[-1]['terminal_exon']
	assert pd.isna(pos_introns.iloc[-1]['end'])
	pos_introns = pos_introns[:-1]
	pos_introns = pos_introns.loc[~pos_introns.terminal_exon]

	# Derive negative-strand introns from exons
	neg_introns = neg_introns.sort_values(['chrom', 'transcript_id', 'num'], ignore_index=True, ascending=[True, True, False])
	neg_introns = neg_introns.rename(columns={'start': 'exon_start', 'end': 'exon_end'})
	neg_introns.num += -1
	neg_introns.transcript_and_num = neg_introns.apply(lambda row: row['transcript_and_num'].split('-')[0] + '-' + str(row['num']), axis=1)
	neg_introns['start'] = neg_introns.exon_end
	neg_introns['end'] = neg_introns.exon_start.shift(-1)
	assert pd.isna(neg_introns.iloc[-1]['end'])
	neg_introns = neg_introns[:-1]
	neg_introns = neg_introns.loc[~(neg_introns.num==0)].reset_index(drop=True)

	# Add plus and negative strand introns
	introns = pd.concat([pos_introns, neg_introns], ignore_index=True)
	introns.end = introns.end.astype(int)
	assert introns[['transcript_id', 'num']].duplicated().sum() == 0

	# Collapse exons that have the same genomic coordinates but different annotated transcripts
	# We'll treat such exons as 1 exon
	introns = introns.groupby(['chrom', 'start', 'end', 'strand'], as_index=False).transcript_and_num.agg('|'.join)
	introns['bed_name'] = introns.apply(lambda row: f"intron;{row['chrom']};{row['start']};{row['end']};{row['strand']};{row['transcript_and_num']}", axis=1)

	# Write introns to file
	print('Writing introns to BED file...')
	introns['zero'] = 0
	introns[['chrom', 'start', 'end', 'bed_name', 'zero', 'strand']].to_csv(out_introns_bed, sep='\t', header=False, index=False)

	print('Exon and intron files done')


def build_fivep(ref_fasta:str, out_dir:str, out_fivep_fasta:str, out_fivep_upstream:str) -> None:
	'''
	
	'''
	print("Processing 5' splice sites...")

	# Load GTF
	fiveps = pd.read_csv('/home/tmooney/Lariat_mapping/testing/new_references/introns.bed', sep='\t', names=['chrom', 'start', 'end', 'name', '_', 'strand'])
	fiveps = fiveps.loc[fiveps.end-fiveps.start>=20]

	# Process info
	fiveps['transcript_nums'] = fiveps.name.transform(lambda name: name.split(';')[-1])
	fiveps['transcript_id'] = fiveps.transcript_nums.transform(lambda transcript_nums: [t_n.split('-')[0] for t_n in transcript_nums.split('|')])
	fiveps['pos'] = fiveps.apply(lambda row: row['start'] if row['strand']=='+' else row['end']-1, axis=1) 
	fiveps['slice_start'] = fiveps.apply(lambda row: row['pos']-5 if row['strand']=='+' else row['pos']-20, axis=1) 
	fiveps['slice_end'] = fiveps.apply(lambda row: row['pos']+19+1 if row['strand']=='+' else row['pos']+4+1, axis=1) 
	fiveps['name'] = fiveps.apply(lambda row: f"{row['chrom']};{row['pos']};{row['strand']};{'|'.join(row['transcript_id'])}", axis=1)
	fiveps['zero'] = 0

	# Write to temp BED file
	temp_bed = f'{out_dir}/temp.bed'
	temp_seqs = f'{out_dir}/temp.tsv'
	fiveps[['chrom', 'slice_start', 'slice_end', 'name', 'zero', 'strand']].to_csv(temp_bed, sep='\t', header=False, index=False)

	# Get sequence from 5nt upstream to 20nt downstream of 5'ss
	command = f"bedtools getfasta -nameOnly -s -tab -fi {ref_fasta} -bed {temp_bed} -fo {temp_seqs}"
	out = subprocess.run(command.split(' '), capture_output=True)
	print(out.stdout.decode())

	# Load 5'ss seqs
	fiveps_seqs = pd.read_csv(temp_seqs, sep='\t', names=['name', 'seq'])
	fiveps_seqs['in_intron'] = fiveps_seqs.seq.str.slice(5)
	fiveps_seqs['in_exon'] = fiveps_seqs.seq.str.slice(0,5)
	fiveps_seqs.name = fiveps_seqs.name.str.slice(0, -3)

	print("Writing 5'ss sequences to FASTA file...")
	# Write 5'ss sequence to FASTA file 
	with open(out_fivep_fasta, 'w') as w:
			for i, row in fiveps_seqs.iterrows():
				w.write(f">{row['name']}\n{row['in_intron']}\n")

	print("Writing 5nt upstream sequences to TSV file...")
	# Write 5'ss 5nt upstream sequence to TSV file
	with open(out_fivep_upstream, 'w') as w:
		for i, row in fiveps_seqs.iterrows():
			w.write(f"{row['name']}\t{row['in_exon']}\n")

	# Delete temp files
	os.remove(temp_bed)
	os.remove(temp_seqs)

	print("5'ss files done")


def build_bowtie2_index(ref_fasta:str, out_index:str) -> None:
	'''
	
	'''
	print('Building bowtie2 index...')
	command = f'bowtie2-build {ref_fasta} {out_index}'.split(' ')
	# command = f'bowtie2-build --threads 4 {ref_fasta} {out_index}'.split(' ')
	out = subprocess.run(command, capture_output=True)
	print(out.stdout.decode())

	print('bowtie2 index done.')


def build_hisat2_index(ref_fasta:str, out_index:str) -> None:
	'''
	
	'''
	print('Building hisat2 index...')
	command = f'hisat2-build {ref_fasta} {out_index}'.split(' ')
	out = subprocess.run(command, capture_output=True)
	print(out.stdout.decode())

	print('hisat2 index done.')



#=============================================================================#
#                                    Main                                     #
#=============================================================================#
def main():
	ref_fasta, ref_gtf, ref_repeatmasker, out_dir = sys.argv[1:]

	# Make dir
	if not os.path.isdir(out_dir):
		os.mkdir(out_dir)

	# # Copy input reference files to dir
	# print('Copying input reference files...')
	# shutil.copyfile(ref_fasta, f'{out_dir}/genome.fa')
	# shutil.copyfile(ref_gtf, f'{out_dir}/annotation.gtf')
	# shutil.copyfile(ref_repeatmasker, f'{out_dir}/repeatmasker.bed')

	# # Create remaining reference files
	# out_exons_bed = f'{out_dir}/exons.bed'
	# out_introns_bed = f'{out_dir}/introns.bed'
	# build_exons_introns(ref_gtf, out_exons_bed, out_introns_bed)

	# out_fivep_fasta = f'{out_dir}/fivep_sites.fa'
	# out_fivep_upstream = f'{out_dir}/fivep_sites_upstream.tsv'
	# build_fivep(ref_gtf, out_dir, out_fivep_fasta, out_fivep_upstream)

	out_index = f'{out_dir}/index'
	build_bowtie2_index(ref_fasta, out_index)

	print('Reference-building complete')





if __name__ == '__main__':
	main()



# input        
	#     ref_fasta=$OPTARG
	#     ref_gtf=$OPTARG 
	#     ref_repeatmasker=$OPTARG 

# output
#    ref_b2index=$OPTARG 
        #     ref_5p_fasta=$OPTARG 
        #     ref_5p_upstream=$OPTARG 
        #     ref_introns=$OPTARG
        #     ref_exons=$OPTARG


# t_info = functions.get_annotations(GENCODE_v44_COMP_ANNOTATIONS_FILE, feature='transcript')
# t_info = {anno['transcript_id']: [anno['gene_id'], anno['gene_name'], anno['gene_type']] for anno in t_info}

# feature_type = 'intron'
# features = load_bed(f'/home/tmooney/Lariat_mapping/testing/references/hg38.gencode.v44.comprehensive.{feature_type}s.bed', make_end_inclusive=False)

# features = features.rename(columns={'feat': 'chrom'})
# features = features[features.chrom.isin(HUMAN_CHROMOSOMES)]
# features['tid'] = features.name.transform(lambda name: name.split('_')[0])
# features['num'] = features.name.transform(lambda name: name.split('_')[2])
# features[['gene_id', 'gene_name', 'gene_type']] = features.tid.map(t_info).to_list()
# features['mod_name'] = features.apply(lambda row: f"{feature_type};{row['tid']};{row['num']};{row['start']};{row['end']};{row['strand']};{row['gene_id']};{row['gene_name']};{row['gene_type']}", axis=1)

# print(features.at[0, 'mod_name'])
# features

# features[['chrom', 'start', 'end', 'mod_name', 'score', 'strand']].to_csv(f'/home/tmooney/Lariat_mapping/testing/references/hg38.gencode.v44.comprehensive.{feature_type}s_mod.bed', sep='\t', index=False, header=False)
# features[['chrom', 'start', 'end', 'mod_name', 'score', 'strand']]