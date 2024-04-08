import sys
import os
import shutil
import subprocess
import time

import pandas as pd



# =============================================================================#
#                                  Constants                                   #
# =============================================================================#
BAD_GENE_TYPES = ('pseudogene', 
				  'processed_pseudogene', 
				  'polymorphic_pseudogene', 
				  'transcribed_processed_pseudogene', 
				  'transcribed_unprocessed_pseudogene', 
				  'transcribed_unitary_pseudogene', 
				  'translated_processed_pseudogene', 
				  'translated_unprocessed_pseudogene', 
				  'unitary_pseudogene',
				  'unprocessed_pseudogene',
				  'artifact',
				  )



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

	return transcripts


def build_exons_introns(ref_gtf:str, out_exons_bed:str, out_introns_bed:str) -> None:
	'''
	
	'''
	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Processing exons and introns...')
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
	exons['gene_type'] = exons.attributes.transform(lambda attrs: attrs['gene_type'])
	exons['transcript_id'] = exons.attributes.transform(lambda attributes: attributes['transcript_id'])
	exons['num'] = exons.attributes.transform(lambda attributes: int(attributes['exon_number']))
	exons['transcript_and_num'] = exons.transcript_id + '-' + exons.num.astype(str)
	exons['terminal_exon'] = exons.transcript_id.map(exons.groupby('transcript_id').num.agg('max'))
	exons['terminal_exon'] = exons.num == exons.terminal_exon

	# Filter out unwanted exons
	exons = exons.loc[~exons.gene_type.isin(BAD_GENE_TYPES)]

	# Copy exons DF at this point to derive introns from it
	pos_introns = exons.loc[exons.strand=='+'].copy()
	neg_introns = exons.loc[exons.strand=='-'].copy()

	# Collapse exons that have the same genomic coordinates but different annotated transcripts
	# We'll treat such exons as 1 exon
	exons = exons.groupby(['chrom', 'start', 'end', 'strand'], as_index=False).transcript_and_num.agg('|'.join)
	exons['bed_name'] = exons.apply(lambda row: f"exon;{row['chrom']};{row['start']};{row['end']};{row['strand']};{row['transcript_and_num']}", axis=1)
	
	# Write exons to file
	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Writing exons to BED file...')
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

	# Collapse introns that have the same genomic coordinates but different annotated transcripts
	# We'll treat such introns as 1 intron
	introns = introns.groupby(['chrom', 'start', 'end', 'strand'], as_index=False).transcript_and_num.agg('|'.join)
	introns['bed_name'] = introns.apply(lambda row: f"intron;{row['chrom']};{row['start']};{row['end']};{row['strand']};{row['transcript_and_num']}", axis=1)

	# Write introns to file
	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Writing introns to BED file...')
	introns['zero'] = 0
	introns[['chrom', 'start', 'end', 'bed_name', 'zero', 'strand']].to_csv(out_introns_bed, sep='\t', header=False, index=False)

	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Exon and intron files done.')


def build_fivep(ref_gtf:str, ref_fasta:str, out_dir:str, out_fivep_fasta:str, out_fivep_upstream:str) -> None:
	'''
	
	'''
	print(time.strftime('%m/%d/%y - %H:%M:%S') + " | Processing 5' splice sites...")

	# Load the introns BED we just built
	fiveps = pd.read_csv('/home/tmooney/Lariat_mapping/testing/new_references/introns.bed', sep='\t', names=['chrom', 'start', 'end', 'name', 'zero', 'strand'])
	fiveps = fiveps.loc[fiveps.end-fiveps.start>=20]
	
	transcripts = parse_transcript_info(ref_gtf)

	# Get gene ids
	fiveps['transcript_and_num'] = fiveps.name.transform(lambda name: name.split(';')[-1])
	fiveps['transcript_id'] = fiveps.transcript_and_num.transform(lambda transcript_and_num: [t_n.split('-')[0] for t_n in transcript_and_num.split('|')])
	fiveps = pd.merge(fiveps.explode('transcript_id'), transcripts[['transcript_id', 'gene_id']], on='transcript_id')
	fiveps = fiveps[['chrom', 'start', 'end', 'strand', 'gene_id']].drop_duplicates()
	fiveps = fiveps.groupby(['chrom', 'start', 'end', 'strand'], as_index=False).gene_id.agg('|'.join)
	# fiveps = fiveps.groupby(['chrom', 'start', 'end', 'strand'], as_index=False).gene_id.agg(','.join)

	# Prep rest of info
	fiveps['pos'] = fiveps.apply(lambda row: row['start'] if row['strand']=='+' else row['end']-1, axis=1) 
	fiveps['slice_start'] = fiveps.apply(lambda row: row['pos']-5 if row['strand']=='+' else row['pos']-19, axis=1) 
	fiveps['slice_end'] = fiveps.apply(lambda row: row['pos']+20 if row['strand']=='+' else row['pos']+6, axis=1) 
	fiveps['fivep_site'] = fiveps.apply(lambda row: f"{row['chrom']};{row['pos']};{row['strand']};{row['gene_id']}", axis=1)
	# fiveps['name'] = fiveps.apply(lambda row: f"{row['chrom']};{row['pos']};{row['strand']}", axis=1)
	fiveps['zero'] = 0

	# Write to temp BED file
	temp_bed = f'{out_dir}/temp.bed'
	temp_seqs = f'{out_dir}/temp.tsv'
	fiveps[['chrom', 'slice_start', 'slice_end', 'fivep_site', 'zero', 'strand']].to_csv(temp_bed, sep='\t', header=False, index=False)
	# fiveps[['chrom', 'slice_start', 'slice_end', 'name', 'zero', 'strand']].to_csv(temp_bed, sep='\t', header=False, index=False)

	# Get sequence from 5nt upstream to 20nt downstream of 5'ss
	command = f"bedtools getfasta -nameOnly -s -tab -fi {ref_fasta} -bed {temp_bed} -fo {temp_seqs}"
	out = subprocess.run(command.split(' '), capture_output=True)
	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | ' + out.stdout.decode() + out.stderr.decode())

	# Load 5'ss seqs
	fiveps_seqs = pd.read_csv(temp_seqs, sep='\t', names=['fivep_site', 'seq']).drop_duplicates()
	# fiveps_seqs = pd.read_csv(temp_seqs, sep='\t', names=['name', 'seq']).drop_duplicates()
	fiveps_seqs['intron_sequence'] = fiveps_seqs.seq.str.slice(5)
	fiveps_seqs['upstream_sequence'] = fiveps_seqs.seq.str.slice(0, 5)
	fiveps_seqs.fivep_site = fiveps_seqs.fivep_site.str.slice(0, -3)
	# fiveps_seqs.name = fiveps_seqs.name.str.slice(0, -3)
	# fiveps_seqs = pd.merge(fiveps_seqs, fiveps[['name', 'chrom', 'pos', 'strand', 'gene_id']], on='name', validate='one_to_one')

	print(time.strftime('%m/%d/%y - %H:%M:%S') + " | Writing 5'ss sequences to FASTA file...")
	# Write 5'ss sequence to FASTA file 
	with open(out_fivep_fasta, 'w') as w:
			for i, row in fiveps_seqs.iterrows():
				w.write(f">{row['fivep_site']}\n{row['intron_sequence']}\n")

	print(time.strftime('%m/%d/%y - %H:%M:%S') + " | Writing 5nt upstream sequences to TSV file...")
	# Write 5'ss 5nt upstream sequence to TSV file
	fiveps_seqs[['fivep_site', 'upstream_sequence']].to_csv(out_fivep_upstream, sep='\t', index=False)
	with open(out_fivep_upstream, 'w') as w:
		for i, row in fiveps_seqs.iterrows():
			w.write(f"{row['fivep_site']}\t{row['upstream_sequence']}\n")
	# print(time.strftime('%m/%d/%y - %H:%M:%S') + " | Writing 5'ss information to TSV file...")
	# # Write 5'ss info to TSV file
	# fiveps_seqs[['name', 'chrom', 'pos', 'strand', 'gene_id', 'upstream_sequence']].to_csv(out_fivep_upstream, sep='\t', index=False)

	# Delete temp files
	os.remove(temp_bed)
	os.remove(temp_seqs)

	print(time.strftime('%m/%d/%y - %H:%M:%S') + " | 5'ss files done.")


def build_bowtie2_index(ref_fasta:str, out_bowtie2_index:str, threads) -> None:
	'''
	
	'''
	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Building bowtie2 index...')
	# command = f'bowtie2-build {ref_fasta} {out_bowtie2_index}'.split(' ')
	command = f'bowtie2-build --threads {threads} {ref_fasta} {out_bowtie2_index}'.split(' ')
	out = subprocess.run(command, capture_output=True)
	print(time.strftime('%m/%d/%y - %H:%M:%S') + out.stdout.decode() + out.stderr.decode())

	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | bowtie2 index done.')


# def build_hisat2_index(ref_fasta:str, out_hisat2_index:str, threads) -> None:
# 	'''
	
# 	'''
# 	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Building hisat2 index...')
# 	command = f'hisat2-build dwdwdwdw--threads {threads} {ref_fasta} {out_hisat2_index}'.split(' ')
# 	out = subprocess.run(command, capture_output=True)
# 	print(time.strftime('%m/%d/%y - %H:%M:%S') + out.stdout.decode() + out.stderr.decode())

# 	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | hisat2 index done.')



#=============================================================================#
#                                    Main                                     #
#=============================================================================#
def main():
	threads, ref_fasta, ref_gtf, ref_repeatmasker, out_dir = sys.argv[1:]

	# Make dir
	if not os.path.isdir(out_dir):
		os.mkdir(out_dir)

	# Copy input reference files to dir
	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Copying input reference files...')
	shutil.copyfile(ref_fasta, f'{out_dir}/genome.fa')
	shutil.copyfile(ref_gtf, f'{out_dir}/annotation.gtf')
	shutil.copyfile(ref_repeatmasker, f'{out_dir}/repeatmasker.bed')
	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Copying done.')

	# Create remaining reference files
	out_exons_bed = f'{out_dir}/exons.bed'
	out_introns_bed = f'{out_dir}/introns.bed'
	build_exons_introns(ref_gtf, out_exons_bed, out_introns_bed)

	out_fivep_fasta = f'{out_dir}/fivep_sites.fa'
	out_fivep_upstream = f'{out_dir}/fivep_sites_upstream.tsv'
	build_fivep(ref_gtf, ref_fasta, out_dir, out_fivep_fasta, out_fivep_upstream)

	out_bowtie2_index = f'{out_dir}/bowtie2_index'
	build_bowtie2_index(ref_fasta, out_bowtie2_index, threads)

	# out_hisat2_index = f'{out_dir}/hisat2_index'
	# build_hisat2_index(ref_fasta, out_hisat2_index, threads)

	print('\n' + time.strftime('%m/%d/%y - %H:%M:%S') + ' | Reference-building complete')





if __name__ == '__main__':
	main()