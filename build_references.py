#!/usr/bin/env python3

import os
import subprocess
import shutil
import time
import gzip
import argparse
import collections
import tempfile


import pandas as pd

from scripts import functions



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
EXON_INTRON_COLUMNS = ('chrom', 'strand', 'start', 'end', 'gene_id')
REF_GENOME_FILE = 'genome.fa'
REF_HISAT2_INDEX = 'hisat2_index'
REF_EXONS_FILE = 'exons.tsv.gz'
REF_INTRONS_FILE = 'introns.tsv.gz'
REF_FIVEP_FILE = 'fivep_sites.fa'
REF_FIVEP_INDEX = 'fivep_sites'
REF_REPEATMASKER_FILE = 'repeatmasker.bed'
HISAT2_EXTENSIONS = ('1.ht2', '2.ht2', '.3.ht2', '.4.ht2','.5.ht2', '.6.ht2','.7.ht2', '.8.ht2')



#=============================================================================#
#                                  Functions                                  #
#=============================================================================#
def process_args(args:argparse.Namespace, parser:argparse.ArgumentParser, log):
	# Confirm that input files exist
	ref_names = ['Genome fasta', 'Reference annotation']
	ref_files = [args.genome_fasta, args.genome_anno]
	for rn, rf in zip(ref_names, ref_files):
		if not os.path.isfile(rf):
			raise FileNotFoundError(f'{rn} file does not exist at {rf}')
	if args.repeatmasker_bed is not None and not os.path.isfile(args.repeatmasker_bed):
		parser.error(f'RepeatMasker file does not exist at {args.repeatmasker_bed}')
		
	# Determine the hisat2 file extensions and confirm the files exist
	if os.path.isfile(f'{args.hisat2_index}.1.ht2'):
		hisat2_extensions = [f'.{i}.ht2' for i in range(1,9)]
	elif os.path.isfile(f'{args.hisat2_index}.1.ht2l'):
		hisat2_extensions = [f'.{i}.ht2l' for i in range(1,9)]
	else:
		raise FileNotFoundError(f'hisat2 index does not exist at {args.hisat2_index}')
	
	# Determine the annotation file format
	prev_ext, last_ext = args.genome_anno.split('.')[-2:]
	if last_ext == 'gz':
		anno_type, gunzip = prev_ext, True
	else:
		anno_type, gunzip = last_ext, False
	if not anno_type in ('gtf', 'gff'):
		parser.error(f'Annotation file must be in .gtf or .gff format, not .{anno_type}')

	# Validate threads arg
	if not args.threads>0:
		parser.error(f'-t/--threads must be a positive integer. Input was "{repr(args.threads)}"')
	
	return hisat2_extensions, anno_type, gunzip


def parse_attributes(attribute_string:str, file_type:str) -> dict:
	if file_type == 'gtf':
		attributes = attribute_string.rstrip('";').split('; ')
		attributes = [attr.split(' ') for attr in attributes]
		tags = [attr_val.strip('"') for attr_name, attr_val in attributes if attr_name=='tag']
		attributes = {attr_name: attr_val.strip('"') for attr_name, attr_val in attributes if attr_name!='tag'}
		attributes['tags'] = tags
	elif file_type == 'gff':
		attributes = [attr.split('=') for attr in attribute_string.split('; ')]
		attributes = [(attr[0].lstrip(), attr[1]) for attr in attributes]
		attributes = dict(attributes)
		if 'tag' in attributes:
			attributes['tags'] = attributes['tag'].split(',')
	else:
		raise ValueError(f'Invalid file type "{file_type}"')

	return attributes


def parse_transcripts(genome_anno:str, anno_type:str, gunzip:bool, transcript_attribute:str, gene_attribute:str, log):
	if gunzip:
		in_file = gzip.open(genome_anno, 'rt')
	else:
		in_file = open(genome_anno)

	transcripts = {}
	for line in in_file:
		if line.startswith('#'):
			continue
		
		# Parse the line, making sure it's an exon feature
		try:
			chrom, _, feature, start, end, _, strand, _, attributes = line.strip().split('\t')
		except ValueError:
			raise ValueError(f'Line: {repr(line)}')
		
		if feature != 'exon':
			continue
		start = int(start) - 1
		end = int(end)
		attributes = parse_attributes(attributes, anno_type)

		if transcript_attribute not in attributes:
			raise ValueError(f'Attribute for transcript id "{transcript_attribute}" not found in attributes of line "{line}"')
		if gene_attribute not in attributes:
			raise ValueError(f'Attribute for gene id "{gene_attribute}" not found in attributes {repr(attributes)} of line "{line}"')
		
		# Add the exon to its transcript's list
		exon = (start, end)
		transcript_id = attributes[transcript_attribute]
		if transcript_id not in transcripts:
			transcripts[transcript_id] = {'chrom': chrom, 'strand': strand, 'gene': attributes[gene_attribute], 'exons': []}
		transcripts[transcript_id]['exons'].append(exon)

	return transcripts


def build_exons_introns(transcripts:dict, out_dir:str, log) -> pd.DataFrame:
	exons = []
	introns = []

	# For each transcripts, create its exons and introns and add them to the lists
	for transcript_id in transcripts:
		chrom = transcripts[transcript_id]['chrom']
		strand = transcripts[transcript_id]['strand']
		gene_id = transcripts[transcript_id]['gene']
		transcript_exons = sorted(transcripts[transcript_id]['exons'])
		n_exons = len(transcript_exons)

		for i in range(n_exons):
			# Add the ith exon to list
			exon = (chrom, strand, transcript_exons[i][0], transcript_exons[i][1], gene_id)
			exons.append(exon)

			# Don't create an intron after the last exon 
			if i == n_exons-1:
				break
			
			# Add the ith intron to list
			intron_start = transcript_exons[i][1]
			intron_end = transcript_exons[i+1][0]
			intron = (chrom, strand, intron_start, intron_end, gene_id)
			introns.append(intron)
	
	# Collapse gene ids
	exons = pd.DataFrame(exons, columns=EXON_INTRON_COLUMNS)
	exons = exons.groupby(['chrom', 'strand', 'start', 'end'], as_index=False).agg({'gene_id': functions.str_join})
	# Sort to make debugging easier
	exons = exons.sort_values(['chrom', 'start', 'end'])
	# Write to file
	exons.to_csv(f'{out_dir}/exons.tsv.gz', sep='\t', index=False, compression='gzip')

	# Collapse gene ids
	introns = pd.DataFrame(introns, columns=EXON_INTRON_COLUMNS)
	introns = introns.groupby(['chrom', 'strand', 'start', 'end'], as_index=False).agg({'gene_id': functions.str_join})
	# Remove introns 20bp or shorter
	log.info(f'{sum(introns.end-introns.start<20):,} of {len(introns):,} introns excluded for being shorter than 20nt')
	introns = introns.loc[introns.end-introns.start>=20]
	# Sort to make debugging easier
	introns = introns.sort_values(['chrom', 'start', 'end'])
	# Write to file
	introns.to_csv(f'{out_dir}/introns.tsv.gz', sep='\t', index=False, compression='gzip')

	return introns


def build_fivep(introns:pd.DataFrame, genome_fasta:str, threads:int, out_dir:str, log) -> None:
	#TODO: Refactor this through the "bedtools_input += fivep_line" line, it's the biggest time-consumer
	introns = [row.to_list() for i, row in introns.iterrows()]
	fivep_sites = set()
	for intron in introns:
		chrom, strand, intron_start, intron_end, gene_ids = intron
		fivep_pos = intron_start if strand=='+' else intron_end-1
		fivep_sites.add((chrom, strand, fivep_pos))
	
	# Prepare fivep coordinates input that will be passed to bedtools getfasta 
	bedtools_input = ''
	for fivep in fivep_sites:
		chrom, strand, fivep_pos = fivep
		if strand == '+':
			window_start, window_end = fivep_pos, fivep_pos+20
		else:
			window_start, window_end = fivep_pos-19, fivep_pos+1
		fivep_site = f'{chrom};{fivep_pos};{strand}'
		fivep_line = f'{chrom}\t{window_start}\t{window_end}\t{fivep_site}\t0\t{strand}\n'
		bedtools_input += fivep_line
	
	# Get fivep sequences using bedtools getfasta
	fivep_seqs = functions.getfasta(genome_fasta, bedtools_input, log=log)
	fivep_seqs = fivep_seqs.rename(columns={'name': 'fivep_site'})

	# Parse output
	fivep_seqs[['chrom', 'pos', 'strand']] = fivep_seqs.fivep_site.str.split(';', expand=True)
	fivep_seqs.pos = fivep_seqs.pos.astype('int')
	# Sort to make debugging easier
	fivep_seqs = fivep_seqs.sort_values(['chrom', 'pos', 'strand'])

	# Write sequences to fasta file
	with open(f'{out_dir}/{REF_FIVEP_FILE}', 'w') as fasta_out:
		for i, row in fivep_seqs.iterrows():
			fasta_out.write(f">{row['fivep_site']}\n{row['seq']}\n")
	
	# log.debug('Building fivep index')
	# build_index_call = f'bowtie2-build --quiet --threads {threads} {out_dir}/{REF_FIVEP_FILE} {out_dir}/{REF_FIVEP_INDEX}'
	# functions.run_command(build_index_call, log=log)



#=============================================================================#
#                                    Main                                     #
#=============================================================================#
if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='build_references',
								  	description='Create a reference directory for running LariatMapper. It will contain symbolic links to the input files and custom-built reference files.')
	
	# Required arguments
	parser.add_argument('-f', '--genome_fasta', required=True, help='Path to reference genome fasta file')
	parser.add_argument('-a', '--genome_anno', required=True, help='Path to reference gene annotation file in GTF or GFF format (may be gzipped with .gz extension)')
	parser.add_argument('-i', '--hisat2_index', required=True, help='Path to base name of hisat2 index of reference genome (i.e. everything before .1.ht2 extension)')
	parser.add_argument('-o', '--out_dir', required=True, help='Path to directory where reference files will be output (will be created if it does not exist)')
	# Optional arguments
	parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use for parallel processing (default = 1)')
	parser.add_argument('-r', '--repeatmasker_bed', help='Path to BED file with RepeatMasker annotation of reference genome')
	parser.add_argument('-c', '--copy', action='store_true', help='Create deep copies of the input files in out_dir (default = create symbolic links)')
	parser.add_argument('-x', '--transcript_attribute', default='transcript_id', help='The attribute in the annotation file that uniquely identifies each transcript. Each exon feature must have this attribute (default=transcript_id)',)
	parser.add_argument('-g', '--gene_attribute', default='gene_id', help='The attribute in the annotation file that uniquely identifies each gene. Each exon feature must have this attribute (default=gene_id)',)
	log_levels = parser.add_mutually_exclusive_group()
	log_levels.add_argument('-q', '--quiet', action='store_true', help="Don't print any status messages")
	log_levels.add_argument('-d', '--debug', action='store_true', help="Print extensive status messages")

	# Parse args
	args = parser.parse_args()

	# Setup logging
	if args.quiet is True:
		log_level = 'ERROR'
	elif args.debug is True:
		log_level = 'DEBUG'
	else:
		log_level = 'INFO'
	log = functions.get_logger(log_level)

	# Report version
	pipeline_dir = os.path.dirname(os.path.realpath(__file__))
	with open(f'{pipeline_dir}/pyproject.toml') as toml:
		for line in toml:
			if line.startswith('version = "'):
				version = line.lstrip('version = "').rstrip('"\n')
				log.info(f'LariatMapper {version}')

	# Print arguments
	arg_message = [f'{key}={val}' for key, val in vars(args).items() if val is not None and val is not False]
	arg_message = '\n\t'.join(arg_message)
	log.info(f'Arguments: \n\t{arg_message}')

	# Validate the args and determine additional variables
	hisat2_extensions, anno_type, gunzip = process_args(args, parser, log)

	genome_fasta, genome_anno, repeatmasker_bed, hisat2_index, out_dir, threads, copy, transcript_attribute, gene_attribute  = args.genome_fasta, args.genome_anno, args.repeatmasker_bed, args.hisat2_index, args.out_dir, args.threads, args.copy, args.transcript_attribute, args.gene_attribute

	# Make dir
	if not os.path.isdir(out_dir):
		os.mkdir(out_dir)

	# Link or copy the neccesary input files
	if copy is True:
		log.info('Copying input files...')
		shutil.copyfile(genome_fasta, f'{out_dir}/{REF_GENOME_FILE}')
		for ext in hisat2_extensions:
			shutil.copyfile(f'{hisat2_index}{ext}', f'{out_dir}/{REF_HISAT2_INDEX}{ext}')
		if repeatmasker_bed is not None:
			shutil.copyfile(repeatmasker_bed, f'{out_dir}/{REF_REPEATMASKER_FILE}')	
	else:
		log.info('Creating links to input files...')
		if not os.path.isfile(f'{out_dir}/{REF_GENOME_FILE}'):
			os.symlink(genome_fasta, f'{out_dir}/{REF_GENOME_FILE}')
		for ext in hisat2_extensions:
			if not os.path.isfile(f'{out_dir}/{REF_HISAT2_INDEX}{ext}'):
				os.symlink(f'{hisat2_index}{ext}', f'{out_dir}/{REF_HISAT2_INDEX}{ext}')
		if repeatmasker_bed is not None and not os.path.isfile(f'{out_dir}/{REF_REPEATMASKER_FILE}'):
			os.symlink(repeatmasker_bed, f'{out_dir}/{REF_REPEATMASKER_FILE}')	

	log.info('Parsing transcripts from annotation file...')
	transcripts = parse_transcripts(genome_anno, anno_type, gunzip, transcript_attribute, gene_attribute, log)
	
	log.info('Processing exons and introns...')
	introns = build_exons_introns(transcripts, out_dir, log)

	log.info("Processing five-prime splice sites...")
	build_fivep(introns, genome_fasta, threads, out_dir, log)

	log.info('Reference building complete.')