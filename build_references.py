#!/usr/bin/env python3
import argparse
import datetime
import gzip
import os
import shutil

import pyfaidx
import pandas as pd

from scripts import functions



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
EXON_INTRON_COLUMNS = ('chrom', 'strand', 'start', 'end', 'gene_id')
REF_GENOME_FILE = 'genome'
REF_HISAT2_INDEX = 'hisat2_index'
REF_EXONS_FILE = 'exons.tsv.gz'
REF_INTRONS_FILE = 'introns.tsv.gz'
REF_FIVEP_FILE = 'fivep_sites'
REF_REPEATMASKER_FILE = 'repeatmasker.bed'
HISAT2_EXTENSIONS = ('1.ht2', '2.ht2', '.3.ht2', '.4.ht2','.5.ht2', '.6.ht2','.7.ht2', '.8.ht2')



#=============================================================================#
#                                  Functions                                  #
#=============================================================================#
def write_metadata(args:argparse.Namespace, out_dir:str):
	now = datetime.datetime.now(datetime.timezone.utc)
	time = now.ctime()

	with open(f'{out_dir}/meta.txt', 'w') as meta_out:
		meta_out.write(
		"----------------------------------------\n"
		"                Metadata                \n"
		"----------------------------------------\n"
		)
		meta_out.write(f'Version: {functions.version()}\n')
		meta_out.write(f'Run time (UTC): {time}\n')

		meta_out.write(
		"\n"
		"----------------------------------------\n"
		"                Settings                \n"
		"----------------------------------------\n"
		)
		for key, val in vars(args).items():
			val = '' if val is None else val
			meta_out.write(f'{key}: {val}\n')


def process_args(args:argparse.Namespace, parser:argparse.ArgumentParser, log):
	# Confirm that input files exist
	ref_names = ['Genome fasta', 'Reference annotation']
	ref_files = [args.genome_fasta, args.genome_anno]
	for rn, rf in zip(ref_names, ref_files):
		if not os.path.isfile(rf):
			raise FileNotFoundError(f'{rn} file does not exist at {rf}')
	if args.repeatmasker_bed is not None and not os.path.isfile(args.repeatmasker_bed):
		parser.error(f'RepeatMasker file does not exist at {args.repeatmasker_bed}')

	# Determine the genome fasta file format
	if args.genome_fasta.endswith('.gz'):
		fasta_ext = '.fa.gz'
	else:
		fasta_ext = '.fa'

	# Determine the annotation file format
	prev_ext, last_ext = args.genome_anno.split('.')[-2:]
	if last_ext == 'gz':
		anno_type, anno_gzip = prev_ext, True
	else:
		anno_type, anno_gzip = last_ext, False
	if not anno_type in ('gtf', 'gff', 'gff3'):
		parser.error(f'Annotation file must be in .gtf or .gff format, not .{anno_type}')

	# Determine the hisat2 file extensions and confirm the files exist
	if os.path.isfile(f'{args.hisat2_index}.1.ht2'):
		hisat2_extensions = [f'.{i}.ht2' for i in range(1,9)]
	elif os.path.isfile(f'{args.hisat2_index}.1.ht2l'):
		hisat2_extensions = [f'.{i}.ht2l' for i in range(1,9)]
	else:
		raise FileNotFoundError(f'hisat2 index does not exist at {args.hisat2_index}')
	
	# Validate threads arg
	if not args.threads>0:
		parser.error(f'-t/--threads must be a positive integer. Input was "{repr(args.threads)}"')
	
	return hisat2_extensions, anno_type, anno_gzip, fasta_ext


def parse_attributes(attribute_string:str, file_type:str) -> dict:
	if file_type == 'gtf':
		attributes = attribute_string.rstrip('";').split('; ')
		attributes = [attr.split(' ') for attr in attributes]
		tags = [attr_val.strip('"') for attr_name, attr_val in attributes if attr_name=='tag']
		attributes = {attr_name: attr_val.strip('"') for attr_name, attr_val in attributes if attr_name!='tag'}
		attributes['tags'] = tags
	
	elif file_type in ('gff', 'gff3'):
		attributes = [attr.split('=') for attr in attribute_string.rstrip(';').split(';')]
		attributes = [(attr[0].lstrip(), attr[1]) for attr in attributes]
		attributes = dict(attributes)
		if 'tag' in attributes:
			attributes['tags'] = attributes['tag'].split(',')

	else:
		raise ValueError(f'Invalid file type "{file_type}"')

	return attributes


def parse_transcripts(genome_anno:str, anno_type:str, anno_gzip:bool, t_attr:str, g_attr:str, log):
	if anno_gzip:
		in_file = gzip.open(genome_anno, 'rt')
	else:
		in_file = open(genome_anno)

	transcripts = {}
	for line in in_file:
		if line=='\n' or line.startswith('#'):
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

		if g_attr not in attributes:
			raise ValueError(f'Attribute for gene id "{g_attr}" not found in attributes {repr(attributes)} of line "{line}"')
		if t_attr not in attributes:
			raise ValueError(f'Attribute for transcript id "{t_attr}" not found in attributes {repr(attributes)}of line "{line}"')
		
		# Add the exon to its transcript's list
		exon = (start, end)
		transcript_id = attributes[t_attr]
		if transcript_id not in transcripts:
			transcripts[transcript_id] = {'chrom': chrom, 'strand': strand, 'gene': attributes[g_attr], 'exons': []}
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
	with open(f'{out_dir}/{REF_FIVEP_FILE}.fa', 'w') as fasta_out:
		for i, row in fivep_seqs.iterrows():
			fasta_out.write(f">{row['fivep_site']}\n{row['seq']}\n")

	# 	print(fasta_ext)
	# Compress
	log.debug('Compressing five-prime splice sites FASTA...')
	functions.run_command(f'bgzip --threads {threads} {out_dir}/{REF_FIVEP_FILE}.fa', log=log)
	


#=============================================================================#
#                                    Main                                     #
#=============================================================================#
if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='build_references',
								  	description='Create a directory with the required reference files for running LariatMapper. A reference directory is specific to one reference genome.')
	parser.add_argument('-v', '--version', action='version', version=f'LariatMapper {functions.version()}', help='print the version id and exit')
	
	# Required arguments
	parser.add_argument('-f', '--genome_fasta', required=True, help='FASTA file of the reference genome')
	parser.add_argument('-a', '--genome_anno', required=True, help='GTF or GFF file of the reference gene annotation. May be gzip-compressed')
	parser.add_argument('-i', '--hisat2_index', required=True, help='HISAT2 index of the reference genome')
	parser.add_argument('-o', '--out_dir', required=True, help='Output reference directory. Will be created if it does not exist at runtime')
	# Optional arguments
	optional_args = parser.add_argument_group(title='Optional arguments')
		# Experimentally-revelant options 
	optional_args.add_argument('-r', '--repeatmasker_bed', help='Path to BED file with RepeatMasker annotation of reference genome')
	optional_args.add_argument('-g', '--g_attr', default='gene_id', help='The attribute in the annotation file that uniquely identifies each gene. Each exon feature must have this attribute. (Default = gene_id)',)
	optional_args.add_argument('-x', '--t_attr', default='transcript_id', help='The attribute in the annotation file that uniquely identifies each transcript. Each exon feature must have this attribute. (Default = transcript_id)',)
		# Output options
	optional_args.add_argument('-c', '--copy', action='store_true', help='Create deep copies of the input files in out_dir (Default = create symbolic links)')
		# Technical options
	optional_args.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use. (Default = 1)')
	log_levels = optional_args.add_mutually_exclusive_group()
	log_levels.add_argument('-q', '--quiet', action='store_true', help="Only print fatal error messages. Mutually exclusive with -w and -d")
	log_levels.add_argument('-w', '--warning', action='store_true', help="Print warning messages and fatal error messages. Mutually exclusive with -q and -d")
	log_levels.add_argument('-d', '--debug', action='store_true', help="Print extensive status messages. Mutually exclusive with -q and -w")
	# Just for development
	# We use this to skip the build_r_refs.R call when testing the python script
	# since the custom-made test annotation files just don't work in the R script
	parser.add_argument('--skip_r', action='store_true', help=argparse.SUPPRESS)

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
	arg_message = '\n'.join(arg_message)
	log.info(f'Arguments: \n{arg_message}')

	# Validate the args and determine additional variables
	hisat2_extensions, anno_type, anno_gzip, fasta_ext = process_args(args, parser, log)

	genome_fasta, genome_anno, repeatmasker_bed, hisat2_index, out_dir, threads, copy, t_attr, g_attr  = args.genome_fasta, args.genome_anno, args.repeatmasker_bed, args.hisat2_index, args.out_dir, args.threads, args.copy, args.t_attr, args.g_attr

	# Make dir
	if not os.path.isdir(out_dir):
		os.mkdir(out_dir)

	# Write metadata file
	log.info('Recording metadata...')
	write_metadata(args, out_dir)

	# Link or copy the neccesary input files
	if copy is True:
		log.info('Copying input files...')
		shutil.copyfile(genome_fasta, f'{out_dir}/{REF_GENOME_FILE}{fasta_ext}', follow_symlinks=False)
		for ext in hisat2_extensions:
			shutil.copyfile(f'{hisat2_index}{ext}', f'{out_dir}/{REF_HISAT2_INDEX}{ext}', follow_symlinks=False)
		if repeatmasker_bed is not None:
			shutil.copyfile(repeatmasker_bed, f'{out_dir}/{REF_REPEATMASKER_FILE}', follow_symlinks=False)	
	else:
		log.info('Creating links to input files...')
		if not os.path.isfile(f'{out_dir}/{REF_GENOME_FILE}{fasta_ext}'):
			os.symlink(genome_fasta, f'{out_dir}/{REF_GENOME_FILE}{fasta_ext}')
		for ext in hisat2_extensions:
			if not os.path.isfile(f'{out_dir}/{REF_HISAT2_INDEX}{ext}'):
				os.symlink(f'{hisat2_index}{ext}', f'{out_dir}/{REF_HISAT2_INDEX}{ext}')
		if repeatmasker_bed is not None and not os.path.isfile(f'{out_dir}/{REF_REPEATMASKER_FILE}'):
			os.symlink(repeatmasker_bed, f'{out_dir}/{REF_REPEATMASKER_FILE}')	

	log.info('Parsing transcripts from annotation file...')
	transcripts = parse_transcripts(genome_anno, anno_type, anno_gzip, t_attr, g_attr, log)
	
	log.info('Processing exons and introns...')
	introns = build_exons_introns(transcripts, out_dir, log)

	log.info("Processing five-prime splice sites...")
	build_fivep(introns, genome_fasta, threads, out_dir, log)

	log.info('Building FASTA indices...')
	functions.run_command(f'samtools faidx {out_dir}/{REF_GENOME_FILE}{fasta_ext}', log=log)
	functions.run_command(f'samtools faidx {out_dir}/{REF_FIVEP_FILE}.fa.gz', log=log)

	if not args.skip_r:
		log.info('Building R objects...')
		cmd = f"Rscript {pipeline_dir}/scripts/build_R_refs.R" +\
				f" --anno {genome_anno} --g_attr {g_attr} --t_attr {t_attr} --output {out_dir}"
		functions.run_command(cmd, log=log)

	log.info('Reference building complete.')
