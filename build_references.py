#!/usr/bin/env python3
import argparse
import datetime
import gzip
import os
import shutil
import pathlib

import pandas as pd

from scripts import utils





# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
EXON_INTRON_COLUMNS = ('chrom', 'strand', 'start', 'end', 'gene_id')
REF_GENOME_FILE = 'genome'	# ends up as genome.fa or genome.fa.gz
REF_FIVEP_FILE = 'fivep_sites' # ends up as fivep_sites.fa.gz, accompanied by the index files fivep_sites.fa.gz.fai and fivep_sites.fa.gz.gzi
REF_BP_FILE = 'bp_sites' # ends up as bp_sites.fa.gz, accompanied by the index files bp_sites.fa.gz.fai and bp_sites.fa.gz.gzi
REF_GENES_FILE = 'genes.tsv.gz'
REF_INTRONS_FILE = 'introns.tsv.gz'
REF_EXONS_FILE = 'exons.tsv.gz'
REF_REPEATMASKER_FILE = 'repeatmasker.bed'
REF_GENOME_UNCOMPRESSED_EXTENSIONS = ('.fa', '.fasta', '.fas', '.fna', '.ffn', '.faa', '.mpfa')
REF_GENOME_COMPRESSED_EXTENSIONS = ('.gz', '.bgzip')



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
		meta_out.write(f'Version: {utils.version()}\n')
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
	ref_names = ['Genome fasta', 'Reference annotation', 'Branchpoint annotation']
	ref_files = [args.genome_fasta, args.genome_anno, args.bp_bed]
	for rn, rf in zip(ref_names, ref_files):
		if not os.path.isfile(rf):
			raise FileNotFoundError(f'{rn} file does not exist at {rf}')
	genome_fasta, genome_anno, bp_bed = args.genome_fasta.resolve(), args.genome_anno.resolve(), args.bp_bed.resolve()

	repeatmasker_bed = args.repeatmasker_bed
	if repeatmasker_bed is not None:
		if not os.path.isfile(repeatmasker_bed):
			parser.error(f'RepeatMasker file does not exist at {repeatmasker_bed}')
		else:
			repeatmasker_bed = repeatmasker_bed.resolve()
	
	# Determine the genome fasta file format
	last_ext = genome_fasta.suffix
	if last_ext in REF_GENOME_UNCOMPRESSED_EXTENSIONS:
		genome_fasta_out = f"{args.out_dir}/{REF_GENOME_FILE}.fa"
	elif last_ext in REF_GENOME_COMPRESSED_EXTENSIONS:
		genome_fasta_out = f"{args.out_dir}/{REF_GENOME_FILE}.fa.gz"
	else:
		parser.error(f'Genome FASTA file extension is invalid ({last_ext}), must be one of' +\
			   f'{REF_GENOME_UNCOMPRESSED_EXTENSIONS}, {REF_GENOME_COMPRESSED_EXTENSIONS}')

	# Determine the annotation file format
	exts = genome_anno.suffixes
	if exts[-1] == '.gz':
		anno_type, anno_gzip = exts[-2], True
	else:
		anno_type, anno_gzip = exts[-1], False
	if not anno_type in ('.gtf', '.gff', '.gff3'):
		parser.error(f'Annotation file must be in GTF or GFF format, not {anno_type}')

	# Determine the hisat2 file extensions and confirm the files exist
	if os.path.isfile(f'{args.ref_hisat2_index}.1.ht2'):
		hisat2_extensions = [f'.{i}.ht2' for i in range(1,9)]
	elif os.path.isfile(f'{args.ref_hisat2_index}.1.ht2l'):
		hisat2_extensions = [f'.{i}.ht2l' for i in range(1,9)]
	else:
		raise FileNotFoundError(f'hisat2 index does not exist at {args.ref_hisat2_index}')
	ref_hisat2_index = args.ref_hisat2_index.resolve()

	# Confirm the STAR index exists
	if os.path.isfile(f'{args.ref_STAR_index}/SA') and os.path.isfile(f'{args.ref_STAR_index}/SAindex'):
		ref_STAR_index = args.ref_STAR_index.resolve()
	else:
		raise FileNotFoundError(f'STAR index does not exist at {args.ref_STAR_index}')
	
	# Validate threads arg
	if not args.threads>0:
		parser.error(f'-t/--threads must be a positive integer. Input was "{repr(args.threads)}"')
		
	return genome_fasta, genome_anno, bp_bed, repeatmasker_bed, ref_STAR_index, ref_hisat2_index, hisat2_extensions, anno_type, anno_gzip, genome_fasta_out


def parse_attributes(attribute_string:str, file_type:str) -> dict:
	if file_type == '.gtf':
		attributes = attribute_string.rstrip('";').split('; ')
		attributes = [attr.split(' ') for attr in attributes]
		tags = [attr_val.strip('"') for attr_name, attr_val in attributes if attr_name=='tag']
		attributes = {attr_name: attr_val.strip('"') for attr_name, attr_val in attributes if attr_name!='tag'}
		attributes['tags'] = tags
	
	elif file_type in ('.gff', '.gff3'):
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

	transcripts, genes = {}, {}
	for line in in_file:
		if line=='\n' or line.startswith('#'):
			continue
		
		# Parse the line, making sure it's an exon feature
		try:
			chrom, _, feature, start, end, _, strand, _, attributes = line.strip().split('\t')
		except ValueError:
			raise ValueError(f'Line: {repr(line)}')
		
		if feature == 'gene':
			attributes = parse_attributes(attributes, anno_type)
			genes[attributes[g_attr]] = (chrom, strand, int(start)-1, int(end))
		elif feature == 'exon':
			start, end = int(start)-1, int(end)
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

	return transcripts, genes

def build_genes(genes:dict, out_dir:str, log) -> None:

	with gzip.open(f'{out_dir}/{REF_GENES_FILE}', 'wt') as out_file:
		out_file.write('gene_id\tchrom\tstrand\tstart\tend\n')
		for gene_id in sorted(genes.keys(), key=lambda g:(genes[g][0], genes[g][2], genes[g][1])):
			chrom, strand, start, end = genes[gene_id]
			out_file.write(f'{gene_id}\t{chrom}\t{strand}\t{start}\t{end}\n')


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
	exons = (exons
		  .groupby(['chrom', 'strand', 'start', 'end'], as_index=False)
		  .agg({'gene_id': lambda gids: utils.str_join(gids, unique=True)})
	)
	# Sort to make debugging easier
	exons = exons.sort_values(['chrom', 'start', 'end'])
	# Write to file
	exons.to_csv(f'{out_dir}/{REF_EXONS_FILE}', sep='\t', index=False, compression='gzip')

	# Collapse gene ids
	introns = pd.DataFrame(introns, columns=EXON_INTRON_COLUMNS)
	introns = (introns
			.groupby(['chrom', 'strand', 'start', 'end'], as_index=False)
			.agg({'gene_id': lambda gids: utils.str_join(gids, unique=True)})
	)
	# Remove introns 20bp or shorter
	log.info(f'{sum(introns.end-introns.start<20):,} of {len(introns):,} introns excluded for being shorter than 20nt')
	introns = introns.loc[introns.end-introns.start>=20]
	# Sort to make debugging easier
	introns = introns.sort_values(['chrom', 'start', 'end'])
	# Write to file
	introns.to_csv(f'{out_dir}/{REF_INTRONS_FILE}', sep='\t', index=False, compression='gzip')

	return introns

def build_bp(bp_bed:pathlib.Path, genome_fasta:str, threads:int, out_dir:str, log) -> None:

	bp_file = gzip.open(bp_bed, 'rt') if str(bp_bed).endswith('gz') else open(bp_bed)
	bedtools_input = ''
	bp_window_size = 20
	for line in bp_file:
		chrom, bp_pos, _, _, _, strand = line.strip().split('\t')[:6]
		bp_site = f'{chrom};{bp_pos};{strand}'
		#bp_window_start = int(bp_pos)-24 if strand=='+' else int(bp_pos)
		bp_window_start = int(bp_pos)-bp_window_size-4 if strand=='+' else int(bp_pos)+5
		bedtools_input += f'{chrom}\t{bp_window_start}\t{bp_window_start+bp_window_size}\t{bp_site}\t0\t{strand}\n'
	bp_file.close()
	
	# Get fivep sequences using bedtools getfasta
	bp_seqs = utils.getfasta(genome_fasta, bedtools_input, log=log)
	bp_seqs = bp_seqs.rename(columns={'name': 'bp_site'})

	# Parse output
	bp_seqs[['chrom', 'pos', 'strand']] = bp_seqs.bp_site.str.split(';', expand=True)
	bp_seqs.pos = bp_seqs.pos.astype('int')
	# Sort to make debugging easier
	bp_seqs = bp_seqs.sort_values(['chrom', 'pos', 'strand'])

	# Write sequences to fasta file
	with open(f'{out_dir}/{REF_BP_FILE}.fa', 'w') as fasta_out:
		for i, row in bp_seqs.iterrows():
			fasta_out.write(f">{row['bp_site']}\n{row['seq']}\n")

	# Compress
	log.debug('Compressing branchpoint upstream sequences FASTA...')
	compressed_out = f"{out_dir}/{REF_BP_FILE}.fa.gz"
	if os.path.isfile(compressed_out):
		log.warning(f'Overwriting existing file at {compressed_out}')
		os.remove(compressed_out)
	utils.run_command(f'bgzip --threads {threads} {out_dir}/{REF_BP_FILE}.fa', log=log)
	
	# Index
	log.debug("Building branchpoint upstream sequences FASTA index...")
	utils.run_command(f'samtools faidx {compressed_out} --fai-idx {compressed_out}.fai --gzi-idx {compressed_out}.gzi', log=log)


def build_fivep(introns:pd.DataFrame, genome_fasta:str, threads:int, out_dir:str, log) -> None:
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
	fivep_seqs = utils.getfasta(genome_fasta, bedtools_input, log=log)
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

	# Compress
	log.debug('Compressing five-prime splice sites FASTA...')
	compressed_out = f"{out_dir}/{REF_FIVEP_FILE}.fa.gz"
	if os.path.isfile(compressed_out):
		log.warning(f'Overwriting existing file at {compressed_out}')
		os.remove(compressed_out)
	utils.run_command(f'bgzip --threads {threads} {out_dir}/{REF_FIVEP_FILE}.fa', log=log)
	
	# Index
	log.debug("Building 5' splice sites FASTA index...")
	utils.run_command(f'samtools faidx {compressed_out} --fai-idx {compressed_out}.fai --gzi-idx {compressed_out}.gzi', log=log)





#=============================================================================#
#                                    Main                                     #
#=============================================================================#
if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='build_references',
								  	description='Create a directory with the required reference files for running LariatMapper. A reference directory is specific to one reference genome.')
	parser.add_argument('-v', '--version', action='version', version=f'LariatMapper {utils.version()}', help='print the version id and exit')
	
	# Required arguments
	parser.add_argument('-f', '--genome_fasta', type=pathlib.Path, required=True, help='FASTA file of the reference genome. May be compressed with bgzip, but not with standard gzip.')
	parser.add_argument('-a', '--genome_anno', type=pathlib.Path, required=True, help='GTF or GFF file of the reference gene annotation. May be gzip-compressed')
	parser.add_argument('-b', '--bp_bed', type=pathlib.Path, required=True, help='BED 6+ file of branchpoint annotation. May be gzip-compressed')
	parser.add_argument('-S', '--ref_STAR_index', type=pathlib.Path, required=True, help='STAR index of the reference genome')
	parser.add_argument('-H', '--ref_hisat2_index', type=pathlib.Path, required=True, help='HISAT2 index of the reference genome')
	parser.add_argument('-o', '--out_dir', type=pathlib.Path, required=True, help='Output reference directory. Will be created if it does not exist at runtime')
	# Optional arguments
	optional_args = parser.add_argument_group(title='Optional arguments')
		# Experimentally-revelant options 
	optional_args.add_argument('-r', '--repeatmasker_bed', type=pathlib.Path, help='Path to BED file with RepeatMasker annotation of reference genome')
	optional_args.add_argument('-g', '--g_attr', default='gene_id', help='The attribute in the annotation file that uniquely identifies each gene. Each exon feature must have this attribute. (Default = gene_id)',)
	optional_args.add_argument('-x', '--t_attr', default='transcript_id', help='The attribute in the annotation file that uniquely identifies each transcript. Each exon feature must have this attribute. (Default = transcript_id)',)
		# Output options
	optional_args.add_argument('-c', '--copy', action='store_true', help='Create copies of the input files in out_dir (Default = create symbolic links)')
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
	log = utils.get_logger(log_level)

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
	genome_fasta, genome_anno, bp_bed, repeatmasker_bed, ref_STAR_index, ref_hisat2_index, hisat2_extensions, anno_type, anno_gzip, genome_fasta_out = process_args(args, parser, log)

	# Unpack args
	out_dir, threads, copy, t_attr, g_attr = args.out_dir, args.threads, args.copy, args.t_attr, args.g_attr

	# Make dir
	os.makedirs(out_dir, exist_ok=True)

	# Write metadata file
	log.info('Recording metadata...')
	write_metadata(args, out_dir)

	# Link or copy the neccesary input files
	# Copying
	ref_hisat2_index_prefix = str(ref_hisat2_index).removeprefix(f'{ref_hisat2_index.parent}/')
	if copy is True:
		log.info('Copying input files...')
		if os.path.isfile(genome_fasta_out):
			log.warning(f'Overwriting existing file at {genome_fasta_out}')
		shutil.copyfile(genome_fasta, genome_fasta_out, follow_symlinks=False)

		ref_hisat2_index_dir = os.path.join(out_dir, 'hisat2_index')
		os.makedirs(ref_hisat2_index_dir, exist_ok=True)
		for ext in hisat2_extensions:
			out_file = os.path.join(ref_hisat2_index_dir, f'genome{ext}')
			if os.path.isfile(out_file):
				log.warning(f'Overwriting existing file at {out_file}')
			shutil.copyfile(f'{ref_hisat2_index}{ext}', out_file, follow_symlinks=False)
		
		ref_STAR_index_dir = os.path.join(out_dir, 'STAR_index')
		if os.path.exists(ref_STAR_index_dir):
			log.warning(f'Overwriting existing STAR index directory at {ref_STAR_index_dir}')
		shutil.copytree(ref_STAR_index, ref_STAR_index_dir, symlinks=False, dirs_exist_ok=True)

		if repeatmasker_bed is not None:
			if os.path.isfile(f'{out_dir}/{REF_REPEATMASKER_FILE}'):
				log.warning(f'Overwriting existing file at {out_dir}/{REF_REPEATMASKER_FILE}')
			shutil.copyfile(repeatmasker_bed, f'{out_dir}/{REF_REPEATMASKER_FILE}', follow_symlinks=False)	
	# Linking
	else:
		log.info('Creating links to input files...')
		if os.path.isfile(genome_fasta_out):
			log.warning(f'Overwriting existing symbolic link at {genome_fasta_out}')
			os.remove(genome_fasta_out) 
		os.symlink(genome_fasta, genome_fasta_out)

		ref_hisat2_index_dir = os.path.join(out_dir, 'hisat2_index')
		os.makedirs(ref_hisat2_index_dir, exist_ok=True)
		for ext in hisat2_extensions:
			out_file = os.path.join(ref_hisat2_index_dir, f'genome{ext}')
			if os.path.isfile(out_file):
				log.warning(f'Overwriting existing symbolic link at {out_file}')
				os.remove(out_file)
			os.symlink(f'{ref_hisat2_index}{ext}', out_file)
		
		ref_STAR_index_dir = os.path.join(out_dir, 'STAR_index')
		if os.path.exists(ref_STAR_index_dir):
			log.warning(f'Overwriting existing symbolic link at {ref_STAR_index_dir}')
			os.remove(ref_STAR_index_dir)
		os.symlink(ref_STAR_index, ref_STAR_index_dir, target_is_directory=True)

		if repeatmasker_bed is not None:
			if os.path.isfile(f'{out_dir}/{REF_REPEATMASKER_FILE}'):
				log.warning(f'Overwriting existing symbolic link at {out_dir}/{REF_REPEATMASKER_FILE}')
				os.remove(f'{out_dir}/{REF_REPEATMASKER_FILE}')
			os.symlink(repeatmasker_bed, f'{out_dir}/{REF_REPEATMASKER_FILE}')	

	log.debug('Building genome FASTA index...')
	command = f'samtools faidx {genome_fasta_out} --fai-idx {genome_fasta_out}.fai'
	if genome_fasta_out.endswith('.gz'):
		command += f" --gzi-idx {genome_fasta_out}.gzi"
	utils.run_command(command, log=log)

	log.info('Parsing transcripts from annotation file...')
	transcripts, genes = parse_transcripts(genome_anno, anno_type, anno_gzip, t_attr, g_attr, log)
	
	log.info('Processing genes...')
	build_genes(genes, out_dir, log)

	log.info('Processing exons and introns...')
	introns = build_exons_introns(transcripts, out_dir, log)

	log.info("Processing 5' splice sites...")
	build_fivep(introns, genome_fasta_out, threads, out_dir, log)

	log.info("Processing branchpoints...")
	build_bp(bp_bed, genome_fasta_out, threads, out_dir, log)

	if not args.skip_r:
		log.info('Building R objects...')
		cmd = f"Rscript {pipeline_dir}/scripts/build_R_refs.R" +\
				f" --anno {genome_anno} --g_attr {g_attr} --t_attr {t_attr} --output {out_dir}"
		utils.run_command(cmd, log=log)

	log.info('Reference building complete.')
