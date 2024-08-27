#!/usr/bin/env python3

import argparse
import json
import os
import time
import multiprocessing
import logging
import subprocess



from scripts import functions



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
FORBIDDEN_CHARS = ('/', '\\', '>', '<', '&', '~', '*', '?', '[', ']', ';', '|', '!', '$', "'", '"')



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def process_args(args:argparse.Namespace, parser:argparse.ArgumentParser, log:logging.Logger):
	# Determine whether input is single-end or paired-end and confirm that the read file(s) exit
	if args.read_file is not None:
		seq_type = 'single'
		if not os.path.isfile(args.read_file):
			parser.error(f'"{args.read_file}" is not an existing file')
	elif args.read_one is not None and args.read_two is not None:
		seq_type = 'paired'
		if not os.path.isfile(args.read_one):
			parser.error(f'"{args.read_one}" is not an existing file')
		if not os.path.isfile(args.read_two):
			parser.error(f'"{args.read_two}" is not an existing file')
	else:
		parser.error('Provide either -f/--read_file (for single-end read) OR -1/--read_one and -2/--read_two (for paired-end reads)')

	# Get reference files from reference dir OR from direct input, and confirm that they exist
	if args.ref_dir is not None:
		ref_files = ['hisat2_index', 'genome.fa', 'fivep_sites.fa', 'exons.tsv.gz', 'introns.tsv.gz', 'repeatmasker.bed']
		ref_h2index, ref_fasta, ref_5p_fasta, ref_exons, ref_introns, ref_repeatmasker = [os.path.join(args.ref_dir, f) for f in ref_files]
	else:
		ref_h2index, ref_fasta, ref_5p_fasta, ref_exons, ref_introns, ref_repeatmasker = args.ref_h2index, args.ref_fasta, args.ref_5p_fasta, args.ref_exons, args.ref_introns, args.ref_repeatmasker

	for file in (ref_fasta, ref_5p_fasta, ref_exons, ref_introns):
		if not os.path.isfile(file):
			parser.error(f'"{file}" is not an existing file')
	if (not os.path.isfile(ref_h2index + '.1.ht2')) and (not os.path.isfile(ref_h2index + '.1.ht2l')):
		parser.error(f'"{ref_h2index}" is not an existing hisat2 index')
	
	# Determine the output_base
	# All output files will be formatted like f"{output_base}file.ext"
	if args.output_prefix is None:
		output_prefix = ''
	else:
		for char in FORBIDDEN_CHARS:
			if char in args.output_prefix:
				parser.error(f'Illegal character in output prefix: {char}')
		output_prefix = args.output_prefix + '_'
	output_base = os.path.join(args.output_dir, output_prefix)

	# Validate threads arg
	if not args.threads>0:
		parser.error(f'-t/--threads must be a positive integer. Input was "{repr(args.threads)}"')

	return args, seq_type, ref_h2index, ref_fasta, ref_5p_fasta, ref_exons, ref_introns, ref_repeatmasker, output_base


def check_up_to_date(pipeline_dir, log):
	# Get repo status, ignoring file mode changes with -c core.fileMode=false
	fetch_command = f'git --git-dir {pipeline_dir}/.git --work-tree {pipeline_dir} -c core.fileMode=false fetch --all'
	functions.run_command(fetch_command, log=log, timeout=60)
	check_command = f'git --git-dir {pipeline_dir}/.git --work-tree {pipeline_dir} -c core.fileMode=false status'
	status = functions.run_command(check_command, log=log, timeout=60)
	log.debug(f'Git status: {status}')
	status = [line for line in status.split('\n') if line!='']

	branch = status[0].split(' ')[-1]
	if branch != 'main':
		log.warning(f'Git: LariatMapper is currently on branch {branch}, NOT the main branch!')
		time.sleep(60)	# Give the user a chance to see the warning

	if not status[1].startswith('Your branch is up-to-date with '):
		log.warning('Git: LariatMapper is not up-to-date!\n\t' + '\n\t'.join(status))
		time.sleep(60)	# Give the user a chance to see the warning
	elif not any(line.startswith('nothing to commit') for line in status):
		log.warning("Git: Working directory is not clean!\n\t" + '\n\t'.join(status))
		time.sleep(60)	# Give the user a chance to see the warning
	
	

# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	# Get path to the pipeline's directory 
	pipeline_dir = os.path.dirname(os.path.realpath(__file__))

	# Get the current version
	with open(f'{pipeline_dir}/pyproject.toml') as toml:
		for line in toml:
			if line.startswith('version = "'):
				version = line[11:].rstrip('"\n')

	# Argument parser
	parser = argparse.ArgumentParser(prog='Lariat mapping', description='Performs annotation-based mapping of lariat-derived RNA-seq reads')
	parser.add_argument('-v', '--version', action='version', version=f'LariatMapper {version}', help='Print the version id and exit')

	# Required arguments
	read_group = parser.add_mutually_exclusive_group(required=True)
	read_group.add_argument('-f', '--read_file', help='Input FASTQ file when processing single-end RNA-seq data. Can be uncompressed or gzip-compressed. Requires either this argument or both -1 and -2')
	read_group.add_argument('-1', '--read_one', help='Read 1 input FASTQ file when processing paired-end RNA-seq data. Can be uncompressed or gzip-compressed. ')
	parser.add_argument('-2', '--read_two', help='Read 2 input FASTQ file when processing paired-end RNA-seq data. Can be uncompressed or gzip-compressed. ')
	parser.add_argument('-o', '--output_dir', required=True, help='Directory for output files (will be created if it does not exist)')
	reference_group = parser.add_mutually_exclusive_group(required=True)
	reference_group.add_argument('-r', '--ref_dir', help='Directory with reference files for lariat mapping. Create by running build_references.py. Requires either this argument or all of -i, -g, -a, -5, -n, and -m')
	reference_group.add_argument('-i', '--ref_h2index', help='hisat2 index of the reference genome')
	parser.add_argument('-g', '--ref_fasta', help='FASTA file of the reference genome')
	parser.add_argument('-5', '--ref_5p_fasta', help='FASTA file with sequences of first 20nt of annotated introns')
	parser.add_argument('-e', '--ref_exons', help='TSV file of all annotated exons')
	parser.add_argument('-n', '--ref_introns', help='TSV file of all annotated introns')
	# Optional arguments
	parser.add_argument('-x', '--ignore_version', action='store_true', help='Don\'t check if LariatMapper is up-to-date with the main branch on GitHub (default=check and warn if not up-to-date)')
	parser.add_argument('-m', '--ref_repeatmasker', help='BED file of repetitive regions annotated by RepeatMasker. Putative lariats that map to a repetitive region will be filtered out as false positives (Optional)')
	parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use for parallel processing (default=1)')
	parser.add_argument('-p', '--output_prefix', help='Add a prefix to output file names (-o OUT -p ABC   ->   OUT/ABC_lariat_reads.tsv)')
	parser.add_argument('-u', '--ucsc_track', action='store_true', help='Add an output file named "lariat_reads.bed" which can be used as a custom track in the UCSC Genome Browser (https://www.genome.ucsc.edu/cgi-bin/hgCustom) to visualize lariat alignments')
	parser.add_argument('-k', '--keep_temp', action='store_true', help='Don\'t delete the temporary files created while running the pipeline (default=delete)')
	log_levels = parser.add_mutually_exclusive_group()
	log_levels.add_argument('-q', '--quiet', action='store_true', help="Only print fatal error messages (sets logging level to ERROR)")
	log_levels.add_argument('-w', '--warning', action='store_true', help="Print warning messages and fatal error messages (sets logging level to WARNING)")
	log_levels.add_argument('-d', '--debug', action='store_true', help="Print extensive status messages (sets logging level to DEBUG)")

	# Parse arguments
	args = parser.parse_args()
	
	# Setup logging
	if args.quiet is True:
		log_level = 'ERROR'
	if args.warning is True:
		log_level = 'WARNING'
	elif args.debug is True:
		log_level = 'DEBUG'
	else:
		log_level = 'INFO'
	log = functions.get_logger(log_level)

	# Report version
	log.info(f'LariatMapper {version}')

	# Check if up-to-date
	if args.ignore_version is False:
		log.debug('Checking if LariatMapper is up-to-date with the main branch on GitHub...')
		try:
			check_up_to_date(pipeline_dir, log)
		except Exception as e:
			log.debug(e)
			log.warning('Could not check if LariatMapper is up-to-date with the main branch on GitHub. Continuing anyway...')

	# Report arguments
	arg_message = [f'{key}={val}' for key, val in vars(args).items() if val is not None and val is not False]
	arg_message = '\n\t'.join(arg_message)
	log.info(f'Arguments: \n\t{arg_message}')

	# Validate the args and determine additional variables
	args, seq_type, ref_h2index, ref_fasta, ref_5p_fasta, ref_exons, ref_introns, ref_repeatmasker, output_base = process_args(args, parser, log)
	log.debug(f'seq_type: {repr(seq_type)}, output_base: {repr(output_base)}')

	# Make output dir
	log.info('Preparing directories...')
	if not os.path.isdir(args.output_dir):
		os.mkdir(args.output_dir)
	# Move to output dir
	os.chdir(args.output_dir)

	# Set start method for multiprocessing in filter_fivep_aligns.py and filter_head_aligns.py
	# Using the default "fork" method causes memory errors when processing bigger RNA-seq inputs
	# This has to be defined in the first python script and only once, or else we get "RuntimeError: context has already been set" 
	multiprocessing.set_start_method('spawn')	

	# Prepare call to map_lariats.sh
	map_lariats_args = [output_base, str(args.threads), ref_h2index, ref_fasta, ref_5p_fasta, ref_exons, ref_introns, ref_repeatmasker, str(args.keep_temp).lower(), str(args.ucsc_track).lower(), pipeline_dir, log_level, seq_type]
	if seq_type == 'single':
		log.debug('Processing single-end read file...')
		map_lariats_args += [args.read_file]
	elif seq_type == 'paired':
		log.debug('Processing paired-end read files...')
		map_lariats_args += [args.read_one, args.read_two]
	log.debug(f'map_lariats args: {map_lariats_args}')

	# Dump arguments to a JSON file so summarize.py can report them
	with open(f'{args.output_dir}/args.json', 'w') as json_file:
		json.dump(map_lariats_args, json_file)

	# Run it
	map_call = f'{os.path.join(pipeline_dir, "scripts", "map_lariats.sh")} {" ".join(map_lariats_args)}'
	subprocess.run(map_call.split(' '))
