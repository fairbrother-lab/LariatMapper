#!/usr/bin/env python3

import argparse
import os
import sys
import time
from subprocess import run
import multiprocessing
import logging



from scripts import functions


# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def process_args(args:argparse.ArgumentParser, log):
	# Determine whether input is single-end or paired-end and confirm that the read file(s) exit
	if args.read_file is not None:
		seq_type = 'single'
		if not os.path.isfile(args.read_file):
			raise ValueError(f'"{args.read_file}" is not an existing file')
	elif args.read_one is not None and args.read_two is not None:
		seq_type = 'paired'
		if not os.path.isfile(args.read_one):
			raise ValueError(f'"{args.read_one}" is not an existing file')
		if not os.path.isfile(args.read_two):
			raise ValueError(f'"{args.read_two}" is not an existing file')
	else:
		parser.error('Provide either -f/--read_file (for single-end read) OR -1/--read_one and -2/--read_two (for paired-end reads)')

	# Get reference files from reference dir OR from direct input, and confirm that they exist
	if args.ref_dir is not None:
		ref_files = ['hisat2_index', 'genome.fa', 'fivep_sites.fa', 'introns.tsv.gz', 'repeatmasker.bed']
		ref_h2index, ref_fasta, ref_5p_fasta, ref_introns, ref_repeatmasker = [os.path.join(args.ref_dir, f) for f in ref_files]
	else:
		ref_h2index, ref_fasta, ref_5p_fasta, ref_introns, ref_repeatmasker = args.ref_h2index, args.ref_fasta, args.ref_5p_fasta, args.ref_introns, args.ref_repeatmasker

	for file in (ref_fasta, ref_5p_fasta, ref_introns):
		if not os.path.isfile(file):
			parser.error(f'"{file}" is not an existing file')
	if (not os.path.isfile(ref_h2index + '.1.ht2')) and (not os.path.isfile(ref_h2index + '.1.ht2l')):
		parser.error(f'"{ref_h2index}" is not an existing hisat2 index')
	
	# Determine the output_base
	# All output files will be formatted like f"{output_base}file.ext"
	if args.output_prefix is None:
		output_prefix = ''
	elif not args.output_prefx.isidentifier():
		parser.error(f'-p/--output_prefix can only contain alphanumerics and underscores. Input was "{repr(output_prefix)}"')
	else:
		output_prefix = args.output_prefix + '_'
	output_base = os.path.join(args.output_dir, output_prefix)

	# Validate threads arg
	if not args.threads>0:
		parser.error(f'-t/--threads must be a positive integer. Input was "{repr(args.threads)}"')

	return args, seq_type, ref_h2index, ref_fasta, ref_5p_fasta, ref_introns, ref_repeatmasker, output_base
		


# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='Lariat mapping', description='Performs annotation-based mapping of lariat-derived RNA-seq reads')
	
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
	parser.add_argument('-n', '--ref_introns', help='BED file of all annotated introns')
	# Optional arguments
	parser.add_argument('-m', '--ref_repeatmasker', help='BED file of repetitive regions annotated by RepeatMasker. Putative lariats that map to a repetitive region will be filtered out as false positives (Optional)')
	parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use for parallel processing (default=1)')
	parser.add_argument('-p', '--output_prefix', help='Add a prefix to output file names (-o OUT -p ABC   ->   OUT/ABC_lariat_reads.tsv)')
	parser.add_argument('-u', '--ucsc_track', action='store_true', help='Add an output file named "lariat_reads.bed" which can be used as a custom track in the UCSC Genome Browser (https://www.genome.ucsc.edu/cgi-bin/hgCustom) to visualize lariat alignments')
	parser.add_argument('-k', '--keep_intermediates', action='store_true', help='Don\'t delete the intermediate files created while running the pipeline (default=delete)')
	log_levels = parser.add_mutually_exclusive_group()
	log_levels.add_argument('-q', '--quiet', action='store_true', help="Don't print any status messages")
	log_levels.add_argument('-d', '--debug', action='store_true', help="Print extensive status messages")

	# Parse
	args = parser.parse_args()

	# Get path to the pipeline's directory so we can call scripts in the "scripts" folder
	pipeline_dir = os.path.dirname(os.path.realpath(__file__))

	# Setup logging
	if args.quiet is True:
		log_level = 'ERROR'
	elif args.debug is True:
		log_level = 'DEBUG'
	else:
		log_level = 'INFO'
	log = functions.get_logger(log_level)

	# Print arguments
	arg_message = [f'{key}={val}' for key, val in vars(args).items() if val is not None and val is not False]
	arg_message = '\n\t'.join(arg_message)
	log.info(f'Arguments: \n\t{arg_message}')

	# Validate the args and determine additional variables
	args, seq_type, ref_h2index, ref_fasta, ref_5p_fasta, ref_introns, ref_repeatmasker, output_base = process_args(args, log)
	log.debug(f'seq_type: {repr(seq_type)}, output_base: {repr(output_base)}')

	# Set start method for multiprocessing in filter_fivep_alignments.py and filter_trimmed_alignments.py
	# Using the default "fork" method causes memory errors when processing bigger RNA-seq inputs
	# This has to be defined in the first python script and only once, or else we get "RuntimeError: context has already been set" 
	multiprocessing.set_start_method('spawn')	

	# Make output dir
	# print(time.strftime('%m/%d/%y - %H:%M:%S | Preparing directories...'), flush=True)
	log.info('Preparing directories...')
	if not os.path.isdir(args.output_dir):
		os.mkdir(args.output_dir)

	# Prepare call to map_lariats.sh
	map_lariats_args = [output_base, str(args.threads), ref_h2index, ref_fasta, ref_5p_fasta, ref_introns, ref_repeatmasker, str(args.keep_intermediates).lower(), str(args.ucsc_track).lower(), pipeline_dir, log_level]
	if seq_type == 'single':
		# print(time.strftime('%m/%d/%y - %H:%M:%S | Processing single-end read file...'), flush=True)
		log.debug('Processing single-end read file...')
		map_lariats_args += [args.read_file]
	else:
		# print(time.strftime('%m/%d/%y - %H:%M:%S | Processing paired-end read files...'), flush=True)
		log.debug('Processing paired-end read files...')
		map_lariats_args += [args.read_one, args.read_two]
	log.debug(f'map_lariats args: {"\n\t".join(map_lariats_args)}')

	# Run it
	map_call = f'{os.path.join(pipeline_dir, "scripts", "map_lariats.sh")} {" ".join(map_lariats_args)}'
	run(map_call.split(' '))
