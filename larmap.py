#!/usr/bin/env python3

import argparse
import os
import time
from subprocess import run



# =============================================================================#
#                                    Main                                      #
# =============================================================================#
def main():

	parser = argparse.ArgumentParser(prog='Lariat mapping', description='Performs annotation-based mapping of lariat-derived RNA-seq reads')
	
	# Required arguments
	read_group = parser.add_mutually_exclusive_group(required=True)
	read_group.add_argument('-f', '--read_file', help='Input FASTQ file when processing single-end RNA-seq data. Requires either this argument or both -1 and -2')
	read_group.add_argument('-1', '--read_one', help='Read 1 input FASTQ file when processing paired-end RNA-seq data')
	parser.add_argument('-2', '--read_two', help='Read 2 input FASTQ file when processing paired-end RNA-seq data')
	parser.add_argument('-o', '--output_dir', help='Directory for output files (will be created if it does not exist)', required=True)
	reference_group = parser.add_mutually_exclusive_group(required=True)
	reference_group.add_argument('-r', '--ref_dir', help='Directory with reference files for lariat mapping. Create by running build_references.py. Requires either this argument or all of -i, -g, -a, -5, -n, and -m')
	reference_group.add_argument('-i', '--ref_h2index', help='hisat2 index of the reference genome')
	parser.add_argument('-g', '--ref_fasta', help='FASTA file of the reference genome')
	parser.add_argument('-a', '--ref_anno', help='Gene annotation of the reference genome in GTF or GFF format (may be gzipped with .gz extension)')
	parser.add_argument('-5', '--ref_5p_fasta', help='FASTA file with sequences of first 20nt of annotated introns')
	parser.add_argument('-n', '--ref_introns', help='BED file of all annotated introns')
	parser.add_argument('-m', '--ref_repeatmasker', help='BED file of repetitive element annotation from RepeatMasker')
	# Optional arguments
	parser.add_argument('-t', '--threads', help='Number of threads to use for parallel processing (default=1)', type=int, default=1)
	parser.add_argument('-p', '--output_prefix', help='Prefix to add to output file names')
	parser.add_argument('-k', '--keep_intermediates', action='store_true', help='Don\'t delete the intermediate files created while running the pipeline')

	args = parser.parse_args()

	# Print arguments recieved
	print('\nArguments received:', flush=True)
	# if seq_type == 'single':
	# 	arg_names, arg_values = ['read_file'], [read_file]
	# else:
	# 	arg_names, arg_values = ['read_one', 'read_two'], [read_one, read_two]
	# arg_names += ['output_dir', 'output_prefix', 'threads', 'ref_h2index', 'ref_fasta', 'ref_anno', 'ref_5p_fasta', 'ref_introns', 'ref_repeatmasker', 'keep_intermediates']
	# arg_values += [output_dir, output_prefix, threads, ref_h2index, ref_fasta, ref_anno, ref_5p_fasta, ref_introns, ref_repeatmasker, keep_intermediates]
	# for n, v in zip(arg_names, arg_values):
	# 	print(f'{n}: {v}', flush=True)
	arg_message = [f'{key}={val}' for key, val in vars(args).items() if val is not None]
	arg_message = '\n'.join(arg_message) + '\n'
	print(arg_message)

	# Determine whether input is single-end or paired-end 
	if args.read_file is not None:
		read_file, seq_type = args.read_file, 'single'
	elif args.read_one is not None and args.read_two is not None:
		read_one, read_two, seq_type = args.read_one, args.read_two, 'paired'
	else:
		parser.error('Provide either -f/--read_file (for single-end read) OR -1/--read_one and -2/--read_two (for paired-end reads)')
	
	output_dir, output_prefix, threads, keep_intermediates = args.output_dir, args.output_prefix, args.threads, args.keep_intermediates

	# Get references from reference dir OR from direct input
	if args.ref_dir is not None:
		ref_files = ['hisat2_index', 'genome.fa', 'fivep_sites.fa', 'introns.bed', 'repeatmasker.bed']
		ref_h2index, ref_fasta, ref_5p_fasta, ref_introns, ref_repeatmasker = [os.path.join(args.ref_dir, f) for f in ref_files]
		ref_anno = os.path.join(args.ref_dir, [f for f in os.listdir(args.ref_dir) if f.split('.')[0]=='annotation'][0])
	else:
		ref_h2index, ref_fasta, ref_anno, ref_5p_fasta, ref_introns, ref_repeatmasker = args.ref_h2index, args.ref_fasta, args.ref_anno, args.ref_5p_fasta, args.ref_introns, args.ref_repeatmasker

	print(time.strftime('\n%m/%d/%y - %H:%M:%S | Starting lariat mapping run...'), flush=True)

	# Make output dir
	print(time.strftime('%m/%d/%y - %H:%M:%S | Preparing directories...'), flush=True)
	pipeline_dir = os.path.dirname(os.path.realpath(__file__))
	if not os.path.isdir(output_dir):
		os.mkdir(output_dir)
	if output_prefix is None:
		output_prefix = ''
	else:
		output_prefix += '_'
	output_base = os.path.join(output_dir, output_prefix)
	
	# Run map_lariats.sh
	map_lariats_args = [output_base, str(threads), ref_h2index, ref_fasta, ref_anno, ref_5p_fasta, ref_introns, ref_repeatmasker, str(keep_intermediates), pipeline_dir]
	if seq_type == 'single':
		print(time.strftime('%m/%d/%y - %H:%M:%S | Processing single-end read file...'), flush=True)
		map_lariats_args += [read_file]
	else:
		print(time.strftime('%m/%d/%y - %H:%M:%S | Processing paired-end read files...'), flush=True)
		map_lariats_args += [read_one, read_two]
	map_call = f'{os.path.join(pipeline_dir, "scripts", "map_lariats.sh")} {" ".join(map_lariats_args)}'
	run(map_call.split(' '))


if __name__ == '__main__':

	main()