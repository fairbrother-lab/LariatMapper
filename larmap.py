#!/usr/bin/env python3

import argparse
import json
import os
import multiprocessing
import subprocess
import dataclasses
import pathlib

# This line is where our third-party imports would go, if we had any

from scripts import utils





# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
FORBIDDEN_CHARS = (" ", '\\', '>', '<', '&', '~', '*', '?', '[', ']', ',', ';', '|', '!', '$', "'", '"')

@dataclasses.dataclass
class Settings:
	REQ_REFS = (('ref_h2index', 'hisat2_index'), 
				('ref_fasta', 'genome.fa'), 
				('ref_5p_fasta', 'fivep_sites.fa'),
				('ref_exons', 'exons.tsv.gz'),
				('ref_introns', 'introns.tsv.gz'))
	PATH_SETTINGS = ('read_file', 'read_one', 'read_two', 'ref_dir', 'ref_h2index', 'ref_fasta', 
					'ref_5p_fasta', 'ref_exons', 'ref_introns', 'output_dir', 'ref_repeatmasker', 
					'model_correction', 'pipeline_dir')
	ARGS_TO_MAP_LARIATS = ('input_reads', 'ref_dir', 'ref_h2index', 'ref_fasta', 'ref_5p_fasta', 
						'ref_exons', 'ref_introns', 'strand', 'temp_switch_filter', 
						'ref_repeatmasker', 'pwm_correction', 'model_correction', 'ucsc_track', 
						'keep_bam', 'keep_classes', 'keep_temp', 'threads', 'seq_type', 
						'output_base', 'log_level', 'pipeline_dir')

	# Supplied argument attributes
	# If an argument is not supplied, it will be None
	read_file: pathlib.Path
	read_one: pathlib.Path
	read_two: pathlib.Path
	ref_dir: pathlib.Path
	ref_h2index: pathlib.Path		# If ref_dir is supplied and this is not, will be set to {ref_dir}/hisat2_index
	ref_fasta: pathlib.Path			# If ref_dir is supplied and this is not, will be set to {ref_dir}/genome.fa
	ref_5p_fasta: pathlib.Path		# If ref_dir is supplied and this is not, will be set to {ref_dir}/fivep_sites.fa
	ref_exons: pathlib.Path			# If ref_dir is supplied and this is not, will be set to {ref_dir}/exons.tsv.gz
	ref_introns: pathlib.Path		# If ref_dir is supplied and this is not, will be set to {ref_dir}/introns.tsv.gz
	output_dir: pathlib.Path
	strand: str
	temp_switch_filter: str
	ref_repeatmasker: pathlib.Path	# If ref_dir is supplied and this is not, will be set to {ref_dir}/repeatmasker.bed
	pwm_correction: str
	model_correction: str
	output_prefix: str
	ucsc_track: bool
	keep_bam: bool
	keep_classes: bool
	keep_temp: bool
	threads: str
	quiet: bool
	warning: bool
	debug: bool
	# Extrapolated attributes
	input_reads: str = dataclasses.field(init=False)
	seq_type: str = dataclasses.field(init=False)
	output_base: pathlib.Path = dataclasses.field(init=False)
	log_level: str = dataclasses.field(init=False)
	# Automatically determined attributes
	pipeline_dir = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))


	def __post_init__(self):
		# Reference files
		for ref_attr_name, ref_file_name in Settings.REQ_REFS:
			if getattr(self, ref_attr_name) is None:
				default_path = pathlib.Path(os.path.join(self.ref_dir, ref_file_name))

				if ref_attr_name=='ref_fasta' and not default_path.is_file():
					default_path = pathlib.Path(os.path.join(self.ref_dir, ref_file_name + '.gz'))
				elif ref_attr_name=='ref_5p_fasta' and not default_path.is_file():
					default_path = pathlib.Path(os.path.join(self.ref_dir, ref_file_name + '.gz'))

				setattr(self, ref_attr_name, default_path)
		
		# If ref_repeatmasker wasn't input, set it to {ref_dir}/repeatmasker.bed
		# If it doesn't exist that's fine, filter_lariats.py will check before trying to use it
		if self.ref_repeatmasker is None:
			setattr(self, 'ref_repeatmasker', pathlib.Path(os.path.join(self.ref_dir, 'repeatmasker.bed')))

		# Make sure all paths are absolute
		for attr in Settings.PATH_SETTINGS:
			if getattr(self, attr) is not None:
				setattr(self, attr, pathlib.Path(getattr(self, attr)).resolve())

		# input_reads and seq_type
		if self.read_file is not None and self.read_one is None and self.read_two is None:
			self.input_reads = self.read_file
			self.seq_type = 'single'
		elif self.read_file is None and self.read_one is not None and self.read_two is not None:
			self.input_reads = f'{self.read_one},{self.read_two}'
			self.seq_type = 'paired'
		else:
			raise ValueError('Provide either -f/--read_file (for single-end read) OR -1/--read_one and -2/--read_two (for paired-end reads)')
		
		# output_base
		# All output files will be formatted like f"{output_base}file.ext"
		if self.output_prefix is None:
			self.output_base = f'{self.output_dir}/'
		else:
			for char in FORBIDDEN_CHARS:
				if char in str(self.output_prefix):
					parser.error(f'Illegal character in output prefix: {char}')
			self.output_base = f'{self.output_dir/self.output_prefix}_'

		# pwm_correction and model_correction
		if self.pwm_correction is None:
			self.pwm_correction = ''
		if self.model_correction is None:
			self.model_correction = ''

		# log_level
		if self.quiet is True:
			self.log_level = 'ERROR'
		elif self.warning is True:
			self.log_level = 'WARNING'
		elif self.debug is True:
			self.log_level = 'DEBUG'
		else:
			self.log_level = 'INFO'


	def validate_args(self):
		# Confirm that the read file(s) exit
		if self.seq_type=='single' and self.read_file.is_file() is False:
			raise ValueError(f'"{self.read_file}" is not an existing file')
		elif self.seq_type=='paired' and self.read_one.is_file() is False:
			raise ValueError(f'"{self.read_one}" is not an existing file')
		elif self.seq_type=='paired' and self.read_two.is_file() is False:
			raise ValueError(f'"{self.read_two}" is not an existing file')

		# Confirm that the reference files exist and don't have forbidden characters
		for file in (self.ref_fasta, self.ref_5p_fasta, self.ref_introns):
			if file.is_file() is False:
				raise ValueError(f'"{file}" is not an existing file')
			for char in FORBIDDEN_CHARS:
				if char in str(file):
					raise ValueError(f'Illegal character in {file}: "{char}"')
		if (not os.path.isfile(f'{self.ref_h2index}.1.ht2')) and (not os.path.isfile(f'{self.ref_h2index}.1.ht2l')):
			raise ValueError(f'"{self.ref_h2index}" is not an existing hisat2 index')
		
		# Confirm that the output directory parent exists and it doesn't have forbidden characters
		if self.output_dir.parent.is_dir() is False:
			raise ValueError(f'"{self.output_dir.parent}" is not an existing directory')
		for char in FORBIDDEN_CHARS:
			if char in str(self.output_dir):
				raise ValueError(f'Illegal character in output directory: "{char}"')
			
		# Confirm the temp_switch_filter arg is formatted correctly, if supplied
		if self.temp_switch_filter != '':
			if self.temp_switch_filter.count(',') != 1 or self.temp_switch_filter.startswith(',') or self.temp_switch_filter.endswith(','):
				raise ValueError(f'--temp_switch_filter must be formatted as "N,M". Input was "{self.temp_switch_filter}"')
			n, m = self.temp_switch_filter.split(',')
			if not n.isdigit() or not m.isdigit():
				raise ValueError(f'--temp_switch_filter must be formatted as "N,M". Input was "{self.temp_switch_filter}"')
			if int(n) < int(m):
				raise ValueError(f'--temp_switch_filter: N must be greater than or equal to M. Input was "{self.temp_switch_filter}"')
			
		# Confirm the pwm correction files exist OR model correction file exists, if supplied
		if self.pwm_correction != '':
			for file in self.pwm_correction.split(','):
				if pathlib.Path(file).is_file() is False:
					raise ValueError(f'"{file}" is not an existing file')
		if self.model_correction != '':
			if pathlib.Path(self.model_correction).is_file() is False:
				raise ValueError(f'"{self.model_correction}" is not an existing file')

		# Validate threads arg
		if not self.threads>0:
			raise ValueError(f'-t/--threads must be a positive integer. Input was "{repr(self.threads)}"')


	@property
	def map_lariats_args(self):
		args = []
		for attr in Settings.ARGS_TO_MAP_LARIATS:
			arg_val = str(getattr(self, attr))
			if arg_val in ('True', 'False'):
				arg_val = arg_val.lower()
			args.append(arg_val)
			
		return args
	

	def json_dump(self, file):
		dict_ = {key: str(val) for key,val in dataclasses.asdict(self).items()}
		with open(file, 'w') as json_file:
			json.dump(dict_, json_file)

	
	


# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	# Argument parser
	parser = argparse.ArgumentParser(prog='larmap.py', description='Extracts lariats and their branchpoint positions from RNA-seq data')
	parser.add_argument('-v', '--version', action='version', version=f'LariatMapper {utils.version()}', help='print the version id and exit')

	# Required arguments
	# We use argument groups to make the help message more readable, but we have to enforce 
	# more complicated mutual exclusion relationships (like with -1 and -2 vs -f) manually later
	# because argparse only supports basic mutual exclusivity 
	read_group = parser.add_argument_group(title='Input read files', 
										description='Provide either two FASTQ files from paired-end RNA-seq or one FASTQ file from single-end RNA-seq. Files can be gzip-compressed')
	read_group.add_argument('-1', '--read_one', type=pathlib.Path, help='FASTQ file of read mate 1. Use for paired-end RNA-seq data. Mutually exclusive with -f, requires -2')
	read_group.add_argument('-2', '--read_two', type=pathlib.Path, help='FASTQ file of read mate 2. Use for paired-end RNA-seq data. Mutually exclusive with -f, requires -1')
	read_group.add_argument('-f', '--read_file', type=pathlib.Path, help='FASTQ file. Use for single-end RNA-seq data. Mutually exclusive with -1 and -2')
	reference_group = parser.add_argument_group(title='Reference data')
	reference_group.add_argument('-r', '--ref_dir', required=True, type=pathlib.Path, help='Directory with reference files created by build_references.py')
	out_group = parser.add_argument_group(title='Output')
	out_group.add_argument('-o', '--output_dir', required=True, type=pathlib.Path, help='Directory for output files. Will be created if it does not exist at runtime')

	# Optional arguments
	optional_args = parser.add_argument_group(title='Optional arguments')
		# Experimentally-revelant options 
	#TODO: Test --strand arg with stranded data and verify that it works correctly
	# Strandedness of the input reads. Choices: Unstranded = Library preparation wasn't strand-specific; 
	# First = READ_ONE/READ_FILE reads match the RNA sequence (i.e. 2nd cDNA synthesis strand); 
	# Second = READ_ONE/READ_FILE reads are reverse-complementary to the RNA sequence (i.e. 1st cDNA synthesis strand) 
	# (Default = Unstranded)")
	optional_args.add_argument('-s', '--strand', choices=('Unstranded', 'First', 'Second'), default='Unstranded', help=argparse.SUPPRESS)
	optional_args.add_argument('-T', '--temp_switch_filter', default='5,5', help='Set the parameters of the template-switching filter in the head-filtering step. Format = "N,M", where N is the number of downstream bases to check, and M is the minimum number of matches required to identify an alignment as template-switching. (Default = 5,5)')
	optional_args.add_argument('-m', '--ref_repeatmasker', type=pathlib.Path, help="BED file of repetitive regions in the genome. Putative lariats that map to a repetitive region will be filtered out as false positives. May be gzip-compressed. (Default = REF_DIR/repeatmasker.bed if it's an existing file, otherwise skip repetitive region filtering)")
	optional_args.add_argument('-H', '--ref_h2index', type=pathlib.Path, help='HISAT2 index of the reference genome. (Default = REF_DIR/hisat2_index)')
	optional_args.add_argument('-g', '--ref_fasta', type=pathlib.Path, help='FASTA file of the reference genome. May be gzip-compressed. (Default = REF_DIR/genome.fa.gz)')
	optional_args.add_argument('-5', '--ref_5p_fasta', type=pathlib.Path, help="FASTA file of 5' splice site sequences, i.e. the first 20nt of all annotated introns. (Default = REF_DIR/fivep_sites.fa)")
	optional_args.add_argument('-e', '--ref_exons', type=pathlib.Path, help='TSV file of all annotated introns. (Default = REF_DIR/exons.tsv.gz)')
	optional_args.add_argument('-i', '--ref_introns', type=pathlib.Path, help='TSV file of all annotated introns. (Default = REF_DIR/introns.tsv.gz)')
	bp_correction = optional_args.add_mutually_exclusive_group()
	bp_correction.add_argument('-P', '--pwm_correction', help='RDS file with a position weight matrix to correct apparent branchpoint positions. Multiple files can be provided in comma-seperated format. Mutually exclusive with --model_correction. See https://doi.org/10.5281/zenodo.14735947 to download prebuilt PWMs. See scripts/pwm_build.R to build a custom matrix (Default = no correction)')
	bp_correction.add_argument('-M', '--model_correction', help='RDS file with predictions from DeepEnsemble, a deep-learning-based branchpoint prediction model. Mutually exclusive with --pwm_correction. See https://doi.org/10.5281/zenodo.14735947 to download predictions for specific reference genomes. (Default = no correction)')
		# Output options
	optional_args.add_argument('-p', '--output_prefix', help='Add a prefix to output file names (-o OUT -p ABC -> OUT/ABC_lariat_reads.tsv). (Default = no prefix)')
	optional_args.add_argument('-u', '--ucsc_track', action='store_true', help='Add an output file named "lariat_reads.bed". This can be used as a custom track in the UCSC Genome Browser to visualize lariat read alignments')
	optional_args.add_argument('-b', '--keep_bam', action='store_true', help='Keep the BAM file produced in the initial linear mapping step (Default = delete)')
	optional_args.add_argument('-c', '--keep_classes', action='store_true', help='Keep a file with per-read classification of non-linearly-aligned reads named "read_classes.tsv.gz" in the output (Default = delete)')
	optional_args.add_argument('-k', '--keep_temp', action='store_true', help='Keep all temporary files created while running the pipeline. Forces -c and -b (Default = delete)')
		# Technical options
	optional_args.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use. (Default = 1)')
	log_levels = optional_args.add_mutually_exclusive_group()
	log_levels.add_argument('-q', '--quiet', action='store_true', help="Only print fatal error messages. Mutually exclusive with -w and -d")
	log_levels.add_argument('-w', '--warning', action='store_true', help="Print warning messages and fatal error messages. Mutually exclusive with -q and -d")
	log_levels.add_argument('-d', '--debug', action='store_true', help="Print extensive status messages. Mutually exclusive with -q and -w")

	# Parse args into a dict
	args = vars(parser.parse_args())
	# Convert to a Settings object
	settings = Settings(**args)

	# Set up logging
	log = utils.get_logger(settings.log_level)

	# Report version
	log.info(f'LariatMapper {utils.version()}')

	# Report arguments
	arg_message = [f'{key}={val}' for key, val in args.items() if val is not None and val is not False]
	arg_message = '\n\t'.join(arg_message)
	log.info(f'Arguments: \n\t{arg_message}')

	# Validate arguments
	settings.validate_args()

	# Make output dir
	log.info('Preparing directories...')
	if not os.path.isdir(settings.output_dir):
		os.mkdir(settings.output_dir)
	# Move to output dir
	os.chdir(settings.output_dir)

	# Dump arguments to a JSON file so summarize.py can report them
	settings.json_dump(f'{settings.output_base}settings.json')

	# Set start method for multiprocessing in filter_fivep_aligns.py and filter_head_aligns.py
	# Using the default "fork" method causes memory errors when processing bigger RNA-seq inputs
	# This has to be defined in the first python script and only once, or else we get "RuntimeError: context has already been set" 
	multiprocessing.set_start_method('spawn')	

	# Call map_lariats.sh
	log.debug(f'map_lariats args: {settings.map_lariats_args}')
	map_call = f"{settings.pipeline_dir/'scripts'/'map_lariats.sh'} {' '.join(settings.map_lariats_args)}"
	subprocess.run(map_call.split(' '))
