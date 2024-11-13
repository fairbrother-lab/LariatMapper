#!/usr/bin/env python3

import argparse
import json
import os
import time
import multiprocessing
import subprocess
import dataclasses
import pathlib

# This line is where our third-party imports go, if we had any

from scripts import functions



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
FORBIDDEN_CHARS = ('\\', '>', '<', '&', '~', '*', '?', '[', ']', ',', ';', '|', '!', '$', "'", '"')

@dataclasses.dataclass
class Settings:
	REQ_REFS = (('ref_h2index', 'hisat2_index'), 
				('ref_fasta', 'genome.fa'), 
				('ref_5p_fasta', 'fivep_sites.fa'),
				('ref_introns', 'introns.tsv.gz'))
	PATH_SETTINGS = ('read_file', 'read_one', 'read_two', 'ref_dir', 'ref_h2index', 'ref_fasta', 
					'ref_5p_fasta', 'ref_introns', 'output_dir', 'ref_repeatmasker',
					'pipeline_dir')
	ARGS_TO_MAP_LARIATS = ('ref_dir', 'ref_h2index', 'ref_fasta', 'ref_5p_fasta', 'ref_introns', 
						'strand', 'ref_repeatmasker', 'ucsc_track', 'keep_classes', 'keep_temp', 
						'threads', 'input_reads', 'seq_type', 'output_base', 'log_level', 'pipeline_dir')

	# Supplied argument attributes
	# If an argument is not supplied, it will be None
	read_file: pathlib.Path
	read_one: pathlib.Path
	read_two: pathlib.Path
	ref_dir: pathlib.Path
	ref_h2index: pathlib.Path		# If ref_dir is supplied and this is not, will be set to {ref_dir}/hisat2_index
	ref_fasta: pathlib.Path			# If ref_dir is supplied and this is not, will be set to {ref_dir}/genome.fa
	ref_5p_fasta: pathlib.Path		# If ref_dir is supplied and this is not, will be set to {ref_dir}/fivep_sites.fa
	ref_introns: pathlib.Path		# If ref_dir is supplied and this is not, will be set to {ref_dir}/introns.tsv.gz
	output_dir: pathlib.Path
	strand: str
	ref_repeatmasker: pathlib.Path	# If ref_dir is supplied and this is not, will be set to {ref_dir}/repeatmasker.bed
	output_prefix: str
	ucsc_track: bool
	keep_classes: bool
	keep_temp: bool
	threads: str
	skip_version_check: bool
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
		# # Reference files
		for ref_attr_name, ref_file_name in Settings.REQ_REFS:
			if getattr(self, ref_attr_name) is None:
				setattr(self, ref_attr_name, pathlib.Path(os.path.join(self.ref_dir, ref_file_name)))
		
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
		
		# log_level
		if self.quiet is True:
			self.log_level = 'ERROR'
		elif self.warning is True:
			self.log_level = 'WARNING'
		elif self.debug is True:
			self.log_level = 'DEBUG'
		else:
			self.log_level = 'INFO'
		
		# output_base
		# All output files will be formatted like f"{output_base}file.ext"
		if self.output_prefix is None:
			self.output_base = f'{self.output_dir}/'
		else:
			for char in FORBIDDEN_CHARS:
				if char in str(self.output_prefix):
					parser.error(f'Illegal character in output prefix: {char}')
			self.output_base = f'{self.output_dir/self.output_prefix}_'


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
					raise ValueError(f'Illegal character in {file}: {char}')
		if (not os.path.isfile(f'{self.ref_h2index}.1.ht2')) and (not os.path.isfile(f'{self.ref_h2index}.1.ht2l')):
			raise ValueError(f'"{self.ref_h2index}" is not an existing hisat2 index')
		
		# Confirm that the output directory parent exists and it doesn't have forbidden characters
		if self.output_dir.parent.is_dir() is False:
			raise ValueError(f'"{self.output_dir.parent}" is not an existing directory')
		for char in FORBIDDEN_CHARS:
			if char in str(self.output_dir):
				raise ValueError(f'Illegal character in output directory: {char}')

		# Validate threads arg
		if not self.threads>0:
			raise ValueError(f'-t/--threads must be a positive integer. Input was "{repr(self.threads)}"')


	@property
	def map_lariats_args(self):
		args = []
		for attr in Settings.ARGS_TO_MAP_LARIATS:
			arg_val = str(getattr(self, attr))
			if attr in ('ucsc_track', 'keep_classes', 'keep_temp'):
				arg_val = arg_val.lower()
			args.append(arg_val)
			
		return args
	

	def json_dump(self, file):
		dict_ = {key: str(val) for key,val in dataclasses.asdict(self).items()}
		with open(file, 'w') as json_file:
			json.dump(dict_, json_file)



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
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
	# Argument parser
	parser = argparse.ArgumentParser(prog='Lariat mapping', description='Performs annotation-based mapping of lariat-derived RNA-seq reads')
	parser.add_argument('-v', '--version', action='version', version=f'LariatMapper {functions.version()}', help='Print the version id and exit')

	# Required arguments
	# We use argument groups to make the help message more readable, but we have to enforce 
	# mutually exclusive arguments later because argparse doesn't support mutually exclusive groups
	read_group = parser.add_argument_group(title='Input read files', 
										description='Provide either two paired-end read files or one single-end read file. Files can be uncompressed or gzip-compressed')
	read_group.add_argument('-1', '--read_one', type=pathlib.Path, help='Read 1 input FASTQ file when processing paired-end RNA-seq data. Mutually exclusive with -f')
	read_group.add_argument('-2', '--read_two', type=pathlib.Path, help='Read 2 input FASTQ file when processing paired-end RNA-seq data. Mutually exclusive with -f')
	read_group.add_argument('-f', '--read_file', type=pathlib.Path, help='Input FASTQ file when processing single-end RNA-seq data. Mutually exclusive with -1 and -2')
	reference_group = parser.add_argument_group(title='Reference data')
	reference_group.add_argument('-r', '--ref_dir', required=True, type=pathlib.Path, help='Directory with reference files created by build_references.py')
	out_group = parser.add_argument_group(title='Output')
	out_group.add_argument('-o', '--output_dir', required=True, type=pathlib.Path, help='Directory for output files. Will be created if it does not exist at runtime')
	# Optional arguments
	optional_args = parser.add_argument_group(title='Optional arguments')
		# Experimentally-revelant options 
	optional_args.add_argument('-s', '--strand', choices=('Unstranded', 'First', 'Second'), default='Unstranded', help="WARNING, EXPERIMENTAL FEATURE STILL IN DEVELOPMENT! Strandedness of the input reads. Choices: Unstranded = Library preparation wasn't strand-specific; First = READ_ONE/READ_FILE reads match the RNA sequence (i.e. 2nd cDNA synthesis strand); Second = READ_ONE/READ_FILE reads are reverse-complementary to the RNA sequence (i.e. 1st cDNA synthesis strand) (Default = Unstranded)")
	optional_args.add_argument('-m', '--ref_repeatmasker', type=pathlib.Path, help="BED file of repetitive regions in the genome. Putative lariats that map to a repetitive region will be filtered out as false positives (Default = REF_DIR/repeatmasker.bed if it's an existing file, otherwise skip repetitive region filtering")
	optional_args.add_argument('-i', '--ref_h2index', type=pathlib.Path, help='hisat2 index of the reference genome (Default = REF_DIR/hisat2_index)')
	optional_args.add_argument('-g', '--ref_fasta', type=pathlib.Path, help='FASTA file of the reference genome (Default = REF_DIR/genome.fa)')
	optional_args.add_argument('-5', '--ref_5p_fasta', type=pathlib.Path, help='FASTA file with sequences of first 20nt of annotated introns (Default = REF_DIR/fivep_sites.fa)')
	optional_args.add_argument('-n', '--ref_introns', type=pathlib.Path, help='TSV file of all annotated introns (Default = REF_DIR/introns.tsv.gz)')
		# Output options
	optional_args.add_argument('-p', '--output_prefix', help='Add a prefix to output file names (-o OUT -p ABC   ->   OUT/ABC_lariat_reads.tsv)')
	optional_args.add_argument('-u', '--ucsc_track', action='store_true', help='Add an output file named "lariat_reads.bed" which can be used as a custom track in the UCSC Genome Browser (https://www.genome.ucsc.edu/cgi-bin/hgCustom) to visualize lariat alignments')
	optional_args.add_argument('-c', '--keep_classes', action='store_true', help='Keep a file with per-read classification named "read_classes.tsv.gz" in the output (Default = delete)')
	optional_args.add_argument('-k', '--keep_temp', action='store_true', help='Keep all temporary files created while running the pipeline. Forces -c/--keep_classes (Default = delete)')
		# Technical options
	optional_args.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use for parallel processing (Default = 1)')
	optional_args.add_argument('-x', '--skip_version_check', action='store_true', help='Don\'t check if LariatMapper is up-to-date with the main branch on GitHub (Default = check and warn if not up-to-date)')
	log_levels = optional_args.add_mutually_exclusive_group()
	log_levels.add_argument('-q', '--quiet', action='store_true', help="Only print fatal error messages (sets logging level to ERROR, Default = INFO). Mutually exclusive with -w and -d")
	log_levels.add_argument('-w', '--warning', action='store_true', help="Print warning messages and fatal error messages (sets logging level to WARNING, Default = INFO). Mutually exclusive with -q and -d")
	log_levels.add_argument('-d', '--debug', action='store_true', help="Print extensive status messages (sets logging level to DEBUG, Default = INFO). Mutually exclusive with -q and -w")

	# Parse args into a dict
	args = vars(parser.parse_args())
	# Convert to a Settings object
	settings = Settings(**args)

	# Set up logging
	log = functions.get_logger(settings.log_level)

	# Report version
	log.info(f'LariatMapper {functions.version()}')

	# Report arguments
	arg_message = [f'{key}={val}' for key, val in args.items() if val is not None and val is not False]
	arg_message = '\n\t'.join(arg_message)
	log.info(f'Arguments: \n\t{arg_message}')

	# Check if up-to-date
	if settings.skip_version_check is False:
		log.debug('Checking if LariatMapper is up-to-date with the main branch on GitHub...')
		try:
			check_up_to_date(settings.pipeline_dir, log)
		except Exception as e:
			log.debug(e)
			log.warning('Could not check if LariatMapper is up-to-date with the main branch on GitHub. Continuing anyway...')
			time.sleep(60)	# Give the user a chance to see the warning

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
