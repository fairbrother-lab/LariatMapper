from logging import Logger
from os.path import join, getsize, isdir, isfile
from os import listdir, mkdir, chmod
import stat
from shutil import move
from subprocess import run
from gzip import open as gzopen
from csv import reader
from dataclasses import dataclass
from logger import setup_logger
from argparse import ArgumentParser

#=============================================================================#
#                                  Constants                                  #
#=============================================================================#
STEPS = ['mkdir','fqsplit', 'map', 'combine', 'merge_filter']

ENV_DIR = "/datasets2/lariat_mapping/larmap_env"

TIMESTAMP_SCRIPT = "/datasets2/lariat_mapping/scripts/timestamp.sh"
MAP_LARIATS_SCRIPT = "/datasets2/lariat_mapping/scripts/map_lariats.sh"
FILTER_5p_SCRIPT = "/datasets2/lariat_mapping/scripts/filter_fivep_alignments.py"
FILTER_3p_SCRIPT = "/datasets2/lariat_mapping/scripts/filter_threep_alignments.py"
COMBINE_SCRIPT = "/datasets2/lariat_mapping/scripts/combine_split_maps.py"

TRANSLATE_TO_IMP = {'Raw Fastq Files Directory:': 'fastq_dir',	
					'Child Scripts Directory:': 'child_dir',		
					'Log Directory:': 'log_dir',						
					'Mapping Out Directory:': 'out_dir',
					'Final Results Path:': 'results_path',						
					'SBATCH ntasks:': 'SBntasks',						
					'SBATCH cpus:': 'SBcpus',						
					'SBATCH memory (in Gb):': 'SBmem',
					'Bowtie2 Genome Index Base:': 'ref_b2index',
					'Genome FASTA File:': 'ref_genome',
					'Genome Annotations File:': 'ref_anno',
					'Fiveprime Splice Sites File:': 'ref_fivepsites',
					'Threeprime Splice Sites Bowtie2 Index Base:': 'ref_threepb2index',
					'Threeprime Lengths File:': 'ref_threeplengths',
					'Introns BED File:': 'ref_introns',		
					'Repeat-Masker BED File:': 'ref_repeatmasker'
					}





#=============================================================================#
#                          Implementation Info Class                          #
#=============================================================================#
@dataclass
class Imp:
	ATTRS = ['fastq_dir', 'child_dir', 'log_dir', 'out_dir', 'results_path', 'SBntasks', 'SBcpus', 'SBmem', 'ref_b2index', 'ref_genome', 'ref_anno', 'ref_fivepsites', 'ref_threeplengths', 'ref_threepb2index', 'ref_introns', 'ref_repeatmasker', 'fq_files', 'fq_categories']

	fastq_dir: str
	child_dir: str
	log_dir: str
	out_dir: str
	results_path: str
	ref_files: str
	SBntasks: int
	SBcpus: int
	SBmem: int
	ref_b2index: str
	ref_genome: str
	ref_anno: str
	ref_fivepsites: str
	ref_threeplengths: str
	ref_threepb2index: str
	ref_introns: str
	ref_repeatmasker: str
	fq_files: list
	fq_categories: list



	def __init__(self, info_file, log) -> None:
		# with open(join(INFO_DIR, f'{info_file}')) as info:
		with open(info_file, 'r') as info:
			info_reader = reader(info, delimiter='\t')

			for row in info_reader:
				log.info(f'Pulled row: {row}')
				row = [item for item in row if item != '']	#Remove any empty values from extra tabs
				
				if row == [] or row[0].startswith('#'):
					continue 
				if row[0] in TRANSLATE_TO_IMP.keys():
					variable, value = row[:2]
					log.info(f'Got {variable}: {value}')
					variable = TRANSLATE_TO_IMP[variable]
					setattr(self, variable, value)
				if row[0] == 'Original File Name':
					break

			log.info(f'The imp, halfway complete: {self}')

			#row 'should' be at the Original File Name row
			self.fq_categories = row
			self.fq_files = []
			log.info(f'Imp set for {len(self.fq_categories)} categories (including 3 names): {self.fq_categories}')
			for row in info_reader:
				log.info(f'Pulled row: {row}')

				fq = {}
				for i in range(len(self.fq_categories)):
					fq[self.fq_categories[i]] = row[i]
				self.fq_files.append(fq)

			for fq in self.fq_files:
				for i in range(len(self.fq_files)):
					if self.fq_files[i]["New Name"] == fq["Mate New Name"]:
						fq["Mate Position"] = i

			self.fq_files, self.fq_categories = tuple(self.fq_files), tuple(self.fq_categories)



	def __repr__(self):
		rep = ""
		for attr, val in vars(self).items():
			if not attr.startswith('__'):
				rep += f'{attr}: {val}\t\t({type(val)})\n\t\t'
		return rep



	def validate(self) -> None:
		'''
		Double-checks the provided info file for common errors
		'''
		#Do all the imp's attributes have values?
		for attr in Imp.ATTRS:
			if vars(self)[attr] in ('', None):
				raise ValueError(f"Invalid value ({vars(self)[attr]}) for the imp's {attr} attribute")
		

		#Does the fastq directory exist?
		if not isdir(self.fastq_dir):
			raise ValueError(f'The directory for the raw fastq files, "{self.fastq_dir}", was not found.')

		#Are the SBATCH settings valid?
		for SB_var, SB_val in (('ntasks', self.SBntasks), 
								('number of cpus', self.SBcpus),
								('memory allocation (in Gb)', self.SBmem)):
			if not SB_val.isdigit() or int(SB_val) < 1:
				raise ValueError(f'The value for the {SB_var} SBATCH setting, {SB_val}, is not a positive integer.')

		#Does the bowtie2 index exist?
		for suffix in ('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2'):
			if not isfile(self.ref_b2index + suffix):
				raise ValueError(f'The Bowtie2 index file "{self.ref_b2index + suffix}" does not exist')

		#Do the reference files exist?
		for ref in (self.ref_b2index + '.1.bt2',
					self.ref_b2index + '.rev.2.bt2',
					self.ref_genome,
					self.ref_anno,
					self.ref_fivepsites,
					self.ref_threepb2index + '.1.bt2',
					self.ref_threepb2index + '.rev.2.bt2',
					self.ref_threeplengths,
					self.ref_introns,
					self.ref_repeatmasker):
			if not isfile(ref):
				raise ValueError(f'The reference file "{ref}" does not exist.')

		# #Are there an even number of fastq files, allowing for opposite read pairs?
		# if len(self.fq_files) % 2 != 0:
		# 	raise ValueError(f'There are an odd number of fastq files ({len(self.fq_files)}, which means at least one is missing an opposite read partner')

		#Are the values for each fq file's variables valid? 
		for fq in self.fq_files:
			filepath = join(self.fastq_dir, fq['Original File Name'])
			if not isfile(filepath):
				raise ValueError(f'The raw fastq file "{filepath}" was not found.')
			for bad_char in ('split', 'map', '$', '%', '#', "'", '"', '/', '\\', '{', '}', ','):
				if bad_char in fq['Original File Name']:
					raise ValueError(f'{fq["Original File Name"]} cannot contain the substring " {bad_char} "')
				if bad_char in fq['New Name']:
					raise ValueError(f'The New Name for the file "{fq["Original File Name"]}", {fq["New Name"]}, cannot contain the substring " {bad_char} "')
			




#=============================================================================#
#                                  Functions                                  #
#=============================================================================#
#TODO: Finish writing function descriptions
#TODO: Change setuplogger timing from miliseconds to something else
def parse(log: Logger) -> dict:
	'''
	Parse the arguments
	'''
	parser = ArgumentParser(description='Top-level script for mapping lariats reads from raw RNA sequences reads')
	#The path to the implementation info file 
	parser.add_argument('info_file',
						help="The name of the information file (including the file extension) "\
							"that provides this script the details necessary for the lariat mapping process. "\
							"The path may be relative or absolute")
	#Optional: When writing the child scripts, skip to the start of the given step in the pipeline
	parser.add_argument('-skto', '--skip-to', 	#skip to [step]
						choices=STEPS,
						help="Skip all the steps in the pipeline that come before the given step. This is achieved by omitting "\
							"the appropriate sections of code when creating each child bash script. "\
							"example: -skto=map")
	#Optional: When writing the child scripts, stop before the start of the given step in the pipeline
	parser.add_argument('-stbe', '--stop-before',	#stop before [step]
						choices=STEPS,
						help="Stop before the given step. This is achieved by omitting "\
							"the appropriate sections of code when creating each child bash script. "\
							"example: -stbe=combine")

	#Collect all the passed arguments into a dict with their keys equal to
	#their name specified in parser.add arguments. vars() makes the default
	#return (a Namespace obj) into a dict.				
	args = vars(parser.parse_args())
	log.info(f'Arguments passed:   {args}')

	#Set things up in accordance with skto and stbe
	if args['skip_to'] and args['stop_before']:
		skip_index = STEPS.index(args['skip_to'])
		stop_index = STEPS.index(args['stop_before'])
		if skip_index >= stop_index:
			raise ValueError('You provided an improper combination of skip-to and stop-before arguments '\
								f'("{args["skip_to"]}" and "{args["stop_before"]}", respectively)).'\
								'The skip-to step must come before the stop-before step.')

	args['steps'] = STEPS
	if args['skip_to']:
		log.info(f'Skipping to {args["skip_to"]}')
		#remove steps from the start until you reach the target step
		while True:
			if args['steps'][0] == args['skip_to']:
				break
			else:
				args['steps'].pop(0)

	if args['stop_before']:
		log.info(f'Stopping before {args["stop_before"]}')
		#remove steps from the end until you reach the target step
		while  args['steps']:
			if args['steps'][-1] == args['stop_before']:
				args['steps'].pop(-1)
				break
			else:
				args['steps'].pop(-1)



	log.info(f'Steps to perform: {args["steps"]}')
	return args



def prep_directories(imp: Imp, log: Logger) -> None:
	'''
	Checks whether or not the child, log, and out directories provided in the info file exist.
	If not, they are created. 
	Also checks for the child_logs and map_logs subdirectories in the log directory
	'''
	for d in [imp.child_dir, imp.log_dir, imp.out_dir]:
		if not isdir(d):
			log.info(f'Did not find the directory "{d}", so I will create it')
			mkdir(d)
	log_contents = listdir(imp.log_dir)
	log.info(f'Log contents: {log_contents}')
	# for sd in ['child_logs', 'map_logs']:
	for sd in ['child_logs']:
		if sd not in log_contents:
			log.info(f'Did not find the directory "{sd}" in {imp.log_dir}, so I will create it')
			mkdir(join(imp.log_dir, sd))
	log.info('Dictionaries have been prepared.')



def get_chromosomes(imp: Imp, log: Logger) -> list:
	'''
	Collects chromosome info from the bowtie2 index
	'''
	log.info(f'Grabbing chromosomes from the bowtie2 index at {imp.ref_b2index}')
	chroms = []
	inspection = run(f'bowtie2-inspect --summary {imp.ref_b2index}'.split(' '),
					capture_output=True)
	out = inspection.stdout.decode()
	err = inspection.stderr.decode()
	returncode = inspection.returncode
	log.info(f'Out: {out}')
	assert returncode == 0, f'Chromosome-grab error. Returncode = {returncode}, Error = {err}'

	seqs=out.split('\n')[5:-1]
	for s in seqs:
		chrom = s.split('\t')[1]
		if chrom == '' or not chrom.startswith('chr'): continue
		chroms.append(chrom)
	
	assert chroms != [], f'Failed to get any chromosomes from the bowtie2 index. Inspect summary: {out}'
	log.info(f'Collected {len(chroms)} chromosomes:\t{chroms}')
	return chroms



def get_nsplits(fq: dict, imp: Imp, log: Logger) -> int:
	'''
	Measures the size of the given fastq file in Gb and returns the size // 3
	This is how many splits fastqsplitter will be ordered to make for that file
	'''
	file_path = join(imp.fastq_dir, fq['Original File Name'])
	size = getsize(file_path) / 1000000000		# Get size in Gb
	nsplits = size // 3							# Splits to make the split files 3-5.99 Gb
	if nsplits == 0: 
		nsplits = 1 
	if nsplits == 1:
		log.info(f'"{fq["Original File Name"]}" is {size} Gb, small enough to process it without splitting. I will proceed with making 1 "split" of it which is just a copy')
	else:
		log.info(f'"{fq["Original File Name"]}" is {size} Gb, and will be split into {nsplits} chunks of ~{size/nsplits} Gb')
	return nsplits



def get_read_len(fq:dict, imp: Imp, log: Logger) -> int:
	'''
	Peeks at the first read in the given fastq file and measures the read length.
	'''
	file_path = join(imp.fastq_dir, fq['Original File Name'])
	# Check whether or not the fastq file is gzipped
	with open(file_path, 'rb') as f:
		is_gzipped = f.read(2) == b'\x1f\x8b'

	if is_gzipped:
		with gzopen(file_path, 'r') as gzf:
			#skip the first line
			gzf.readline()
			#Get the line with an RNA read and take its length
			read = gzf.readline().strip(b'\n')
	else:
		with open(file_path, 'r') as f:
			#skip the first line
			f.readline()
			#Get the line with an RNA read and take its length
			read = f.readline().strip('\n')
			
	log.info(f'Grabbed "{read}" ({int(len(read))} characters) from {fq["Original File Name"]}')
	return int(len(read))



def write_mkdir_calls(fq:dict, nsplits:int, imp, log) -> str:
	'''
	make a filewise results dir in the full results dir to hold all split files from the same og file
	'''
	call = f'mkdir {imp.out_dir}/{fq["New Name"]} \n'

	for ns in range(1, nsplits+1):
		call += f'mkdir {imp.out_dir}/{fq["New Name"]}/{fq["New Name"]}_split{ns} \n'

	call += '\n\n\n'
	log.info(f'mkdir calls: \n{call}')
	return call



def write_fastqsplit_call(fq:dir, nsplits:int, imp, log) -> str:
	'''
	The call to fastqsplitter to split the input fastq file into n fastq file "chunks"
	-i==input file, -o=output split files
	'''
	#echo start
	call = f'timestamp "Beginning fastsplit..." 1 \n\n'
	#start of call + in-file
	call += f'fastqsplitter -i {imp.fastq_dir}/{fq["Original File Name"]} \\\n'

	# iterate for n outfiles
	for n in range(1, nsplits+1):
		split_name = f'{fq["New Name"]}_split{n}'
		call +=f'-o {imp.out_dir}/{fq["New Name"]}/{split_name}/{split_name}.fq.gz \\\n'

	#chop of last line's \\\n
	call = call[:-3]
	#echo end
	call += f'\n\ntimestamp "Fastqsplit complete" 1 \n\n\n\n'
	log.info(f'fastqsplitter call: \n{call}')
	return call



def write_map_calls(fq:dict, read_len:int, nsplits:int, chroms:list, imp, log) -> str:
	'''
	n map_lariats calls to map each split file, one after the other
	'''
	# echo start
	calls = f'timestamp "Beginning lariat mapping..." 1 \n\n'

	# convert the chromsomes list into a comma-delimited string for the map_lariats script so that it can pass them into the filter_threep_alignment.py file, 
	# which needs the set of the sample organism's chromosomes 
	chroms = ','.join(chroms)

	# iterate for n outfiles
	for n in range(1, nsplits+1):
		split_name = f'{fq["New Name"]}_split{n}'

		calls += f'{MAP_LARIATS_SCRIPT} \\\n'\
				f'{TIMESTAMP_SCRIPT} \\\n'\
				f'{imp.out_dir}/{fq["New Name"]}/{split_name}/{split_name}.fq.gz \\\n'\
				f'{read_len} \\\n'\
				f'{split_name} \\\n'\
				f'{imp.ref_b2index} \\\n'\
				f'{imp.ref_genome} \\\n'\
				f'{imp.ref_fivepsites} \\\n'\
				f'{imp.ref_threepb2index} \\\n'\
				f'{imp.ref_threeplengths} \\\n'\
				f'{FILTER_5p_SCRIPT} \\\n'\
				f'{FILTER_3p_SCRIPT} \\\n'\
				f'{imp.ref_anno} \\\n'\
				f'{chroms} \\\n'\
				f'{imp.SBcpus} \\\n'\
				f'{imp.out_dir}/{fq["New Name"]}/{split_name} \n\n'
				
		# echo for the call
		calls += f'timestamp "Lariat mapping of split {n} complete" 1 \n\n'

	# echo end
	calls += f'timestamp "All maps complete" 3 \n\n\n\n'
	log.info(f'map_lariats calls: \n{calls}')
	return calls



def write_combine_call(fq:dict, imp, nsplits:int, log) -> str:
		call = f'timestamp "Combining mapping results..." 1 \n\n'
		call += f'python {COMBINE_SCRIPT} {imp.out_dir} $LOGFILE_PATH {fq["New Name"]} {nsplits}\n\n'
		call += f'timestamp "Combination complete" 1 \n\n\n\n'
		log.info(f'combine_split_maps call: \n{call}')
		return call



def write_child(fq: dict, nsplits: int, read_len: int, chroms:list, args:dict, imp, log) -> None:
	'''
	
	'''
	log.info(f'Creating child for {fq["Original File Name"]}, AKA {fq["New Name"]}')
	nsplits = int(nsplits)	#For some reason this has to be done even though type(nsplits) is 'int' in main() >:|

	child_name = f'larmap_{fq["New Name"]}'
	logfile_path = f'"{imp.log_dir}/child_logs/{child_name}-$SLURM_JOB_ID.out"'

	with open(f'{imp.child_dir}/{child_name}.sh', 'w') as child:
		# Write the sbatch settings for the script
		# child.write('#!/bin/bash\n\n'\
		# Running the slurm job in login mode by adding -l  to the first line allows the conda activate lariat_mapping-specific env process
		# from https://stackoverflow.com/questions/68094835/how-to-load-anaconda-virtual-environment-from-slurm
		child.write('#!/bin/bash -l\n\n'\
					f'#SBATCH --ntasks={imp.SBntasks}\n'\
					f'#SBATCH --cpus-per-task={imp.SBcpus}\n'\
					f'#SBATCH --mem={imp.SBmem}G\n'\
					f'#SBATCH --job-name={child_name}\n'\
					f'#SBATCH --output={imp.log_dir}/child_logs/{child_name}-%j.out\n'\
					f'#SBATCH --error={imp.log_dir}/child_logs/{child_name}-%j.out\n\n')

		# SECONDS environmental variable to track time elapsed and its subsequent referals
		child.write('SECONDS=0 \n')
		# SLURM_JOB_ID environmental variable to check its own log
		child.write('JOB_ID="$SLURM_JOB_ID" \n')
		#NSPLITS variable to pass in the number of split files made
		child.write(f'NSPLITS={nsplits}\n')
		# timestamp function to timestamp log messages
		child.write(f'source {TIMESTAMP_SCRIPT} \n')
		# Infer the child's log file path from the job id
		child.write(f'LOGFILE_PATH={logfile_path} \n')
		child.write('timestamp "My job id is $SLURM_JOB_ID, and my output will be written to $LOGFILE_PATH." 1 \n')
		# source your conda installation
		child.write('source ~/anaconda3/etc/profile.d/conda.sh \n')
		# activate the conda env
		child.write(f'if [[ "$(conda config --show envs_dirs)" != *"- {ENV_DIR}"* ]]; then \n'\
						f'\tconda config --append envs_dirs {ENV_DIR} \n'\
					'fi \n')
		child.write(f'conda activate {ENV_DIR} \n')
		child.write(f'timestamp "The conda environment is now $CONDA_PREFIX" \n\n\n\n')

		# make a filewise results dir in the full results dir to hold all split files from the same og file
		mkdir_calls = write_mkdir_calls(fq, nsplits, imp, log)
		# split the input fastq file into n fastq file "chunks"
		fastqsplit_call = write_fastqsplit_call(fq, nsplits, imp, log)
		# map_calls = write_map_calls(fq, read_len, nsplits, time_elapsed, time_now, imp)
		map_calls = write_map_calls(fq, read_len, nsplits, chroms, imp, log)
		# combine the results of the mapping for all splits
		# combine_call = write_combine_call(fq, imp, logfile_path, nsplits, log)
		combine_call = write_combine_call(fq, imp, nsplits, log)

		#for each call, check that the step it represents is in the arg['steps'] list
		#if it IS, write it to the child. if it IS NOT, omit it.
		stepcalls = (('mkdir', mkdir_calls),
					('fqsplit', fastqsplit_call), 
					('map', map_calls),
					('combine', combine_call))
		for step, call in stepcalls:
			if step in args['steps']: child.write(call)

	#Set the child's permissions so everyone group can rwx 
	chmod(f'{imp.child_dir}/{child_name}.sh', stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)

	log.info(f'The child for {fq["New Name"]} is complete')



def write_sbatch_all(imp: Imp, log: Logger) -> None:
	'''
	Write a bash script in imp.child_dir that sbatch's all children
	'''
	with open(f'{imp.child_dir}/_sbatch_all.sh', 'w') as f:
		f.write('!#/bin/bash\n\n\n\n')
		for fq in imp.fq_files:
			f.write(f'sbatch {imp.child_dir}/larmap_{fq["New Name"]}.sh\n')

	#Set the child's permissions so everyone group can rwx 
	chmod(f'{imp.child_dir}/_sbatch_all.sh', stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)

	log.info(f'"{imp.child_dir}/_sbatch_all.sh" has been created\n\n')


	


#=============================================================================#
#                                    Main                                     #
#=============================================================================#
def main():
	log = setup_logger('top_log.out')

	args = parse(log)

	imp = Imp(args['info_file'], log)
	log.info(f'The imp: {imp}')
	imp.validate()						#make sure the imp's attributes have appropriate values
	log.info('The imp has passed inspection')

	#make sure all required directories and subdirectories exist
	prep_directories(imp, log)

	#move the log file to the imp.log_dir directory
	move('top_log.out', f'{imp.log_dir}/top_log.out')

	#collect a list of the organism's chromsomes
	chroms = get_chromosomes(imp, log)

	for fq in imp.fq_files:
		nsplits = get_nsplits(fq, imp, log)
		read_len = get_read_len(fq, imp, log)
		#Create the child bash script that will direct the lariat mapping of fq
		write_child(fq, nsplits, read_len, chroms, args, imp, log)

	write_sbatch_all(imp, log)

	log.info('Program complete.')
			


if __name__ == '__main__':
	main()