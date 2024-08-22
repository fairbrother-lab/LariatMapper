import os
import logging
import logging.config
import json
import subprocess
import tempfile

import pandas as pd
import pyfaidx



# =============================================================================#
#                                  Classes                                     #
# =============================================================================#
class RunCommandError(Exception):
	'''
	Custom Exception class raised when subprocess.run(command) fails, 
	allowing us to print the stdout and stderr from the command so we know why 
	it failed
	'''

	def __init__(self, process:subprocess.CompletedProcess):
		super().__init__()
		self.process = process
		self.response = process.stdout + process.stderr

	def __str__(self):
		return f'Command returned non-zero exit status {self.process.returncode}. \n{self.response}'



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
COMP_NTS = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
def reverse_complement(seq:str):
	'''
	Return the reverse-complement of an input DNA strand.
	<seq> must be a string with only the characters "A", "C", "G", "T", and "N" 
	reverse_complement("ACGTN") = "NACGT"
	'''
	return ''.join([COMP_NTS[seq[i]] for i in range(len(seq)-1,-1,-1)])


def str_join(items:list|tuple|pd.Series, join_string:str=',', unique:bool=False) -> str:
	'''
	If <items> contains only 1 unique item, return str(that item)
	Else, return a <join_string>-delimited string of <items>
		e.g.) 	[toy1, toy2, toy2, toy3]			  	  -> "toy1,toy2,toy2,toy3"
			 	[toy1, toy2, toy2, toy3], join_string=";" -> "toy1;toy2;toy2;toy3"
			 	[toy1, toy2, toy2, toy3], unique=True     -> "toy1,toy2,toy3"
	'''
	if len(set(items)) == 1:
		if isinstance(items, pd.Series):
			return str(items.iloc[0])
		else:
			return str(items[0])
	
	if unique is True:
		out = []
		for item in items:
			if item not in out:
				out.append(str(item))
		return join_string.join(out)
	
	else:
		return join_string.join([str(item) for item in items])
	


def align_is_reverse(flag:int) -> bool:
	'''
	Return True if the reverse flag bit 0x10 is 1, else False
	'''
	bit_flags = bin(int(flag))
	is_reverse = True if len(bit_flags)>=7 and bit_flags[-5]=='1' else False
	return is_reverse


def get_logger(level:str) -> logging.Logger:
	'''
	Configure logging settings based on <pipeline_dir>/resources/log_config.json file, and return a logging.Logger object for the specified <level>
	<level> must be "DEBUG", "INFO", "WARNING", or "ERROR"
	'''
	assert level in ('DEBUG', 'INFO', 'WARNING', 'ERROR')
	config_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../resources/log_config.json')
	with open(config_file) as file_in:
		config = json.load(file_in)
	logging.config.dictConfig(config)

	log = logging.getLogger(level.lower())
	return log


def run_command(command:str, log:logging.Logger=None, input:str=None, timeout:int=None) -> str:
	'''
	Wrapper for subprocess.run(call.split(' '), input=input, capture_output=True, text=True) for handling errors and extracting stdout, when appropriate
	'''
	if log is not None:
		log.debug(f'Running command: {repr(command)}')
		if input is not None: 
			if len(input) <= 1_000:
				log.debug(f'Input: \n{input}')
			else:
				log.debug(f'Input (TRUNCATED): \n{input[:1_000]}')

	response = subprocess.run(command.split(' '), 
						   capture_output=True, 
						   text=True, 
						   input=input, 
						   timeout=timeout)
	if response.returncode != 0:
		raise RunCommandError(process=response)
	
	return response.stderr.strip() + response.stdout.strip()


def getfasta(genome_fasta:str, bedtools_input:str, log:logging.Logger=None) -> pd.DataFrame:
	'''
	Get the sequences from <genome_fasta> for the regions specified in <bedtools_input>
	Uses the bedtools command "getfasta" with the args -s -tab -nameOnly
	<bedtools_input> must be a string in BED format with a unique name, 
		e.g. "chr1	100	200	feature1	0	+\nchr1	150	250	feature2	0	-\n"
	Returns a pandas DataFrame with columns ['name', 'seq'], where 'name' is the feature name
	We can't parse the standard output for the sequences because warnings will be included 
	in some lines in a non-deterministic pattern 
	'''
	with tempfile.NamedTemporaryFile() as tmp:
		bedtools_call = f'bedtools getfasta -s -tab -nameOnly -fi {genome_fasta} -fo {tmp.name} -bed -'
		bedtools_output = run_command(bedtools_call, input=bedtools_input, log=log)
		# Check for warnings in bedtools_output
		if bedtools_output != '' and log is not None:
			bedtools_output = bedtools_output.split('\n')
			for warning in bedtools_output:
				log.warning(warning)
		seqs = pd.read_csv(tmp, sep='\t', header=None, names=['name', 'seq'], na_filter=False)
	
	# Remove the strand suffix from the names
	seqs.name = seqs.name.str.slice(0,-3)
	
	return seqs


def get_seq(genome_fasta:str, chrom:str, start:int, end:int, rev_comp:bool) -> str:
	"""
	Retrieves a sequence from a genome FASTA file.
	Args:
		genome_fasta (str): The path to the genome FASTA file.
		chrom (str): The chromosome or sequence identifier.
		start (int): The start position of the sequence (0-based inclusive)
		end (int): The end position of the sequence (0-based exclusive)
		rev_comp (bool): Flag indicating whether to retrieve the reverse complement of the sequence.
	"""
	return pyfaidx.Fasta(genome_fasta, sequence_always_upper=True, as_raw=True).get_seq(chrom, start+1, end+1, rev_comp)


def version() -> str:
	'''
	Return the current version of LariatMapper (format Major.Minor.Patch)
	'''
	toml_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../pyproject.toml')
	with open(toml_file) as file_in:
		for line in file_in:
			if line.startswith('version = "'):
				v = line.lstrip('version = "').rstrip('"\n')
				return v


def linecount(file:str) -> int:
	'''
	Return the number of lines in the input file
	'''
	with open(file) as file_in:
		count = sum(1 for line in file_in)
	
	return count