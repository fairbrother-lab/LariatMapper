import os
import logging
import logging.config
import json
import subprocess

import pandas as pd

from scripts import exceptions




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


def run_command(command:str, input:str=None, log:logging.Logger=None) -> str:
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

	response = subprocess.run(command.split(' '), input=input, capture_output=True, text=True)
	if response.returncode != 0:
		raise exceptions.RunCommandError(process=response)
	
	return response.stderr.strip() + response.stdout.strip()