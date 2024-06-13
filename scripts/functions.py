import os
import logging
import logging.config
import json




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


def comma_join(items) -> str:
	'''
	Return ','.join(set(<items>))
	'''
	return ','.join(set(items))


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
	<level> must be "DEBUG", "INFO", or "ERROR"
	'''
	assert level in ('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL')
	config_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../resources/log_config.json')
	with open(config_file) as x:
		config = json.load(x)
	logging.config.dictConfig(config)

	log = logging.getLogger(level.lower())
	return log
	