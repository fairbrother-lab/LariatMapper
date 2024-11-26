import os
import pathlib
import shutil
import gzip
import subprocess

import pandas as pd
import pytest



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
PACKAGE_DIR = pathlib.Path(__file__).parent.parent.resolve()
BUILD_REFERENCES_DIR = PACKAGE_DIR/'tests'/'build_references'
FILTER_FIVEP_ALIGNS_DIR = PACKAGE_DIR/'tests'/'filter_fivep_aligns'
FILTER_HEAD_ALIGNS_DIR = PACKAGE_DIR/'tests'/'filter_head_aligns'
FILTER_LARIATS_DIR = PACKAGE_DIR/'tests'/'filter_lariats'
INPUTS_DIR = PACKAGE_DIR/'tests'/'inputs'
OUTPUTS_DIR = PACKAGE_DIR/'tests'/'outputs'



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def load_file_lines(file_a:str, file_b:str, sort:bool=True) -> tuple[list, list]:
	"""
	Load lines from two files and optionally sort them.
	Parameters:
	- file_a (str): Path to the first file.
	- file_b (str): Path to the second file.
	- sort (bool, optional): Whether to sort the lines. Defaults to True.
	Returns:
	- lines_a (list): List of lines from file_a.
	- lines_b (list): List of lines from file_b.
	"""
	# Load
	if str(file_a).endswith('.gz'):
		with gzip.open(file_a, 'rt') as in_file:
			lines_a = in_file.readlines()
	else:
		with open(file_a) as in_file:
			lines_a = in_file.readlines()
			
	if str(file_b).endswith('.gz'):
		with gzip.open(file_b, 'rt') as in_file:
			lines_b = in_file.readlines()
	else:
		with open(file_b) as in_file:
			lines_b = in_file.readlines()

	# Sort
	if sort is True:
		lines_a.sort()
		lines_b.sort()

	# Compare
	return lines_a, lines_b


def vscode_available() -> bool:
	"""
	Check if the Visual Studio Code executable is available by running shutil.which('code')
	Returns:
	- bool: True if available, False otherwise.
	"""
	return shutil.which('code') is not None
	

def vscode_compare(file_a, file_b):
	"""
	Compare two files using Visual Studio Code.
	Parameters:
	- file_a (str): Path to the first file.
	- file_b (str): Path to the second file.
	"""
	subprocess.run(['code', '--diff', file_a, file_b])


def check_read_bp_pos(file:pathlib.Path, circularized_introns:bool=False):
	"""
	Checks the base pair position in a read sequence from a given file.
	This function reads a tab-separated values (TSV) file into a DataFrame and checks if the base pair 
	position in the read sequence matches the expected nucleotide. If the `circularized_introns` flag 
	is set to True, it checks the 'read_head_end_nt' and 'read_head_end_pos' columns; otherwise, it 
	checks the 'read_bp_nt' and 'read_bp_pos' columns. If any mismatches are found, the function 
	raises a pytest failure with details of the mismatched rows.
	Args:
		file (pathlib.Path): The path to the TSV file containing the read sequences and positions.
		circularized_introns (bool, optional): Flag indicating whether to check for circularized introns. 
												Defaults to False.
	Raises:
		pytest.fail: If any mismatches are found between the expected and actual base pair positions.
	"""
	df = pd.read_csv(file, sep='\t')
	if circularized_introns is True:
		nt_name = 'read_head_end_nt'
		pos_name = 'read_head_end_pos'
	else:
		nt_name = 'read_bp_nt'
		pos_name = 'read_bp_pos'

	read_bp_wrong = df.apply(lambda row: row[nt_name] != row['read_seq_forward'][row[pos_name]], 
								axis=1)
	if read_bp_wrong.any():
		pytest.fail(f'Wrong read_bp_nt in {file}: {df.loc[read_bp_wrong].index}')
