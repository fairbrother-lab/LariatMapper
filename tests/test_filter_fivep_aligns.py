import pathlib
import subprocess
import os
import shutil

import pytest

from . import functions



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
PACKAGE_DIR = pathlib.Path(__file__).parent.parent.resolve()
SCRIPT = PACKAGE_DIR/'scripts'/'filter_fivep_aligns.py'
INPUTS_DIR = PACKAGE_DIR/'tests'/'inputs'
TEMPS_DIR = PACKAGE_DIR/'tests'/'temps'



# =============================================================================#
#                                   Tests                                      #
# =============================================================================#
@pytest.mark.parametrize('threads',
						['1', '2', '3', '4', '5'])
@pytest.mark.parametrize('prefix',
						['', 'p'])
def test_run(threads, prefix, tmp_path):
	# Link files needed in working dir
	os.symlink(TEMPS_DIR/'unmapped_reads.fa', tmp_path/f'{prefix}unmapped_reads.fa')
	os.symlink(TEMPS_DIR/'unmapped_reads.fa.fai', tmp_path/f'{prefix}unmapped_reads.fa.fai')
	os.symlink(TEMPS_DIR/'fivep_to_reads.sam', tmp_path/f'{prefix}fivep_to_reads.sam')
	os.symlink(TEMPS_DIR/'fivep_sites.fa', tmp_path/f'{prefix}fivep_sites.fa')

	# Run script
	command = f"python {SCRIPT} {threads} {INPUTS_DIR/'genome.fa'} {INPUTS_DIR/'fivep_sites.fa'} {tmp_path}/{prefix} DEBUG"
	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	print(response_text)
	response_code = response.returncode
	assert response_code  == 0, response_text

	# Check output
	functions.files_equal(TEMPS_DIR/'tails.tsv', tmp_path/f'{prefix}tails.tsv')
	functions.files_equal(TEMPS_DIR/'failed_fivep_alignments.tsv', tmp_path/f'{prefix}failed_fivep_alignments.tsv')

