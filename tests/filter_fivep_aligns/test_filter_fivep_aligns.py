import os
import pathlib
import subprocess
import sys

import pytest

sys.path.append(str(pathlib.Path(__file__).parent.parent.resolve()))
import test_utils



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
PACKAGE_DIR = pathlib.Path(__file__).parent.parent.parent.resolve()
SCRIPTS_DIR = PACKAGE_DIR/'LariatMapper'/'scripts'
TEST_DIR = PACKAGE_DIR/'tests'/'filter_fivep_aligns'



# =============================================================================#
#                                   Tests                                      #
# =============================================================================#
@pytest.mark.parametrize('threads',
						['1', 
						'3', 
						'8'])
@pytest.mark.parametrize('prefix',
						['', 
						'prefix_'])
@pytest.mark.parametrize('strand',
						 ['Unstranded'])
@pytest.mark.parametrize('verbosity',
						 ['DEBUG', 
						'ERROR'])
def test_filter_fivep_aligns(threads, prefix, strand, verbosity, tmp_path):
	# Link files needed in working dir
	for file in ('unmapped_reads.fa', 'unmapped_reads.fa.fai', 'fivep_to_reads.sam', 'fivep_sites.fa'):
		os.symlink(TEST_DIR/'inputs'/file, tmp_path/f'{prefix}{file}')

	# Run script
	command = f"python {SCRIPTS_DIR/'filter_fivep_aligns.py'} {tmp_path}/{prefix} {verbosity}" \
			  f" {TEST_DIR/'inputs'/'genome.fa'} {TEST_DIR/'inputs'/'fivep_sites.fa'} {strand} {threads}"
	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	if response.returncode != 0:
		pytest.fail(response_text)
		
	# Check output
	ref_and_out_files = (
		(TEST_DIR/'outputs'/'failed_fivep_alignments.tsv', tmp_path/f'{prefix}failed_fivep_alignments.tsv'),
		(TEST_DIR/'outputs'/'tails.tsv', tmp_path/f'{prefix}tails.tsv'),
		(TEST_DIR/'outputs'/'heads.fa', tmp_path/f'{prefix}heads.fa')
	)
	test_utils.confirm_outputs_match_references(ref_and_out_files, response_text)
