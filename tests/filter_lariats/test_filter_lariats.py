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
SCRIPTS_DIR = PACKAGE_DIR/'lariatmapper'/'scripts'
TEST_DIR = PACKAGE_DIR/'tests'/'filter_lariats'



# =============================================================================#
#                                   Tests                                      #
# =============================================================================#
@pytest.mark.parametrize('prefix',
						['', 
	   					'prefix_'])
@pytest.mark.parametrize('seq_type', 
						 ['paired'])
@pytest.mark.parametrize('verbosity',
						 ['DEBUG', 
						'ERROR'])
def test_filter_lariats(prefix, seq_type, verbosity, tmp_path):
	# Link files needed in working dir
	for file in ('output.bam', 'template_switching_reads.tsv', 'circularized_intron_reads.tsv', 'putative_lariats.tsv'):
		os.symlink(TEST_DIR/'inputs'/file, tmp_path/f'{prefix}{file}')

	# Run script
	command = f"python {SCRIPTS_DIR/'filter_lariats.py'} {tmp_path}/{prefix} {verbosity}" \
			  f" {seq_type} {TEST_DIR/'inputs'/'genome.fa'} {TEST_DIR/'inputs'/'repeatmasker.bed'}" 
	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	if response.returncode != 0:
		pytest.fail(response_text)

	# Check output
	ref_and_out_files = (
		(TEST_DIR/'outputs'/'failed_lariat_alignments.tsv', tmp_path/f'{prefix}failed_lariat_alignments.tsv'),
		(TEST_DIR/'outputs'/'lariat_reads.tsv', tmp_path/f'{prefix}lariat_reads.tsv')
	)
	test_utils.confirm_outputs_match_references(ref_and_out_files, response_text)