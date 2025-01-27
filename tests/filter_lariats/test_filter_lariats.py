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
	for file in ('output.bam', 'circularized_intron_reads.tsv', 'putative_lariats.tsv'):
		os.symlink(TEST_DIR/'inputs'/file, tmp_path/f'{prefix}{file}')

	# Run script
	command = f"python {PACKAGE_DIR/'scripts'/'filter_lariats.py'} {tmp_path}/{prefix} {verbosity}" \
			  f" {seq_type} {TEST_DIR/'inputs'/'genome.fa'} {TEST_DIR/'inputs'/'repeatmasker.bed'}" 
	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	if response.returncode != 0:
		pytest.fail(response_text)

	# Check output
	for ref, out in ((TEST_DIR/'outputs'/'failed_lariat_alignments.tsv', tmp_path/f'{prefix}failed_lariat_alignments.tsv'),
					(TEST_DIR/'outputs'/'lariat_reads.tsv', tmp_path/f'{prefix}lariat_reads.tsv')):
		
		test_utils.check_read_bp_pos(ref)

		ref_lines, out_lines = test_utils.load_file_lines(ref, out)
		if ref_lines == out_lines:
			continue

		### If the lines don't match, report it
		# Print the response text
		print(response_text)

		# If in vscode, open the files for comparison
		if test_utils.vscode_available():
			test_utils.vscode_compare_sorted(ref, out)
		
		# Trigger pytest fail
		assert ref_lines == out_lines