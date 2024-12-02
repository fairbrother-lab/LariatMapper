import gzip
import os
import pathlib
import shutil
import subprocess
import sys

import pandas as pd
import pytest

sys.path.append(str(pathlib.Path(__file__).parent.parent.resolve()))
import test_utils



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
PACKAGE_DIR = pathlib.Path(__file__).parent.parent.parent.resolve()
FILTER_LARIATS_DIR = PACKAGE_DIR/'tests'/'filter_lariats'
INPUTS_DIR = PACKAGE_DIR/'tests'/'inputs'
OUTPUTS_DIR = PACKAGE_DIR/'tests'/'outputs'



# =============================================================================#
#                                   Tests                                      #
# =============================================================================#
@pytest.mark.parametrize('prefix',
						['', 'prefix_'])
@pytest.mark.parametrize('seq_type', 
						 ['paired'])
def test_filter_lariats(prefix, seq_type, tmp_path):
	# Link files needed in working dir
	for file in ('output.bam', 'circularized_intron_reads.tsv', 'putative_lariats.tsv'):
		os.symlink(FILTER_LARIATS_DIR/'inputs'/file, tmp_path/f'{prefix}{file}')

	# Run script
	command = f"python {PACKAGE_DIR/'scripts'/'filter_lariats.py'} {tmp_path}/{prefix} DEBUG" \
			  f" {seq_type} {FILTER_LARIATS_DIR/'inputs'/'genome.fa'} {FILTER_LARIATS_DIR/'inputs'/'repeatmasker.bed'}" 
	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	assert response.returncode == 0, response_text

	# Check output
	for ref, out in ((FILTER_LARIATS_DIR/'outputs'/'failed_lariat_alignments.tsv', tmp_path/f'{prefix}failed_lariat_alignments.tsv'),
					(FILTER_LARIATS_DIR/'outputs'/'lariat_reads.tsv', tmp_path/f'{prefix}lariat_reads.tsv')):
		
		test_utils.check_read_bp_pos(ref)

		ref_lines, out_lines = test_utils.load_file_lines(ref, out)
		if ref_lines == out_lines:
			continue

		# File contents differ, decide how to report
		if test_utils.vscode_available():
			test_utils.vscode_compare(ref, out)
			print(response_text)
			pytest.fail(f'Output file differs from expected output: {out.name}')
		else:
			assert ref_lines == out_lines, response_text
