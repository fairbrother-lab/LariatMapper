import gzip
import os
import pathlib
import shutil
import subprocess
import sys

import pytest

sys.path.append(str(pathlib.Path(__file__).parent.parent.resolve()))
import testing_tools as test_funcs



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
PACKAGE_DIR = pathlib.Path(__file__).parent.parent.parent.resolve()
FILTER_FIVEP_ALIGNS_DIR = PACKAGE_DIR/'tests'/'filter_fivep_aligns'
INPUTS_DIR = PACKAGE_DIR/'tests'/'inputs'
OUTPUTS_DIR = PACKAGE_DIR/'tests'/'outputs'



# =============================================================================#
#                                   Tests                                      #
# =============================================================================#
@pytest.mark.parametrize('threads',
						['1', '2', '3', '4', '5'])
@pytest.mark.parametrize('prefix',
						['', 'prefix_'])
@pytest.mark.parametrize('strand',
						 ['Unstranded'])
def test_filter_fivep_aligns(threads, prefix, strand, tmp_path):
	# Link files needed in working dir
	for file in ('unmapped_reads.fa', 'unmapped_reads.fa.fai', 'fivep_to_reads.sam', 'fivep_sites.fa'):
		os.symlink(FILTER_FIVEP_ALIGNS_DIR/'inputs'/file, tmp_path/f'{prefix}{file}')

	# Run script
	command = f"python {PACKAGE_DIR/'scripts'/'filter_fivep_aligns.py'} {tmp_path}/{prefix} DEBUG" \
			  f" {FILTER_FIVEP_ALIGNS_DIR/'inputs'/'genome.fa'} {FILTER_FIVEP_ALIGNS_DIR/'inputs'/'fivep_sites.fa'} {strand} {threads}"
	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	assert response.returncode == 0, response_text

	# Check output
	for ref, out in ((FILTER_FIVEP_ALIGNS_DIR/'outputs'/'failed_fivep_alignments.tsv', tmp_path/f'{prefix}failed_fivep_alignments.tsv'),
				  	(FILTER_FIVEP_ALIGNS_DIR/'outputs'/'tails.tsv', tmp_path/f'{prefix}tails.tsv'),
				  	(FILTER_FIVEP_ALIGNS_DIR/'outputs'/'heads.fa', tmp_path/f'{prefix}heads.fa')):
		ref_lines, out_lines = test_funcs.load_file_lines(ref, out)
		if ref_lines == out_lines:
			continue

		# If file contents differ, decide how to report
		if test_funcs.vscode_available():
			test_funcs.vscode_compare(ref, out)
			pytest.fail(f'Output file differs from expected output: {out.name}')
		else:
			assert ref_lines == out_lines