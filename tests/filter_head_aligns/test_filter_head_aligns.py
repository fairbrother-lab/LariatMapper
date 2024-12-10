import gzip
import os
import pathlib
import shutil
import subprocess
import sys

import pytest

sys.path.append(str(pathlib.Path(__file__).parent.parent.resolve()))
import test_utils



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
PACKAGE_DIR = pathlib.Path(__file__).parent.parent.parent.resolve()
TEST_DIR = PACKAGE_DIR/'tests'/'filter_head_aligns'



# =============================================================================#
#                                   Tests                                      #
# =============================================================================#
@pytest.mark.parametrize('threads',
						['1', '2', '3', '4', '5'])
@pytest.mark.parametrize('prefix',
						['', 'prefix_'])
def test_filter_head_aligns(threads, prefix, tmp_path):
	# Link files needed in working dir
	for file in ('heads_to_genome.sam', 'tails.tsv'):
		os.symlink(TEST_DIR/'inputs'/file, tmp_path/f'{prefix}{file}')

	# Run script
	command = f"python {PACKAGE_DIR/'scripts'/'filter_head_aligns.py'} {threads}" \
			  f" {TEST_DIR/'inputs'/'introns.tsv'} {TEST_DIR/'inputs'/'genome.fa'} {tmp_path}/{prefix} DEBUG"  
	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	if response.returncode != 0:
		pytest.fail(response_text)

	# Check output
	for ref, out in ((TEST_DIR/'outputs'/'failed_head_alignments.tsv', tmp_path/f'{prefix}failed_head_alignments.tsv'),
					(TEST_DIR/'outputs'/'template_switching_reads.tsv', tmp_path/f'{prefix}template_switching_reads.tsv'),
					(TEST_DIR/'outputs'/'circularized_intron_reads.tsv', tmp_path/f'{prefix}circularized_intron_reads.tsv'),
					(TEST_DIR/'outputs'/'putative_lariats.tsv', tmp_path/f'{prefix}putative_lariats.tsv')):
		
		if ref.name == 'putative_lariats.tsv':
			test_utils.check_read_bp_pos(ref)
		elif ref.name == 'circularized_intron_reads.tsv':
			test_utils.check_read_bp_pos(ref, True)
		
		ref_lines, out_lines = test_utils.load_file_lines(ref, out)
		if ref_lines == out_lines:
			continue

		# File contents differ, decide how to report
		print(response_text)
		if test_utils.vscode_available():
			test_utils.vscode_compare(ref, out)
			pytest.fail(f'Output file differs from expected output: {out.name}')
		else:
			assert ref_lines == out_lines