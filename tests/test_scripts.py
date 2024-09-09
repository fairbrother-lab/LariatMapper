import os
import pathlib
import shutil
import gzip
import subprocess

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



# =============================================================================#
#                                   Tests                                      #
# =============================================================================#
# Required args + the transcript attribute and gene attribute args since their values depend on the genome gtf
@pytest.mark.parametrize('base',
					[f"-f {BUILD_REFERENCES_DIR/'inputs'/'genome.fa'} -i {BUILD_REFERENCES_DIR/'inputs'/'hisat2_index'} -x tid -g gid", 
					f"--genome_fasta {BUILD_REFERENCES_DIR/'inputs'/'genome.fa'}  --hisat2_index {BUILD_REFERENCES_DIR/'inputs'/'hisat2_index'} --transcript_attribute tid --gene_attribute gid",],)
@pytest.mark.parametrize('anno',
					[f"-a {BUILD_REFERENCES_DIR/'inputs'/'anno.gtf'}", 
					f"-a {BUILD_REFERENCES_DIR/'inputs'/'anno.gff'}", 
					f"--genome_anno {BUILD_REFERENCES_DIR/'inputs'/'anno.gtf'}"])
# Optional args
@pytest.mark.parametrize('repeatmasker_bed',
					[None, 
					f"-r {BUILD_REFERENCES_DIR/'inputs'/'repeatmasker.bed'}", 
					f"--repeatmasker_bed {BUILD_REFERENCES_DIR/'inputs'/'repeatmasker.bed'}"])
@pytest.mark.parametrize('threads',
					[None, '-t 4', '--threads 4'])
@pytest.mark.parametrize('copy',
					[None, '--copy'])
def test_build_references(base, anno, repeatmasker_bed, threads, copy, tmp_path):
	command = f"python {PACKAGE_DIR/'build_references.py'} {base} {anno} -o {tmp_path}"
	for optional_arg in (repeatmasker_bed, threads, copy):
		if optional_arg is not None:
			command = f'{command} {optional_arg}'

	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	assert response.returncode == 0, response_text

	# Check output
	for ref, out in ((BUILD_REFERENCES_DIR/'outputs'/'exons.tsv', tmp_path/'exons.tsv.gz'),
					(BUILD_REFERENCES_DIR/'outputs'/'introns.tsv', tmp_path/'introns.tsv.gz'),
					(BUILD_REFERENCES_DIR/'outputs'/'fivep_sites.fa', tmp_path/'fivep_sites.fa')):
		ref_lines, out_lines = load_file_lines(ref, out)
		if ref_lines == out_lines:
			continue

		# If file contents differ, decide how to report
		if vscode_available():
			vscode_compare(ref, out)
			pytest.fail(f'Output file differs from expected output: {out.name}')
		else:
			assert ref_lines == out_lines


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
		ref_lines, out_lines = load_file_lines(ref, out)
		if ref_lines == out_lines:
			continue

		# If file contents differ, decide how to report
		if vscode_available():
			vscode_compare(ref, out)
			pytest.fail(f'Output file differs from expected output: {out.name}')
		else:
			assert ref_lines == out_lines


@pytest.mark.parametrize('threads',
						['1', '2', '3', '4', '5'])
@pytest.mark.parametrize('prefix',
						['', 'prefix_'])
def test_filter_head_aligns(threads, prefix, tmp_path):
	# Link files needed in working dir
	for file in ('heads_to_genome.sam', 'tails.tsv'):
		os.symlink(FILTER_HEAD_ALIGNS_DIR/'inputs'/file, tmp_path/f'{prefix}{file}')

	# Run script
	command = f"python {PACKAGE_DIR/'scripts'/'filter_head_aligns.py'} {threads}" \
			  f" {FILTER_HEAD_ALIGNS_DIR/'inputs'/'introns.tsv'} {FILTER_HEAD_ALIGNS_DIR/'inputs'/'genome.fa'} {tmp_path}/{prefix} DEBUG"  
	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	assert response.returncode == 0, response_text

	# Check output
	for ref, out in ((FILTER_HEAD_ALIGNS_DIR/'outputs'/'failed_head_alignments.tsv', tmp_path/f'{prefix}failed_head_alignments.tsv'),
					(FILTER_HEAD_ALIGNS_DIR/'outputs'/'template_switching_reads.tsv', tmp_path/f'{prefix}template_switching_reads.tsv'),
					(FILTER_HEAD_ALIGNS_DIR/'outputs'/'circularized_intron_reads.tsv', tmp_path/f'{prefix}circularized_intron_reads.tsv'),
					(FILTER_HEAD_ALIGNS_DIR/'outputs'/'putative_lariats.tsv', tmp_path/f'{prefix}putative_lariats.tsv')):
		ref_lines, out_lines = load_file_lines(ref, out)
		if ref_lines == out_lines:
			continue

		# File contents differ, decide how to report
		print(response_text)
		if vscode_available():
			vscode_compare(ref, out)
			pytest.fail(f'Output file differs from expected output: {out.name}')
		else:
			assert ref_lines == out_lines


@pytest.mark.parametrize('prefix',
						['', 'prefix_'])
@pytest.mark.parametrize('seq_type', 
						 ['paired'])
def test_filter_lariats(prefix, seq_type, tmp_path):
	# Link files needed in working dir
	for file in ('output.bam', 'template_switching_reads.tsv', 'circularized_intron_reads.tsv', 'putative_lariats.tsv'):
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
		ref_lines, out_lines = load_file_lines(ref, out)
		if ref_lines == out_lines:
			continue

		# File contents differ, decide how to report
		if vscode_available():
			vscode_compare(ref, out)
			print(response_text)
			pytest.fail(f'Output file differs from expected output: {out.name}')
		else:
			assert ref_lines == out_lines, response_text
