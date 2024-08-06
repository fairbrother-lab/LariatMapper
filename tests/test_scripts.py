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
	- sort (bool, optional): Whether to sort the lines. Defaults to False.
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
	which_check = shutil.which('code')
	if which_check is None:
		return False
	else:
		return True
	

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
					[f"-f {INPUTS_DIR/'genome.fa'} -i {INPUTS_DIR/'hisat2_index'} -x tid -g gid", 
					f"--genome_fasta {INPUTS_DIR/'genome.fa'}  --hisat2_index {INPUTS_DIR/'hisat2_index'} --transcript_attribute tid --gene_attribute gid",],)
@pytest.mark.parametrize('anno',
					[f"-a {INPUTS_DIR/'anno.gtf'}", 
					f"-a {INPUTS_DIR/'anno.gff'}", 
					f"--genome_anno {INPUTS_DIR/'anno.gtf'}"])
# Optional args
@pytest.mark.parametrize('repeatmasker_bed',
					[None, 
					f"-r {INPUTS_DIR/'repeatmasker.bed'}", 
					f"--repeatmasker_bed {INPUTS_DIR/'repeatmasker.bed'}"])
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
	for ref, out in ((INPUTS_DIR/'exons.tsv', tmp_path/'exons.tsv.gz'),
					(INPUTS_DIR/'introns.tsv', tmp_path/'introns.tsv.gz'),
					(INPUTS_DIR/'fivep_sites.fa', tmp_path/'fivep_sites.fa')):
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
@pytest.mark.parametrize('seq_type', 
						 ['paired'])
def test_filter_fivep_aligns(threads, prefix, seq_type, tmp_path):
	# Link files needed in working dir
	for file in ('unmapped_reads.fa', 'unmapped_reads.fa.fai', 'fivep_to_reads.sam', 'fivep_sites.fa'):
		os.symlink(OUTPUTS_DIR/seq_type/file, tmp_path/f'{prefix}{file}')

	# Run script
	command = f"python {PACKAGE_DIR/'scripts'/'filter_fivep_aligns.py'} {threads}" \
			  f" {INPUTS_DIR/'genome.fa'} {INPUTS_DIR/'fivep_sites.fa'} {tmp_path}/{prefix} DEBUG"
	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	assert response.returncode == 0, response_text

	# Check output
	for ref, out in ((OUTPUTS_DIR/seq_type/'failed_fivep_alignments.tsv', tmp_path/f'{prefix}failed_fivep_alignments.tsv'),
				  	(OUTPUTS_DIR/seq_type/'tails.tsv', tmp_path/f'{prefix}tails.tsv')):
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
@pytest.mark.parametrize('seq_type', 
						 ['paired'])
def test_filter_head_aligns(threads, prefix, seq_type, tmp_path):
	# Link files needed in working dir
	for file in ('heads_to_genome.sam', 'tails.tsv'):
		os.symlink(OUTPUTS_DIR/seq_type/file, tmp_path/f'{prefix}{file}')

	# Run script
	command = f"python {PACKAGE_DIR/'scripts'/'filter_head_aligns.py'} {threads}" \
			  f" {INPUTS_DIR/'introns.tsv'} {INPUTS_DIR/'genome.fa'} {tmp_path}/{prefix} DEBUG"  
	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	assert response.returncode == 0, response_text

	# Check output
	for ref, out in ((OUTPUTS_DIR/seq_type/'failed_head_alignments.tsv', tmp_path/f'{prefix}failed_head_alignments.tsv'),
					(OUTPUTS_DIR/seq_type/'template_switching_reads.tsv', tmp_path/f'{prefix}template_switching_reads.tsv'),
					(OUTPUTS_DIR/seq_type/'circularized_intron_reads.tsv', tmp_path/f'{prefix}circularized_intron_reads.tsv'),
					(OUTPUTS_DIR/seq_type/'putative_lariats.tsv', tmp_path/f'{prefix}putative_lariats.tsv')):
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
		os.symlink(OUTPUTS_DIR/seq_type/file, tmp_path/f'{prefix}{file}')

	# Run script
	command = f"python {PACKAGE_DIR/'scripts'/'filter_lariats.py'} {tmp_path}/{prefix} DEBUG " \
			  f" {seq_type} {INPUTS_DIR/'genome.fa'} {INPUTS_DIR/'repeatmasker.bed'}" 
	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	assert response.returncode == 0, response_text

	# Check output
	for ref, out in ((OUTPUTS_DIR/seq_type/'failed_lariat_alignments.tsv', tmp_path/f'{prefix}failed_lariat_alignments.tsv'),
					(OUTPUTS_DIR/seq_type/'lariat_reads.tsv', tmp_path/f'{prefix}lariat_reads.tsv')):
		ref_lines, out_lines = load_file_lines(ref, out)
		if ref_lines == out_lines:
			continue

		# File contents differ, decide how to report
		if vscode_available():
			vscode_compare(ref, out)
			pytest.fail(f'Output file differs from expected output: {out.name}')
		else:
			assert ref_lines == out_lines
