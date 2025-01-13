import hashlib
import pathlib
import subprocess
import sys
import tempfile

import pytest

sys.path.append(str(pathlib.Path(__file__).parent.parent.resolve()))
import test_utils


# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
PACKAGE_DIR = pathlib.Path(__file__).parent.parent.parent.resolve()
TEST_DIR = PACKAGE_DIR/'tests'/'build_references'

# =============================================================================#
#                                   Tests                                      #
# =============================================================================#
# Required anno arg + the transcript attribute and gene attribute args since their values depend on the anno file
@pytest.mark.parametrize('anno',
					[f"--genome_anno {TEST_DIR/'inputs'/'anno.gtf'} --t_attr tid --g_attr gid", 
					f"--genome_anno {TEST_DIR/'inputs'/'anno.gff'} --t_attr tid --g_attr gid"])
# Optional args
@pytest.mark.parametrize('repeatmasker_bed',
					[None, 
	  				f"--repeatmasker_bed {TEST_DIR/'inputs'/'repeatmasker.bed'}"])
@pytest.mark.parametrize('threads',
					[None, 
	  				'--threads 1',
	  				'--threads 4'])
@pytest.mark.parametrize('copy',
					[None, 
	  				'--copy'])
@pytest.mark.parametrize('verbosity',
					[None, 
	  				'--quiet',
					'--debug'])
def test_build_references(anno, repeatmasker_bed, threads, copy, verbosity, tmp_path):
	command = f"python {PACKAGE_DIR/'build_references.py'} --skip_r" +\
			f" --genome_fasta {TEST_DIR/'inputs'/'genome.fa'}" +\
			f" --hisat2_index {TEST_DIR/'inputs'/'hisat2_index'}" +\
			f" {anno} -o {tmp_path}"
	for optional_arg in (repeatmasker_bed, threads, copy):
		if optional_arg is not None:
			command += ' ' + optional_arg

	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	assert response.returncode == 0, response_text

	# Check output
	for ref, out in (
					(TEST_DIR/'outputs'/'introns.tsv', tmp_path/'introns.tsv.gz'),
					(TEST_DIR/'outputs'/'fivep_sites.fa', tmp_path/'fivep_sites.fa')
					):
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



# Optional args
@pytest.mark.parametrize('verbosity',
					[None, 
	  				'--quiet',
					'--debug'])
def test_build_R_refs(verbosity, tmp_path):
	command = f"Rscript {PACKAGE_DIR}/scripts/build_R_refs.R" +\
			f" --anno {TEST_DIR/'inputs'/'hg38.gencode.v44.sample.gtf'}" +\
			f" --g_attr gene_id --t_attr transcript_id --output {tmp_path}"
	for optional_arg in (verbosity,):
		if optional_arg is not None:
			command += ' ' + optional_arg
	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	if response.returncode != 0:
		pytest.fail(response_text)

	# Check if objects match ref objects in R
	# with open('./tempy.R', 'w') as tmp:
	with tempfile.NamedTemporaryFile(suffix='.R', mode='w') as tmp:

		for feature in ('gene', 'exon', 'intron'):
			tmp.write(f"""
				ref_obj = readRDS('{TEST_DIR/'outputs'}/{feature}_gr.rds')
				test_obj = readRDS('{tmp_path}/{feature}_gr.rds')
				if (!identical(ref_obj, test_obj)){{
					stop('The {feature} object did not match the reference {feature} object')
				}}
			""")

		# response = subprocess.run(['Rscript', './tempy.R'], capture_output=True, text=True)
		response = subprocess.run(['Rscript', tmp.name], capture_output=True, text=True)
		response_text = '\n' + response.stdout + response.stderr
		if response.returncode != 0:
			pytest.fail(response_text)