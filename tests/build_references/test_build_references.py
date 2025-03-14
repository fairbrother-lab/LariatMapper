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
# Required args
@pytest.mark.parametrize('fasta',
					# [TEST_DIR/'inputs'/'genome.fa',] 
					[TEST_DIR/'inputs'/'genome.fa', 
					TEST_DIR/'inputs'/'genome.fa.gz']
					)
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
	  				'--threads 4'])
@pytest.mark.parametrize('copy',
					[None, 
	  				'--copy'])
@pytest.mark.parametrize('verbosity',
					[None, 
					'--debug'])
def test_build_references(fasta, anno, repeatmasker_bed, threads, copy, verbosity, tmp_path):
	command = f"python {PACKAGE_DIR/'build_references.py'} --skip_r" +\
			f" --hisat2_index {TEST_DIR/'inputs'/'hisat2_index'}" +\
			f" --genome_fasta {fasta}" +\
			f" {anno} -o {tmp_path}"
	for optional_arg in (repeatmasker_bed, threads, copy, verbosity):
		if optional_arg is not None:
			command += ' ' + optional_arg

	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	assert response.returncode == 0, response_text

	fasta_ext = '.fa' if fasta.suffix == '.fa' else '.fa.gz'

	# Check output
	ref_and_out_files = (
					(TEST_DIR/'outputs'/f'genome{fasta_ext}.fai', tmp_path/f'genome{fasta_ext}.fai'),
					(TEST_DIR/'outputs'/'exons.tsv.gz', tmp_path/'exons.tsv.gz'),
					(TEST_DIR/'outputs'/'introns.tsv.gz', tmp_path/'introns.tsv.gz'),
					(TEST_DIR/'outputs'/'fivep_sites.fa.gz', tmp_path/'fivep_sites.fa.gz'),
					(TEST_DIR/'outputs'/'fivep_sites.fa.gz.fai', tmp_path/'fivep_sites.fa.gz.fai'),
	)
	test_utils.confirm_outputs_match_references(ref_and_out_files, response_text)


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