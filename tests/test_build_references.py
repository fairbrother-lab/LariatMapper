import pathlib
import pytest
import subprocess



from . import functions



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
PACKAGE_DIR = pathlib.Path(__file__).parent.parent.resolve()
SCRIPT = PACKAGE_DIR / 'build_references.py'
INPUTS_DIR = PACKAGE_DIR / 'tests' / 'inputs'



# =============================================================================#
#                                   Tests                                      #
# =============================================================================#
# Required args + the transcript attribute and gene attribute args since their values depend on the genome gtf
@pytest.mark.parametrize('req_args',
						[f"-f {INPUTS_DIR / 'genome.fa'} -a {INPUTS_DIR / 'anno.gtf'} -i {INPUTS_DIR / 'genome'} -x tid -g gid", 
						f"--genome_fasta {INPUTS_DIR / 'genome.fa'} --genome_anno {INPUTS_DIR / 'anno.gtf'} --hisat2_index {INPUTS_DIR / 'genome'} --transcript_attribute tid --gene_attribute gid",],
						 )
# Optional args
@pytest.mark.parametrize('repeatmasker_bed',
						 [None, f"-r {INPUTS_DIR / 'repeatmasker.bed'}", f"--repeatmasker_bed {INPUTS_DIR / 'repeatmasker.bed'}",
						])
@pytest.mark.parametrize('threads',
						 [None, '-t 4', '--threads 4'])
def test_run(req_args, repeatmasker_bed, threads, tmp_path):
	command = f"python {SCRIPT} {req_args} -o {tmp_path}"
	for optional_arg in (repeatmasker_bed, threads,):
		if optional_arg is not None:
			command = f'{command} {optional_arg}'

	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	response_code = response.returncode
	assert response_code  == 0, response_text

	functions.files_equal(INPUTS_DIR / 'exons.tsv', tmp_path / 'exons.tsv.gz')
	functions.files_equal(INPUTS_DIR / 'introns.tsv', tmp_path / 'introns.tsv.gz')
	functions.files_equal(INPUTS_DIR / 'fivep.fa', tmp_path / 'fivep_sites.fa')
