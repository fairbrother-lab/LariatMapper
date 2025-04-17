import collections
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
TEST_DIR = PACKAGE_DIR/'tests'/'bp_correction'
Method_Combo = collections.namedtuple('Method_Combo', ['method_arg', 'path_arg', 'output'])

# =============================================================================#
#                                   Tests                                      #
# =============================================================================#
# Required args
@pytest.mark.parametrize('method_combo',
					[Method_Combo(method_arg='--method Model-based', 
									path_arg=f"--model_path {TEST_DIR/'inputs'/'gencode_v44_U2_U12_BP_pred_prob_tx.rds'}", 
									output=pathlib.Path(f"{TEST_DIR}/outputs/lariat_reads_model.tsv")),
					Method_Combo(method_arg='--method PWM',
									path_arg=f"--PWM_path {TEST_DIR/'inputs'/'gencode_v44_U2_pwm.rds'}",
									output=pathlib.Path(f"{TEST_DIR}/outputs/lariat_reads_pwm_U2.tsv")),
					Method_Combo(method_arg='--method PWM',
									path_arg=f"--PWM_path {TEST_DIR/'inputs'/'gencode_v44_U2_pwm.rds'},{TEST_DIR/'inputs'/'gencode_v44_U12_pwm.rds'}",
									output=pathlib.Path(f"{TEST_DIR}/outputs/lariat_reads_pwm_U2_U12.tsv")),
					]
)
# Optional args
@pytest.mark.parametrize('log_level',
					[None, 
	  				'--log_level DEBUG',
					'--log_level ERROR',]
)
def test_bp_correction(method_combo, log_level, tmp_path):	
	subprocess.run(f"gunzip --keep --stdout {TEST_DIR/'inputs'/'hg38.chr22.fa.gz'} > {tmp_path}/hg38.chr22.fa", shell=True)
	subprocess.run(f"cp {TEST_DIR/'inputs'/'hg38.chr22.fa.fai'} {tmp_path}/hg38.chr22.fa.fai", shell=True)

	command = f"Rscript {SCRIPTS_DIR/'bp_correction_wrapper.R'}" +\
				f" --input {TEST_DIR/'inputs'/'lariat_reads.tsv'}" +\
				f" --ref_fasta {tmp_path/'hg38.chr22.fa'}" +\
				f" --file {SCRIPTS_DIR/'bp_correction.R'}" +\
				' ' + method_combo.method_arg +\
				' ' + method_combo.path_arg +\
				f" --output_base {tmp_path}/"
	for optional_arg in (log_level,):
		if optional_arg is not None:
			command += ' ' + optional_arg
	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	assert response.returncode == 0, response_text

	# Check output
	ref_and_out_files = (
		(method_combo.output, tmp_path/'lariat_reads.tsv'),
	)
	test_utils.confirm_outputs_match_references(ref_and_out_files, response_text)
