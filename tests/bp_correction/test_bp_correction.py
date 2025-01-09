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
TEST_DIR = PACKAGE_DIR/'tests'/'bp_correction'
DATA_DIR = PACKAGE_DIR/'data'
Method_Combo = collections.namedtuple('Method_Combo', ['method_arg', 'path_arg', 'output'])

# =============================================================================#
#                                   Tests                                      #
# =============================================================================#
# Required args
@pytest.mark.parametrize('method_combo',
					[Method_Combo(method_arg='--method Model-based', 
									path_arg=f'--model_path {DATA_DIR}/gencode_v44_U2_U12_BP_pred_prob_tx.rds', 
									output=f'{TEST_DIR}/outputs/lariat_reads_model.tsv'),
					Method_Combo(method_arg='--method PWM',
									path_arg=f'--PWM_path {DATA_DIR}/U2_pwm.rds',
									output=f'{TEST_DIR}/outputs/lariat_reads_pwm_U2.tsv'),
					Method_Combo(method_arg='--method PWM',
									path_arg=f'--PWM_path {DATA_DIR}/U2_pwm.rds,{DATA_DIR}/U12_pwm.rds',
									output=f'{TEST_DIR}/outputs/lariat_reads_pwm_U2_U12.tsv'),
					]
)
# Optional args
@pytest.mark.parametrize('log_level',
					[None, 
	  				'--log_level DEBUG',
					'--log_level ERROR',]
)
def test_bp_correction(method_combo, log_level, tmp_path):
# def test_bp_correction(method_combo, verbosity, tmp_path):
	command = f"Rscript {PACKAGE_DIR/'scripts'/'bp_correction_wrapper.R'}" +\
				f" --input {TEST_DIR/'inputs'/'lariat_reads.tsv'}" +\
				f" --ref_fasta {TEST_DIR/'inputs'/'hg38.demo.fa'}" +\
				f" --file {PACKAGE_DIR/'scripts'/'bp_correction.R'}" +\
				' ' + method_combo.method_arg +\
				' ' + method_combo.path_arg +\
				f" --output_base {tmp_path}/"
	for optional_arg in (log_level,):
		if optional_arg is not None:
			command += ' ' + optional_arg
	print(command)
	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response_text = '\n' + response.stdout + response.stderr
	assert response.returncode == 0, response_text

	# Check output
	ref = method_combo.output
	out = tmp_path/'lariat_reads.tsv'
	ref_lines, out_lines = test_utils.load_file_lines(ref, out)
	if ref_lines == out_lines:
		return
	
	# If the lines don't match, report the differences
	# If in vscode, open the files for comparison
	if test_utils.vscode_available():
		test_utils.vscode_compare_sorted(ref, out)
		
	assert ref_lines == out_lines
