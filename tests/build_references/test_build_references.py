# import gzip
# import os
# import pathlib
# import shutil
# import subprocess
# import sys

# import pytest

# sys.path.append(str(pathlib.Path(__file__).parent.parent.resolve()))
# import testing_tools as test_funcs


# # =============================================================================#
# #                                  Globals                                     #
# # =============================================================================#
# PACKAGE_DIR = pathlib.Path(__file__).parent.parent.parent.resolve()
# BUILD_REFERENCES_DIR = PACKAGE_DIR/'tests'/'build_references'
# INPUTS_DIR = PACKAGE_DIR/'tests'/'inputs'
# OUTPUTS_DIR = PACKAGE_DIR/'tests'/'outputs'



# # =============================================================================#
# #                                   Tests                                      #
# # =============================================================================#
# # Required args + the transcript attribute and gene attribute args since their values depend on the genome gtf
# @pytest.mark.parametrize('base',
# 					[f"-f {BUILD_REFERENCES_DIR/'inputs'/'genome.fa'} -i {BUILD_REFERENCES_DIR/'inputs'/'hisat2_index'} -x tid -g gid", 
# 					f"--genome_fasta {BUILD_REFERENCES_DIR/'inputs'/'genome.fa'}  --hisat2_index {BUILD_REFERENCES_DIR/'inputs'/'hisat2_index'} --transcript_attribute tid --gene_attribute gid",],)
# @pytest.mark.parametrize('anno',
# 					[f"-a {BUILD_REFERENCES_DIR/'inputs'/'anno.gtf'}", 
# 					f"-a {BUILD_REFERENCES_DIR/'inputs'/'anno.gff'}", 
# 					f"--genome_anno {BUILD_REFERENCES_DIR/'inputs'/'anno.gtf'}"])
# # Optional args
# @pytest.mark.parametrize('repeatmasker_bed',
# 					[None, 
# 					f"-r {BUILD_REFERENCES_DIR/'inputs'/'repeatmasker.bed'}", 
# 					f"--repeatmasker_bed {BUILD_REFERENCES_DIR/'inputs'/'repeatmasker.bed'}"])
# @pytest.mark.parametrize('threads',
# 					[None, '-t 4', '--threads 4'])
# @pytest.mark.parametrize('copy',
# 					[None, '--copy'])
# def test_build_references(base, anno, repeatmasker_bed, threads, copy, tmp_path):
# 	command = f"python {PACKAGE_DIR/'build_references.py'} {base} {anno} -o {tmp_path}"
# 	for optional_arg in (repeatmasker_bed, threads, copy):
# 		if optional_arg is not None:
# 			command = f'{command} {optional_arg}'

# 	response = subprocess.run(command, shell=True, capture_output=True, text=True)
# 	response_text = '\n' + response.stdout + response.stderr
# 	assert response.returncode == 0, response_text

# 	# Check output
# 	for ref, out in (
# 					# (BUILD_REFERENCES_DIR/'outputs'/'exons.tsv', tmp_path/'exons.tsv.gz'),
# 					(BUILD_REFERENCES_DIR/'outputs'/'introns.tsv', tmp_path/'introns.tsv.gz'),
# 					(BUILD_REFERENCES_DIR/'outputs'/'fivep_sites.fa', tmp_path/'fivep_sites.fa')
# 					):
# 		ref_lines, out_lines = test_funcs.load_file_lines(ref, out)
# 		if ref_lines == out_lines:
# 			continue

# 		# If file contents differ, decide how to report
# 		if test_funcs.vscode_available():
# 			test_funcs.vscode_compare(ref, out)
# 			pytest.fail(f'Output file differs from expected output: {out.name}')
# 		else:
# 			assert ref_lines == out_lines
# 		pass