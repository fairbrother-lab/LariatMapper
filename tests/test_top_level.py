import pathlib
import pytest
import subprocess
import os
import shutil

# hisat2-build /Users/trumanmooney/Library/CloudStorage/OneDrive-BrownUniversity/Documents/Projects/Lariat_mapping/LariatMapper/tests/inputs/genome_1.fa /Users/trumanmooney/Library/CloudStorage/OneDrive-BrownUniversity/Documents/Projects/Lariat_mapping/LariatMapper/tests/inputs/genome_1  
PACKAGE_DIR = pathlib.Path(__file__).parent.parent.resolve()
INPUTS_DIR = PACKAGE_DIR / 'tests' / 'inputs'
BUILD_REF_SCRIPT = PACKAGE_DIR / 'build_references.py'



@pytest.mark.datafiles(INPUTS_DIR / 'genome_1.fa',
					INPUTS_DIR / 'anno_1.gtf',
					INPUTS_DIR / 'genome_1.1.ht2',
					INPUTS_DIR / 'genome_1.2.ht2',
					INPUTS_DIR / 'genome_1.3.ht2',
					INPUTS_DIR / 'genome_1.4.ht2',
					INPUTS_DIR / 'genome_1.5.ht2',
					INPUTS_DIR / 'genome_1.6.ht2',
					INPUTS_DIR / 'genome_1.7.ht2',
					INPUTS_DIR / 'genome_1.8.ht2',
					   )
def test_build_references(datafiles, tmp_path):
	genome_fasta = datafiles / 'genome_1.fa'
	genome_anno = datafiles / 'anno_1.gtf'
	hisat2_index = datafiles / 'genome_1'
	out_dir = tmp_path / 'out'
	command = f'python -u {BUILD_REF_SCRIPT} -f {genome_fasta} -a {genome_anno} -i {hisat2_index} -o {out_dir} -x tid -g gid --debug'
	response = subprocess.run(command, shell=True, capture_output=True, text=True)
	response = response.stdout + response.stderr
	shutil.copytree(out_dir, '/Users/trumanmooney/Library/CloudStorage/OneDrive-BrownUniversity/Documents/Projects/Lariat_mapping/LariatMapper/tests/out')
	assert response == ''
	# assert pathlib.Path.is_file(datafiles / 'genome_1.fa')
	assert pathlib.Path.is_file(out_dir / 'genome_1.fa')
	assert pathlib.Path.is_file(out_dir / 'genome_1.1.ht2')


