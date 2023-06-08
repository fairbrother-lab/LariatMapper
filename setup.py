from setuptools import setup, find_packages
import subprocess

setup(
    name='larmap',
    version='0.0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'intervaltree',
        'pyfaidx'
    ],
    python_requires='>=3.7.6',
    author='nealneal',
    description='lariat mapping (test package)',
    url='https://github.com/nyin01/lariat_mapping_test',
)

subprocess.run(["python", "set_alias.py"])