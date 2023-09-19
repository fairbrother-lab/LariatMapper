
from sys import argv
from logging import Logger
from os.path import join, getsize, isdir, isfile
from os import listdir, mkdir, chmod
import stat
from shutil import move
from subprocess import run
from csv import reader
from dataclasses import dataclass
from logger import setup_logger
from argparse import ArgumentParser

# =====================================================================#
#                          RunData Class                              #
# =====================================================================#
@dataclass
class RunData:

    ATTRS = ['fastq_dir', 'scripts_dir', 'log_dir', 'output_dir', 'top_output_dir', 'results_path', 'num_cpus',
             'ref_b2index', 'ref_fasta', 'ref_gtf', 'ref_5p_fasta', 'ref_3p_b2index', 'ref_3p_lengths',
             'ref_introns', 'ref_repeatmasker', 'sample_info', 'sample_categories']

    fastq_dir: str
    read_one_file: str
    read_two_file: str
    output_base_name: str
    results_path: str
    num_cpus: int
    ref_b2index: str
    ref_fasta: str
    ############################
    # can provide these 5 in package
    ref_gtf: str
    ref_5p_fasta: str
    ref_3p_b2index: str
    ref_3p_lengths: str
    ref_introns: str
    ############################
    ref_repeatmasker: str

    top_output_dir: str     # default: larmap_out
    scripts_dir: str        # default: {top_output_dir}/scripts
    log_dir: str            # default: {top_output_dir}/logs
    output_dir: str         # default: {top_output_dir}/output
    sample_categories: list # default: [read_1_file, read_2_file, output_base_name]
    sample_info: list       # default: dict {read_1_file: , read_2_file: , output_base_name: }


    def __init__(self, 
        fastq_dir: str, 
        read_one_file: str,
        read_two_file: str,
        output_base_name: str,
        results_path: str,
        num_cpus: int,
        ref_b2index: str,
        ref_fasta: str,
        ref_gtf: str,
        ref_5p_fasta: str,
        ref_3p_b2index: str,
        ref_3p_lengths: str,
        ref_introns: str,
        ref_repeatmasker: str,
        log: Logger) -> None:

        # custom params
        setattr(self, 'fastq_dir', fastq_dir)
        setattr(self, 'read_one_file', read_one_file)
        setattr(self, 'read_two_file', read_two_file)
        setattr(self, 'output_base_name', output_base_name)
        setattr(self, 'num_cpus', num_cpus)
        setattr(self, 'ref_b2index', ref_b2index)
        setattr(self, 'ref_fasta', ref_fasta)
        setattr(self, 'ref_gtf', ref_gtf)
        setattr(self, 'ref_5p_fasta', ref_5p_fasta)
        setattr(self, 'ref_3p_b2index', ref_3p_b2index)
        setattr(self, 'ref_3p_lengths', ref_3p_lengths)
        setattr(self, 'ref_introns', ref_introns)
        setattr(self, 'ref_repeatmasker', ref_repeatmasker)
        # default params
        setattr(self, 'results_path', 'larmap_out/' + results_path)
        setattr(self, 'top_output_dir', 'larmap_out')
        setattr(self, 'scripts_dir', 'larmap_out' + '/scripts')
        setattr(self, 'log_dir', 'larmap_out' + '/logs')
        setattr(self, 'output_dir', 'larmap_out' + '/output')
        setattr(self, 'sample_categories', ['read_one_file', 'read_two_file', 'output_base_name'])
        setattr(self, 'sample_info', [{'read_one_file': read_one_file, 'read_two_file': read_two_file, 'output_base_name': output_base_name}])

        log.info(f'{self.sample_info[0]}')

    def __repr__(self):
        rep = ""
        for attr, val in vars(self).items():
            if not attr.startswith('__'):
                rep += f'{attr}: {val}\t\t({type(val)})\n\t\t'
        return rep

    def validate(self) -> None:
        '''
        Double-checks the provided info file for common errors
        '''
        # Do all the attributes have values?
        for attr in RunData.ATTRS:
            if vars(self)[attr] in ('', None):
                raise ValueError(
                    f"Invalid value ({vars(self)[attr]}) for {attr} attribute")

        # Does the fastq directory exist?
        if not isdir(self.fastq_dir):
            raise ValueError(
                f'The directory for the raw fastq files, "{self.fastq_dir}", was not found.')

        # Is the num_cpus setting valid?
        if not self.num_cpus.isdigit() or int(self.num_cpus) < 1:
            raise ValueError(
                f'num_cpus value {self.num_cpus} is not a positive integer')

        # Does the bowtie2 index exist?
        for suffix in ('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2'):
            if not isfile(self.ref_b2index + suffix):
                raise ValueError(
                    f'The Bowtie2 index file "{self.ref_b2index + suffix}" does not exist')

        # Do the reference files exist?
        for ref in (self.ref_b2index + '.1.bt2',
                    self.ref_b2index + '.rev.2.bt2',
                    self.ref_fasta,
                    self.ref_gtf,
                    self.ref_5p_fasta,
                    self.ref_3p_b2index + '.1.bt2',
                    self.ref_3p_b2index + '.rev.2.bt2',
                    self.ref_3p_lengths,
                    self.ref_introns,
                    self.ref_repeatmasker):
            if not isfile(ref):
                raise ValueError(f'The reference file "{ref}" does not exist.')

        # Are the values for each sample file's variables valid?
        for sample in self.sample_info:
            read_one_path = join(self.fastq_dir, sample['read_one_file'])
            read_two_path = join(self.fastq_dir, sample['read_two_file'])
            if not isfile(read_one_path):
                raise ValueError(
                    f'The read one fastq file "{read_one_path}" was not found.')
            if not isfile(read_two_path):
                raise ValueError(
                    f'The read two fastq file "{read_two_path}" was not found.')
            for bad_char in ('split', 'map', '$', '%', '#', "'", '"', '/', '\\', '{', '}', ','):
                if bad_char in sample['read_one_file']:
                    raise ValueError(
                        f'{sample["read_one_file"]} cannot contain the character " {bad_char} "')
                if bad_char in sample['read_two_file']:
                    raise ValueError(
                        f'{sample["read_two_file"]} cannot contain the character " {bad_char} "')
                if bad_char in sample['output_base_name']:
                    raise ValueError(
                        f'The output base name for the file {sample["output_base_name"]} cannot contain the character " {bad_char} "')
