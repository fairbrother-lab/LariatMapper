
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

# =============================================================================#
#                                  Constants                                  #
# =============================================================================#

ENV_DIR = 'larmap_env'
SCRIPTS_DIR = 'scripts'

TIMESTAMP_SCRIPT = join(SCRIPTS_DIR, 'timestamp.sh')
MAP_LARIATS_SCRIPT = join(SCRIPTS_DIR, 'map_lariats.sh')
FILTER_5p_SCRIPT = join(SCRIPTS_DIR, 'filter_fivep_alignments.py')
FILTER_3p_SCRIPT = join(SCRIPTS_DIR, 'filter_threep_alignments.py')

TRANSLATE_TO_LOG = {'fastq_dir': 'Fastq Files Directory',
                    'scripts_dir': 'Scripts Directory',
                    'log_dir': 'Log Directory',
                    'output_dir': 'Mapping Output Directory',
                    'results_path': 'Final Results Path',
                    'num_cpus': 'Number of CPUs',
                    'ref_b2index': 'Bowtie2 Genome Index Base',
                    'ref_fasta': 'Genome FASTA File',
                    'ref_gtf': 'Genome GTF Annotation File',
                    'ref_5p_fasta': '5p Splice Sites File',
                    'ref_3p_b2index': '3p Splice Sites Bowtie2 Index Base',
                    'ref_3p_lengths': '3p Sequence Lengths File',
                    'ref_introns': 'Introns BED File',
                    'ref_repeatmasker': 'RepeatMasker BED File'
                    }
RMAP = {'one': 'R1', 'two': 'R2'}


# =====================================================================#
#                          RunData Class                              #
# =====================================================================#
@dataclass
class RunData:
    ATTRS = ['fastq_dir', 'scripts_dir', 'log_dir', 'output_dir', 'results_path', 'num_cpus',
             'ref_b2index', 'ref_fasta', 'ref_gtf', 'ref_5p_fasta', 'ref_3p_b2index', 'ref_3p_lengths',
             'ref_introns', 'ref_repeatmasker', 'sample_info', 'sample_categories']

    fastq_dir: str
    scripts_dir: str
    log_dir: str
    output_dir: str
    results_path: str
    num_cpus: int
    ref_b2index: str
    ref_fasta: str
    ref_gtf: str
    ref_5p_fasta: str
    ref_3p_b2index: str
    ref_3p_lengths: str
    ref_introns: str
    ref_repeatmasker: str
    sample_info: list
    sample_categories: list

    def __init__(self, info_file, log) -> None:
        with open(info_file, 'r') as info:
            info_reader = reader(info, delimiter='\t')

            log.info('Run parameters from info file:')
            for row in info_reader:

                # Remove any empty values from extra tabs
                row = [item for item in row if item != '']

                if row == [] or row[0].startswith('#'):
                    continue
                if row[0] in TRANSLATE_TO_LOG.keys():
                    variable, value = row[:2]
                    log.info(f'{TRANSLATE_TO_LOG[variable]}: {value}')
                    setattr(self, variable, value)

                if row[0] == 'read_one_file':
                    break

            self.sample_categories = row
            self.sample_info = []
            log.info(
                f'{len(self.sample_categories)} categories pulled from info file: {self.sample_categories}')

            for row in info_reader:
                print(row)
                if row[0].startswith('#'):
                    continue
                log.info(f'Pulled row: {row}')

                sample = {}
                for i in range(len(self.sample_categories)):
                    sample[self.sample_categories[i]] = row[i]
                self.sample_info.append(sample)

            self.sample_info, self.sample_categories = tuple(
                self.sample_info), tuple(self.sample_categories)

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


# =============================================================================#
#                                  Functions                                  #
# =============================================================================#

def write_mapping_scripts(sample: dict, args: dict, run_data: RunData, log: Logger) -> None:
    '''
    For each sample write two lariat mapping scripts which process the sample's read one and read two files
    '''
    for read_num in ['one', 'two']:

        script_name = f'larmap_{sample["output_base_name"]}_{RMAP[read_num]}'
        logfile_path = join(run_data.log_dir, 'child_logs',
                            f'{script_name}.out')
        script_path = join(run_data.scripts_dir, f'{script_name}.sh')

        read_file_name = sample[f'read_{read_num}_file']
        read_file_path = join(run_data.fastq_dir, read_file_name)
        read_output_name = f'{sample["output_base_name"]}_{RMAP[read_num]}'
        read_output_dir = join(run_data.output_dir, read_output_name)
        read_output_base = join(read_output_dir, read_output_name)

        log.info(f'Creating mapping script for {read_file_name}')

        with open(script_path, 'w') as script_file:
            script_file.write('#!/bin/bash -l\n\n')
            # SECONDS environmental variable to track time elapsed and its subsequent referals
            script_file.write('SECONDS=0 \n')
            # Infer the scripts log file path from the job id
            script_file.write(f'LOGFILE_PATH={logfile_path} \n')
            script_file.write('printf "Log file path = $LOGFILE_PATH\\n"\n')
            # Source your conda installation
            script_file.write('source ~/anaconda3/etc/profile.d/conda.sh \n')
            # Activate the conda env
            script_file.write(f'if [[ "$(conda config --show envs_dirs)" != *"- {ENV_DIR}"* ]]; then \n'
                              f'\tconda config --append envs_dirs {ENV_DIR} \n'
                              'fi \n')
            script_file.write(f'conda activate {ENV_DIR} \n')
            script_file.write(
                f'printf "conda environment = $CONDA_PREFIX\\n"\n\n')
            # Create a results directory for the read file in the top results directory
            script_file.write(f'mkdir -p {read_output_dir} \n\n')
            # Call the map_lariat.sh script to perform the lariat mapping on the read file
            script_file.write(
                f'printf "Beginning {sample["output_base_name"]} read {read_num} lariat mapping...\\n"\n')
            mapping_call = f'{MAP_LARIATS_SCRIPT} \\\n'\
                f'{read_file_path} \\\n'\
                f'{read_output_base} \\\n'\
                f'{run_data.num_cpus} \\\n'\
                f'{run_data.ref_b2index} \\\n'\
                f'{run_data.ref_fasta} \\\n'\
                f'{run_data.ref_gtf} \\\n'\
                f'{run_data.ref_5p_fasta} \\\n'\
                f'{run_data.ref_3p_b2index} \\\n'\
                f'{run_data.ref_3p_lengths} \\\n'\
                f'{FILTER_5p_SCRIPT} \\\n'\
                f'{FILTER_3p_SCRIPT} \\\n'\
                f'{TIMESTAMP_SCRIPT} \n'
            script_file.write(mapping_call)
            script_file.write(
                f'printf "Mapping of {sample["output_base_name"]} read {read_num} file complete\\n"\n\n')

        # Set the scripts permissions so everyone group can rwx
        chmod(script_path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)


def write_bash_all(run_data: RunData, log: Logger) -> None:
    '''
    Write a bash_all.sh script in run_data.scripts_dir that sequentially processes all the read lariat mapping scripts
    '''
    bash_all_path = join(run_data.scripts_dir, 'bash_all.sh')
    with open(bash_all_path, 'w') as script_file:
        script_file.write('#!/bin/bash\n\n')
        for sample in run_data.sample_info:
            for read_suffix in ['R1', 'R2']:
                script_name = f'larmap_{sample["output_base_name"]}_{read_suffix}'
                script_path = join(run_data.scripts_dir, f'{script_name}.sh')
                logfile_path = join(
                    run_data.log_dir, 'child_logs', f'{script_name}.out')
                script_file.write(f'bash {script_path} &> {logfile_path}\n')

    # Set the child's permissions so everyone group can rwx
    chmod(bash_all_path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)

    log.info(f'Creating top level mapping script at {bash_all_path}')


# =============================================================================#
#                                    Main                                     #
# =============================================================================#
def main():

    log = setup_logger('top_log.out')

    # Parse info file argument and load information into RunData object
    parser = ArgumentParser(
        description='Top-level script for mapping lariats reads from RNA-seq data')
    parser.add_argument('info_file',
                        help='The path to the information file that contains the settings and fastq file '
                        'details needed for the lariat mapping process')
    args = vars(parser.parse_args())
    run_data = RunData(args['info_file'], log)
    # make sure the run attributes have appropriate values
    run_data.validate()
    log.info('Run data has passed inspection')

    # Prepare script, log and output directories needed for the mapping run
    dir_types = ['script', 'log', 'child logs', 'output']
    dir_list = [run_data.scripts_dir, run_data.log_dir, join(
        run_data.log_dir, 'child_logs'), run_data.output_dir]
    for dtype, d in zip(dir_types, dir_list):
        if not isdir(d):
            log.info(f'Creating {dtype} directory at {d}')
            mkdir(d)
    log.info('Directories have been prepared.')
    # move the log file to the run_data.log_dir directory
    move('top_log.out', join(run_data.log_dir, 'top_log.out'))

    print(run_data.sample_info)

    # For each sample (one line in fastq files table) write two scripts which perform the lariat mapping
    # of the sample's read one and read two files
    for sample in run_data.sample_info:
        write_mapping_scripts(sample, args, run_data, log)

    write_bash_all(run_data, log)


if __name__ == '__main__':
    main()
