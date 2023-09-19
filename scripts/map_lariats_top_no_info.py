
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
from rundata import RunData

# =============================================================================#
#                                   Constants                                  #
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

# =============================================================================#
#                                  Functions                                  #
# =============================================================================#


def write_mapping_scripts(sample: dict, run_data: RunData, log: Logger) -> None:
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
        # add parallel
        script_file.write("(trap 'kill 0' SIGINT; ")

        proc_length = len(run_data.sample_info) * 2
        count = 0

        for sample in run_data.sample_info:
            for read_suffix in ['R1', 'R2']:
                script_name = f'larmap_{sample["output_base_name"]}_{read_suffix}'
                script_path = join(run_data.scripts_dir, f'{script_name}.sh')
                logfile_path = join(
                    run_data.log_dir, 'child_logs', f'{script_name}.out')
                # script_file.write(f'bash {script_path} &> {logfile_path}\n')
                count += 1
                if (count < proc_length):
                    script_file.write(
                        f'bash {script_path} &> {logfile_path} & ')
                else:
                    script_file.write(
                        f'bash {script_path} &> {logfile_path} )')

    # Set the child's permissions so everyone group can rwx
    chmod(bash_all_path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)

    log.info(f'Creating top level mapping script at {bash_all_path}')


# =============================================================================#
#                                    Main                                     #
# =============================================================================#
def main():

    log = setup_logger('top_log.out')
    run_data = RunData(*argv[1:], log)
    # make sure the run attributes have appropriate values
    run_data.validate()
    log.info('Run data has passed inspection')

    # NEAL: create top output dir if doesn't exist
    if not isdir(run_data.top_output_dir):
        mkdir(run_data.top_output_dir)

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
        write_mapping_scripts(sample, run_data, log)

    write_bash_all(run_data, log)


if __name__ == '__main__':
    main()
