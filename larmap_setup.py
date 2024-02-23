
import sys
from os.path import join

if __name__ == '__main__':

	run_info_file, sample_info_file = sys.argv[1:]

	run_info = {}
	with open(run_info_file) as in_file:
		for line in in_file:
			info = line.strip()
			if info:
				run_key, run_value = info.split('\t')
				run_info[run_key] = run_value

	with open(sample_info_file) as in_file:
		next(in_file)
		for line in in_file:
			sinfo = line.strip().split('\t')
			sample_output_base = '_'.join(sinfo[:-1])
			read_file = join(run_info["fastq_dir"], sinfo[-1])
			script_path = sample_output_base + '_larmap.sh'
			with open(script_path, 'w') as out_file:
				out_file.write('#!/bin/bash\n\n')
				out_file.write(f'{run_info["scripts_dir"]}/larmap_run.sh \\\n')
				out_file.write(f'-r {read_file} \\\n')
				out_file.write(f'-o {run_info["output_dir"]} \\\n')
				out_file.write(f'-e {sample_output_base} \\\n')
				out_file.write(f'-c {run_info["num_cpus"]} \\\n')
				out_file.write(f'-f {run_info["ref_fasta"]} \\\n')
				out_file.write(f'-i {run_info["ref_b2index"]} \\\n')
				out_file.write(f'-g {run_info["ref_gtf"]} \\\n')
				out_file.write(f'-n {run_info["ref_introns"]} \\\n')
				out_file.write(f'-5 {run_info["ref_5p_fasta"]} \\\n')
				out_file.write(f'-u {run_info["ref_5p_upstream"]} \\\n')
				out_file.write(f'-3 {run_info["ref_3p_b2index"]} \\\n')
				out_file.write(f'-l {run_info["ref_3p_lengths"]} \\\n')
				out_file.write(f'-m {run_info["ref_repeatmasker"]}\n')

