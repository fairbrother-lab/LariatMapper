
import sys
from os.path import join

if __name__ == '__main__':

	sample_info_file, merged_out_file = sys.argv[1:]

	with open(sample_info_file) as in_file, open(merged_out_file, 'w') as out_file:
		sample_categories = '\t'.join(in_file.readline().strip().split('\t')[:-2])
		print_header = True
		for line in in_file:
			sinfo = line.strip().split('\t')
			sample_data = '\t'.join(sinfo[:-2])
			sample_output_base = '_'.join(sinfo[:-2])
			read_one_file, read_two_file = sinfo[-2:]
			sample_dir = f'{sample_output_base}_lariat_mapping'
			sample_file = f'{sample_output_base}_lariat_reads.txt'
			with open(join(sample_dir, sample_file)) as sample_in:
				header_line = sample_in.readline()
				if print_header:
					out_file.write(f'{sample_categories}\t{header_line}')
					print_header = False
				for lar_line in sample_in:
					out_file.write(f'{sample_data}\t{lar_line}')


