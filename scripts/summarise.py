import os
import sys
import json

import pandas as pd
import pyfaidx
import pysam

import functions



#=============================================================================#
#                                   Globals                                   #
#=============================================================================#
# In files
ARGS_FILE = "{}args.json"
OUTPUT_BAM_FILE = "{}output.bam"
BAM_COUNTS_FILE = "{}output.bam_count.tsv"
READ_CLASSES_FILE = "{}read_classes.tsv.gz"
# Out files
SUMMARY_FILE = "{}summary.txt"
READ_COUNTS_FILE = "{}read_counts.tsv"

SUMMARY_TEMPLATE = (
					"----------------------------------------\n"
					"                Metadata                \n"
					"----------------------------------------\n"
					"Version:\t{version}\n"
					"Threads:\t{threads}\n"
					#TODO: Time and resources
					"\n"
					"----------------------------------------\n"
					"                Settings                \n"
					"----------------------------------------\n"
					"Input data type:\t{seq_type}\n"
					"Input reads:\t{input_reads}\n"
					"Reference HISAT2 index:\t{hisat2_index}\n"
					"Reference genome FASTA:\t{genome_fasta}\n"
					"Reference 5'ss FASTA:\t{fivep_fasta}\n"
					"Reference exons:\t{exons_tsv}\n"
					"Reference introns:\t{introns_tsv}\n"
					"Reference RepeatMasker:\t{repeats_bed}\n"
					"Output:\t{output_base}\n"
					"Log level:\t{log_level}\n"
					"Keep temporary files:\t{keep_temp}\n"
					"Make UCSC track:\t{ucsc_track}\n"
					"\n"
					"----------------------------------------\n"
					"              Read classes              \n"
					"----------------------------------------\n"
					"Linear:\t{Linear}\n"
					"Unmapped:\t{Unmapped}\n"
					"Unmapped with 5'ss alignment:\t{Unmapped_with_5ss_alignment}\n"
					"Template-switching:\t{Template_switching}\n"
					"Circularized intron:\t{Circularized_intron}\n"
					"In repetitive region:\t{In_repetitive_region}\n"
					"Lariat:\t{Lariat}\n"
					"\n"
					"----------------------------------------\n"
					"      Read count after each stage       \n"
					"----------------------------------------\n"
					"Input:\t{input_count}\n"
					"Linear mapping:\t{not_linear}\n"
					"5'ss mapping:\t{5ss_mapping}\n"
					"5'ss alignment filtering:\t{5ss_alignment_filtering}\n"
					"Head mapping:\t{Head_mapping}\n"
					"Head alignment filtering:\t{Head_alignment_filtering}\n"
					"Lariat filtering:\t{Lariat_filtering}\n"
)
READ_COUNTS_TEMPLATE = (
						"Category\tSubcategory\tReads\n"
						"Input\tTotal\t{input_count}\n"
						"Linear\tTotal\t{Linear}\n"
						"Not linear\tTotal\t{not_linear}\n"
						"Not linear\tUnmapped\t{Unmapped}\n"
						"Not linear\tUnmapped with 5'ss alignment\t{Unmapped_with_5ss_alignment}\n"
						"Not linear\tTemplate-switching\t{Template_switching}\n"
						"Not linear\tCircularized intron\t{Circularized_intron}\n"
						"Not linear\tIn repetitive region\t{In_repetitive_region}\n"
						"Not linear\tLariat\t{Lariat}"
						
)

SETTINGS_VARS = ('output_base', 'threads', 'hisat2_index', 'genome_fasta', 'fivep_fasta', 
				 'exons_tsv', 'introns_tsv', 'repeats_bed', 'keep_temp', 'ucsc_track', 
				 'pipeline_dir', 'log_level', 'seq_type')
READ_CLASSES = ("Linear", "Unmapped", "Unmapped_with_5ss_alignment", 'In_repetitive_region', 
				'Template_switching', 'Circularized_intron', 'Lariat')
STAGES = ('Linear_mapping', "5ss_mapping", "5ss_alignment_filtering", 
		  'Head_mapping', 'Head_alignment_filtering', 'Lariat_filtering', 'To_the_end')


	

#=============================================================================#
#                                  Functions                                  #
#=============================================================================#
def run_settings(output_base:str) -> dict:
	with open(ARGS_FILE.format(output_base), 'r') as json_file:
		settings = json.load(json_file)
	
	# Get input reads and convert settings to dict
	input_reads = ','.join(settings[13:])
	settings = {key: val for key, val in zip(SETTINGS_VARS, settings[:13])}
	settings['input_reads'] = input_reads

	# Remove args file, as it is no longer needed
	os.remove(f'{output_base}args.json')

	return settings



#=============================================================================#
#                                    Main                                     #
#=============================================================================#
if __name__ == '__main__':
	# Get args
	output_base, log_level, seq_type = sys.argv[1:]

	# Get logger
	log = functions.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}')

	# Initialise stats dict
	stats = {}

	# Run information
	stats['version'] = functions.version()
	settings = run_settings(output_base)
	stats.update(settings)

	# Add input read count
	r1 = pysam.FastxFile(settings['input_reads'].split(',')[0])
	stats['input_count'] = sum(1 for read in r1)

	# Add read class counts
	read_classes = pd.read_csv(READ_CLASSES_FILE.format(output_base), sep='\t', na_filter=False)
	read_classes.read_class = (read_classes.read_class
								.str.replace(' ', '_')
								.str.replace('-', '_')
								.str.replace("'", ''))
	classes = (pd.Categorical(read_classes.read_class, categories=READ_CLASSES, ordered=True)
				.value_counts()
				.to_dict())
	log.debug(f'Read classes: {classes}')
	stats.update(classes)
	stats['not_linear'] = stats['input_count'] - stats['Linear']

	# Add counts after each stage reached
	read_classes.stage_reached = (read_classes.stage_reached
									.str.replace(' ', '_')
									.str.replace('-', '_')
									.str.replace("'", ''))
	stages = (pd.Categorical(read_classes.stage_reached, categories=STAGES, ordered=True)
				.value_counts()
				.to_dict())
	log.debug(f'Stages reached: {stages}')
	reads_left = stats['input_count']
	for stage in STAGES:
		reads_left += -stages[stage]
		stats[stage] = reads_left

	# Write summary info to file
	log.debug(f'Summary stats: {stats}')
	with open(SUMMARY_FILE.format(output_base), 'w') as w:
		w.write(SUMMARY_TEMPLATE.format(**stats))

	# Write read counts to file
	read_count_stats = {key: stats[key] for key in stats if '{'+key+'}' in READ_COUNTS_TEMPLATE}
	log.debug(f'Read count stats: {read_count_stats}')
	with open(READ_COUNTS_FILE.format(output_base), 'w') as w:
		w.write(READ_COUNTS_TEMPLATE.format(**read_count_stats))

	log.debug('End of script')

