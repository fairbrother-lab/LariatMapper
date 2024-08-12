import os
import sys

import pandas as pd
import pyfaidx
import pysam

import functions



#=============================================================================#
#                                   Globals                                   #
#=============================================================================#
OUT_TEMPLATE = (
				"----------------------------------------\n"
				"Run\n"
				"----------------------------------------\n"
				"Version: {version}\n"
				#TODO: Time and resources
				#TODO: Settings

				"\n----------------------------------------\n"
				"Read class counts\n"
				"----------------------------------------\n"
				"mRNA: {mRNA}\n" 
				"pre-mRNA: {pre_mRNA}\n"
				"Unmapped: {unmapped}\n"
				"Unmapped with 5'ss alignment: {unmapped_5p}\n"
				"Template-switching: {temp_switch}\n"
				"Circularized intron: {circ_intron}\n"
				"In repetitive region: {repetitive}\n"
				"Lariat: {lariat}\n"

				"\n----------------------------------------\n"
				"Pipeline read counts\n"
				"----------------------------------------\n"
				"Input: {input}\n"
				"Linear mapping: {linear_map}\n"
				"5'ss mapping: {fivep_map}\n"
				"5'ss alignment filtering: {fivep_filter}\n"
				"Head mapping: {head_map}\n"
				"Head alignment filtering: {head_filter}\n"
				"Lariat filtering: {lariat_filter}\n"
)



#=============================================================================#
#                                  Functions                                  #
#=============================================================================#
# def (:) -> :

	# return



# def (:) -> :

	# return





#=============================================================================#
#                                    Main                                     #
#=============================================================================#
if __name__ == '__main__':
	# Get args
	output_base, log_level = sys.argv[1:]

	# Get logger
	log = functions.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}')

	# Initialise stats dict
	stats = {}

	# Run information
	stats['version'] = functions.version()

	# Read classification stats
	read_classes = pd.read_csv(f'{output_base}read_classes.tsv', sep='\t')
	stats.update(read_classes.read_class.value_counts().to_dict())
	stats['spliced'] = read_classes.spliced.sum()

	# Pipline stats
	stats['linear_map'] = pysam.AlignmentFile(f'{output_base}output.bam').mapped
	stats['linear_unmap'] = pysam.AlignmentFile(f'{output_base}output.bam').unmapped

	fivep_map_rids = set()
	with open(f'{output_base}fivep_to_reads.sam') as r:
		for line in r:
			rid = line.split('\t')[2]
			rid = rid[:-2]
			fivep_map_rids.add(rid)
	stats['fivep_map'] = len(fivep_map_rids)
	fivep_filter_rids = pd.read_csv(f'{output_base}tails.tsv', sep='\t', usecols=['read_id'])['read_id']
	stats['fivep_filter'] = fivep_filter_rids.str.slice(0, -6).nunique()

	head_map_rids = set()
	with open(f'{output_base}heads_to_genome.sam') as r:
		for line in r:
			rid = line.split('\t')[0]
			rid = rid[:-6]
			head_map_rids.add(rid)
	stats['head_map'] = len(head_map_rids)
	head_filter_rids = pd.read_csv(f'{output_base}putative_lariats.tsv', sep='\t', usecols=['read_id'])['read_id']
	stats['head_filter'] = head_filter_rids.str.slice(0, -6).nunique()

	lariat_filter_rids = pd.read_csv(f'{output_base}lariats.tsv', sep='\t', usecols=['read_id'])['read_id']
	stats['lariat_filter'] = lariat_filter_rids.nunique()

	# Write stats to file
	with open(f'{output_base}summary_stats.tsv', 'w') as w:
		w.write(OUT_TEMPLATE.format(**stats))

	log.debug('End of script')

