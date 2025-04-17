
# This script creates a table of all input reads that didn't map linearly, assigning each read
# a read class and recording the furthest stage they reached in the pipeline before being filtered 
# out (or "To the end", for lariat reads). See OUT_COLS for the full list of columns in the table.
# The read classes and stages are determined by working backwards through the pipeline, starting with
# the lariat reads and ending with the FASTA of unmapped reads from the initial mapping step,
# classifying each read by the last output file it appears in. 

import sys
import os

import numpy as np
import pandas as pd

import utils





# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
STAGES = ('Linear mapping', "Fivep mapping", "Fivep alignment filtering", 
		  'Head mapping', 'Head alignment filtering', 'Lariat filtering', 'To the end')
READ_CLASSES = ("Linear", "No alignment", "Fivep alignment", 'Template-switching', 
				'Circularized intron', 'Lariat',)

OUT_COLS = ['read_id', 			# Unique
			'read_class', 		# One of READ_CLASSES
			'stage_reached',	# One of STAGES
			'filter_failed',	# Comma-delimited list of filters the read failed, or "N/A" 
			'gene_id',			# Comma-delimited list of Gene ID(s) associated the read, or "N/A"
			]





# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def exclude_classed_reads(df:pd.DataFrame, read_classes) -> np.array:
	classed_rids = [row[0] for row in read_classes]
	df = df.loc[~df.read_id.isin(classed_rids)]

	return df


def add_reads(file:str, class_:str, stage:str, read_classes:list, read_id_process=None) -> list:
	if not os.path.isfile(file):
		log.warning(f'{file} not found')
		return read_classes
	
	df = pd.read_csv(file, sep='\t', low_memory=False, na_filter=False)

	if read_id_process is not None:
		df.read_id = df.read_id.transform(read_id_process)
	if 'filter_failed' not in df:
		df['filter_failed'] = ''
	if 'gene_id' not in df:
		df['gene_id'] = ''
	df = (df
	   		.groupby('read_id', as_index=False)
			.agg({'filter_failed': lambda ffs: utils.str_join(ffs, unique=True), 
		 			'gene_id': lambda gids: utils.str_join(gids, unique=True)})
	)

	df = exclude_classed_reads(df, read_classes)

	reads = [[read_id, class_, stage, filter_failed, gene_id] for read_id, filter_failed, gene_id in df.values] 
	read_classes.extend(reads)

	return read_classes


def correct_read_class(stage) -> str:
	if ',' not in stage:
		return stage
	
	class_a, class_b = stage.split(',')
	# Return latest class
	if STAGES.index(class_a) < STAGES.index(class_b):
		return class_b
	else:
		return class_a


def correct_stage_reached(stage) -> str:
	if ',' not in stage:
		return stage
	
	stage_a, stage_b = stage.split(',')
	# Return latest stage
	if STAGES.index(stage_a) < STAGES.index(stage_b):
		return stage_b
	else:
		return stage_a
	




# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	# Get args
	output_base, seq_type, log_level = sys.argv[1:]

	# Get logger
	log = utils.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}') 

	# Moving backwards through the pipeline, collect all reads that didn't map linearly to the genome
	# while recording their class and the furthest they got in the pipeline before being filtered out 
	read_classes = []

	read_classes = add_reads(f'{output_base}lariat_reads.tsv', 
						  'Lariat', 
						  'To the end', 
						  read_classes)

	read_classes = add_reads(f'{output_base}template_switching_reads.tsv', 
						  'Template-switching', 
						  'Head alignment filtering', 
						  read_classes)

	read_classes = add_reads(f'{output_base}circularized_intron_reads.tsv', 
						  'Circularized intron', 
						  'Head alignment filtering', 
						  read_classes)

	read_classes = add_reads(f'{output_base}failed_head_alignments.tsv', 
						  "Fivep alignment", 
						  'Head alignment filtering', 
						  read_classes,
						  lambda read_id: read_id[:-6])

	read_classes = add_reads(f'{output_base}tails.tsv', 
						  "Fivep alignment", 
						  'Head mapping', 
						  read_classes,
						  lambda read_id: read_id[:-6])

	read_classes = add_reads(f'{output_base}failed_fivep_alignments.tsv', 
						  "Fivep alignment", 
						  "Fivep alignment filtering", 
						  read_classes,
						  lambda read_id: read_id[:-2])
	
	if os.path.isfile(f'{output_base}unmapped_reads.fa'):
		# Now get all the reads that didn't get past the first stage of 5'ss mapping
		unmapped_reads = set()
		with open(f'{output_base}unmapped_reads.fa') as r:
			for line in r:
				rid = line[1:-3]
				unmapped_reads.add(rid)
				r.readline()

		unmapped_reads = pd.DataFrame({'read_id': list(unmapped_reads)})
		unmapped_reads = exclude_classed_reads(unmapped_reads, read_classes)
		unmapped_reads = [[read_id, "No alignment", "Fivep mapping", '', '',] for read_id in unmapped_reads.read_id]
		read_classes.extend(unmapped_reads)

	if len(read_classes) == 0:
		log.warning('0 nonlinear read alignments')

	# Convert to DataFrame
	read_classes = pd.DataFrame(read_classes, columns=OUT_COLS)
	log.debug(f'read class counts: {read_classes.read_class.astype("str").value_counts().sort_index().to_dict()}')

	# Do this to prevent ('id_a', 'id_a,id_b') -> 'id_a,id_a,id_b' since the gene_id col may already be comma-delimited
	read_classes.gene_id = read_classes.gene_id.transform(lambda gids: utils.str_join(gids.split(','), unique=True))

	# Write to file
	read_classes.to_csv(f'{output_base}read_classes.tsv.gz', sep='\t', index=False, na_rep='N/A')

	log.debug('End of script')