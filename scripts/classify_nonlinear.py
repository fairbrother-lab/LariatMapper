import sys
import os

import numpy as np
import pandas as pd

import functions



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
CLASS_AND_STEP = (
				('Lariat', 'To the end'),
				('In repetitive region', 'Lariat filtering'),
				('Circularized intron', 'Head alignment filtering'),
				('Template-switching', 'Head alignment filtering'),
				("Unmapped with 5'ss alignment", 'Head alignment filtering'),
				("Unmapped with 5'ss alignment", 'Head mapping'),
				("Unmapped with 5'ss alignment", "5'ss alignment filtering"),
				("Unmapped", "5'ss mapping")
)

OUT_COLS = ['read_id',
			'read_class',
			'stage_reached',
			'filter_failed',
			'spliced',
			'gene_id',
			]



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def exclude_classed_reads(df:pd.DataFrame, read_classes) -> np.array:
	classed_rids = [row[0] for row in read_classes]
	# keep = np.isin(read_ids, classed_rids, assume_unique=True, invert=True)
	# read_ids = read_ids[keep]
	df = df.loc[~df.read_id.isin(classed_rids)]

	return df


def add_reads(file:str, class_:str, stage:str, read_classes:list, read_id_process=None) -> list:
	if not os.path.isfile(file):
		log.warning(f'{file} not found')
		return read_classes
	
	df = pd.read_csv(file, sep='\t', low_memory=False)

	if read_id_process is not None:
		df.read_id = df.read_id.transform(read_id_process)
	if 'filter_failed' not in df:
		df['filter_failed'] = ''
	if 'gene_id' not in df:
		df['gene_id'] = ''
	df = df.groupby('read_id', as_index=False).agg({'filter_failed': lambda ffs: functions.str_join(ffs, unique=True), 'gene_id': lambda gids: functions.str_join(gids, unique=True)})

	df = exclude_classed_reads(df, read_classes)

	reads = [[read_id, class_, stage, filter_failed, gene_id] for read_id, filter_failed, gene_id in df.values] 
	read_classes.extend(reads)

	return read_classes



# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	# Get args
	output_base, seq_type, log_level = sys.argv[1:]

	# Get logger
	log = functions.get_logger(log_level)
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

	if os.path.isfile(f'{output_base}failed_lariat_alignments.tsv'):
		lariat_failed = pd.read_csv(f'{output_base}failed_lariat_alignments.tsv', sep='\t')
		lariat_failed = lariat_failed.loc[lariat_failed.filter_failed=='in_repeat']
		lariat_failed = lariat_failed.groupby('read_id', as_index=False).agg({'gene_id': lambda gids: functions.str_join(gids, unique=True)})
		lariat_failed = exclude_classed_reads(lariat_failed, read_classes)
		lariat_failed = [[read_id, 'In repetitive region', 'Lariat filtering', 'in_repeat', gene_id] for read_id, gene_id in lariat_failed.values] 
		read_classes.extend(lariat_failed) 

	read_classes = add_reads(f'{output_base}failed_head_alignments.tsv', 
						  "Unmapped with 5'ss alignment", 
						  'Head alignment filtering', 
						  read_classes,
						  lambda read_id: read_id[:-6])

	read_classes = add_reads(f'{output_base}tails.tsv', 
						  "Unmapped with 5'ss alignment", 
						  'Head mapping', 
						  read_classes,
						  lambda read_id: read_id[:-6])

	read_classes = add_reads(f'{output_base}failed_fivep_alignments.tsv', 
						  "Unmapped with 5'ss alignment", 
						  "5'ss alignment filtering", 
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
		unmapped_reads = [[read_id, "Unmapped", "5'ss mapping", '', '',] for read_id in unmapped_reads.read_id]
		read_classes.extend(unmapped_reads)

	if len(read_classes) == 0:
		log.warning('0 nonlinear read alignments')

	# Convert to DataFrame
	read_classes = pd.DataFrame(read_classes, columns=['read_id', 'read_class', 'stage_reached', 'filter_failed', 'gene_id'])
	log.debug(f'read class counts: {read_classes.read_class.astype("str").value_counts().sort_index().to_dict()}')

	# Do this to prevent ('id_a', 'id_a,id_b') -> 'id_a,id_a,id_b' since the gene_id col may already be comma-delimited
	read_classes.gene_id = read_classes.gene_id.transform(lambda gids: functions.str_join(gids.split(','), unique=True))

	read_classes['spliced'] = pd.Series(np.nan, dtype='object')
	
	# Concat the linearly-aligned reads 
	linear_classes = pd.read_csv(f'{output_base}linear_classes.tsv', sep='\t', na_filter=False)
	read_classes = pd.concat([linear_classes.drop(columns=['mate']), read_classes])

	# Collapse paired-end reads where one mate got linearly aligned and one mate didn't
	if seq_type == 'paired':
		read_classes = read_classes.groupby('read_id', as_index=False).agg({col: functions.str_join for col in read_classes.columns if col != 'read_id'})

	# Write to file
	read_classes.to_csv(f'{output_base}read_classes.tsv', sep='\t', index=False, na_rep='N/A')

	log.debug('End of script')