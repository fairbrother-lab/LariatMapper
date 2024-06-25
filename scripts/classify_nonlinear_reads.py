import sys
import os

import numpy as np
import pandas as pd

import functions



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
CLASS_AND_STEP = (
				('Lariat', 'Lariat filtering'),
				('Circularized intron', 'Lariat filtering'),
				('Template-switching', 'Lariat filtering'),
				('In repetitive region', 'Lariat filtering'),
				("Unmapped, with 5'ss alignment", 'Trimmed alignment filtering'),
				("Unmapped, with 5'ss alignment", 'Trimmed read mapping'),
				("Unmapped, with 5'ss alignment", "5'ss alignment filtering"),
				("Unmapped", "5'ss mapping")
				  )



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def exclude_classed_reads(read_ids:np.array, read_classes) -> np.array:
	classed_rids = [row[0] for row in read_classes]
	keep = np.isin(read_ids, classed_rids, assume_unique=True, invert=True)
	read_ids = read_ids[keep]

	return read_ids



# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	# Get args
	output_base, log_level = sys.argv[1:]

	# Get logger
	log = functions.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}') 

	# Moving backwards through the pipeline, collect all reads that didn't map linearly to the genome
	# while recording their class and the furthest they got in the pipeline before being filtered out 
	read_classes = []

	lariat = pd.read_csv(f'{output_base}lariat_reads.tsv', sep='\t', usecols=[0]).read_id.unique()
	lariat = [[read_id, 'Lariat', 'To the end'] for read_id in lariat] 
	read_classes.extend(lariat) 

	circular = pd.read_csv(f'{output_base}circularized_intron_reads.tsv', sep='\t', usecols=[0]).read_id.unique()
	circular = [[read_id, 'Circularized intron', 'To the end'] for read_id in circular] 
	read_classes.extend(circular) 

	temp_switch = pd.read_csv(f'{output_base}template_switching_reads.tsv', sep='\t', usecols=[0]).read_id.unique()
	temp_switch = [[read_id, 'Template-switching', 'To the end'] for read_id in temp_switch] 
	read_classes.extend(temp_switch) 

	lariat_failed = pd.read_csv(f'{output_base}failed_lariat_alignments.tsv', sep='\t')
	lariat_failed = lariat_failed.loc[lariat_failed.filter_failed=='in_repeat', 'read_id'].unique()
	lariat_failed = exclude_classed_reads(lariat_failed, read_classes)
	lariat_failed = [[read_id, 'In repetitive region', 'Lariat filtering'] for read_id in lariat_failed] 
	read_classes.extend(lariat_failed) 

	trim_failed = pd.read_csv(f'{output_base}failed_trimmed_alignments.tsv', sep='\t', usecols=[0]).read_id
	trim_failed = trim_failed.str.slice(0,-6).unique()
	trim_failed = exclude_classed_reads(trim_failed, read_classes)
	trim_failed = [[read_id, "Unmapped, with 5'ss alignment", 'Trimmed alignment filtering'] for read_id in trim_failed] 
	read_classes.extend(trim_failed) 

	fivep_passed = pd.read_csv(f'{output_base}fivep_info_table.tsv', sep='\t', usecols=[0]).read_id
	fivep_passed = fivep_passed.str.slice(0,-6).unique()
	fivep_passed = exclude_classed_reads(fivep_passed, read_classes)
	fivep_passed = [[read_id, "Unmapped, with 5'ss alignment", 'Trimmed read mapping'] for read_id in fivep_passed]
	read_classes.extend(fivep_passed) 

	fivep_failed = pd.read_csv(f'{output_base}failed_fivep_alignments.tsv', sep='\t', usecols=[0]).read_id
	fivep_failed = fivep_failed.str.slice(0,-2).unique()
	fivep_failed = exclude_classed_reads(fivep_failed, read_classes)
	fivep_failed = [[read_id, "Unmapped, with 5'ss alignment", "5'ss alignment filtering"] for read_id in fivep_failed]
	read_classes.extend(fivep_failed)

	# Now get all the reads that didn't get past the first stage of 5'ss mapping
	unmapped_reads = set()
	with open(f'{output_base}/unmapped_reads.fa') as r:
		for line in r:
			rid = line[1:-3]
			unmapped_reads.add(rid)
			r.readline()

	unmapped_reads = np.asarray(list(unmapped_reads))
	unmapped_reads = exclude_classed_reads(unmapped_reads, read_classes)
	unmapped_reads = [[read_id, "Unmapped", "5'ss mapping"] for read_id in unmapped_reads]
	read_classes.extend(unmapped_reads)

	# Convert to a dataframe 
	read_classes = pd.DataFrame(read_classes, columns=['read_id', 'read_class', 'stage_reached'])
	read_classes['spliced'] = np.nan

	# Write to file
	read_classes.to_csv(f'{output_base}read_classes.tsv', sep='\t', mode='a', index=False, header=False, na_rep='N/A')

	log.debug('End of script')