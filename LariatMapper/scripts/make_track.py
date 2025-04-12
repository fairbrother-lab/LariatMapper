import sys

import pandas as pd

import utils





#=============================================================================#
#                                  Globals                                    #
#=============================================================================#
COLOR='66,66,66'
COMP_NTS = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}





# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def reverse_complement(seq):
	return ''.join([COMP_NTS[seq[i]] for i in range(len(seq)-1,-1,-1)])





#=============================================================================#
#                                    Main                                     #
#=============================================================================#
if __name__ == '__main__':
	# Get args
	output_base, log_level = sys.argv[1:]

	# Get logger
	log = utils.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}')

	# Load lariat reads table
	lariat_reads = pd.read_csv(f'{output_base}lariat_reads.tsv', sep='\t', na_filter=False)
	if lariat_reads.empty:
		with open(f'{output_base}lariat_reads.bed', 'w') as w:
			pass
		sys.exit(0)

	lariat_reads['head_len'] = lariat_reads.read_bp_pos + 1
	
	pos = lariat_reads.strand=='+'
	lariat_reads.loc[pos, 'head_start'] = lariat_reads.loc[pos].apply(lambda row: row['bp_pos']-(row['head_len']-1), axis=1)
	lariat_reads.loc[pos, 'head_end'] = lariat_reads.loc[pos, 'bp_pos'] + 1  # 1-based
	lariat_reads.loc[pos, 'start'] = lariat_reads.loc[pos, 'fivep_pos']
	lariat_reads.loc[pos, 'end'] = lariat_reads.loc[pos, 'head_end'] # 1-based
	lariat_reads.loc[pos, 'block_sizes'] = lariat_reads.loc[pos].apply(lambda row: f"20,{row['head_len']}", axis=1)
	lariat_reads.loc[pos, 'block_starts'] = lariat_reads.loc[pos].apply(lambda row: f"0,{int(row['head_start']-row['fivep_pos'])}", axis=1)

	neg = lariat_reads.strand=='-'
	lariat_reads.loc[neg, 'head_start'] = lariat_reads.loc[neg, 'bp_pos']
	lariat_reads.loc[neg, 'head_end'] = lariat_reads.loc[neg].apply(lambda row: row['bp_pos']+row['head_len']+1, axis=1) # 1-based
	lariat_reads.loc[neg, 'start'] = lariat_reads.loc[neg, 'head_start']
	lariat_reads.loc[neg, 'end'] = lariat_reads.loc[neg, 'fivep_pos'] + 1 # 1-based
	lariat_reads.loc[neg, 'block_sizes'] = lariat_reads.loc[neg].apply(lambda row: f"{row['head_len']},20", axis=1)
	lariat_reads.loc[neg, 'block_starts'] = lariat_reads.loc[neg].apply(lambda row: f"0,{row['fivep_pos']-(20-1)-row['bp_pos']}", axis=1)
	lariat_reads.loc[neg, 'thick_start'] = lariat_reads.loc[neg, 'bp_pos']

 	# start and end are currently float type, must correct so that 12345.0 -> 12345
	lariat_reads.start = lariat_reads.start.astype(int)
	lariat_reads.end = lariat_reads.end.astype(int)
	
	# Set the remaining column values 
	lariat_reads['score'] = 0
	lariat_reads['thick_start'] = lariat_reads.bp_pos
	lariat_reads['thick_end'] = lariat_reads.thick_start + 1
	lariat_reads['rgb'] = COLOR
	lariat_reads['n_blocks'] = 2

	# Write bed to file
	lariat_reads = lariat_reads[['chrom', 'start', 'end', 'read_id', 'score', 'strand', 'thick_start', 'thick_end', 'rgb', 'n_blocks', 'block_sizes', 'block_starts']]
	with open(f'{output_base}lariat_reads.bed', 'w') as w:
		w.write(f'track name=Lariats itemRgb=On visibility="full"\n')
	lariat_reads.to_csv(f'{output_base}lariat_reads.bed', mode='a', sep='\t', index=False, header=False)

	log.debug('End of script')

