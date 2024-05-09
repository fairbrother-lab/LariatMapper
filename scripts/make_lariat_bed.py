from dataclasses import dataclass, field
import os, sys, subprocess
import copy
import statistics
import math
import time
import datetime as dt
import itertools as it

import numpy as np
import pandas as pd




#=============================================================================#
#                                  Constants                                  #
#=============================================================================#
COLOR='84,39,143'





#=============================================================================#
#                                    Main                                     #
#=============================================================================#
def main():
	output_base = sys.argv[1]

	# Load lariat reads table
	lariat_bed = pd.read_csv(f'{output_base}lariat_reads.tsv', sep='\t')

	lariat_bed['head_len'] = lariat_bed.head_end - lariat_bed.head_start

	lariat_bed.loc[lariat_bed.strand=='+', 'start'] = lariat_bed.loc[lariat_bed.strand=='+', 'fivep_pos']
	lariat_bed.loc[lariat_bed.strand=='+', 'end'] = lariat_bed.loc[lariat_bed.strand=='+', 'head_end']
	lariat_bed.loc[lariat_bed.strand=='+', 'block_sizes'] = lariat_bed.loc[lariat_bed.strand=='+'].apply(lambda row: f"20,{row['head_len']}", axis=1)
	lariat_bed.loc[lariat_bed.strand=='+', 'block_starts'] = lariat_bed.loc[lariat_bed.strand=='+'].apply(lambda row: f"0,{row['head_start']-row['fivep_pos']}", axis=1)

	lariat_bed.loc[lariat_bed.strand=='-', 'start'] = lariat_bed.loc[lariat_bed.strand=='-', 'head_start']
	lariat_bed.loc[lariat_bed.strand=='-', 'end'] = lariat_bed.loc[lariat_bed.strand=='-', 'fivep_pos'] + 1
	lariat_bed.loc[lariat_bed.strand=='-', 'block_sizes'] = lariat_bed.loc[lariat_bed.strand=='-'].apply(lambda row: f"{row['head_len']},20", axis=1)
	lariat_bed.loc[lariat_bed.strand=='-', 'block_starts'] = lariat_bed.loc[lariat_bed.strand=='-'].apply(lambda row: f"0,{row['fivep_pos']-19-row['head_start']}", axis=1)

	# Set the remaining column values 
 	# start and end are currently float type, must correct so that 12345.0 -> 12345
	lariat_bed.start = lariat_bed.start.astype(int)
	lariat_bed.end = lariat_bed.end.astype(int)
	lariat_bed['score'] = 0
	lariat_bed['rgb'] = COLOR
	lariat_bed['n_blocks'] = 2

	# Write bed to file
	lariat_bed = lariat_bed[['chrom', 'start', 'end', 'read_id', 'score', 'strand', 'start', 'end', 'rgb', 'n_blocks', 'block_sizes', 'block_starts']]
	with open(f'{output_base}lariat_reads.bed', 'w') as w:
		w.write(f'track name=Lariats itemRgb=On visibility="squish"\n')
	lariat_bed.to_csv(f'{output_base}lariat_reads.bed', mode='a', sep='\t', index=False, header=False)



if __name__ == '__main__':
	main()