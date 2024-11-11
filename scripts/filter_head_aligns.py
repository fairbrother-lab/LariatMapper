import dataclasses
import itertools as it
import multiprocessing as mp
import os
import sys

from intervaltree import Interval, IntervalTree
import numpy as np
import pandas as pd
import pyfaidx
import pysam

import functions


# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
MAX_MISMATCHES = 5
MAX_MISMATCH_PERCENT = 0.1
MAX_GAP_LENGTH = 3
BP_CONTEXT_LENGTH = 8
TEMP_SWITCH_BASES = 5

TEMPLATE_SWITCHING_COLS = ['read_id',
						'read_bp_pos',
						'read_seq', 
						'fivep_seq',
						'fivep_sites',
						'genomic_bp_context',
						'temp_switch_sites',
						]
CIRCULARS_COLS = ['read_id',
				'chrom',
				'strand',
				'gene_id',
				'fivep_pos',
				'head_end_pos',
				'threep_pos',
				'head_dist_to_threep',
				'read_is_reverse',
				'read_seq',
				'read_head_end_pos',
				'read_head_end_nt',
				'genomic_head_end_nt',
				'genomic_head_end_context',
				]
PUTATITVE_LARIATS_COLS = ['read_id', 
						'read_is_reverse', 
						'chrom', 
						'strand', 
						'gene_id', 
						'fivep_pos', 
						'bp_pos', 
						'threep_pos', 
						'bp_dist_to_threep',
						'read_seq', 
						'read_bp_pos',
						'read_bp_nt', 
						'genomic_bp_nt', 
						'genomic_bp_context', 
						]

output_base = sys.argv[-2]
# In files
HEADS_TO_GENOME_FILE = f"{output_base}heads_to_genome.sam"
TAILS_FILE = f"{output_base}tails.tsv"
# Out files
FAILED_HEADS_FILE = f"{output_base}failed_head_alignments.tsv"
TEMP_SWITCH_FILE = f"{output_base}template_switching_reads.tsv"
CIRCULARS_FILE = f"{output_base}circularized_intron_reads.tsv"
PUTATITVE_LARIATS_FILE = f"{output_base}putative_lariats.tsv"

failed_out_lock = mp.Lock()
temp_switch_lock = mp.Lock()
circle_lock = mp.Lock()
filtered_out_lock = mp.Lock()

# =============================================================================#
#                                  Classes                                     #
# =============================================================================#
@dataclasses.dataclass
class FivepSite():
	# Assigned at creation
	chrom: str
	pos: int
	strand: str
	gene_ids: set[str]


	def __str__(self):
		return f'{self.chrom};{self.pos};{self.strand}'
	

	def from_compact_str(compact_str:str, fivep_genes:dict) -> list:
		out = []
		for site_str in compact_str.split(','):
			chrom, pos, strand = site_str.split(';')
			gene_ids = fivep_genes[site_str]
			fivep_site = FivepSite(chrom, pos, strand, gene_ids)
			out.append(fivep_site)

		return out
	

@dataclasses.dataclass
class ReadTail():
	# Assigned at creation
	read_id: str
	read_is_reverse: bool
	read_seq: str
	fivep_seq: str
	fivep_sites: list[FivepSite]
	read_bp_pos: int

	# Derived at creation
	read_bp_nt: str = None


	def __post_init__(self):
		if self.read_is_reverse is True:
			self.read_bp_nt = functions.reverse_complement(self.read_seq[self.read_bp_pos])
		else:
			self.read_bp_nt = self.read_seq[self.read_bp_pos]


	def from_row(row:pd.Series):
		return ReadTail(read_id=row['read_id'], 
						read_is_reverse=row['read_is_reverse'], 
						read_seq=row['read_seq'],
						fivep_seq=row['fivep_seq'], 
				  		fivep_sites=row['fivep_sites'], 
						read_bp_pos=row['read_bp_pos'])


@dataclasses.dataclass
class ReadHeadAlignment():
	"""
	"""
	TAIL_INFO_ATTRS = ('read_is_reverse', 'read_seq', 'fivep_seq', 'fivep_sites', 'read_bp_pos', 'read_bp_nt')
	
	# Assigned at creation
	read_id: str
	chrom: str
	align_start: int
	align_end: int
	align_is_reverse: bool
	mismatches: int
	mismatches_p: float
	gaps: list[int]
	quality: int
	
	# Copied from the ReadTail object 
	read_is_reverse: bool
	read_seq: str = None
	fivep_seq: str = None
	fivep_sites: list[FivepSite] = None
	read_bp_pos: int = None
	read_bp_nt: str = None

	# Determined during filtering
	strand: str = None
	bp_pos: int	= None
	genomic_bp_context: str	= None
	genomic_bp_nt: str	= None
	introns: list = None
	threep_pos: int	= None
	bp_dist_to_threep: int = None
	gene_id: str = None
	fivep_pos: int = None


	def from_pysam(pysam_align:pysam.AlignmentSegment):
		"""
		"""
		mismatch_p = pysam_align.get_tag('NM')/pysam_align.query_length
		gaps = [length for op, length in pysam_align.cigartuples if op in (1,2)]
		read_is_reverse = pysam_align.query_name.endswith('_rev')
		return ReadHeadAlignment(
					read_id = pysam_align.query_name,
					read_is_reverse = read_is_reverse,
					chrom = pysam_align.reference_name,
					align_start = pysam_align.reference_start,
					align_end = pysam_align.reference_end,
					align_is_reverse = pysam_align.is_reverse,
					mismatches = pysam_align.get_tag('NM'),
					mismatches_p = mismatch_p,
					gaps = gaps,
					quality = pysam_align.mapping_quality
				)


	def fill_tail_info(self, tail:ReadTail):
		for attr in self.TAIL_INFO_ATTRS:
			setattr(self, attr, getattr(tail, attr))


	def write_failed_out(self, filter_failed:str):
		line = ''
		for attr in PUTATITVE_LARIATS_COLS:
			val = getattr(self, attr)
			val = '' if val is None else val
			if attr in ('gene_id', ):
				val = ','.join(val)
			line += f'\t{val}'
		line += f'\t{filter_failed}'
		
		with failed_out_lock:
			with open(FAILED_HEADS_FILE, 'a') as a:
				a.write(line + '\n')


	def write_temp_switch_out(self):
		line = self.read_id[:-6] # Remove the '/X_XXX' suffix
		for attr in TEMPLATE_SWITCHING_COLS[1:-1]:
			val = getattr(self, attr)
			if attr in ('fivep_sites', ):
				val = ','.join(str(val))
			line += f'\t{val}'
		temp_switch_site = self.bp_pos + ';' + self.chrom
		line += f'\t{temp_switch_site}'

		with temp_switch_lock:
			with open(TEMP_SWITCH_FILE, 'a') as a:
				a.write(line + '\n')


	def write_circle_out(self):
		line = self.read_id[:-6] # Remove the '/X_XXX' suffix
		for attr in CIRCULARS_COLS[1:]:
			attr = attr.replace('head_end', 'bp')
			val = getattr(self, attr)
			if attr in ('gene_id',):
				val = functions.str_join(val)

		with circle_lock:
			with open(CIRCULARS_FILE, 'a') as a:
				a.write(line + '\n')


	def write_lariat_out(self):
		line = ''
		for attr in PUTATITVE_LARIATS_COLS:
			val = getattr(self, attr)
			if attr in ('gene_id',):
				val = functions.str_join(val)
			line += f'\t{val}'

		with filtered_out_lock:
			with open(PUTATITVE_LARIATS_FILE, 'a') as a:
				a.write(line + '\n')
	




# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def parse_intron_info(ref_introns):
	introns_df = pd.read_csv(ref_introns, sep='\t', na_filter=False)
	introns_df['gene_id_dict'] = introns_df.gene_id.transform(lambda gene_id: {'gene_id': set(gene_id.split(','))})
	
	introns = {}
	for chrom in introns_df.chrom.unique():
		if chrom not in introns:
			introns[chrom] = {'+': IntervalTree(), '-': IntervalTree()}
		pos_subset = introns_df.loc[(introns_df.chrom==chrom) & (introns_df.strand=='+'), ['start', 'end', 'gene_id_dict']]
		introns[chrom]['+'] = IntervalTree.from_tuples([row for i, row in pos_subset.iterrows()])
		neg_subset = introns_df.loc[(introns_df.chrom==chrom) & (introns_df.strand=='-'), ['start', 'end', 'gene_id_dict']]
		introns[chrom]['-'] = IntervalTree.from_tuples([row for i, row in neg_subset.iterrows()])

	fivep_genes = introns_df.copy()
	fivep_genes['fivep_site'] = fivep_genes.apply(lambda row: f"{row['chrom']};{row['start']};{row['strand']}" if row['strand']=='+' else f"{row['chrom']};{row['end']-1};{row['strand']}", axis=1)
	fivep_genes = fivep_genes.groupby('fivep_site').agg({'gene_id': set}, as_index=False).gene_id.to_dict()

	return introns, fivep_genes	


def parse_tails(fivep_genes:dict) -> dict:
	tails = pd.read_csv(TAILS_FILE, 
						sep='\t', 
						na_filter=False)
	'''
	Parses the tails.tsv file produced by filter_fivep_aligns.py into a dictionary of ReadTail objects
	Output: { str(read ID): ReadTail, ...}
	'''
	# Check for duplicates
	dups = tails[['read_id', 'read_is_reverse']].duplicated().sum()
	assert dups==0, 'Multiple rows in tails.tsv have the same read_id and read_is_reverse value, '\
					'something went wrong in filter_fivep_aligns.py'

	# Parse fivep_sites into lists of FivepSite objects
	tails.fivep_sites = tails.fivep_sites.apply(FivepSite.from_compact_str, fivep_genes=fivep_genes)

	# Make ReadTail objects
	tails['readtail'] = tails.apply(ReadTail.from_row, axis=1)

	# Convert to dictionary of { str(read ID): ReadTail, ...}
	tails = tails.set_index(['read_id'], drop=True, verify_integrity=True)
	tails = tails['readtail'].to_dict()

	return tails


def yield_read_aligns(chunk_start:int, chunk_end:int, n_aligns:int):
	align_file = pysam.AlignmentFile(HEADS_TO_GENOME_FILE)

	# Run to the align directly before the starting alignment
	for _ in range(chunk_start-2):
		next(align_file)

	if chunk_start == 1:
		start_align = ReadHeadAlignment.from_pysam(next(align_file))
	# Check the line before the starting line to see if the chunk starts in the middle of a read's collection of alignments  
	# If it is, move the starting line up until it reaches the next read's alignments
	else:
		prev_align = ReadHeadAlignment.from_pysam(next(align_file))
		start_align = ReadHeadAlignment.from_pysam(next(align_file))
		while prev_align.read_id == start_align.read_id:
			start_align = ReadHeadAlignment.from_pysam(next(align_file))
			chunk_start += 1

	read_aligns = [start_align]
	current_read_id = start_align.read_id
	align_num = chunk_start
	while align_num < n_aligns:
		# Get the next alignment
		align_num += 1
		align = ReadHeadAlignment.from_pysam(next(align_file))

		# If we're still in the same read's collection of alignments, add the info to fivep_sites
		if align.read_id == current_read_id:
			read_aligns.append(align)
		# If we've reached the first alignment for a new read...
		else:
			# Yield the current read's alignments
			yield current_read_id, read_aligns

			# Set to processing next read's alignments
			current_read_id = align.read_id
			read_aligns = [align]

			# If we're at or have passed the end of the assigned chunk, we're done
			# We don't do this check until we know we got all of the last read's alignments, 
			# so we might process a few lines after the assigned chunk_end 
			if align_num >= chunk_end:
				break

	# Yield the last read's alignments
	yield current_read_id, read_aligns


def enveloping_introns(align:ReadHeadAlignment, introns:dict) -> list[Interval]:
	if align.chrom not in introns:
		return []
	
	overlaps = introns[align.chrom][align.strand].overlap(align.align_start, align.align_end)
	if len(overlaps) == 0:
		return []
	envelops = [intron for intron in overlaps if intron.begin<=align.align_start and intron.end>=align.align_end]

	return envelops


def is_template_switch(align:ReadHeadAlignment) -> bool:
	bp_adj_seq = align.genomic_bp_context[BP_CONTEXT_LENGTH+1:]

	base_matches = 0
	for i in range(TEMP_SWITCH_BASES):
		if bp_adj_seq[i]==align.fivep_seq[i]:
			base_matches += 1
	
	return base_matches == TEMP_SWITCH_BASES


def nearest_threep_pos(align:ReadHeadAlignment) -> int:
	if len(align.introns) == 0:
		return np.nan
	
	if align.strand == '+':
		threep_pos = min(align.introns, key=lambda i: i.end - align.bp_pos).end - 1
	else:
		threep_pos = min(align.introns, key=lambda i: align.bp_pos - i.begin).begin

	return threep_pos


def match_introns_to_fivep(align:ReadHeadAlignment) -> tuple[list[Interval], list[FivepSite]]:
	intron_matches = {}
	fivep_matches = {}
	for intron, fivep in it.product(align.introns, align.fivep_sites):
		if len(intron.data['gene_id'].intersection(fivep.gene_ids))>0:
			intron_matches.add(intron)
			fivep_matches.add(fivep)
	
	return list(intron_matches), list(fivep_matches)


def is_intron_circle(align:ReadHeadAlignment) -> bool:
	return align.bp_dist_to_threep in (-2, -1, 0)


def filter_head_alignment(align:ReadHeadAlignment, 
						genome_fasta:str,
						introns:dict
						) -> None:
	# Omit fivep sites that are not on the strand we're investigating
	align.fivep_sites = [fivep for fivep in align.fivep_sites if fivep.strand==align.strand]
	# Skip cases where the read doesn't have any fivep sites with  
	# the given fivep_strand, since they aren't relevant, 
	# and there's no useful information to be gained from writing to failed out
	if len(align.fivep_sites) == 0:
		return

	# Infer bp position in genome
	if align.strand == '+':
		align.bp_pos = align.align_end - 1
	else:
		align.bp_pos = align.align_start

	# Filter out if the bp window would extend outside of the chromosome bounds
	if align.bp_pos - BP_CONTEXT_LENGTH//2 < 0:
		align.write_failed_out('bp_window_less_than_0')
		return
	chrom_len = pyfaidx.Faidx(genome_fasta).index['chr1'].rlen
	if align.bp_pos + BP_CONTEXT_LENGTH//2 > chrom_len:
		align.write_failed_out('bp_window_greater_than_chrom')
		return

	# Add more info
	align.genomic_bp_context = functions.get_seq(
											genome_fasta = genome_fasta, 
											chrom = align.chrom,
											start = align.bp_pos-BP_CONTEXT_LENGTH//2,
											end = align.bp_pos+BP_CONTEXT_LENGTH//2,
											rev_comp = align.strand=='-'
	)
	align.genomic_bp_nt = align.genomic_bp_context[BP_CONTEXT_LENGTH]
	align.introns = enveloping_introns(align, introns)
	align.threep_pos = nearest_threep_pos(align)
	align.bp_dist_to_threep = -abs(align.bp_pos - align.threep_pos)

	# Filter out if low-quality
	if align.mismatches > MAX_MISMATCHES:
		align.write_failed_out('mismatches')
		return
	if align.mismatches_p > MAX_MISMATCH_PERCENT:
		align.write_failed_out('mismatches')
		return
	if len(align.gaps) > 1:
		align.write_failed_out('gaps')
		return
	if len(align.gaps) == 1 and align.gaps[0] > MAX_GAP_LENGTH:
		align.write_failed_out('gaps')
		return

	# Write out if template-switching
	if is_template_switch(align) is True:
		align.write_temp_switch_out()
		return

	# Filter out if no overlaping introns of the given strand
	if len(align.introns) == 0:
		align.write_failed_out('overlap_introns')
		return

	# Filter out if it's a bad alignment orientation combination, 
	# which does NOT leave the branchpoint adjacent to the 5'ss 
	# in the read as is expected of lariats
	if align.read_is_reverse is True:
		if align.align_is_reverse is True and align.strand=='-':
			align.write_failed_out('wrong_orient')
			return
		if align.align_is_reverse is False and align.strand=='+':
			align.write_failed_out('wrong_orient')
			return
	if align.read_is_reverse is False:
		if align.align_is_reverse is True and align.strand=='+':
			align.write_failed_out('wrong_orient')
			return
		if align.align_is_reverse is False and align.strand=='-':
			align.write_failed_out('wrong_orient')
			return

	# Match introns to 5'ss
	align.introns, align.fivep_sites = match_introns_to_fiveps(align)
	# Filter out if no 5'ss-intron matches
	if align.introns == 0:
		align.write_failed_out('fivep_intron_match')
		return

	# Write out if intron circle
	if is_intron_circle(align) is True:
		align.gene_id = [functions.str_join(fivep.gene_ids, ';') for fivep in align.fivep_sites]
		align.write_circle_out()
		return
	
	# If an alignment reaches this point it probably has only 1 5'ss matched to an intron
	# but there are cases where multiple valid 5'ss remain, 
	# in which case we treat each possible 5'ss as a separate putative lariat  
	for fivep in align.fivep_sites:
		align.fivep_pos = fivep.pos
		align.gene_id = fivep.gene_ids

		# Filter out if the 5'ss is at or downstream of the tail's start
		if align.strand == '+' and align.fivep_pos > align.align_start:
			align.write_failed_out('5p_bp_order')
			return
		if align.strand == '-' and align.fivep_pos < align.align_end-1:
			align.write_failed_out('5p_bp_order')
			return

		# Write out putative lariat
		# It survived all the filters
		align.write_lariat_out()


def filter_alignments_chunk(genome_fasta:str, 
						tails:dict, 
						introns:dict, 
						chunk_start:int, 
						chunk_end:int, 
						n_aligns:int, 
						log_level:str):
	# We have to set the log level in each process because the children don't inherit the log level from their parent,
	# even if you pass the log object itself
	log = functions.get_logger(log_level)
	log.debug(f'Process {os.getpid()}: Born and assigned lines {chunk_start:,}-{chunk_end:,}')

	for read_id, read_aligns in yield_read_aligns(chunk_start, chunk_end, n_aligns):
		read_tail = tails[read_id]

		# Go through each alignment for the read
		# filling in information and filtering out bad alignments
		for align, fivep_strand in it.product(read_aligns, ('+', '-')):
			align.strand = fivep_strand
			align.fill_tail_info(read_tail)

			filter_head_alignment(align, genome_fasta, introns)

	log.debug(f'Process {os.getpid()}: Finished')
		

def post_processing():
	"""
	Perform post-processing on template-switching and circularized reads data.
	This function performs the following steps:
	1. Load data tables from specified files.
	2. Remove template-switching reads from the circularized reads table.
	3. Collapse the template-switching and circularized reads tables to one row per read.
	4. Overwrite the original files with the corrected tables.
	"""
	# Load data tables
	temp_switches = pd.read_csv(TEMP_SWITCH_FILE, sep='\t', na_filter=False)
	circulars = pd.read_csv(CIRCULARS_FILE, sep='\t', na_filter=False)
	
	# Remove template-switching reads from circularized reads
	circulars = circulars.loc[~circulars.read_id.isin(temp_switches.read_id.values)]

	# Collapse temp switch and circular tables to one row per read
	temp_switches = (temp_switches
					.groupby('read_id', as_index=False)
					.agg({col: functions.str_join for col in TEMPLATE_SWITCHING_COLS[1:]})
	)
	circulars = (circulars
					.groupby('read_id', as_index=False)
					.agg({col: functions.str_join for col in CIRCULARS_COLS[1:]})
	)

	# Overwrite files with the corrected tables
	temp_switches.to_csv(TEMP_SWITCH_FILE, sep='\t', index=False)
	circulars.to_csv(CIRCULARS_FILE, sep='\t', index=False)




			





# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	# Get args
	threads, ref_introns, genome_fasta, output_base, log_level = sys.argv[1:]

	# Get logger
	log = functions.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}')

	threads = int(threads)

	# Count the number of head alignments
	n_aligns = pysam.AlignmentFile(HEADS_TO_GENOME_FILE).count()
	log.debug(f'{n_aligns:,} head alignments')
	# If there are no reads left, end the run early 
	if n_aligns == 0:
		sys.exit(4)

	# Decide how to divide the collection of alignments into chunks
	# for parallel processing, so each thread gets a roughly equal number
	chunk_ranges = functions.decide_chunk_ranges(n_aligns, threads)
	log.debug(f'chunk_ranges: {chunk_ranges}')

	# Load reference data for processing alignments
	introns, fivep_genes = parse_intron_info(ref_introns)
	tails = parse_tails(output_base, fivep_genes)

	# Write headers for the outfiles
	# The rows will be written one by one as they are processed
	with open(FAILED_HEADS_FILE, 'w') as w:
		w.write('\t'.join(PUTATITVE_LARIATS_COLS + ['filter_failed']) + '\n')
	with open(TEMP_SWITCH_FILE, 'w') as w:
		w.write('\t'.join(TEMPLATE_SWITCHING_COLS) + '\n')
	with open(CIRCULARS_FILE, 'w') as w:
		w.write('\t'.join(CIRCULARS_COLS) + '\n')
	with open(PUTATITVE_LARIATS_FILE, 'w') as w:
		w.write('\t'.join(PUTATITVE_LARIATS_COLS) + '\n')

	# Start parallel processes, leaving the first chunk for the main process
	log.debug(f'Parallel processing {len(chunk_ranges):,} chunks...')
	processes = []
	for chunk_start, chunk_end in chunk_ranges[1:]:
		process = mp.Process(target=filter_alignments_chunk, 
							args=(genome_fasta, tails, introns, 
								chunk_start, chunk_end, n_aligns, log_level,))
		process.start()
		processes.append(process)

	# Process the first chunk in the main process
	filter_alignments_chunk(genome_fasta, tails, introns, 
						 	chunk_start, chunk_end, n_aligns, log_level,)

	# Check if any processes hit an error
	for process in processes:
		process.join()
		if process.exitcode != 0:
			raise RuntimeError()

	# Post processing
	post_processing()



	log.debug('End of script')