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

TEMP_SWITCH_COLS = ['read_id',
					'read_orient_to_gene',
					'read_seq_forward', 
					'read_head_end_pos',
					'fivep_seq',
					'fivep_sites',
					'genomic_head_end_context',
					'temp_switch_sites',
					]
CIRCULARS_COLS = ['read_id',
				'chrom',
				'strand',
				'gene_id',
				'fivep_pos',
				'head_end_pos',
				'threep_pos',
				'head_end_dist_to_threep',
				'read_orient_to_gene',
				'read_seq_forward',
				'read_head_end_pos',
				'read_head_end_nt',
				'genomic_head_end_nt',
				'genomic_head_end_context',
				]
HEAD_END_BP_RENAMES = {'head_end_pos': 'bp_pos',
						'head_end_dist_to_threep': 'bp_dist_to_threep',
						'read_head_end_pos': 'read_bp_pos',
						'read_head_end_nt': 'read_bp_nt',
						'genomic_head_end_nt': 'genomic_bp_nt',
						'genomic_head_end_context': 'genomic_bp_context',}
PUTATITVE_LARIATS_COLS = ['read_id', 
						'read_orient_to_gene', 
						'chrom', 
						'strand', 
						'gene_id', 
						'fivep_pos', 
						'bp_pos', 
						'threep_pos', 
						'bp_dist_to_threep',
						'read_seq_forward', 
						'read_bp_pos',
						'read_bp_nt', 
						'genomic_bp_nt', 
						'genomic_bp_context', 
						'head_align_quality',
						]
COMMA_JOIN_COLS = ('fivep_sites', 'gene_id', 'gaps', 'fivep_sites', 'introns', 'gene_id')

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
@dataclasses.dataclass(frozen=True)
class FivepSite():
	# Assigned at creation
	chrom: str
	pos: int
	strand: str
	gene_ids: frozenset[str] = dataclasses.field(default_factory=frozenset)


	def __str__(self):
		return f'{self.chrom};{self.pos};{self.strand}'
	

	def from_compact_str(compact_str:str, fivep_genes:dict=None) -> list:
		out = []
		for site_str in compact_str.split(','):
			chrom, pos, strand = site_str.split(';')
			pos = int(pos)
			if fivep_genes is None:
				gene_ids = frozenset()
			else:
				gene_ids = fivep_genes[site_str]

			fivep_site = FivepSite(chrom, pos, strand, gene_ids)
			out.append(fivep_site)

		return out
	

@dataclasses.dataclass(frozen=True)
class ReadTail:
	# Assigned at creation
	read_id: str
	read_orient_to_gene: str
	read_seq: str
	fivep_seq: str
	fivep_sites: list[FivepSite]
	read_bp_pos: int


	@property
	def read_bp_nt(self):
		if self.read_orient_to_gene == 'Reverse':
			return functions.reverse_complement(self.read_seq[self.read_bp_pos])
		else:
			return self.read_seq[self.read_bp_pos]


	def from_row(row:pd.Series):
		return ReadTail(read_id=row['read_id'], 
						read_orient_to_gene=row['read_orient_to_gene'], 
						read_seq=row['read_seq'],
						fivep_seq=row['fivep_seq'], 
				  		fivep_sites=row['fivep_sites'], 
						read_bp_pos=row['read_bp_pos'])


@dataclasses.dataclass
class ReadHeadAlignment():
	"""
	"""
	TAIL_INFO_ATTRS = ('read_orient_to_gene', 'read_seq', 'fivep_seq', 'fivep_sites',
						'read_bp_pos', 'read_bp_nt')
	
	# Assigned at creation
	read_id: str					# Included in output
	chrom: str						# Included in output
	align_start: int
	align_end: int
	align_is_reverse: bool
	mismatches: int
	mismatches_p: float
	gaps: list[int]
	head_align_quality: int			# Included in output
	
	# Copied from the ReadTail object 
	read_orient_to_gene: str = None		# Included in output
	read_seq: str = None	# Included in output
	fivep_seq: str = None
	fivep_sites: list[FivepSite] = None
	read_bp_pos: int = None			# Included in output
	read_bp_nt: str = None			# Included in output

	# Determined during filtering
	strand: str = None				# Included in output
	bp_pos: int	= None				# Included in output
	genomic_bp_context: str	= None	# Included in output
	genomic_bp_nt: str	= None		# Included in output
	introns: list = None
	threep_pos: int	= None			# Included in output
	bp_dist_to_threep: int = None	# Included in output
	gene_id: str = None				# Included in output
	fivep_pos: int = None			# Included in output

	# Properties
	@property
	def read_seq_forward(self):
		if self.read_seq is None:
			return None
		if self.read_orient_to_gene == 'Reverse':
			return functions.reverse_complement(self.read_seq) 
		elif self.read_orient_to_gene == 'Forward':
			return self.read_seq
		else:
			raise ValueError(f'Invalid read_orient_to_gene value: {self.read_orient_to_gene}')
		

	# Methods
	@classmethod
	def fail_out_fields(cls):
		return tuple(cls.__annotations__.keys()) + ('read_seq_forward', 'filter_failed',)


	def from_pysam(pysam_align:pysam.AlignedSegment):
		"""
		"""
		mismatch_p = pysam_align.get_tag('NM')/pysam_align.query_length
		gaps = [length for op, length in pysam_align.cigartuples if op in (1,2)]
		return ReadHeadAlignment(
					read_id = pysam_align.query_name,
					chrom = pysam_align.reference_name,
					align_start = pysam_align.reference_start,
					align_end = pysam_align.reference_end,
					align_is_reverse = pysam_align.is_reverse,
					mismatches = pysam_align.get_tag('NM'),
					mismatches_p = mismatch_p,
					gaps = gaps,
					head_align_quality = pysam_align.mapping_quality
		)


	def from_row(row:pd.Series):
		attrs = tuple(ReadHeadAlignment.__annotations__.keys())
		attr_dict = {}
		for attr in attrs:
			if attr in row:
				attr_dict[attr] = row[attr]
			elif attr in HEAD_END_BP_RENAMES.values() and attr.replace('_bp_', '_head_end_') in row:
				attr_dict[attr] = row[attr.replace('_bp_', '_head_end_')]
			else:
				attr_dict[attr] = None

		if 'fivep_sites' in row:
			attr_dict['fivep_sites'] = FivepSite.from_compact_str(row['fivep_sites'])

		# If the row is missing 'read_seq' but has 'read_seq_forward' and 'read_orient_to_gene'
		# we can infer 'read_seq' from the other two
		if 'read_seq_forward' in row and 'read_orient_to_gene' in row and 'read_seq' not in row:
			if row['read_orient_to_gene'] == 'Forward':
				attr_dict['read_seq'] = row['read_seq_forward']
			elif row['read_orient_to_gene'] == 'Reverse':
				attr_dict['read_seq'] = functions.reverse_complement(row['read_seq_forward'])
			else:
				raise ValueError(f'Invalid read_orient_to_gene value: {row["read_orient_to_gene"]}')

		return ReadHeadAlignment(**attr_dict)


	def fill_tail_info(self, tail:ReadTail):
		for attr in self.TAIL_INFO_ATTRS:
			setattr(self, attr, getattr(tail, attr))


	def to_line(self, columns:list):
		line = self.read_id
		for col in columns[1:]:
			if col in HEAD_END_BP_RENAMES.keys():
				col = HEAD_END_BP_RENAMES[col]

			val = getattr(self, col)
			val = '' if val is None else val
			# Val fixes
			if col in COMMA_JOIN_COLS and isinstance(val, (list, tuple, set, frozenset)):
				val = functions.str_join(val)
			elif col == 'read_bp_pos' and self.read_seq is not None:
				val = len(self.read_seq) - val - 1 if self.read_orient_to_gene == 'Reverse' else val
			elif col == 'fivep_pos' and val == '' and self.fivep_sites is not None:
				val = [site.pos for site in self.fivep_sites]
				val = functions.str_join(val)

			line += f'\t{val}'
			
		return line


	def write_failed_out(self, filter_failed:str):
		line = self.to_line(self.fail_out_fields()[:-1])

		line += f'\t{filter_failed}'
		
		with failed_out_lock:
			with open(FAILED_HEADS_FILE, 'a') as a:
				a.write(line + '\n')


	def write_temp_switch_out(self):
		line = self.to_line(TEMP_SWITCH_COLS[:-1])

		temp_switch_site = self.chrom + ';' + str(self.bp_pos) + ';' + self.strand
		line += f'\t{temp_switch_site}'

		with temp_switch_lock:
			with open(TEMP_SWITCH_FILE, 'a') as a:
				a.write(line + '\n')


	def write_circle_out(self):
		line = self.to_line(CIRCULARS_COLS)

		with circle_lock:
			with open(CIRCULARS_FILE, 'a') as a:
				a.write(line + '\n')


	def write_lariat_out(self):
		line = self.to_line(PUTATITVE_LARIATS_COLS)

		with filtered_out_lock:
			with open(PUTATITVE_LARIATS_FILE, 'a') as a:
				a.write(line + '\n')
	




# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def parse_intron_info(ref_introns):
	introns_df = pd.read_csv(ref_introns, sep='\t', na_filter=False)
	introns_df['gene_id_dict'] = introns_df.gene_id.transform(lambda gene_id: {'gene_id': frozenset(gene_id.split(','))})
	
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
	fivep_genes = fivep_genes.groupby('fivep_site').agg({'gene_id': frozenset}, as_index=False).gene_id.to_dict()

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
	dups = tails[['read_id', 'read_orient_to_gene']].duplicated().sum()
	assert dups==0, 'Multiple rows in tails.tsv have the same read_id and read_orient_to_gene value, '\
					'something went wrong in filter_fivep_aligns.py'

	# Parse fivep_sites into lists of FivepSite objects
	tails.fivep_sites = tails.fivep_sites.apply(FivepSite.from_compact_str, fivep_genes=fivep_genes)

	# Make ReadTail objects
	tails['readtail'] = tails.apply(ReadTail.from_row, axis=1)

	# Convert to dictionary of { str(read ID): ReadTail, ...}
	tails = tails.set_index(['read_id'], drop=True, verify_integrity=True)
	tails = tails['readtail'].to_dict()

	return tails


def yield_read_aligns(chunk_start:int, chunk_end:int):
	n_aligns = pysam.AlignmentFile(HEADS_TO_GENOME_FILE).count()
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
			if align_num > chunk_end:
				break

	# Yield the last read's alignments
	if align_num == n_aligns:
		yield current_read_id, read_aligns


def enveloping_introns(align:ReadHeadAlignment, introns:dict) -> list[Interval]:
	if align.chrom not in introns:
		return []
	
	overlaps = introns[align.chrom][align.strand].overlap(align.align_start, align.align_end)
	if len(overlaps) == 0:
		return []
	envelops = [intron for intron in overlaps if intron.begin<=align.align_start and intron.end>=align.align_end]

	return envelops


def nearest_threep_pos(align:ReadHeadAlignment) -> int:
	if len(align.introns) == 0:
		return np.nan
	
	if align.strand == '+':
		threep_pos = min(align.introns, key=lambda i: i.end - align.bp_pos).end - 1
	else:
		threep_pos = min(align.introns, key=lambda i: align.bp_pos - i.begin).begin

	return threep_pos


def match_introns_to_fiveps(align:ReadHeadAlignment) -> tuple[list[Interval], list[FivepSite]]:
	intron_matches = set()
	fivep_matches = set()
	for intron, fivep in it.product(align.introns, align.fivep_sites):
		if len(intron.data['gene_id'].intersection(fivep.gene_ids))>0:
			intron_matches.add(intron)
			fivep_matches.add(fivep)
	
	return list(intron_matches), list(fivep_matches)


def is_template_switch(align:ReadHeadAlignment) -> bool:
	bp_adj_seq = align.genomic_bp_context[BP_CONTEXT_LENGTH+1:]

	base_matches = 0
	for i in range(TEMP_SWITCH_BASES):
		if bp_adj_seq[i]==align.fivep_seq[i]:
			base_matches += 1
	
	return base_matches == TEMP_SWITCH_BASES


def is_intron_circle(align:ReadHeadAlignment) -> bool:
	return align.bp_dist_to_threep in (-2, -1, 0)


def filter_head_alignment(align:ReadHeadAlignment, 
						genome_fasta:str,
						introns:dict
						) -> None:
	# Skip cases where the read doesn't have any fivep sites with  
	# the given fivep_strand, since they aren't relevant, 
	# and there's no useful information to be gained from writing to failed out
	if sum(1 for fp in align.fivep_sites if fp.strand==align.strand) == 0:
		return

	# Infer bp position in genome
	if align.strand == '+':
		align.bp_pos = align.align_end - 1
	else:
		align.bp_pos = align.align_start

	# Filter out if the bp window would extend outside of the chromosome bounds
	if align.bp_pos - BP_CONTEXT_LENGTH < 0:
		align.write_failed_out('bp_window_less_than_0')
		return
	chrom_len = functions.get_chrom_length(f'{genome_fasta}.fai', align.chrom)
	if align.bp_pos + BP_CONTEXT_LENGTH > chrom_len:
		align.write_failed_out('bp_window_past_chrom_end')
		return
	
	# Add more info
	align.genomic_bp_context = functions.get_seq(
											genome_fasta = genome_fasta, 
											chrom = align.chrom,
											start = align.bp_pos-BP_CONTEXT_LENGTH,
											end = align.bp_pos+BP_CONTEXT_LENGTH+1,
											rev_comp = align.strand=='-')
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
	if align.read_orient_to_gene == "Forward" :
		if align.align_is_reverse is True and align.strand=='+':
			align.write_failed_out('wrong_orient')
			return
		if align.align_is_reverse is False and align.strand=='-':
			align.write_failed_out('wrong_orient')
			return
	if align.read_orient_to_gene == "Reverse":
		if align.align_is_reverse is True and align.strand=='-':
			align.write_failed_out('wrong_orient')
			return
		if align.align_is_reverse is False and align.strand=='+':
			align.write_failed_out('wrong_orient')
			return

	# Match introns to 5'ss
	align.introns, align.fivep_sites = match_introns_to_fiveps(align)
	# Filter out if no 5'ss-intron matches
	if len(align.introns) == 0:
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
						log_level:str):
	# We have to set the log level in each process because the children don't inherit the log level from their parent,
	# even if you pass the log object itself
	log = functions.get_logger(log_level)
	log.debug(f'Process {os.getpid()}: Born and assigned lines {chunk_start:,}-{chunk_end:,}')

	for read_id, read_aligns in yield_read_aligns(chunk_start, chunk_end):
		read_tail = tails[read_id]

		# Go through each alignment for the read
		# filling in information and filtering out bad alignments
		for align, fivep_strand in it.product(read_aligns, ('+', '-')):
			align.strand = fivep_strand
			align.fill_tail_info(read_tail)

			filter_head_alignment(align, genome_fasta, introns)

	log.debug(f'Process {os.getpid()}: Finished')
		

def post_processing(log):
	"""
	Perform post-processing on template-switching and circularized reads data.
	This function performs the following steps:
	1. Load data tables from specified files.
	2. Move template-switching reads that are in the circularized reads table to failed output.
	3. Collapse the template-switching reads table to one row per read.
	4. Choose one alignment for circularized intron reads that have multiple alignments.
	5. Overwrite the original files with the corrected tables.
	"""
	log.debug('Post-processing...')

	# Load data tables
	temp_switches = pd.read_csv(TEMP_SWITCH_FILE, sep='\t', na_filter=False)
	circulars = pd.read_csv(CIRCULARS_FILE, sep='\t', na_filter=False)

	# Load putative lariat read ids
	lariat_rids = pd.read_csv(PUTATITVE_LARIATS_FILE, sep='\t', na_filter=False, usecols=['read_id'])
	lariat_rids = set(lariat_rids.read_id.str.slice(0,-6))
		
	# Add read_id_base columns to both tables
	temp_switches['read_id_base'] = temp_switches.read_id.str.slice(0,-6)
	circulars['read_id_base'] = circulars.read_id.str.slice(0,-6)

	# Mark template-switching reads that are in the circularized reads table
	temp_switches['in_circulars'] = temp_switches.read_id_base.isin(circulars.read_id_base)
	# Write template-switching reads that are in the circularized reads table to failed output
	# We have to re-calculate the read_seq for write_failed_out()
	if temp_switches.in_circulars.sum() > 0:
		temp_in_circs = temp_switches[temp_switches.in_circulars].copy()
		temp_in_circs = [ReadHeadAlignment.from_row(row) for i, row in temp_in_circs.iterrows()]
		for align in temp_in_circs:
			align.write_failed_out('temp_switch_but_also_circular')

	# Mark template-switching reads that are in the putative lariat reads table
	temp_switches['in_lariats'] = temp_switches.read_id_base.isin(lariat_rids)
	# Write template-switching reads that are in the putative lariat reads table to failed output
	# We have to re-calculate the read_seq for write_failed_out()
	if temp_switches.in_lariats.sum() > 0:
		temp_in_lariats = temp_switches[temp_switches.in_lariats].copy()
		temp_in_lariats = [ReadHeadAlignment.from_row(row) for i, row in temp_in_lariats.iterrows()]
		for align in temp_in_lariats:
			align.write_failed_out('temp_switch_but_also_lariat')

	# Drop template-switiching reads that are in either other table
	temp_switches = temp_switches[(~temp_switches.in_circulars) & (~temp_switches.in_lariats)]

	# If temp_switch reads remain, collapse the table to one row per read
	if len(temp_switches) > 0:
			temp_switches = (
				temp_switches
				.assign(read_id = temp_switches.read_id_base)
				.groupby('read_id', as_index=False)
				.agg({col: lambda c: functions.str_join(sorted(c)) for col in TEMP_SWITCH_COLS[1:]})
			)

	# Choose one alignment for circularized intron reads that have multiple alignments
	# and write the other alignments to failed output
	circulars = circulars.sort_values(['read_id_base', 'read_orient_to_gene'])
	circulars['read_dup'] = circulars.duplicated('read_id_base', keep='first')
	dup_circs = circulars[circulars.read_dup].copy()
	if len(dup_circs) > 0:
		dup_circs = [ReadHeadAlignment.from_row(row) for i, row in dup_circs.iterrows()]
		for align in dup_circs:
			align.write_failed_out('not_chosen')

	# Drop duplicate circulars
	circulars = circulars[~circulars.read_dup]

	# Prepare circulars table for output
	circulars = (
		circulars
			.assign(read_id = circulars.read_id_base)
			.drop(columns=['read_dup', 'read_id_base'])
	)

	# Overwrite files with the corrected tables
	temp_switches.to_csv(TEMP_SWITCH_FILE, sep='\t', index=False)
	circulars.to_csv(CIRCULARS_FILE, sep='\t', index=False)

	log.debug('Post-processing complete')



			





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
	log.debug('Loading reference data...')
	introns, fivep_genes = parse_intron_info(ref_introns)
	log.debug('Loading read tails...')
	tails = parse_tails(fivep_genes)

	# Write headers for the outfiles
	# The rows will be written one by one as they are processed
	with open(FAILED_HEADS_FILE, 'w') as w:
		w.write('\t'.join(ReadHeadAlignment.fail_out_fields()) + '\n')
	with open(TEMP_SWITCH_FILE, 'w') as w:
		w.write('\t'.join(TEMP_SWITCH_COLS) + '\n')
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
								chunk_start, chunk_end, log_level,))
		process.start()
		processes.append(process)

	# Process the first chunk in the main process
	chunk_start, chunk_end = chunk_ranges[0]
	filter_alignments_chunk(genome_fasta, tails, introns, 
						 	chunk_start, chunk_end, log_level,)

	# Check if any processes hit an error
	for process in processes:
		process.join()
		if process.exitcode != 0:
			raise RuntimeError()

	# Post processing
	post_processing(log)



	log.debug('End of script')