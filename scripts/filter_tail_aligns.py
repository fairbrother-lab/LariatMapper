import bisect
import dataclasses
import itertools as it
import multiprocessing as mp
import os
import sys

from intervaltree import Interval, IntervalTree
import numpy as np
import pandas as pd
import pysam

import utils





# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
MAX_MISMATCHES = 5
MAX_MISMATCH_PERCENT = 0.1
MAX_GAP_LENGTH = 3
BP_CONTEXT_LENGTH = 8
FIVEP_CONTEXT_UP = 2
FIVEP_CONTEXT_DOWN = 7

# Out file columns
TEMP_SWITCH_COLS = ['read_id',
					'read_orient_to_gene',
					'read_seq_forward', 
					'read_head_end_pos',
					'fivep_seq',
					'fivep_sites',
					'head_chrom',
					'head_strand',
					'head_gene_id',
					'head_first_pos',
					'head_last_pos',
					'head_end_pos',
					'threep_pos',
					'head_end_dist_to_threep',
					'genomic_head_end_context',
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
PUTATITVE_LARIATS_COLS = ['read_id', 
						'gene_id', 
						'chrom', 
						'strand', 
						'fivep_pos', 
						'bp_pos',
						'bp_align_pos',
						'threep_pos',
						'fivep_nearest_anno_pos',
						'bp_dist_to_threep',
						'fivep_dist_to_anno',
						'read_orient_to_gene', 
						'read_seq_forward', 
						'read_bp_pos',
						'read_fivep_pos',
						'read_bp_nt', 
						'genomic_bp_align_nt', 
						'bp_mismatch',
						'genomic_fivep_context', 
						'tail_align_quality',
						'tail_cigar_str',
						'tail_align_start_pos',
						'tail_align_end_pos',
						'tail_mismatches',
						'tail_gaps'
						]
TEMP_SWITCH_COLS = PUTATITVE_LARIATS_COLS + ['genomic_bp_downstream']


COL_NAME_TO_ATTR_NAME = {'read_tail_end_pos': 'read_bp_pos',
						 'align_start': 'tail_align_start_pos',
						 'align_end': 'tail_align_end_pos',
						 'mismatches': 'tail_mismatches',
						 'gaps': 'tail_gaps',
						 'head_chrom': 'chrom',
						 'head_strand': 'strand',
						 'head_gene_id': 'gene_id',
						 'head_first_pos': 'align_start',
						 'head_last_pos': 'align_end',
						 'head_end_pos': 'bp_pos',
						 'head_end_dist_to_threep': 'bp_dist_to_threep',
						 'read_head_end_nt': 'read_bp_nt',
						 'genomic_head_end_nt': 'genomic_bp_nt',
						 'genomic_head_end_context': 'genomic_bp_context',
						}
COMMA_JOIN_COLS = ('gaps', 'exons', 'introns', 'gene_id', 'fivep_sites')
ROW_JOIN_CHAR = "|"

output_base = sys.argv[-2]
# In files
TAILS_TO_GENOME_FILE = f"{output_base}tails_to_genome.sam"
HEADS_FILE = f"{output_base}heads.tsv"
# Out files
FAILED_TAILS_FILE = f"{output_base}failed_tail_alignments.tsv"
TEMP_SWITCH_FILE = f"{output_base}template_switching_reads.tsv"
CIRCULARS_FILE = f"{output_base}circularized_intron_reads.tsv"
PUTATITVE_LARIATS_FILE = f"{output_base}putative_lariats.tsv"

failed_out_lock = mp.Lock()
temp_switch_lock = mp.Lock()
circulars_lock = mp.Lock()
filtered_out_lock = mp.Lock()





# =============================================================================#
#                                  Classes                                     #
# =============================================================================#
@dataclasses.dataclass(frozen=True)
class BPSite():
	# Assigned at creation
	chrom: str
	pos: int
	strand: str
	bp_align_pos: int
	gene_ids: frozenset[str] = dataclasses.field(default_factory=frozenset)


	def __str__(self):
		gene_ids = utils.str_join(self.gene_ids, '&')
		return f"{self.chrom};{self.pos};{self.strand};{gene_ids}"
	

	def from_compact_str(compact_str:str, genes:dict=None) -> list:
		out = []
		for site_str in compact_str.split(','):
			chrom, pos, strand, bp_rel_pos = site_str.split(';')[:4]
			pos, bp_rel_pos = int(pos), int(bp_rel_pos)
			if genes is not None:
				gene_ids = frozenset([g.data['gene_id'] for g in genes[chrom][strand].at(pos)])
			else:
				gene_ids = frozenset()
			bp_align_pos = pos-bp_rel_pos if strand=='+' else pos+bp_rel_pos
			bp_site = BPSite(chrom, pos, strand, bp_align_pos, gene_ids)
			out.append(bp_site)

		return out
	

@dataclasses.dataclass(frozen=True)
class ReadHead:
	# Assigned at creation
	read_id: str
	read_orient_to_gene: str
	read_seq: str
	bp_seq: str
	bp_sites: list[BPSite]
	read_bp_pos: int
	read_bp_nt: str
	read_head_len: int

	'''
	@property
	def read_bp_nt(self):
		if self.read_orient_to_gene == 'Reverse':
			return utils.reverse_complement(self.read_seq[self.read_bp_pos])
		else:
			return self.read_seq[self.read_bp_pos]
	'''

	def from_row(row:pd.Series):
		read_head_len = row['read_bp_end'] if row['read_orient_to_gene']=='Forward' else len(row['read_seq'])-row['read_bp_start']
		return ReadHead(read_id=row['read_id'], 
						read_orient_to_gene=row['read_orient_to_gene'],
						read_seq=row['read_seq'],
						bp_seq=row['bp_seq'], 
				  		bp_sites=row['bp_sites'],
						read_bp_pos=row['read_bp_pos'],
						read_bp_nt=row['bp_seq'][-1],
						read_head_len=read_head_len)
	


@dataclasses.dataclass
class ReadTailAlignment():
	"""
	"""
	HEAD_INFO_ATTRS = ('read_orient_to_gene', 'read_seq', 'bp_seq', 'bp_sites',
						'read_bp_pos', 'read_bp_nt', 'read_head_len')
	PROPERTIES = ('read_seq_forward', 'bp_mismatch')
	
	# Assigned at creation
	read_id: str
	chrom: str
	align_start: int
	align_end: int
	align_is_reverse: bool
	mismatches: int
	mismatches_p: float
	mismatch_str: str
	gaps: list[int]
	tail_align_quality: int
	tail_align_score: int
	tail_cigar_str: str
	
	# Copied from the ReadTail object 
	read_orient_to_gene: str = None
	read_seq: str = None
	bp_seq: str = None
	bp_sites: list[BPSite] = None
	read_bp_pos: int = None
	read_bp_nt: str = None

	# Determined during filtering
	max_quality: int = None
	strand: str = None
	bp_pos: int	= None
	genomic_fivep_context: str	= None
	genomic_fivep_dnt: str	= None
	genomic_bp_align_nt: str = None
	exons: set[Interval] = None
	introns: set[Interval] = None
	threep_pos: int	= None
	bp_dist_to_threep: int = None
	gene_id: str = None		
	fivep_pos: int = None
	read_fivep_pos: int = None
	fivep_nearest_anno_pos: int = None
	fivep_dist_to_anno: int = None


	# Properties
	@property
	def read_seq_forward(self):
		if self.read_seq is None:
			return None
		if self.read_orient_to_gene == 'Reverse':
			return utils.reverse_complement(self.read_seq) 
		elif self.read_orient_to_gene == 'Forward':
			return self.read_seq
		else:
			raise ValueError(f'Invalid read_orient_to_gene value: {self.read_orient_to_gene}')

	@property
	def bp_mismatch(self):
		return self.read_bp_nt != self.genomic_bp_align_nt
	
	### Methods
	# Class methods
	@classmethod
	def all_fields(cls):
		return tuple(cls.__annotations__.keys()) + cls.PROPERTIES

	# Construction methods
	def from_pysam(pysam_align:pysam.AlignedSegment):
		"""
		For building ReadTailAlignment objects from the tails_to_genome output file
		"""
		mismatch_count = pysam_align.get_tag('NM')
		mismatch_p = mismatch_count/pysam_align.query_length
		mismatch_str = pysam_align.get_tag('MD')
		gaps = [length for op, length in pysam_align.cigartuples if op in (1,2)]
		return ReadTailAlignment(
					read_id = pysam_align.query_name,
					chrom = pysam_align.reference_name,
					align_start = pysam_align.reference_start,
					align_end = pysam_align.reference_end,
					align_is_reverse = pysam_align.is_reverse,
					mismatches = mismatch_count,
					mismatches_p = mismatch_p,
					mismatch_str=mismatch_str,
					gaps = gaps,
					tail_align_quality = pysam_align.mapping_quality,
					tail_align_score = pysam_align.get_tag('AS'),
					tail_cigar_str = pysam_align.cigarstring
		)
	
	'''
	def from_row(row:pd.Series):
		"""
		For rebuilding ReadTailAlignment objects from rows in post_processing
		"""
		attrs = tuple(ReadTailAlignment.__annotations__.keys())
		attr_dict = {}
		for attr in attrs:
			# Set attribute value to None by default
			attr_dict[attr] = None

			# Now look for the attribute in the row
			if attr in row:
				attr_dict[attr] = row[attr]
			# If neccesary, get the altered column name that corresponds to the attribute and then
			# get the value from the row, if present
			elif attr in COL_NAME_TO_ATTR_NAME.values():
				col = [k for k,v in COL_NAME_TO_ATTR_NAME.items() if v==attr][0]
				if col in row:
					attr_dict[attr] = row[col]

		if 'fivep_sites' in row:
			attr_dict['fivep_sites'] = FivepSite.from_compact_str(row['fivep_sites'])

		# If the row is missing 'read_seq' but has 'read_seq_forward' and 'read_orient_to_gene'
		# we can infer it
		if 'read_seq' not in row and 'read_seq_forward' in row and 'read_orient_to_gene' in row:
			if row['read_orient_to_gene'] == 'Forward':
				attr_dict['read_seq'] = row['read_seq_forward']
			elif row['read_orient_to_gene'] == 'Reverse':
				attr_dict['read_seq'] = utils.reverse_complement(row['read_seq_forward'])
			else:
				raise ValueError(f'Invalid read_orient_to_gene value: {row["read_orient_to_gene"]}')
		
		# If the row is missing 'read_bp_nt' but has 'read_bp_pos' and 'read_seq_forward'
		# we can infer it
		if 'read_bp_nt' not in row and 'read_bp_pos' in row and 'read_seq_forward' in row:
			attr_dict['read_bp_nt'] = row['read_seq_forward'][row['read_bp_pos']]
		
		# If the row is missing 'genomic_bp_nt' but has 'genomic_bp_context'
		# we can infer it
		if 'genomic_bp_context' in row and 'genomic_bp_nt' not in row:
			attr_dict['genomic_bp_nt'] = row['genomic_bp_context'][BP_CONTEXT_LENGTH]
			
		return ReadHeadAlignment(**attr_dict)
	'''

	def fill_head_info(self, head:ReadHead):
		for attr in self.HEAD_INFO_ATTRS:
			setattr(self, attr, getattr(head, attr))

	# Inferring info and filtering methods
	def get_overlapping_features(self, exons:dict, introns:dict) -> tuple[set[Interval], set[Interval]]:
		if self.chrom not in exons:
			over_exons = set()
		else:
			over_exons = exons[self.chrom][self.strand].overlap(self.align_start, self.align_end)

		if self.chrom not in introns:
			over_introns = set()
		else:
			over_introns = introns[self.chrom][self.strand].overlap(self.align_start, self.align_end)
		
		return over_exons, over_introns
	
	def get_nearest_threep_pos(self, threep_sites:dict):
		# If no threep sites for the given chrom and strand, return np.nan
		if self.chrom not in threep_sites or len(threep_sites[self.chrom][self.strand]) == 0:
			return np.nan

		# If + strand, find the intron closest to the bp on the right-hand side
		# and return the end pos of that intron (-1 to convert 0-based exlusive end to 0-based inclusive pos)
		if self.strand == '+':
			# Get the index of the nearest threep site to the right of bp_pos (inclusive of bp_pos)
			ind = bisect.bisect_left(threep_sites[self.chrom][self.strand], self.bp_pos)
			# If bp_pos is further right than the last threep site, return np.nan
			if ind == len(threep_sites[self.chrom][self.strand]):
				return np.nan

		# If - strand, find the intron closest to the bp on the left-hand side
		# and return the start pos of that intron
		else:
			# Get the index of the nearest threep site to the left of bp_pos (inclusive of bp_pos)
			ind = bisect.bisect_right(threep_sites[self.chrom][self.strand], self.bp_pos) - 1 
			# If bp_pos is further left than the first threep site, return np.nan
			if ind == -1:
				return np.nan

		return threep_sites[self.chrom][self.strand][ind]
	
	def head_is_downstream_of_fivep(self, fivep_pos:int) -> bool:
		if self.strand == '+':
			if fivep_pos < self.align_start:
				return True
			else:
				return False
		if self.strand == '-':
			if fivep_pos > self.align_end-1:
				return True
			else:
				return False

	# Output methods
	def to_line(self, columns:list):
		line = self.read_id
		for col in columns[1:]:
			# Rename columns to attribute names to get the right attribute value, if necessary
			if col in COL_NAME_TO_ATTR_NAME.keys():
				col = COL_NAME_TO_ATTR_NAME[col]

			val = getattr(self, col)
			val = '' if val is None else val
			# Val fixes for writing to output file
			if col in COMMA_JOIN_COLS and isinstance(val, (list, tuple, set, frozenset, BPSite)):
				val = utils.str_join(val, join_string=",")
			elif col == 'read_bp_pos' and self.read_seq is not None:
				val = len(self.read_seq) - val - 1 if self.read_orient_to_gene == 'Reverse' else val
			elif col == 'read_fivep_pos' and self.read_fivep_pos is not None and self.read_seq is not None:
				val = len(self.read_seq) - val - 1 if self.read_orient_to_gene == 'Reverse' else val
			elif col == 'bp_pos' and val == '' and self.bp_sites is not None:
				val = [site.pos for site in self.bp_sites]
				val = utils.str_join(val)

			line += f'\t{val}'
			
		return line

	def write_failed_out(self, filter_failed:str):
		line = self.to_line(self.all_fields())

		line += f'\t{filter_failed}'
		
		with failed_out_lock:
			with open(FAILED_TAILS_FILE, 'a') as a:
				a.write(line + '\n')

	def write_temp_switch_out(self):
		'''
		head_gene_ids = set()
		for feat in it.chain(self.exons, self.introns):
			for gene_id in feat.data['gene_id']:
				head_gene_ids.add(gene_id)
		self.gene_id = utils.str_join(head_gene_ids, join_string=",")
		'''

		line = self.to_line(TEMP_SWITCH_COLS)

		with temp_switch_lock:
			with open(TEMP_SWITCH_FILE, 'a') as a:
				a.write(line + '\n')

	def write_circular_out(self):
		self.gene_id = [utils.str_join(intron.data['gene_id']) for intron in self.introns]

		line = self.to_line(CIRCULARS_COLS)

		with circulars_lock:
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
# Loading data

def parse_genes(ref_genes:str, log) -> dict[str, dict[str, IntervalTree]]:
	genes_df, genes = pd.read_csv(ref_genes, sep='\t'), {}
	for _, gene_info in genes_df.iterrows():
		gene_id, chrom, strand, start, end = [gene_info[c] for c in ['gene_id', 'chrom', 'strand', 'start', 'end']]
		if chrom not in genes:
			genes[chrom] = {'+': IntervalTree(), '-': IntervalTree()}
		genes[chrom][strand].add(Interval(start, end, {'gene_id':gene_id}))
	return genes


def parse_exons_introns(ref_exons, ref_introns, log) -> tuple[dict[str, dict[str, IntervalTree]], 
															  dict[str, dict[str, IntervalTree]], 
															  dict[str, dict[str, tuple]], 
															  dict[str, frozenset]]:
	exons_df = pd.read_csv(ref_exons, sep='\t', na_filter=False)
	exons_df['gene_id'] = exons_df.gene_id.transform(lambda gene_id: frozenset(gene_id.split(',')))
	introns_df = pd.read_csv(ref_introns, sep='\t', na_filter=False)
	introns_df['gene_id'] = introns_df.gene_id.transform(lambda gene_id: frozenset(gene_id.split(',')))

	chromosomes = set(exons_df.chrom.unique()).union(introns_df.chrom.unique())

	# Initialize IntervalTrees for exons and introns
	exons, introns = {}, {}
	fivep_sites, threep_sites = {}, {}
	for chrom in chromosomes:
		exons[chrom] = {'+': IntervalTree(), '-': IntervalTree()}
		introns[chrom] = {'+': IntervalTree(), '-': IntervalTree()}
		threep_sites[chrom] = {'+': tuple(), '-': tuple()}
		fivep_sites[chrom] = {'+': tuple(), '-': tuple()}

	# Fill IntervalTrees with exons and introns
	for chrom, strand in it.product(chromosomes, ('+', '-')):
		for df, dict_, feat_type in ((exons_df, exons, 'exon'), (introns_df, introns, 'intron')):
			subset = df.loc[(df.chrom==chrom) & (df.strand==strand)]
			subset = [Interval(row['start'], row['end'], {'gene_id': row['gene_id']}) for _,row in subset.iterrows()]
			if len(subset) == 0:
				log.warning(f'No ({strand})-strand {feat_type}s found in chromosome "{chrom}"')
				continue
			dict_[chrom][strand].update(subset)
	
	# Make a sorted, de-duplicated tuple of fivep_sites and threep_sites for quick lookup
	
	for chrom in introns:
		fivep_sites[chrom]['+'] = tuple(sorted(set(interval.begin for interval in introns[chrom]['+'])))
		fivep_sites[chrom]['-'] = tuple(sorted(set(interval.end-1 for interval in introns[chrom]['-'])))
		threep_sites[chrom]['+'] = tuple(sorted(set(interval.end-1 for interval in introns[chrom]['+'])))
		threep_sites[chrom]['-'] = tuple(sorted(set(interval.begin for interval in introns[chrom]['-'])))

	'''
	# Make a dictionary of {fivep_site: gene_ids} for quick lookup
	fivep_genes = introns_df.copy()
	fivep_genes['fivep_site'] = fivep_genes.apply(lambda row: f"{row['chrom']};{row['start']};{row['strand']}" if row['strand']=='+' else f"{row['chrom']};{row['end']-1};{row['strand']}", axis=1)
	fivep_genes = (fivep_genes
					.groupby('fivep_site')
					.gene_id
					.agg(lambda gid_cells: frozenset(gid for gid_set in gid_cells for gid in gid_set))
					.to_dict()
	)
	'''

	return exons, introns, threep_sites, fivep_sites	


def parse_heads(genes:dict) -> dict:
	heads = pd.read_csv(HEADS_FILE, 
						sep='\t', 
						na_filter=False)
	'''
	Parses the heads.tsv file produced by filter_bp_aligns.py into a dictionary of ReadHead objects
	Output: { str(read ID): ReadHead, ...}
	'''
	# Check for duplicates
	dups = heads[['read_id', 'read_orient_to_gene']].duplicated().sum()
	assert dups==0, 'Multiple rows in heads.tsv have the same read_id and read_orient_to_gene value, '\
					'something went wrong in filter_bp_aligns.py'

	# Parse bp_sites into lists of BPSite objects
	heads.bp_sites = heads.bp_sites.apply(BPSite.from_compact_str, genes=genes)

	# Make ReadHead objects
	heads['readhead'] = heads.apply(ReadHead.from_row, axis=1)

	# Convert to dictionary of { str(read ID): ReadHead, ...}
	heads = heads.set_index(['read_id'], drop=True, verify_integrity=True)
	heads = heads['readhead'].to_dict()

	return heads


def yield_read_aligns(chunk_start:int, chunk_end:int):
	n_aligns = pysam.AlignmentFile(TAILS_TO_GENOME_FILE).count()
	align_file = pysam.AlignmentFile(TAILS_TO_GENOME_FILE)

	# Run to the align directly before the starting alignment
	for _ in range(chunk_start-2):
		next(align_file)

	if chunk_start == 1:
		start_align = ReadTailAlignment.from_pysam(next(align_file))
	# Check the line before the starting line to see if the chunk starts in the middle of a read's collection of alignments  
	# If it is, move the starting line up until it reaches the next read's alignments
	else:
		prev_align = ReadTailAlignment.from_pysam(next(align_file))
		start_align = ReadTailAlignment.from_pysam(next(align_file))
		while prev_align.read_id == start_align.read_id:
			start_align = ReadTailAlignment.from_pysam(next(align_file))
			chunk_start += 1

	read_aligns = [start_align]
	current_read_id = start_align.read_id
	align_num = chunk_start
	while align_num < n_aligns:
		# Get the next alignment
		align_num += 1
		align = ReadTailAlignment.from_pysam(next(align_file))

		# If we're still in the same read's collection of alignments, add the info to bp_sites
		if align.read_id == current_read_id:
			read_aligns.append(align)
		# If we've reached the first alignment for a new read...
		else:
			# Find the max alignment quality among all of the current read's alignments
			max_quality = max(align.tail_align_score for align in read_aligns)
			# Yield the current read's alignments
			yield current_read_id, read_aligns, max_quality

			# Set to processing next read's alignments
			current_read_id = align.read_id
			read_aligns = [align]

			# If we're at or have passed the end of the assigned chunk, we're done
			# We don't do this check until we know we got all of the last read's alignments, 
			# so we might process a few lines after the assigned chunk_end 
			if align_num > chunk_end:
				break

	# If this is the last set of read alignments in the file, yield it
	if align_num == n_aligns:
		# Find the max alignment quality among all of the current read's alignments
		max_quality = max(align.tail_align_score for align in read_aligns)
		# Yield the current read's alignments
		yield current_read_id, read_aligns, max_quality


# Inferring info and filtering
def match_introns_to_bp_sites(introns:set[Interval], bp_sites:list[BPSite]) -> tuple[set[Interval], set[BPSite]]:
	# Return early if no introns or no 5'ss
	if len(introns)==0 or len(bp_sites)==0:
		return set(),set()
	
	# Match introns to BP by gene id
	intron_matches = set()
	bp_matches = set()
	for intron, bp in it.product(introns, bp_sites):
		if len(intron.data['gene_id'].intersection(bp.gene_ids))>0:
			intron_matches.add(intron)
			bp_matches.add(bp)
	
	return intron_matches, bp_matches

def get_nearest_splice_sites(chrom, strand, fp_pos, bp_pos, threep_sites, fivep_sites):

	# For + strand, return closest downstream (right) 3'ss from BP and closest upstream (left) 5'ss from fivep
	# For - strand, return closest upstream (left) 3'ss from BP and closest downstream (right) 5'ss from fivep
	# Allow for overlaps between 3'ss-BP and 5'ss-fivep
	if chrom not in fivep_sites or len(fivep_sites[chrom][strand]) == 0:
		fivep = np.nan
	else:
		fivep_index = bisect.bisect_right(fivep_sites[chrom][strand], fp_pos) - 1 if strand=='+' else \
					  bisect.bisect_left(fivep_sites[chrom][strand], fp_pos)
		fivep = fivep_sites[chrom][strand][fivep_index] if fivep_index not in [-1, len(fivep_sites[chrom][strand])] else \
				np.na
	if chrom not in threep_sites or len(threep_sites[chrom][strand]) == 0:
		threep = np.nan
	else:
		threep_index = bisect.bisect_left(threep_sites[chrom][strand], bp_pos) if strand=='+' else \
					   bisect.bisect_right(threep_sites[chrom][strand], bp_pos) - 1
		threep = threep_sites[chrom][strand][threep_index] if threep_index not in [-1, len(threep_sites[chrom][strand])] else \
				 np.na
	
	return fivep, threep

def is_template_switch(bp_downstream_seq:str, fivep_seq:str, temp_switch_start:int, temp_switch_range:int, temp_switch_min_matches:int) -> bool:

	base_matches = 0
	for i in range(temp_switch_start, temp_switch_start+temp_switch_range):
		if bp_downstream_seq[i]==fivep_seq[i]:
			base_matches += 1
	
	return base_matches >= temp_switch_min_matches

def adjust_tail_alignment_mismatch(strand, mismatch_str):
	
	if strand == '+' and mismatch_str[0] == '0':
		match_found, i = False, 0
		while not match_found:
			if mismatch_str[i]=='0' and mismatch_str[i+1].isalpha():
				i += 2
			else:
				match_found = True
		starting_mismatch_count = i/2
	elif strand == '-' and mismatch_str[-1] == '0':
		match_found, i = False, len(mismatch_str)-1
		while not match_found:
			if mismatch_str[i]=='0' and mismatch_str[i-1].isalpha():
				i -= 2
			else:
				match_found = True
		starting_mismatch_count = (len(mismatch_str)-1-i)/2
	else:
		starting_mismatch_count = 0

	return starting_mismatch_count



# The info-inferring and filtering process
def filter_tail_alignment(align:ReadTailAlignment, 
						  genome_fasta:str,
						  genes:dict,
						  exons:dict,
						  introns:dict,
						  threep_sites:dict,
						  fivep_sites:dict,
						  chrom_lengths:dict,
						  temp_switch_range:int,
						  temp_switch_min_matches:int,
						  log
						  ) -> None:
	
	# Filter out BPs where strand doesn't agree with alignment orientation
	filtered_bp_sites = []
	for bp in align.bp_sites:
		if bp.strand == '+' and (align.read_orient_to_gene == 'Forward') != align.align_is_reverse:
			filtered_bp_sites.append(bp)
		if bp.strand == '-' and (align.read_orient_to_gene == 'Forward') == align.align_is_reverse:
			filtered_bp_sites.append(bp)
	if len(filtered_bp_sites) == 0:
		align.write_failed_out('no_bp_matching_alignment_orientation')
		return
	
	# Only keep BPs from the same gene as the tail alignment
	if align.chrom not in genes:
		align.write_failed_out('no_genes_on_chrom')
		return
	align_genes = {s:set([g.data['gene_id'] for g in genes[align.chrom][s].overlap(align.align_start, align.align_end)]) for s in ['+', '-']}
	gene_match_bps = []
	for bp in filtered_bp_sites:
		matching_genes = bp.gene_ids.intersection(align_genes[bp.strand])
		if len(matching_genes) > 0:
			gene_match_bps.append((bp, matching_genes))
	if len(gene_match_bps) == 0:
		align.write_failed_out('no_bp_in_same_gene')
		return
	
	# Filter out alignments that aren't in inverted order on the transcript
	proper_order_bps = []
	for bp_info in gene_match_bps:
		if (bp_info[0].strand == '+' and align.align_end < bp_info[0].pos) or (bp_info[0].strand == '-' and align.align_start > bp_info[0].pos):
			proper_order_bps.append(bp_info)
	if len(proper_order_bps) == 0:
		align.write_failed_out('5p_bp_out_of_order')
		return
	elif len(proper_order_bps) > 1:
		align.write_failed_out('multi_bp_in_same_gene')
		return
	
	bp, gene_set = proper_order_bps[0]
	align.bp_pos, align.strand, align.bp_align_pos = bp.pos, bp.strand, bp.bp_align_pos
	align.gene_id = gene_set
	
	# Infer 5'ss position in genome
	starting_mismatch_count = adjust_tail_alignment_mismatch(align.strand, align.mismatch_str)
	if starting_mismatch_count > 2:
		align.write_failed_out('fivep_starting_mismatch_too_long')
		return
	
	align.read_fivep_pos = align.read_bp_pos+1+starting_mismatch_count if align.read_orient_to_gene=='Forward' else align.read_bp_pos-1-starting_mismatch_count
	if align.strand == '+':
		align.fivep_pos = align.align_start + starting_mismatch_count
		fivep_window_start, fivep_window_end = align.fivep_pos-FIVEP_CONTEXT_UP, align.fivep_pos+FIVEP_CONTEXT_DOWN
		bp_window_start, bp_window_end = align.bp_align_pos, align.bp_align_pos+5
	else:
		align.fivep_pos = align.align_end - 1 - starting_mismatch_count
		fivep_window_start, fivep_window_end = align.fivep_pos-FIVEP_CONTEXT_DOWN, align.fivep_pos+FIVEP_CONTEXT_UP
		bp_window_start, bp_window_end = align.bp_align_pos-5, align.bp_align_pos

	# Filter out if the 5'ss window would extend outside of the chromosome bounds
	if fivep_window_start < 0:
		align.write_failed_out('fivep_window_less_than_0')
		return

	if fivep_window_end > chrom_lengths[align.chrom]:
		align.write_failed_out('fivep_window_past_chrom_end')
		return
	
	# Add more info
	align.genomic_fivep_context = utils.get_seq(genome_fasta = genome_fasta, 
												chrom = align.chrom,
												start = fivep_window_start,
												end = fivep_window_end+1,
												rev_comp = align.strand=='-')
	align.genomic_fivep_dnt = align.genomic_fivep_context[FIVEP_CONTEXT_UP:FIVEP_CONTEXT_UP+2]
	genomic_bp_window = utils.get_seq(genome_fasta = genome_fasta,
										  	  chrom = align.chrom,
											  start = bp_window_start,
											  end = bp_window_end+1,
											  rev_comp = align.strand=='-')
	align.genomic_bp_align_nt, bp_downstream_seq = genomic_bp_window[0], genomic_bp_window[1:]
	align.fivep_nearest_anno_pos, align.threep_pos = get_nearest_splice_sites(align.chrom, align.strand, align.fivep_pos, align.bp_pos, threep_sites, fivep_sites)
	align.fivep_dist_to_anno = abs(align.fivep_pos - align.fivep_nearest_anno_pos) if not pd.isna(align.fivep_nearest_anno_pos) else np.nan
	align.bp_dist_to_threep = -abs(align.bp_pos - align.threep_pos) if not pd.isna(align.threep_pos) else np.nan


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
	
	if not pd.isna(align.bp_dist_to_threep) and align.bp_dist_to_threep in (-2, -1, 0):
		align.write_failed_out('circular_intron')
		return
	
	# Write out if template-switching
	if is_template_switch(bp_downstream_seq, align.genomic_fivep_context[FIVEP_CONTEXT_UP:], 0, 5, 5):
		align.genomic_bp_downstream = bp_downstream_seq
		align.write_temp_switch_out()
		return

	'''
	# Select features that envelop the head alignment
	align.exons, align.introns = align.get_overlapping_features(exons, introns)
	enveloping_exons = set(feat for feat in align.exons if feat.begin<=align.align_start and feat.end>=align.align_end)
	enveloping_introns = set(feat for feat in align.introns if feat.begin<=align.align_start and feat.end>=align.align_end)

	# Match introns to 5'ss
	matched_introns, matched_fivep_sites = match_introns_to_fivep_sites(enveloping_introns, align.fivep_sites)

	# Write out if circularized intron
	upstream_fivep_sites = [fivep for fivep in matched_fivep_sites if align.head_is_downstream_of_fivep(fivep.pos)]
	if align.is_circularized_intron(matched_introns, upstream_fivep_sites) is True:
		align.introns = matched_introns
		align.fivep_sites = upstream_fivep_sites
		align.write_circular_out()
		return
	
	# Filter out if no overlapping introns 
	if len(align.introns) == 0:
		align.write_failed_out('no_overlapping_introns')
		return

	# Filter out if no enveloping introns 
	if len(enveloping_introns) == 0:
		align.write_failed_out('no_enveloping_introns')
		return
	
	# Reduce exons and introns to enveloping subset
	align.exons = enveloping_exons
	align.introns = enveloping_introns
	
	# Filter out if no matching intron-5'ss matches
	if len(matched_introns) == 0:
		align.write_failed_out('no_matching_introns')
		return
	
	# Reduce introns and 5'ss to matching subset
	align.introns = matched_introns
	align.fivep_sites = matched_fivep_sites
	'''

	# If an alignment reaches this point it probably has only 1 5'ss matched to an intron
	# but there are cases where multiple valid 5'ss remain, 
	# in which case we treat each possible 5'ss as a separate putative lariat  
	#for fivep in align.fivep_sites:

	# Write out putative lariat
	# It survived all the filters
	#log.debug(f'Tail alignment: {align}')
	align.write_lariat_out()


def filter_alignments_chunk(genome_fasta:str, 
						heads:dict, 
						genes: dict,
						exons:dict, 
						introns: dict,
						threep_sites:dict,
						fivep_sites:dict,
						chrom_lengths:dict,
						temp_switch_range:int, 
						temp_switch_min_matches:int,
						chunk_start:int, 
						chunk_end:int, 
						log_level:str):
	# We have to set the log level in each process because the children don't inherit the log level from their parent,
	# even if you pass the log object itself
	log = utils.get_logger(log_level)
	log.debug(f'Process {os.getpid()}: Born and assigned lines {chunk_start:,}-{chunk_end:,}')

	for read_id, read_aligns, max_quality in yield_read_aligns(chunk_start, chunk_end):
		read_head = heads[read_id]

		# Go through each alignment for the read
		# filling in information and filtering out bad alignments
		for align in read_aligns:
			align.max_quality = max_quality
			align.fill_head_info(read_head)

			filter_tail_alignment(align, genome_fasta, genes, exons, introns, threep_sites, fivep_sites, chrom_lengths, temp_switch_range, temp_switch_min_matches, log)
	
	log.debug(f'Process {os.getpid()}: Finished')


'''
def filter_out_shared_reads(class_df, class_name, higher_priority_class_df, higher_priority_class_name):
	# Check for reads in both classes
	shared_reads = class_df.loc[class_df.read_id_base.isin(higher_priority_class_df.read_id_base)]

	# If there are any shared reads, filter out the alignments from those reads class_df
	if shared_reads.shape[0] > 0:
		# Make label for failing alignments
		fail_label = f'{class_name}_but_also_{higher_priority_class_name}'
		# Convert shared_reads rows to ReadHeadAlignment objs
		shared_reads = [ReadHeadAlignment.from_row(row) for i, row in shared_reads.iterrows()]
		# Write shared_reads alignments to failed out
		for align in shared_reads:
			align.write_failed_out(fail_label)
			
		# Filter shared reads out of class_df
		class_df = class_df.loc[~class_df.read_id_base.isin(higher_priority_class_df.read_id_base)]

	return class_df, higher_priority_class_df
'''

def post_processing(log):
	"""

	"""
	log.debug('Post-processing...')

	# Load data tables
	temp_switches = pd.read_csv(TEMP_SWITCH_FILE, sep='\t', na_filter=False)
	circulars = pd.read_csv(CIRCULARS_FILE, sep='\t', na_filter=False)

	# Add read_id_base columns to tables
	temp_switches['read_id_base'] = temp_switches.read_id.str.slice(0,-6)
	circulars['read_id_base'] = circulars.read_id.str.slice(0,-6)

	# If template-switching reads remain, collapse the table to one row per read
	if len(temp_switches) > 0:
			temp_switches = (
				temp_switches
				.assign(read_id = temp_switches.read_id_base)
				.groupby('read_id', as_index=False)
				.agg({col: lambda c: utils.str_join(c, ROW_JOIN_CHAR) for col in TEMP_SWITCH_COLS[1:]})
			)

	# Choose one alignment for circularized intron reads that have multiple alignments
	# and write the other alignments to failed output
	circulars = circulars.sort_values(['read_id_base', 'read_orient_to_gene'])
	circulars['read_dup'] = circulars.duplicated('read_id_base', keep='first')
	dup_circs = circulars[circulars.read_dup].copy()
	if len(dup_circs) > 0:
		dup_circs = [ReadTailAlignment.from_row(row) for i, row in dup_circs.iterrows()]
		for align in dup_circs:
			align.write_failed_out('circular_not_chosen')
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
	threads, ref_genes, ref_exons, ref_introns, genome_fasta, temp_switch_filter, output_base, log_level = sys.argv[1:]

	# Get logger
	log = utils.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}')

	threads = int(threads)
	temp_switch_range = int(temp_switch_filter.split(',')[0])
	temp_switch_min_matches = int(temp_switch_filter.split(',')[1])

	# Count the number of tail alignments
	n_aligns = pysam.AlignmentFile(TAILS_TO_GENOME_FILE).count()
	log.debug(f'{n_aligns:,} head alignments')
	# If there are no reads left, end the run early 
	if n_aligns == 0:
		sys.exit(4)

	# Decide how to divide the collection of alignments into chunks
	# for parallel processing, so each thread gets a roughly equal number
	chunk_ranges = utils.decide_chunk_ranges(n_aligns, threads)
	log.debug(f'chunk_ranges: {chunk_ranges}')

	# Load reference data for processing alignments
	log.debug('Loading reference data...')
	genes = parse_genes(ref_genes, log)
	exons, introns, threep_sites, fivep_sites = parse_exons_introns(ref_exons, ref_introns, log)
	chrom_lengths = utils.get_chrom_lengths(f'{genome_fasta}.fai')
	log.debug('Loading read tails...')
	heads = parse_heads(genes)

	# Write headers for the outfiles
	# The rows will be written one by one as they are processed
	with open(FAILED_TAILS_FILE, 'w') as w:
		w.write('\t'.join(ReadTailAlignment.all_fields()) + '\tfilter_failed\n')
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
							args=(genome_fasta, heads, genes, exons, introns, threep_sites, fivep_sites, chrom_lengths,
			 					temp_switch_range, temp_switch_min_matches, chunk_start, chunk_end, 
								log_level,))
		process.start()
		processes.append(process)

	# Process the first chunk in the main process
	chunk_start, chunk_end = chunk_ranges[0]
	filter_alignments_chunk(genome_fasta, heads, genes, exons, introns, threep_sites, fivep_sites, chrom_lengths,
						 	temp_switch_range, temp_switch_min_matches, chunk_start, chunk_end, 
							log_level,)

	# Check if any processes hit an error
	for process in processes:
		process.join()
		if process.exitcode != 0:
			raise mp.ProcessError(f'Error in process {process.pid}')

	# Post processing
	#post_processing(log)

	log.debug('End of script')