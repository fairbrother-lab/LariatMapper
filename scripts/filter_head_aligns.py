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

TEMP_SWITCH_COLS = ['read_id',
					'read_orient_to_gene',
					'read_seq_forward', 
					'read_head_end_pos',
					'fivep_seq',
					'fivep_sites',
					'head_alignment',
					'head_end_pos',
					'threep_pos',
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
HEAD_END_BP_RENAMES = {'head_end_pos': 'bp_pos',
						'head_end_dist_to_threep': 'bp_dist_to_threep',
						'read_head_end_pos': 'read_bp_pos',
						'read_head_end_nt': 'read_bp_nt',
						'genomic_head_end_nt': 'genomic_bp_nt',
						'genomic_head_end_context': 'genomic_bp_context',}
TRANS_SPLICED_COLS = ['read_id',
					'chrom', 
					'strand', 
					'gene_id', 
					'fivep_sites',
					'fivep_seq',
					'bp_pos',
					'threep_pos',
					'bp_dist_to_threep',
					'read_orient_to_gene',
					'read_seq_forward',
					'read_bp_pos',
					'read_bp_nt',
					'genomic_bp_nt',
					'bp_mismatch',
					'genomic_bp_context',
				]
PUTATITVE_LARIATS_COLS = ['read_id', 
						'gene_id', 
						'chrom', 
						'strand', 
						'fivep_pos', 
						'bp_pos', 
						'threep_pos', 
						'bp_dist_to_threep',
						'read_orient_to_gene', 
						'read_seq_forward', 
						'read_bp_pos',
						'read_bp_nt', 
						'genomic_bp_nt', 
						'bp_mismatch',
						'genomic_bp_context', 
						'head_align_quality',
						'max_quality',
						]
COMMA_JOIN_COLS = ('gaps', 'exons', 'introns', 'gene_id', 'fivep_sites')

output_base = sys.argv[-2]
# In files
HEADS_TO_GENOME_FILE = f"{output_base}heads_to_genome.sam"
TAILS_FILE = f"{output_base}tails.tsv"
# Out files
FAILED_HEADS_FILE = f"{output_base}failed_head_alignments.tsv"
TEMP_SWITCH_FILE = f"{output_base}template_switching_reads.tsv"
TRANS_SPLICED_FILE = f"{output_base}trans_spliced_lariat_reads.tsv"
CIRCULARS_FILE = f"{output_base}circularized_intron_reads.tsv"
PUTATITVE_LARIATS_FILE = f"{output_base}putative_lariats.tsv"

failed_out_lock = mp.Lock()
temp_switch_lock = mp.Lock()
trans_spliced_lock = mp.Lock()
circulars_lock = mp.Lock()
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
		gene_ids = utils.str_join(self.gene_ids, '&')
		return f"{self.chrom};{self.pos};{self.strand};{gene_ids}"
	

	def from_compact_str(compact_str:str, fivep_genes:dict=None) -> list:
		out = []
		for site_str in compact_str.split(','):
			chrom, pos, strand = site_str.split(';')[:3]
			pos = int(pos)
			gene_ids_included = len(site_str.split(';'))==4
			if gene_ids_included is True:
				gene_ids = site_str.split(';')[3]
				gene_ids = frozenset(gene_ids.split('&'))
			elif fivep_genes is not None:
				gene_ids = fivep_genes[site_str]
			else:
				gene_ids = frozenset()

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
			return utils.reverse_complement(self.read_seq[self.read_bp_pos])
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
	PROPERTIES = ('read_seq_forward', 'temp_switch_pos')
	
	# Assigned at creation
	read_id: str
	chrom: str
	align_start: int
	align_end: int
	align_is_reverse: bool
	mismatches: int
	mismatches_p: float
	gaps: list[int]
	head_align_quality: int
	
	# Copied from the ReadTail object 
	read_orient_to_gene: str = None
	read_seq: str = None
	fivep_seq: str = None
	fivep_sites: list[FivepSite] = None
	read_bp_pos: int = None
	read_bp_nt: str = None

	# Determined during filtering
	max_quality: int = None
	strand: str = None
	bp_pos: int	= None
	genomic_bp_context: str	= None
	genomic_bp_nt: str	= None
	exons: list = None
	introns: list = None
	threep_pos: int	= None	
	bp_dist_to_threep: int = None
	gene_id: str = None		
	fivep_pos: int = None	
	head_alignment: str = None	

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
	def temp_switch_pos(self):
		return self.bp_pos
	
	@property
	def bp_mismatch(self):
		return self.read_bp_nt != self.genomic_bp_nt
	
	# Methods
	@classmethod
	def all_fields(cls):
		return tuple(cls.__annotations__.keys()) + cls.PROPERTIES


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
				attr_dict['read_seq'] = utils.reverse_complement(row['read_seq_forward'])
			else:
				raise ValueError(f'Invalid read_orient_to_gene value: {row["read_orient_to_gene"]}')
			
		if 'read_bp_pos' in row and 'read_seq_forward' in row and 'read_bp_nt' not in row:
			attr_dict['read_bp_nt'] = row['read_seq_forward'][row['read_bp_pos']]
		
		if 'genomic_bp_context' in row and 'genomic_bp_nt' not in row:
			attr_dict['genomic_bp_nt'] = row['genomic_bp_context'][BP_CONTEXT_LENGTH]
			
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
			if col in COMMA_JOIN_COLS and isinstance(val, (list, tuple, set, frozenset, FivepSite)):
				val = utils.str_join(val)
			elif col == 'read_bp_pos' and self.read_seq is not None:
				val = len(self.read_seq) - val - 1 if self.read_orient_to_gene == 'Reverse' else val
			elif col == 'fivep_pos' and val == '' and self.fivep_sites is not None:
				val = [site.pos for site in self.fivep_sites]
				val = utils.str_join(val)

			line += f'\t{val}'
			
		return line

	def write_failed_out(self, filter_failed:str):
		line = self.to_line(self.all_fields())

		line += f'\t{filter_failed}'
		
		with failed_out_lock:
			with open(FAILED_HEADS_FILE, 'a') as a:
				a.write(line + '\n')

	def write_temp_switch_out(self):
		head_gene_ids = []
		for feat in self.exons.union(self.introns):
			for gene_id in feat.data['gene_id']:
				head_gene_ids.append(gene_id)
		head_gene_ids = utils.str_join(head_gene_ids, '&')
		orient = 'Forward' if self.strand == '+' else 'Reverse'
		self.head_alignment = f"{self.chrom};{self.align_start};{self.align_end};{orient};{self.strand};{head_gene_ids}"
		line = self.to_line(TEMP_SWITCH_COLS)

		with temp_switch_lock:
			with open(TEMP_SWITCH_FILE, 'a') as a:
				a.write(line + '\n')

	def write_trans_spliced_out(self):
		head_gene_ids = []
		for intron in self.introns:
			for gene_id in intron.data['gene_id']:
				head_gene_ids.append(gene_id)
		self.gene_id = utils.str_join(head_gene_ids, '&')

		line = self.to_line(TRANS_SPLICED_COLS)

		with trans_spliced_lock:
			with open(TRANS_SPLICED_FILE, 'a') as a:
				a.write(line + '\n')

	def write_circulars_out(self):
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
	exons = {}
	introns = {}
	threep_sites = {}
	for chrom in chromosomes:
		exons[chrom] = {'+': IntervalTree(), '-': IntervalTree()}
		introns[chrom] = {'+': IntervalTree(), '-': IntervalTree()}
		threep_sites[chrom] = {'+': tuple(), '-': tuple()}

	# Fill IntervalTrees with exons and introns
	for chrom, strand in it.product(chromosomes, ('+', '-')):
		for df, dict_, feat_type in ((exons_df, exons, 'exon'), (introns_df, introns, 'intron')):
			subset = df.loc[(df.chrom==chrom) & (df.strand==strand)]
			subset = [Interval(row['start'], row['end'], {'gene_id': row['gene_id']}) for _,row in subset.iterrows()]
			if len(subset) == 0:
				log.warning(f'No ({strand})-strand {feat_type}s found in chromosome "{chrom}"')
				continue
			dict_[chrom][strand].update(subset)
	
	# Make a sorted, de-duplicated tuple of threep_sites for quick lookup
	for chrom in chromosomes:
		threep_sites[chrom]['+'] = tuple(sorted(set(interval.end-1 for interval in introns[chrom]['+'])))
		threep_sites[chrom]['-'] = tuple(sorted(set(interval.begin for interval in introns[chrom]['-'])))

	# Make a dictionary of {fivep_site: gene_ids} for quick lookup
	fivep_genes = introns_df.copy()
	fivep_genes['fivep_site'] = fivep_genes.apply(lambda row: f"{row['chrom']};{row['start']};{row['strand']}" if row['strand']=='+' else f"{row['chrom']};{row['end']-1};{row['strand']}", axis=1)
	fivep_genes = (fivep_genes
					.groupby('fivep_site')
					.gene_id
					.agg(lambda gid_cells: frozenset(gid for gid_set in gid_cells for gid in gid_set))
					.to_dict()
	)

	return exons, introns, threep_sites, fivep_genes	


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
			# Find the max alignment quality among all of the current read's alignments
			max_quality = max(align.head_align_quality for align in read_aligns)
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
		max_quality = max(align.head_align_quality for align in read_aligns)
		# Yield the current read's alignments
		yield current_read_id, read_aligns, max_quality


def overlapping_features(align:ReadHeadAlignment, exons:dict, introns:dict) -> tuple[set[Interval], set[Interval]]:
	if align.chrom not in exons:
		over_exons = set()
	else:
		over_exons = exons[align.chrom][align.strand].overlap(align.align_start, align.align_end)

	if align.chrom not in introns:
		over_introns = set()
	else:
		over_introns = introns[align.chrom][align.strand].overlap(align.align_start, align.align_end)
	
	return over_exons, over_introns


def nearest_threep_pos(chrom:str, strand:str, bp_pos:int, threep_sites:dict):
	# If no threep sites for the given chrom and strand, return np.nan
	if len(threep_sites[chrom][strand]) == 0:
		return np.nan

	# If + strand, find the intron closest to the bp on the right-hand side
	# and return the end pos of that intron (-1 to convert 0-based exlusive end to 0-based inclusive pos)
	if strand == '+':
		# Get the index of the nearest threep site to the right of bp_pos (inclusive of bp_pos)
		ind = bisect.bisect_left(threep_sites[chrom][strand], bp_pos)
		# If bp_pos is further right than the last threep site, return np.nan
		if ind == len(threep_sites[chrom][strand]):
			return np.nan

	# If - strand, find the intron closest to the bp on the left-hand side
	# and return the start pos of that intron
	else:
		# Get the index of the nearest threep site to the left of bp_pos (inclusive of bp_pos)
		ind = bisect.bisect_right(threep_sites[chrom][strand], bp_pos) - 1 
		# If bp_pos is further left than the first threep site, return np.nan
		if ind == -1:
			return np.nan

	return threep_sites[chrom][strand][ind]


def match_enveloping_introns_to_fiveps(align:ReadHeadAlignment) -> tuple[set[Interval], set[FivepSite]]:
	# Exclude introns that overlap the head alignment but don't envelop it
	env_introns = set(feat for feat in align.introns if feat.begin<=align.align_start and feat.end>=align.align_end)
	
	# Return early if no introns or no 5'ss
	if len(env_introns)==0 or len(align.fivep_sites)==0:
		return set(),set()
	
	# Match introns to 5'ss by gene id
	intron_matches = set()
	fivep_matches = set()
	for intron, fivep in it.product(env_introns, align.fivep_sites):
		if len(intron.data['gene_id'].intersection(fivep.gene_ids))>0:
			intron_matches.add(intron)
			fivep_matches.add(fivep)
	
	return intron_matches, fivep_matches


def is_template_switch(align:ReadHeadAlignment, temp_switch_range:int, temp_switch_min_matches:int) -> bool:
	bp_adj_seq = align.genomic_bp_context[BP_CONTEXT_LENGTH+1:]

	base_matches = 0
	for i in range(temp_switch_range):
		if bp_adj_seq[i]==align.fivep_seq[i]:
			base_matches += 1
	
	return base_matches >= temp_switch_min_matches


def is_circular(align:ReadHeadAlignment) -> bool:
	return align.bp_dist_to_threep in (-2, -1, 0)


def filter_head_alignment(align:ReadHeadAlignment, 
						genome_fasta:str,
						exons:dict,
						introns:dict,
						threep_sites:dict,
						temp_switch_range:int,
						temp_switch_min_matches:int
						) -> None:
	# Infer bp position in genome
	if align.strand == '+':
		align.bp_pos = align.align_end - 1
	else:
		align.bp_pos = align.align_start

	# Filter out if the bp window would extend outside of the chromosome bounds
	if align.bp_pos - BP_CONTEXT_LENGTH < 0:
		align.write_failed_out('bp_window_less_than_0')
		return
	chrom_len = utils.get_chrom_length(f'{genome_fasta}.fai', align.chrom)
	if align.bp_pos + BP_CONTEXT_LENGTH > chrom_len:
		align.write_failed_out('bp_window_past_chrom_end')
		return
	
	# Add more info
	align.genomic_bp_context = utils.get_seq(genome_fasta = genome_fasta, 
											chrom = align.chrom,
											start = align.bp_pos-BP_CONTEXT_LENGTH,
											end = align.bp_pos+BP_CONTEXT_LENGTH+1,
											rev_comp = align.strand=='-')
	align.genomic_bp_nt = align.genomic_bp_context[BP_CONTEXT_LENGTH]
	align.exons, align.introns = overlapping_features(align, exons, introns)
	align.threep_pos = nearest_threep_pos(align.chrom, align.strand, align.bp_pos, threep_sites)
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
	
	# Write out if template-switching
	if is_template_switch(align, temp_switch_range, temp_switch_min_matches) is True:
		align.write_temp_switch_out()
		return

	# Filter out if no overlaping introns of the given strand
	if len(align.introns) == 0:
		align.write_failed_out('overlap_introns')
		return

	# Match introns to 5'ss
	matched_introns, matched_fivep_sites = match_enveloping_introns_to_fiveps(align)
	# Write out as trans-spliced if no matching introns
	if len(matched_introns) == 0:
		# Make sure the introns include the bp
		align.introns = set(feat for feat in align.introns if feat.overlaps(align.bp_pos))
		align.write_trans_spliced_out()
		return
	# Filter out non-matching introns and 5'ss
	align.introns, align.fivep_sites = matched_introns, matched_fivep_sites

	# Write out if intron circle
	if is_circular(align) is True:
		align.write_circulars_out()
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
						exons:dict, 
						introns: dict,
						threep_sites:dict,
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
		read_tail = tails[read_id]

		# Go through each alignment for the read
		# filling in information and filtering out bad alignments
		for align in read_aligns:
			align.max_quality = max_quality
			align.fill_tail_info(read_tail)

			# Infer strand from the alignment orientation
			# By this point we can assume that the read came from a lariat, whether or not reverse
			# transcriptase read through the branchpoint. In this case, the head cDNA has to have 
			# been built upstream in the 3' to 5' direction (AKA upstream direction) regardless of 
			# the RNA template it was built on. If it aligned in the forward orientation that means
			# the upstream direction is toward the left in the reference genome, which means the 
			# RNA was from a (+) strand gene. If it aligned in the reverse that means the 
			# upstream direction is toward the right in the reference genome, which the RNA was from
			# a (-) strand gene.
			if align.read_orient_to_gene == 'Forward' and align.align_is_reverse is False:
				align.strand = '+'
			elif align.read_orient_to_gene == 'Forward' and align.align_is_reverse is True:
				align.strand = '-'
			elif align.read_orient_to_gene == 'Reverse' and align.align_is_reverse is False:
				align.strand = '-'
			elif align.read_orient_to_gene == 'Reverse' and align.align_is_reverse is True:
				align.strand = '+'

			filter_head_alignment(align, genome_fasta, exons, introns, threep_sites, temp_switch_range, temp_switch_min_matches)
	
	log.debug(f'Process {os.getpid()}: Finished')


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


def post_processing(log):
	"""

	"""
	log.debug('Post-processing...')

	# Load data tables
	temp_switches = pd.read_csv(TEMP_SWITCH_FILE, sep='\t', na_filter=False)
	circulars = pd.read_csv(CIRCULARS_FILE, sep='\t', na_filter=False)
	trans_splices = pd.read_csv(TRANS_SPLICED_FILE, sep='\t', na_filter=False)
	lariats = pd.read_csv(PUTATITVE_LARIATS_FILE, sep='\t', na_filter=False)

	# Add read_id_base columns to tables
	temp_switches['read_id_base'] = temp_switches.read_id.str.slice(0,-6)
	circulars['read_id_base'] = circulars.read_id.str.slice(0,-6)
	trans_splices['read_id_base'] = trans_splices.read_id.str.slice(0,-6)
	lariats['read_id_base'] = lariats.read_id.str.slice(0,-6)

	# Filter out trans-spliced lariat read alignments that also have alignments in a higher-priority table
	trans_splices, temp_switches = filter_out_shared_reads(trans_splices, 'trans_spliced', temp_switches, 'template_switching')
	trans_splices, circulars = filter_out_shared_reads(trans_splices, 'trans_spliced', circulars, 'circular')
	trans_splices, lariats = filter_out_shared_reads(trans_splices, 'trans_spliced', lariats, 'lariat')

	# Filter out template-switching read alignments that also have alignments in a higher-priority table
	temp_switches, circulars = filter_out_shared_reads(temp_switches, 'template_switching', circulars, 'circular')
	temp_switches, lariats = filter_out_shared_reads(temp_switches, 'template_switching', lariats, 'lariat')

	# If template-switching reads remain, collapse the table to one row per read
	if len(temp_switches) > 0:
			temp_switches = (
				temp_switches
				.assign(read_id = temp_switches.read_id_base)
				.groupby('read_id', as_index=False)
				.agg({col: lambda c: utils.str_join(c) for col in TEMP_SWITCH_COLS[1:]})
			)

	# If trans-spliced lariat reads remain, collapse the table to one row per read
	if len(trans_splices) > 0:
			trans_splices = (
				trans_splices
				.assign(read_id = trans_splices.read_id_base)
				.groupby('read_id', as_index=False)
				.agg({col: lambda c: utils.str_join(c) for col in TRANS_SPLICED_COLS[1:]})
			)

	# Choose one alignment for circularized intron reads that have multiple alignments
	# and write the other alignments to failed output
	circulars = circulars.sort_values(['read_id_base', 'read_orient_to_gene'])
	circulars['read_dup'] = circulars.duplicated('read_id_base', keep='first')
	dup_circs = circulars[circulars.read_dup].copy()
	if len(dup_circs) > 0:
		dup_circs = [ReadHeadAlignment.from_row(row) for i, row in dup_circs.iterrows()]
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
	trans_splices.to_csv(TRANS_SPLICED_FILE, sep='\t', index=False)

	log.debug('Post-processing complete')





# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	# Get args
	threads, ref_exons, ref_introns, genome_fasta, temp_switch_filter, output_base, log_level = sys.argv[1:]

	# Get logger
	log = utils.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}')

	threads = int(threads)
	temp_switch_range = int(temp_switch_filter.split(',')[0])
	temp_switch_min_matches = int(temp_switch_filter.split(',')[1])

	# Count the number of head alignments
	n_aligns = pysam.AlignmentFile(HEADS_TO_GENOME_FILE).count()
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
	exons, introns, threep_sites, fivep_genes = parse_exons_introns(ref_exons, ref_introns, log)
	log.debug('Loading read tails...')
	tails = parse_tails(fivep_genes)

	# Write headers for the outfiles
	# The rows will be written one by one as they are processed
	with open(FAILED_HEADS_FILE, 'w') as w:
		w.write('\t'.join(ReadHeadAlignment.all_fields()) + '\tfilter_failed\n')
	with open(TEMP_SWITCH_FILE, 'w') as w:
		w.write('\t'.join(TEMP_SWITCH_COLS) + '\n')
	with open(TRANS_SPLICED_FILE, 'w') as w:
		w.write('\t'.join(TRANS_SPLICED_COLS) + '\n')
	with open(CIRCULARS_FILE, 'w') as w:
		w.write('\t'.join(CIRCULARS_COLS) + '\n')
	with open(PUTATITVE_LARIATS_FILE, 'w') as w:
		w.write('\t'.join(PUTATITVE_LARIATS_COLS) + '\n')

	# Start parallel processes, leaving the first chunk for the main process
	log.debug(f'Parallel processing {len(chunk_ranges):,} chunks...')
	processes = []
	for chunk_start, chunk_end in chunk_ranges[1:]:
		process = mp.Process(target=filter_alignments_chunk, 
							args=(genome_fasta, tails, exons, introns, threep_sites,
			 					temp_switch_range, temp_switch_min_matches, chunk_start, chunk_end, 
								log_level,))
		process.start()
		processes.append(process)

	# Process the first chunk in the main process
	chunk_start, chunk_end = chunk_ranges[0]
	filter_alignments_chunk(genome_fasta, tails, exons, introns, threep_sites, 
						 	temp_switch_range, temp_switch_min_matches, chunk_start, chunk_end, 
							log_level,)

	# Check if any processes hit an error
	for process in processes:
		process.join()
		if process.exitcode != 0:
			raise mp.ProcessError(f'Error in process {process.pid}')

	# Post processing
	post_processing(log)



	log.debug('End of script')