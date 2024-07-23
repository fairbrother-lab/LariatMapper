import sys
import itertools as it

import pysam
import pandas as pd
from intervaltree import Interval, IntervalTree

import functions



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
LINEAR_CLASSES = ('mRNA', 
				'pre-mRNA', 
				'Exonic', 
				'Intronic', 
				"Starts at 5'ss", 
				'Intergenic',
				'Ambiguous',
				# 'Multimap',
				)

INITIAL_COLS = ['read_id', 
				'mate',
				'chrom', 
				'align_start',
				'align_end',
				'align_is_reverse',
				'cigar',
				# 'n_maps',
				'mixed',
				]

OUT_COLS = ['read_id',
			'mate',
			'read_class',
			'stage_reached',
			'filter_failed',
			'spliced',
			'gene_id',
			# 'n_maps',
			]

CIGARTUPLE_CODES = {0: 'M',
					1: 'I',
					2: 'D',
					3: 'N',
					4: 'S',
					5: 'H',
					6: 'P',
					7: '=',
					8: 'X',
					9: 'B',
					}

# MAX_ALIGNS_PER_READ = 5



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def infer_segments(row):
	segs = []

	seg_start = row['align_start']
	seg_end = seg_start
	for op, length in row['cigar']:
		if op in ('M', 'I', 'D'):
			seg_end += length
		elif op == 'N':
			segs.append((seg_start, seg_end))
			seg_end += length
			seg_start = seg_end 
		else:
			raise RuntimeError(f"{row['blocks']}\n{row['cigar']}\n{op}, {length}")
		
	segs.append((seg_start, seg_end))
		
	return tuple(segs)


def infer_common_genes(df):
	gene_ids = []
	for i, seg_row in df.iterrows():
		features = seg_row['exons'].union(seg_row['introns'])
		seg_genes = set()
		for feat in features:
			seg_genes.update(feat.data['gene_id'])
		gene_ids.append(seg_genes)

	common_genes = set(gid for gid in gene_ids[0] if all(gid in set_genes for set_genes in gene_ids))
	return common_genes


def tree_covers_interval(tree:IntervalTree, interval:Interval) -> bool:
	total_coverage = False
	merged_tree = tree.copy()
	merged_tree.merge_overlaps(strict=False)
	for merged_interval in merged_tree:
		if merged_interval.contains_interval(interval):
			total_coverage = True
	
	return total_coverage


def parse_linear_alignments(output_base:str) -> pd.DataFrame:
	# Run through mapped reads file, loading read alignments
	linear_reads = []
	for align in pysam.AlignmentFile(f'{output_base}mapped_reads.bam', 'rb'):
		linear_reads.append([
						align.query_name, 
						align.is_read1,
						align.reference_name, 
						align.reference_start,
						align.reference_end,
						align.is_reverse,
						align.cigartuples,
						# align.get_tag('NH'),
						align.get_tag('YT'),
						])

	linear_reads = pd.DataFrame(linear_reads, columns=INITIAL_COLS)

	if len(linear_reads) == 0:
		log.warning('0 linear read alignments')
		exit()

	# Fix starting columns
	linear_reads.read_id = linear_reads.read_id.astype('string')
	linear_reads.mate = linear_reads.mate.map({True: '1', False: '2'}).astype('category')
	linear_reads.chrom = linear_reads.chrom.astype('category')
	linear_reads.cigar = linear_reads.cigar.transform(lambda cigar: tuple((CIGARTUPLE_CODES[op],length) for op, length in cigar))
	linear_reads.mixed = linear_reads.mixed.transform(lambda m: True if m=='UP' else False)

	return linear_reads


def classify_seg(row):
	if row['Intergenic'] is True:
		return 'Intergenic'
	
	if len(row['exons'])==0 and len(row['introns'])==0:
		return 'Ambiguous I'
	
	if len(row['exons'])==0:
		if row['align_is_reverse'] is False and row['align_start']==row['seg'].begin:
			if any(row['seg'].begin==intron.begin for intron in row['introns'] if intron.data['strand']=='+'):
				return "Starts at 5'ss"
		elif row['align_is_reverse'] is True and row['align_end']==row['seg'].end:
			if any(row['seg'].end==intron.end for intron in row['introns'] if intron.data['strand']=='-'):
				return "Starts at 5'ss"
		else:
			return 'Intronic'
	
	# Check all exon-intron combinations to see if any form an exon-intron junction
	# If there are no introns, Python just skips the loop and continues
	for exon, intron in it.product(row['exons'], row['introns']):
		if len(exon.data['gene_id'].intersection(intron.data['gene_id']))==0:
			continue
		
		# Check if we have an exon-to-intron junction
		if exon.end == intron.begin:
			junc_spot = intron.begin
			exon_5bp = junc_spot-5
			intron_5bp = junc_spot+4

			# Disregard if last 5bp of exon overlaps any introns
			# Disregard if first 5bp of intron overlaps any exons
			# This accounts for alternative splice sites
			if len(row['introns'].overlap(exon_5bp, junc_spot)) > 0:
				continue
			if len(row['exons'].overlap(junc_spot, intron_5bp+1)) > 0:
				continue

			return 'pre-mRNA'
		
		# Check if we have an intron-to-exon junction
		elif intron.end == exon.begin:
			junc_spot = exon.begin
			intron_5bp = junc_spot-5
			exon_5bp = junc_spot+4

			# Disregard if last 5bp of exon overlaps any introns
			# Disregard if first 5bp of intron overlaps any exons
			# This accounts for alternative splice sites
			if len(row['exons'].overlap(intron_5bp, junc_spot)) > 0:
				continue
			if len(row['introns'].overlap(junc_spot, exon_5bp+1)) > 0:
				continue

			return 'pre-mRNA'
		
	if any([row['seg'].begin==exon.begin for exon in row['exons']]) or any([row['seg'].end==exon.end for exon in row['exons']]):
		return 'Splice junction'

	if len(row['introns'])==0:
		return 'Exonic'
		
	return 'Ambiguous EI'
	

# def classify_alignment(seg_rows:pd.DataFrame) -> str:
def classify_read(seg_rows:pd.DataFrame) -> str:
	classes_set = set(seg_rows['seg_class'])
	if len(classes_set)==1:
		if classes_set == {'Splice junction'}:
			if seg_rows.spliced.any():
				return 'mRNA'
			else:
				return 'Exonic'
		else:
			return classes_set.pop()
		
	# Now we know it's 2 classes at minimum
	if 'Intergenic' in classes_set:
		return 'Ambiguous IS'
	if "Starts at 5'ss" in classes_set:
		return "Starts at 5'ss"
	if 'pre-mRNA' in classes_set:
		return 'pre-mRNA'
	if 'Splice junction' in classes_set:
		return 'mRNA'
	
	return tuple(classes_set)


# def classify_read(align_rows:pd.DataFrame) -> str:



# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	# Get args
	output_base, ref_exons, ref_introns, seq_type, log_level = sys.argv[1:]

	# Get logger
	log = functions.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}')

	# Load exons into dict of IntervalTrees
	exons = pd.read_csv(ref_exons, sep='\t')
	exons['gene_id'] = exons.gene_id.transform(lambda gid: set(gid.split(',')))
	exons['interval'] = exons.apply(lambda row: Interval(row['start'], row['end'], {'gene_id': row['gene_id'], 'strand': row['strand']}), axis=1)
	exons = exons.groupby('chrom').interval.agg(IntervalTree).to_dict()

	# Load introns into dict of IntervalTrees
	introns = pd.read_csv(ref_introns, sep='\t')
	introns['gene_id'] = introns.gene_id.transform(lambda gid: set(gid.split(',')))
	introns['interval'] = introns.apply(lambda row: Interval(row['start'], row['end'], {'gene_id': row['gene_id'], 'strand': row['strand']}), axis=1)
	introns = introns.groupby('chrom').interval.agg(IntervalTree).to_dict()

	# Load linear read alignments
	linear_reads = parse_linear_alignments(output_base)
	log.debug(f"{linear_reads.read_id.nunique():,} reads with linear alignments. {linear_reads.loc[linear_reads.mixed==True, 'read_id'].nunique():,} have mixed mate alignments")

	# Infer spliced/unspliced and explode read alignments into alignment segments
	# linear_reads['align_identifier'] = linear_reads.apply(lambda row: row['read_id'] if row['mixed'] is True else row['read_id']+row['mate'], axis=1)
	# linear_reads['align_id'] = linear_reads['read_id'] + linear_reads.groupby(['align_identifier']).cumcount().astype('str')
	linear_reads.loc[linear_reads.mixed, 'read_id'] = linear_reads.loc[linear_reads.mixed].apply(lambda row: f"{row['read_id']}/{row['mate']}", axis=1)
	linear_reads['segs'] = linear_reads.apply(infer_segments, axis=1, result_type='reduce')
	linear_reads['spliced'] = linear_reads.segs.transform(lambda segs: len(segs)>1)
	linear_reads = linear_reads.explode('segs', ignore_index=True)
	linear_reads['seg'] = linear_reads.segs.transform(lambda segs: Interval(*segs))

	# Chromosome by chromosome, add all exons and introns 
	for chrom in linear_reads.chrom.unique():
		if chrom not in exons.keys():
			log.warning(f'No exons in {chrom}')
			continue
		chrom_exons = exons[chrom]
		linear_reads.loc[linear_reads.chrom==chrom, 'exons'] = linear_reads.loc[linear_reads.chrom==chrom, 'seg'].transform(chrom_exons.overlap)
		
		if chrom not in introns.keys():
			log.warning(f'No introns in {chrom}')
			continue
		chrom_introns = introns[chrom]
		linear_reads.loc[linear_reads.chrom==chrom, 'introns'] = linear_reads.loc[linear_reads.chrom==chrom, 'seg'].transform(chrom_introns.overlap)

	linear_reads.exons = linear_reads.exons.fillna('').transform(set)
	linear_reads.introns = linear_reads.introns.fillna('').transform(set)

	# Infer intergenic
	linear_reads['Intergenic'] = (linear_reads.exons.transform(len)==0) & (linear_reads.introns.transform(len)==0)

	# For each segment, filter out exons and introns whose genes (yes, unfortunately they can have multiple) don't cover all segments of the read
	linear_reads['common_genes'] = linear_reads.read_id.map(linear_reads.groupby('read_id').apply(infer_common_genes, include_groups=False))
	linear_reads['exons'] = linear_reads.apply(lambda row: IntervalTree(exon for exon in row['exons'] if len(exon.data['gene_id'].intersection(row['common_genes']))>0), axis=1)
	linear_reads['introns'] = linear_reads.apply(lambda row: IntervalTree(intron for intron in row['introns'] if len(intron.data['gene_id'].intersection(row['common_genes']))>0), axis=1)

	# Classify segments
	linear_reads['seg_class'] = linear_reads.apply(classify_seg, axis=1)
	log.debug(f'segment class counts: {linear_reads.seg_class.value_counts().sort_index().to_dict()}')
	linear_reads.seg_class = linear_reads.seg_class.transform(lambda sc: 'Ambiguous' if sc in ('Ambiguous I', 'Ambiguous EI') else sc)

	# Classify reads based on their segment(s) class(es)
	# linear_reads['align_class'] = linear_reads.align_id.map(linear_reads.groupby('align_id').apply(classify_alignment, include_groups=False))
	linear_reads['read_class'] = linear_reads.read_id.map(linear_reads.groupby('read_id').apply(classify_read, include_groups=False))
	log.debug(f'read class counts: {linear_reads.read_class.astype("str").value_counts().sort_index().to_dict()}')
	if log_level=='DEBUG':
		linear_reads.to_csv(f'{output_base}linear_classes_raw.tsv', sep='\t', index=False)
 
	# Designate any unclassifiable combo of seg classes as 'Ambiguous'
	linear_reads.read_class = linear_reads.read_class.transform(lambda rc: 'Ambiguous' if rc in ('Ambiguous IS',) or isinstance(rc, tuple) else rc)

	# Collapse segments back into one row per read and prepare to write to file 
	linear_reads['gene_id'] = linear_reads.common_genes.transform(lambda cg: functions.str_join(tuple(cg)))
	linear_reads = (linear_reads
						.groupby(['read_id', 'read_class', 'gene_id', 'mixed'], as_index=False, observed=True)
						.agg({'spliced': any, 'mate': lambda mate: functions.str_join(mate.values, unique=True)})
	)
	if seq_type == 'paired':
		linear_reads.loc[linear_reads.mixed, 'read_id'] = linear_reads.loc[linear_reads.mixed, 'read_id'].str.slice(0,-2)
		linear_reads = (linear_reads
							.groupby(['read_id'], as_index=False, observed=True)
							.agg({col: functions.str_join for col in linear_reads.columns if col!='read_id'})
		)
	linear_reads['stage_reached'] = 'Linear mapping'
	linear_reads['filter_failed'] = ''
	linear_reads = linear_reads[OUT_COLS]

	# Write to file
	linear_reads.to_csv(f'{output_base}linear_classes.tsv', sep='\t', index=False)

	log.debug('End of script')