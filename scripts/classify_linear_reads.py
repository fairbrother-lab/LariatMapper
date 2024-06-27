import sys
import itertools as it

import pysam
import pandas as pd
from intervaltree import Interval, IntervalTree

import functions



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
COLUMNS = ('read_id', 
			'chrom', 
			'align_is_reverse',
			'blocks', 
			'cigar',
			'read',
			)
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



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def infer_segments(row):
	if len(row['blocks']) == 1:
		return row['blocks']
	
	segs = []
	current_block = row['blocks'][0]
	for blocks_index in range(1, len(row['blocks'])):
		next_block = row['blocks'][blocks_index]
		cigar_index = 2*blocks_index - 1
		gap_type = row['cigar'][cigar_index][0]

		if gap_type == 'N':
			segs.append(current_block)
			current_block = next_block
		elif gap_type in ('I', 'D'):
			current_block = (current_block[0], next_block[1])
		else:
			raise RuntimeError(f"{row['blocks']}\n{row['cigar']}\n{blocks_index}")
	
	segs.append(current_block)
		
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
						align.reference_name, 
						align.is_reverse,
						align.get_blocks(),
						align.cigartuples,
						align.is_read1,
						])

	linear_reads = pd.DataFrame(linear_reads, columns=COLUMNS)

	# Fix starting columns
	linear_reads.read_id = linear_reads.read_id.astype('string')
	linear_reads.chrom = linear_reads.chrom.astype('category')
	linear_reads.cigar = linear_reads.cigar.transform(lambda cigar: tuple((CIGARTUPLE_CODES[op],length) for op, length in cigar))
	linear_reads.read = linear_reads.read.map({True: '1', False: '2'}).astype('category')

	return linear_reads


def classify_seg(row):
	if row['Intergenic'] is True:
		return 'Intergenic'
	
	if len(row['exons'])==0 and len(row['introns'])==0:
		return 'Ambiguous'
	
	if len(row['exons'])==0:
		if row['align_is_reverse'] is False and any(row['seg'].begin==intron.begin for intron in row['introns'] if intron.data['strand']=='+'):
			return "Starts at 5'ss"
		elif row['align_is_reverse'] is True and any(row['seg'].end==intron.end for intron in row['introns'] if intron.data['strand']=='-'):
			return "Starts at 5'ss"
		else:
			return 'Intronic'
	
	#TODO: Maybe move this inside next if and add after pre-mRNA check since segment could cross intron-exon junction AND get spliced at another exon junction
	if any([row['seg'].begin==exon.begin for exon in row['exons']]) or any([row['seg'].end==exon.end for exon in row['exons']]):
		return 'Exon junction'

	if len(row['introns'])==0:
		return 'Exonic'
	
	# At this point there's at least 1 filtered exon and 1 filtered intron
	# We need to see if it's pre-mRNA
	for exon, intron in it.product(row['exons'], row['introns']):
		if len(exon.data['gene_id'].intersection(intron.data['gene_id']))==0:
			continue
		
		# Check if we have an exon-to-intron junction
		if exon.end-1 == intron.begin:
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
		elif intron.end-1 == exon.begin:
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
		
	return 'Ambiguous'
	

def classify_read(seg_rows:pd.DataFrame) -> str:
	classes_set = set(seg_rows['seg_class'])
	if len(classes_set)==1:
		if classes_set == {'Exon junction'}:
			if seg_rows.spliced.any():
				return 'mRNA'
			else:
				return 'Exonic'
		else:
			return classes_set.pop()
		
	# Now we know it's 2 classes at minimum
	if 'Intergenic' in classes_set:
		return 'Ambiguous'
	
	if classes_set in (
					{'Exon junction', 'Exonic',},
					{'Exon junction', 'Exonic', 'Intronic'},
					{'Exon junction', 'Ambiguous'},
					{'Exon junction', 'Exonic', 'Ambiguous'},
					{'Exon junction', 'Exonic', 'Intronic', 'Ambiguous'},
					):
		return 'mRNA'
	
	if classes_set in (
					{'pre-mRNA', 'Exonic'},
					{'pre-mRNA', 'Intronic'},
					{'pre-mRNA', 'Exonic', 'Intronic'},
					{'pre-mRNA', 'Ambiguous'},
					{'pre-mRNA', 'Exonic', 'Ambiguous'},
					{'pre-mRNA', 'Intronic', 'Ambiguous'},
					{'pre-mRNA', 'Exonic', 'Intronic', 'Ambiguous'},
					):
		return 'pre-mRNA'
	
	if classes_set in (
					{"Starts at 5'ss", 'Intronic'},
					{"Starts at 5'ss", 'Ambiguous'},
					{"Starts at 5'ss", 'Intronic', 'Ambiguous'},
					):
		return "Starts at 5'ss"
	
	return tuple(classes_set)



# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	# Get args
	ref_exons, ref_introns, output_base, log_level = sys.argv[1:]

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

	# Infer spliced/unspliced and explode read alignments into alignment segments
	linear_reads['spliced'] = linear_reads.blocks.transform(lambda blocks: len(blocks)>1)
	linear_reads['segs'] = linear_reads.apply(infer_segments, axis=1, result_type='reduce')
	linear_reads = linear_reads.explode('segs', ignore_index=True)
	linear_reads['seg'] = linear_reads.segs.transform(lambda segs: Interval(*segs))

	# Chromosome by chromosome, add all exons and introns 
	for chrom in linear_reads.chrom.unique():
		if chrom not in exons.keys():
			log.debug(f'No exons in {chrom}')
			continue
		chrom_exons = exons[chrom]
		linear_reads.loc[linear_reads.chrom==chrom, 'exons'] = linear_reads.loc[linear_reads.chrom==chrom, 'seg'].transform(chrom_exons.overlap)
		
		if chrom not in introns.keys():
			log.debug(f'No introns in {chrom}')
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

	# Classify reads based on their segment(s) class(es)
	log.debug(f'segment class counts: {linear_reads.seg_class.value_counts().sort_index().to_dict()}')
	linear_reads['read_class'] = linear_reads.read_id.map(linear_reads.groupby('read_id').apply(classify_read, include_groups=False))

	# Designate any unclassifiable combo of seg classes as 'Ambiguous'
	log.debug(f'read class counts: {linear_reads.read_class.astype("str").value_counts().sort_index().to_dict()}')
	linear_reads.read_class = linear_reads.read_class.transform(lambda rc: 'Ambiguous' if isinstance(rc, tuple) else rc)

	# Collapse segments back into one row per read and prepare to write to file 
	linear_reads = linear_reads.groupby(['read_id', 'read_class'], as_index=False).agg({'spliced': any})
	linear_reads['stage_reached'] = 'Linear mapping'
	linear_reads = linear_reads[['read_id', 'read_class', 'stage_reached', 'spliced']]

	# Write to file
	linear_reads.to_csv(f'{output_base}read_classes.tsv', sep='\t', index=False)

	log.debug('End of script')