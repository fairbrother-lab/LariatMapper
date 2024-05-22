
# coding: utf-8

# In[34]:

from dataclasses import dataclass, field
import os, sys, subprocess
import copy
import statistics
import math
import time
import datetime as dt
import itertools as it
import random
import string
import gzip

import pysam
import numpy as np
import pandas as pd
import pyfaidx
from intervaltree import Interval, IntervalTree


# In[42]:

##### Constants #####

LINEAR_COLUMNS = ('read_id', 'chrom', 'start', 'end', 'blocks', 'read',)
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


# In[36]:

##### Functions #####
def parse_attributes(attribute_string:str, anno_type:str) -> dict:
	if anno_type == 'gtf':
		attributes = attribute_string.rstrip('";').split('; ')
		attributes = [attr.split(' ') for attr in attributes]
		tags = [attr_val.strip('"') for attr_name, attr_val in attributes if attr_name=='tag']
		attributes = {attr_name: attr_val.strip('"') for attr_name, attr_val in attributes if attr_name!='tag'}
		attributes['tags'] = tags
	else:
		attributes = [attr.split('=') for attr in attribute_string.split(';')]
		attributes = [(attr[0].lstrip(), attr[1]) for attr in attributes]
		attributes = dict(attributes)
		if 'tag' in attributes:
			attributes['tags'] = attributes['tag'].split(',')

	return attributes


def parse_gene_info(ref_anno:str) -> dict:
	prev_ext, last_ext = ref_anno.split('.')[-2:]
	if last_ext == 'gz':
		in_file, anno_type = gzip.open(ref_anno, 'rt'), prev_ext
	else:
		in_file, anno_type = open(ref_anno), last_ext

	genes = {}
	for line in in_file:
		if line[0] != '#':
			_, _, feature, _, _, _, _, _, attributes = line.strip().split('\t')
			if feature == 'transcript':
				attributes = parse_attributes(attributes, anno_type)
				genes[attributes['gene_id']] = [attributes[e] for e in ['gene_name', 'gene_type']]
	in_file.close()
	
	return genes


def parse_introns(ref_introns:str) -> tuple:
	'''
	Returns a dict formatted as follows:
	{Chromosome: {Strand(+ or -): Intervaltree(StartPosition(int), EndPosition(int), {"gene_id": GeneID})}}
	'''

	if ref_introns.split('.')[-1] == 'gz':
		intron_file = gzip.open(ref_introns, 'rt')
	else:
		intron_file = open(ref_introns)

	introns, introns_done = {}, set()
	for line in intron_file:
		chrom, start, end, intron_info, _, strand = line.strip().split('\t')
		if chrom not in introns:
			introns[chrom] = IntervalTree()
		intron_id = '{}_{}_{}_{}'.format(chrom, strand, start, end)
		if intron_id not in introns_done:
			start, end = int(start), int(end)
			if end-start > 20:
				intron_genes = [i.split('-')[0] for i in intron_info.split(';')[-1].split('|')]
				introns[chrom].add(Interval(start, end, {'strand': strand, 'gene_id':intron_genes}))
				introns_done.add(intron_id)

	intron_file.close()
	return introns


def parse_exons(ref_exons:str) -> tuple:
	'''
	Returns a dict formatted as follows:
	{Chromosome: {Strand(+ or -): Intervaltree(StartPosition(int), EndPosition(int), {"gene_id": GeneID})}}
	'''
	if ref_exons.split('.')[-1] == 'gz':
		exon_file = gzip.open(ref_exons, 'rt')
	else:
		exon_file = open(ref_exons)

	exons, exons_done = {}, set()
	for line in exon_file:
		chrom, start, end, exon_info, _, strand = line.strip().split('\t')
		if chrom not in exons:
			exons[chrom] = IntervalTree()
		exon_id = '{}_{}_{}_{}'.format(chrom, strand, start, end)
		if exon_id not in exons_done:
			start, end = int(start), int(end)
			if end-start > 20:
				exon_genes = [i.split('-')[0] for i in exon_info.split(';')[-1].split('|')]
				exons[chrom].add(Interval(start, end, {'strand': strand, 'gene_id':exon_genes}))
				exons_done.add(exon_id)

	exon_file.close()
	return exons


def tree_covers_interval(tree:IntervalTree, interval:Interval) -> bool:
	total_coverage = False
	merged_tree = tree.copy()
	merged_tree.merge_overlaps(strict=False)
	for merged_interval in merged_tree:
		if merged_interval.contains_interval(interval):
			total_coverage = True
	
	return total_coverage


# In[37]:

##### Settings #####
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_rows', 100)

ref_anno='/datasets2/lariat_mapping/reference_files/human/ref_dir_a/annotation.gtf'
ref_introns='/datasets2/lariat_mapping/reference_files/human/ref_dir_a/introns.bed'
ref_exons='/datasets2/lariat_mapping/reference_files/human/ref_dir_a/exons.bed'
# output_base='/datasets2/lariat_mapping/testing/output/C22-1/'
# output_base='/datasets2/lariat_mapping/testing/output/HEK293T-1/'
# output_base='/datasets2/lariat_mapping/testing/output/C22-1_10m/'
output_base='/datasets2/lariat_mapping/testing/output/HEK293T-1_10m/'
# output_base='/datasets2/lariat_mapping/testing/100k_truncs/C22-1/'
# output_base='/datasets2/lariat_mapping/testing/100k_truncs/HEK293T-1/'
print(output_base)


# In[38]:

genes = parse_gene_info(ref_anno)
print(list(genes.items())[0])
introns = parse_introns(ref_introns)
print(list(introns['chr1'])[:3])
exons = parse_exons(ref_exons)
print(list(exons['chr1'])[:3])


# In[39]:

# Add gene types to introns
for chrom in introns.keys():
	for intron in introns[chrom]:
		gene_info = [genes[gene_id] for gene_id in intron.data['gene_id']]
		gene_types = set([type_ for name, type_ in gene_info])
		intron.data['gene_type'] = gene_types
print(list(introns['chr1'])[:3])

# Add gene types to exons
for chrom in exons.keys():
	for exon in exons[chrom]:
		gene_info = [genes[gene_id] for gene_id in exon.data['gene_id']]
		gene_types = set([type_ for name, type_ in gene_info])
		exon.data['gene_type'] = gene_types
print(list(exons['chr1'])[:3])


# In[65]:

linear_reads = []
for align in pysam.AlignmentFile(f'{output_base}/mapped_reads.bam', 'rb'):
	if align.is_unmapped:
		continue
	linear_reads.append([
					align.query_name, 
					align.reference_name, 
					align.reference_start,
					align.reference_end, 
					align.get_blocks(),
					# align.cigartuples,
					# align.get_tags(),
					align.is_read1,
					])

linear_reads = pd.DataFrame(linear_reads, columns=LINEAR_COLUMNS)
# linear_reads.cigar = linear_reads.cigar.transform(lambda cigar: tuple((CIGARTUPLE_CODES[op],length) for op, length in cigar))
# linear_reads.tags = linear_reads.tags.transform(lambda tags: {tag: val for tag, val in tags})
linear_reads.read = linear_reads.read.map({True: '1', False: '2'})
linear_reads = linear_reads.sort_values(['read_id'])
linear_reads


# In[44]:

# #TODO: finish
# def infer_segments(row):
# 	if row['indels'] == 0:
# 		return row['blocks']
	
# 	segs = []
# 	for i in range(len(row['blocks'])):
# 		row['blocks'][i]

# 	return segs
# linear_reads['segs'] = linear_reads.apply(infer_segments, axis=1, result_type='reduce')

linear_reads['segs'] = linear_reads.blocks


# In[45]:

linear_reads['read_id_f'] = linear_reads.read_id + '/' + linear_reads.read
assert linear_reads.read_id_f.is_unique
# linear_reads['indels'] = linear_reads.cigar.transform(lambda cigar: sum(length for op, length in cigar if op in ('I', 'D')))
# linear_reads['mismatches'] = linear_reads.tags.transform(lambda tags: tags['XM'])
# linear_reads['edit_dist'] = linear_reads.tags.transform(lambda tags: tags['NM'])
# print(linear_reads.indels.value_counts())
# print(linear_reads.mismatches.value_counts())
# print(linear_reads.edit_dist.value_counts())
linear_reads['n_segs'] = linear_reads.segs.transform(len)
linear_reads['spliced'] = linear_reads.n_segs>1
print(linear_reads.n_segs.value_counts())
print(linear_reads.spliced.value_counts())
linear_reads


# In[46]:

linear_reads = linear_reads.explode('segs').rename(columns={'segs': 'seg'}).reset_index(drop=True)
linear_reads.seg = linear_reads.seg.transform(lambda seg: Interval(*seg))
linear_reads


# In[47]:

def add_exons(row, exons):
	overlap_exons = exons[row['chrom']].overlap(row['seg'].begin, row['seg'].end)
	# if any('protein_coding' in exon.data['gene_type'] for exon in overlap_exons):
	overlap_exons = {exon for exon in overlap_exons if 'protein_coding' in exon.data['gene_type']}

	return IntervalTree(overlap_exons)


def add_introns(row, introns):
	if row['chrom'] not in introns.keys():
		return IntervalTree()

	overlap_introns = introns[row['chrom']].overlap(row['seg'].begin, row['seg'].end)
	# if any('protein_coding' in intron.data['gene_type'] for intron in overlap_introns):
	overlap_introns = {intron for intron in overlap_introns if 'protein_coding' in intron.data['gene_type']}

	return IntervalTree(overlap_introns)

linear_reads['exons'] = linear_reads.apply(add_exons, exons=exons, axis=1)
linear_reads['introns'] = linear_reads.apply(add_introns, introns=introns, axis=1)
linear_reads


# In[48]:

# def pull_gene_types(row):
# 	out = set()
# 	for exon in row['exons']:
# 		out = out.union(exon.data['gene_type'])
# 	for intron in row['introns']:
# 		out = out.union(intron.data['gene_type'])
		
# 	return out

# linear_reads['gene_types'] = linear_reads.apply(pull_gene_types, axis=1) 
# print(linear_reads.gene_types.value_counts())


# In[49]:

def classify_seg(row):
	if any(exon.begin==row['seg'].begin for exon in row['exons']) or any(exon.end==row['seg'].end for exon in row['exons']):
		return 'Exon junction'
		
	if len(row['exons'])==0 and len(row['introns'])==0:
		return 'Inter-genic'

	if len(row['exons'])==0 and len(row['introns'])>0:
		#TODO: Check for start at 5'ss
		return 'Intronic'

	if len(row['exons'])>0 and len(row['introns'])==0:
		return 'Exonic'

	# We can deduce that there's at least 1 exon and 1 intron at this point
	if tree_covers_interval(IntervalTree(row['exons']), row['seg']) is True:
		return 'Ambiguous'
	if tree_covers_interval(IntervalTree(row['introns']), row['seg']) is True:
		return 'Ambiguous'

	return 'pre-mRNA'

linear_reads['seg_class'] = linear_reads.apply(classify_seg, axis=1)
print(linear_reads.seg_class.value_counts())


# In[50]:

def classify_read(seg_classes):
	seg_classes_set = set(seg_classes)

	if seg_classes_set == {'Ambiguous', }:
		return 'Ambiguous'

	if len(seg_classes_set)>1 and 'Inter-genic' in seg_classes_set:
		return 'Weird mix'

	seg_classes_set.discard('Ambiguous')

	if seg_classes_set in ({'Exon junction', },
							{'Exon junction', 'Exonic'}
							):
		return 'mRNA'
	
	if seg_classes_set in (
						{'Intronic', 'Exonic'},
						{'Intronic', 'Exon junction'},
						{'Intronic', 'Exonic', 'Exon junction'},
						{'pre-mRNA', 'Intronic'},
						{'pre-mRNA', 'Exonic'},
						{'pre-mRNA', 'Exon junction'},
						{'pre-mRNA', 'Exon junction', 'Intronic'},
						{'pre-mRNA', 'Exon junction', 'Exonic'},
						{'pre-mRNA', 'Exonic', 'Intronic'},
						{'pre-mRNA', 'Exon junction', 'Intronic', 'Exonic'},
						):
		return 'pre-mRNA'

	if len(seg_classes_set) == 1:
		return seg_classes_set.pop()

	# if 'pre-mRNA' in seg_classes:
	# 	return 'pre-mRNA'

	# if 'Intronic' in seg_classes and ('Exonic' in seg_classes or 'Exon_junction' in seg_classes):
	# 	return 'pre-mRNA'

	return seg_classes_set


linear_reads['class_'] = linear_reads.read_id_f.map(linear_reads.groupby('read_id_f').seg_class.agg(classify_read))
print(linear_reads.class_.value_counts())


# In[51]:

def classify_pe(classes):
	classes_set = set(classes)

	if len(classes_set)==1:
		return classes_set.pop()

	if classes_set in (
					{'mRNA', 'Exonic'},
					{'mRNA', 'Ambiguous'},
					):
		return 'mRNA'

	if classes_set in (
					{'pre-mRNA', 'Exonic'},
					{'pre-mRNA', 'Intronic'},
					{'pre-mRNA', 'mRNA'},
					{'pre-mRNA', 'Ambiguous'},
					{'mRNA', 'Intronic'},
					{'Exonic', 'Intronic'},
					):
		return 'pre-mRNA'

	if 'Weird mix' in classes_set:
		return 'Weird mix'

	return tuple(classes_set)

print(linear_reads.groupby('read_id').class_.agg(classify_pe).value_counts().sum())
print(linear_reads.groupby('read_id').class_.agg(classify_pe).value_counts())


# In[52]:

print(linear_reads.read_id.nunique())


# In[ ]:




# In[ ]:




# In[53]:

lariat_rids = pd.read_csv(f'{output_base}lariat_reads.tsv', sep='\t').read_id
lariat_rids = set(lariat_rids)
print(len(lariat_rids))
lariat_rids


# In[54]:

circular_rids = pd.read_csv(f'{output_base}circularized_intron_reads.tsv', sep='\t').read_id
circular_rids = set(circular_rids)
print(len(circular_rids))
circular_rids


# In[55]:

temp_switch_rids = pd.read_csv(f'{output_base}template_switching_alignments.tsv', sep='\t').read_id
temp_switch_rids = set(temp_switch_rids)
print(len(temp_switch_rids))
temp_switch_rids


# In[56]:

assert len(lariat_rids.intersection(circular_rids)) == 0, lariat_rids.intersection(circular_rids)
assert len(lariat_rids.intersection(temp_switch_rids)) == 0, lariat_rids.intersection(temp_switch_rids)
assert len(temp_switch_rids.intersection(circular_rids)) == 0, temp_switch_rids.intersection(circular_rids)


# In[57]:

lariat_failed = pd.read_csv(f'{output_base}failed_lariat_alignments.tsv', sep='\t')
lariat_failed = lariat_failed[['read_id', 'filter_failed']]
repeat_rids = lariat_failed.loc[lariat_failed.filter_failed.isin(('in_repeat', 'ubiquitin_gene'))]
print(len(repeat_rids))
lariat_failed


# In[58]:

trim_failed = pd.read_csv(f'{output_base}failed_trimmed_alignments.tsv', sep='\t')
trim_failed.read_id = trim_failed.read_id.transform(lambda rid: rid[:-4].split('/')[0])
trim_failed = trim_failed[['read_id', 'filter_failed']].drop_duplicates(ignore_index=True)
trim_failed


# In[59]:

fivep_passed_rids = pd.read_csv(f'{output_base}fivep_info_table.tsv', sep='\t').read_id
fivep_passed_rids = set(fivep_passed_rids.transform(lambda rid: rid[:-4].split('/')[0]))
print(len(fivep_passed_rids))
fivep_passed_rids


# In[60]:

genome_unmapped_rids = pyfaidx.Fasta(f'{output_base}unmapped_reads.fa')
genome_unmapped_rids = set([read_id.split('/')[0] for read_id in genome_unmapped_rids.keys()])
print(len(genome_unmapped_rids))


# In[61]:

read_class = []
for read_id in genome_unmapped_rids:
	if read_id in lariat_rids:
		read_class.append((read_id, 'Lariat'))
	elif read_id in repeat_rids:
		read_class.append((read_id, 'Repeat region'))
	elif read_id in circular_rids:
		read_class.append((read_id, 'Circularized intron'))
	elif read_id in temp_switch_rids:
		read_class.append((read_id, 'Template-switching'))
	elif read_id in fivep_passed_rids:
		read_class.append((read_id, "Has 5'ss alignment"))
	else:
		read_class.append((read_id, 'Unmapped'))

read_class = pd.DataFrame(read_class, columns=('read_id', 'class_'))
print(read_class.class_.value_counts())
read_class


# In[62]:

read_class_counts = read_class.class_.value_counts().to_dict()
print(read_class_counts)


# In[ ]:




# In[ ]:




# In[ ]:




# ##### Code snippets

# In[63]:

# plot = (p9.ggplot(, p9.aes(x="", y="")) +
# 		# p9.geom_() +
# 		# p9.geom_() +
# 		# p9.scale_x_() +
# 		# p9.scale_y_() +
# 		# p9.scale_color_brewer(type='', palette=1, direction=1) +
# 		# p9.scale_fill_brewer(type='', palette=1, direction=1) +
# 		p9.labs(title="") +
# 		p9.theme()
# )
# plot.show()
# # plot.save(PLOTS_DIR + '/.png', dpi=500)


# In[ ]:




# ##### Discarded code

# In[64]:

get_ipython().run_cell_magic(u'script', u'false --no-raise-error', u'')

