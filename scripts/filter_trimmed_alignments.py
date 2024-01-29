import sys
import os
import gzip
from intervaltree import Interval, IntervalTree
from collections import Counter
from subprocess import run
from random import sample
from os.path import join
import pysam

import pandas as pd



# =============================================================================#
#                                  Constants                                   #
# =============================================================================#
MAX_MISMATCHES = 5
MAX_MISMATCH_PERCENT = 0.1



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def parse_gtf_line(line:str) -> dict:
	'''
	Returns {'Chromosome': str, 'Source': str, 'Feature': str, 'Start': int, 'End': int, 'Score': str, 'Strand': str, 'Frame': str,
			[attribute 1]: [attribute 1 value], [attribute 2]: [attribute 2 value]...}
	'''
	info = {}
	chrom, source, feature, start, end, score, strand, frame, attributes = line.strip().split('\t')

	for key, val in (('Chromosome', chrom),
					('Source', source),
					('Feature', feature),
					('Start', int(start)),
					('End', int(end)),
					('Score', score),
					('Strand', strand),
					('Frame', frame)
					):
		info[key] = val

	attributes = attributes.split('; ')
	attribute_dir = {attr.split(' ')[0]: attr.split(' ')[1].strip('"') for attr in attributes if not attr.startswith('tag "')}
	tags = [attr.lstrip('tag "').rstrip('"') for attr in attributes if attr.startswith('tag "')]
	attribute_dir['tags'] = tuple(tags)

	info.update(attribute_dir)

	return info



def get_annotations(gtf_path:str, number = 'all', chromosome : str = None, feature : str = None, start : str = None, end : str = None, strand : str = None, none_except : bool = True) -> list:
	'''
	number = 'all' or an int >= 1
	'''
	assert isinstance(number, (str, int)), f'{number}, {type(number)}'
	assert number == 'all' or number >= 1, number
	if [chromosome, feature, start, end, strand] == [None]*5:
		raise RuntimeError(f'You must specify at least one argument from chromosome, start, end, strand, or feature')

	annotations = []
	with open(gtf_path, 'r') as gtf:
		for line in gtf:
			if line.startswith('##'):
				continue

			line_info = parse_gtf_line(line)
			
			if feature is not None and line_info['Feature'] != feature:
				continue
			if start is not None and line_info['Start'] != start:
				continue
			if end is not None and line_info['Snd'] != end:
				continue
			if strand is not None and line_info['Strand'] != strand:
				continue

			annotations.append(line_info)
			if number != 'all' and len(annotations) == number:
				break
	
	if annotations == [] and none_except:
		raise RuntimeError(f'Annotations not found in {gtf_path} matching args "{chromosome} {feature} {start} {end} {strand}"')

	return annotations



def load_gene_ranges(gtf_file):
	'''
	Read through the GTF annotation file, return a dict with all annotated gene coordinates grouped by chromosome, then strand
	Output format = { chromosome: {strand: IntervalTree(gene start position, gene end position, gene name) } }
	'''
	# open file
	if gtf_file[-2:] == 'gz':
		in_file = gzip.open(gtf_file, 'rt')
	else:
		in_file = open(gtf_file)

	gene_ranges = {}		
	for line in in_file:
		if line[0] != '#':
			chrom, _, feat, start, end, _, strand, _, attributes = line.strip().split('\t')
			if feat == 'gene':
				# Convert attributes string to dict of {tag name: value}
				attributes = attributes[:-1].split('; ')
				attributes = {a.split(' ')[0]:a.split(' ')[1].replace('\"', '') for a in attributes}
				if chrom not in gene_ranges:
					gene_ranges[chrom] = {s:IntervalTree() for s in ('-', '+')}
				gene_ranges[chrom][strand].add(Interval(int(start)-1, int(end), {'gene_name':attributes['gene_name']}))

	in_file.close()
	return gene_ranges



def load_fivep_info(fivep_info_table:str, annotation_gtf:str):
	'''
	Read through TSV file of filtered read alignments, return a dict
	Output format = { read id: [read sequence, 5'ss sequence, coordinates of all aligned 5'ss's that passed filtering, sequence alignment is reverse-complementary, start of alignment in read, end of alignment in read, read passed 5'ss filtering, 5'ss filter that was failed] }
	'''
	gene_ranges = load_gene_ranges(annotation_gtf)

	fivep_info = {}
	with open(fivep_info_table) as in_file:
		in_file.readline()
		for line in in_file:
			rid, read_seq, fivep_seq, fivep_sites, read_is_reverse, fivep_start, fivep_end = line.strip().split('\t')
			read_is_reverse = True if read_is_reverse == 'True' else False

			# Unpack 5'ss sites and add gene ids
			fivep_sites = fivep_sites.split(',')
			for i, site in enumerate(fivep_sites):
				chrom, start, end, strand = site.split(';')
				start = int(start)
				end = int(end)
				gene_ids = [g.data['ensembl_id'] for g in gene_ranges[chrom][strand].overlap(start, start+1)]
				fivep_sites[i] = (chrom, start, end, strand, gene_ids)
			
			fivep_info[rid] = (read_seq, fivep_seq, fivep_sites, read_is_reverse, int(fivep_start), int(fivep_end))

	return fivep_info



def load_transcriptome_alignments(trimmed_to_transcriptome_bam) -> dict:
	'''
	{read id: [(trimmed read sequence, transcript id, alignment start position, alignment end position, alignment is reverse-complementary, alignment quality, number of mismatches, [gap length,...]),...]}
	'''
	alignments = {}
	for align in pysam.AlignmentFile(trimmed_to_transcriptome_bam, 'rb'):	

		for tag, val in align.get_tags():
			if tag == 'XM':
				mismatches = val
				break

		gap_lengths = []
		for tag_num, length in align.cigartuples:
			if tag_num in (1, 2):
				gap_lengths.append(length)
			assert tag_num in (0, 1, 2), align

		if align.query_name not in alignments:
			alignments[align.query_name] = []
		alignments[align.query_name].append([align.query_seq, align.reference_name, align.reference_start, align.reference_end, align.is_reverse, align.mapping_quality, mismatches, gap_lengths])
	
	return alignments



def load_exon_overlaps(exon_overlaps_bed) -> list:
	exon_overlaps = []
	with open(exon_overlaps_bed) as r:
		r.readline()
		for line in r:
			transcript_id, align_start, align_end, rid, mapping_quality, align_orient, _, _, _, _, _, _, _, exon_start, exon_end = line.strip().split('\t')
			exon_overlaps.append([rid, transcript_id, int(align_start), int(align_end), int(exon_start), int(exon_end)])

	return exon_overlaps

			

def add_info_to_alignments(alignments:dict, exon_overlaps:list, transcript_info:dict) -> dict:
	'''
	Add annotation info and exons overlapping the trimmed read alignment
	'''
	out_alignments = {}
	for rid in alignments:
		out_alignments[rid] = []

		for alignment in alignments[rid]:
			_, transcript_id, align_start, align_end, _, _, _, _ = alignment
			
			overlaps = [x for x in exon_overlaps if x[0]==rid and x[1]==transcript_id and x[2]==align_start and x[3]==align_end]
			transcript = transcript_info[transcript_id]

			if transcript['Strand']=='+':
				align_start_genomic = transcript['Start'] + align_start
				align_end_genomic = transcript['Start'] + align_end
			else:
				align_start_genomic = transcript['End'] - align_start
				align_end_genomic = transcript['End'] - align_end
			
			out = [*alignment, 
					align_start_genomic, 
					align_end_genomic, 
					transcript['Chromosome'],  
					transcript['Strand'], 
					transcript['Start'], 
					transcript['End'], 
					transcript['gene_id'], 
					transcript['gene_name'],
					transcript['gene_type'],
					overlaps]
			out_alignments[rid].append(out)

	return out_alignments



def filter_alignments(alignments:dict, fivep_info:dict, output_base:str) -> tuple:
	'''
	
	'''
	filtered_alignments = []
	failed_alignments = []
	for rid in alignments:
		read_seq, fivep_seq, fivep_sites, read_is_reverse, fivep_start, fivep_end = fivep_info[rid]

		# Get best mapping quality score
		max_score = max([align[4] for align in alignments[rid]])
	
		for align in alignments[rid]:
			trimmed_seq, transcript_id, align_start, align_end, align_is_reverse, mapping_quality, mismatches, gap_lengths, align_start_genomic, align_end_genomic, transcript_chrom, transcript_strand, transcript_start, transcript_end, transcript_gene_id, transcript_gene_name, transcript_gene_type, exon_overlaps = align
			passed_filtering = True
			fail_reason = None

			# Check if mapping quality is not max
			if mapping_quality < max_score:
				passed_filtering = False
				fail_reason = 'suboptimal_map_quality' if fail_reason is None else fail_reason

			# Check if too many mismatches
			if mismatches > MAX_MISMATCHES or mismatches/len(trimmed_seq) > MAX_MISMATCH_PERCENT:
				passed_filtering = False
				fail_reason = 'mismatches' if fail_reason is None else fail_reason

			# Check if more than one insertion/deletion OR one insertion/deletion that's longer than 3bp
			if len(gap_lengths) > 1 or (len(gap_lengths) == 1 and gap_lengths[0] > 3):
				passed_filtering = False
				fail_reason = 'indel' if fail_reason is None else fail_reason

			# Check if the 5'ss alignment orientation doesn't match the 3'ss's orientation
			if read_is_reverse != align_is_reverse:
				passed_filtering = False
				fail_reason = 'same_orientation' if fail_reason is None else fail_reason

			# Check if the alignment overlaps any exons
			if exon_overlaps != []:
				passed_filtering = False
				fail_reason = 'exon_overlap' if fail_reason is None else fail_reason
			
			# Pair the transcript alignment to a 5'ss alignment
			matching_fivep_sites = []
			for chrom, start, end, strand, fivep_gene_ids in fivep_sites:
				if transcript_gene_id in fivep_gene_ids:
					matching_fivep_sites.append((chrom, start, end, strand))

			if matching_fivep_sites == []:
				passed_filtering = False
				fail_reason = 'no_fivep_matches' if fail_reason is None else fail_reason
				failed_alignments.append([*align, fail_reason])
				continue
			elif len(matching_fivep_sites) > 1:
				passed_filtering = False
				fail_reason = 'multiple_fivep_matches' if fail_reason is None else fail_reason
				failed_alignments.append([*align, fail_reason])
				continue

			fivep_site = matching_fivep_sites[0][1]

			# Infer branchpoint coordinate
			if transcript_strand == '+':
				bp_site = transcript_start + align_end
			else:
				bp_site = transcript_end - align_end

			# Check if the 5'ss is at the same coordinate or downstream of the branchpoint
			if (transcript_strand == '+' and fivep_site >= bp_site) or (transcript_strand == '-' and fivep_site <= bp_site):
				passed_filtering = False
				fail_reason = '5p_bp_order' if fail_reason is None else fail_reason

			# Add alignment to dict
			if passed_filtering is True: 
				filtered_alignments.append([rid, read_seq, transcript_chrom, transcript_strand, transcript_gene_name, transcript_gene_id, transcript_gene_type, fivep_site, read_is_reverse, fivep_start, fivep_end, bp_site, trimmed_seq[-1]])
			else:
				failed_alignments.append((rid, read_seq, transcript_chrom, transcript_strand, transcript_gene_name, transcript_gene_id, transcript_gene_type, fivep_site, read_is_reverse, fivep_start, fivep_end, bp_site, trimmed_seq[-1], fail_reason))



def get_intron_info(ref_introns) -> tuple:
    '''
    Returns a dict formatted as follows:
    {Chromosome: {Strand(+ or -): Intervaltree(StartPosition(int), EndPosition(int))}}
    for 3' splice sites (+/- 2 bases), 5' splice sites (+/- 2 bases), and the introns they come from (start to end)
    '''
    threep_sites, fivep_sites, introns = {}, {}, {}

    if ref_introns[-2:] == 'gz':
        intron_file = gzip.open(ref_introns, 'rt')
    else:
        intron_file = open(ref_introns)

    introns_done = set()
    for line in intron_file:
        chrom, start, end, _, _, strand = line.strip().split()
        if chrom not in threep_sites:
            threep_sites[chrom] = {s: IntervalTree() for s in ['+', '-']}
            fivep_sites[chrom] = {s: IntervalTree() for s in ['+', '-']}
            introns[chrom] = {s: IntervalTree() for s in ['+', '-']}
        intron_id = '{}_{}_{}_{}'.format(chrom, strand, start, end)
        if intron_id not in introns_done:
            start, end = int(start), int(end)
            if strand == '+':
                threep_sites[chrom][strand].add(Interval(end-2, end+2))
                fivep_sites[chrom][strand].add(Interval(start-2, start+2))
            else:
                threep_sites[chrom][strand].add(Interval(start-2, start+2))
                fivep_sites[chrom][strand].add(Interval(end-2, end+2))
            introns[chrom][strand].add(Interval(start, end))
            introns_done.add(intron_id)

    intron_file.close()

    return threep_sites, fivep_sites, introns



def add_more_info_to_alignments(filtered_alignments:list, introns:dict, genome_fasta:str, output_base:str):
	out_alignments = []
	for alignment in filtered_alignments:
		rid, read_seq, transcript_chrom, transcript_strand, fivep_site, read_is_reverse, fivep_start, fivep_end, bp_site, read_bp_nt = alignment

		# Get the closest 3'ss in an intron that overlaps the bp
		overlap_introns = list(introns[transcript_chrom][transcript_strand].overlap(bp_site, bp_site+1))
		if transcript_strand == '+':
			threep_site = min(overlap_introns, key=lambda s: s.end-bp_site).end
		else:
			threep_site = min(overlap_introns, key=lambda s: bp_site-s.begin).begin
		bp_dist_to_threep = bp_site-threep_site if transcript_strand == '+' else threep_site-bp_site

		# Get BP sequence
		temp_bp_bed, temp_bp_seq = output_base+'_temp_bp_seqs.bed', output_base+'_temp_bp_seqs.txt'
		temp_file = open(temp_bp_bed, 'w')
		if transcript_strand == '+':
			bp_start, bp_end = bp_site-4, bp_site+6
		else:
			bp_start, bp_end = bp_site-5, bp_site+5
		temp_file.write(f'{transcript_chrom}\t{bp_start}\t{bp_end}\t{transcript_chrom};{bp_site};{transcript_strand}\t0\t{transcript_strand}\n')
		temp_file.close()
		run(f'bedtools getfasta -fi {genome_fasta} -bed {temp_bp_bed} -fo {temp_bp_seq} -nameOnly -s -tab'.split(' '))

		temp_file = open(temp_bp_seq)
		name, genomic_bp_window = temp_file.readline().strip().split()
		genomic_bp_window = genomic_bp_window.upper()
		genomic_bp_nt = genomic_bp_window[4]
		temp_file.close()

		out_alignments.append([*alignment, genomic_bp_nt, genomic_bp_window, threep_site, bp_dist_to_threep])

	os.remove(temp_bp_bed)
	os.remove(temp_bp_seq)
	return out_alignments



# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	trimmed_to_transcriptome_bam, fivep_info_table, exon_overlaps_bed, annotation_gtf, ref_introns, genome_fasta, output_base = sys.argv[1:]

	# Load data
	fivep_info = load_fivep_info(fivep_info_table, annotation_gtf)
	alignments = load_transcriptome_alignments(trimmed_to_transcriptome_bam)
	exon_overlaps = load_exon_overlaps(exon_overlaps_bed)
	transcript_info = get_annotations(annotation_gtf, feature='transcript')
	transcript_info = {info['transcript_id']: info for info in transcript_info}
	threep_sites, fivep_sites, introns = get_intron_info(ref_introns)

	# 
	mapped_rids = set([x[:-4] for x in alignments.keys()])
	with open(f'{output_base}_run_data.tsv', 'a') as a:
		a.write(f'transcriptome_mapped_reads\t{len(mapped_rids)}\n')

	alignments = add_info_to_alignments(alignments, exon_overlaps, transcript_info)

	# Filter trimmed read alignments
	filtered_alignments, failed_alignments = filter_alignments(alignments, fivep_info, output_base)

	# Write failed alignments to file
	with open(f'{output_base}_failed_transcriptome_alignments.tsv', 'w') as w:
		w.write('read_id\tread_seq\tchrom\tstrand\tfivep_site\tread_is_reverse\tfivep_start\tfivep_end\tbp_site\tread_bp_nt\n')
		for alignment in failed_alignments:
			out_row = '\t'.join([str(x) for x in alignments])
			w.write(out_row + '\n')

	# Add BP sequence info and a 3'ss 
	filtered_alignments = add_more_info_to_alignments(filtered_alignments, introns, genome_fasta, output_base)

	# Write filtered alignments out
	with open(f'{output_base}_final_info_table.tsv', 'w') as w:
		w.write('read_id\tread_seq\tchrom\tstrand\tgene_name\tgene_id\tgene_type\tfivep_site\tread_is_reverse\tfivep_start\tfivep_end\tbp_site\tread_bp_nt\tgenomic_bp_nt\tgenomic_bp_window\tthreep_site\tbp_dist_to_threep\n')
		for alignment in filtered_alignments:
			out_row = '\t'.join([str(x) for x in alignments])
			w.write(out_row + '\n')