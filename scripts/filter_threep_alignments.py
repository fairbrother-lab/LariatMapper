import sys
import gzip
from intervaltree import Interval, IntervalTree
from collections import Counter
from subprocess import run
from random import sample
from os.path import join

MAX_MISMATCHES = 5
MAX_MISMATCH_PERCENT = 0.1



def parse_gene_ranges(gtf_file):
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



def parse_threep_lengths(threep_lengths):
	'''
	Read through TSV file of 3'ss coordinates and lengths (max 250), return a dict 
	Output format = { (chromosome, strand, 3'ss position): length }
	'''
	lengths = {}
	with open(threep_lengths) as in_file:
		for line in in_file:
			chrom, threep_coord, strand, length = line.strip().split('\t')
			lengths[(chrom, strand, int(threep_coord))] = int(length)

	return lengths



def parse_read_info(read_info_table):
	'''
	Read through TSV file of filtered read alignments, return a dict
	Output format = { read id: [read sequence, 5'ss sequence, coordinates of all aligned 5'ss's that passed filtering, sequence alignment is reverse-complementary, start of alignment in read, end of alignment in read] }
	'''
	read_info = {}
	with open(read_info_table) as in_file:
		for line in in_file:
			items = line.strip().split('\t')
			if items[0] != 'read_id':
				read_info[items[0]] = items[1:]

	return read_info



def parse_cigar(cig_str):
	''' 
	Split cigar string into tuples of (bases, operator)
	'''
	cig_tuples, curr_str = [], ''
	for c in cig_str:
		if not c.isalpha():
			curr_str = curr_str + c
		else:
			cig_tuples.append((int(curr_str), c))
			curr_str = ''
	return cig_tuples



def filter_threep_reads(trimmed_reads_to_threep, threep_lengths, read_info, gene_ranges, genome_fasta, output_base):
	'''
	
	'''
	threep_info = {}	# { read id + "_for" or "_rev" : (chromosome, strand, 3'ss coordinate, alignment start in 3'ss sequence, alignment end in 3'ss sequence, alignment to 3'ss is reverse-complementary, read sequence, mapping quality) }
	with open(trimmed_reads_to_threep) as read_file:
		# Loop through lines in SAM file
		for line in read_file:
			alignment_info = line.strip().split('\t')
			rid, flag, threep_site, alignment_start, mapping_quality, read_cig, _, _, _, trimmed_read_seq = alignment_info[:10]
			
			# Parse the number of mismatches
			for alignment_tag in alignment_info[11:]:
				if alignment_tag[:2] == 'XM':
					num_mismatch = int(alignment_tag.split(':')[-1])
			
			read_cig = parse_cigar(read_cig)
			read_len = float(sum(c[0] for c in read_cig if c[1] in ('S', 'I', 'M')))

			# If too many mismatches, skip alignment
			if num_mismatch > MAX_MISMATCHES or num_mismatch/read_len > MAX_MISMATCH_PERCENT:
				continue

			# If more than one insertion/deletion OR one insertion/deletion that's longer than 3bp, skip alignment
			indels = [c[0] for c in read_cig if c[1] in ('D', 'I')]
			if len(indels) > 1 or (len(indels) == 1 and indels[0] > 3):
				continue

			# Prep information for addition to dict
			chrom, start, end, strand = threep_site[:-3].split(';')
			threep_coord = end if strand=='+' else start
			alignment_start = int(alignment_start)-1
			alignment_len = sum(c[0] for c in read_cig if c[1] in ('M', 'D', 'N'))
			alignment_end = alignment_start+alignment_len
			bit_flags = bin(int(flag))
			threep_is_reverse = True if len(bit_flags)>=7 and bit_flags[-5]=='1' else False

			# Add alignment to dict
			if rid not in threep_info:
				threep_info[rid] = []
			threep_info[rid].append((chrom, strand, threep_coord, alignment_start, alignment_end, threep_is_reverse, trimmed_read_seq, int(mapping_quality)))
	
	# Loop through reads with at least 1 3'ss alignment
	potential_alignments = {}			# { read id: [read sequence, chromosome, strand, 5'ss coordinate, alignment to genome is reverse-complementary, 5'ss start in read, 5'ss end in read, 3'ss coordinate, branchpoint coordinate, branchpoint base] }
	for rid in threep_info:
		read_seq, fivep_seq, fivep_sites, fivep_is_reverse, fivep_start, fivep_end = read_info[rid]
		fivep_is_reverse = True if fivep_is_reverse=='True' else False

		# Get the alignment(s) with the highest mapping quality score
		max_score = max(s[-1] for s in threep_info[rid])
		top_alignments = [s for s in threep_info[rid] if s[-1]==max_score]

		# Make a dict with all 5'ss's aligned to the read
		fivep_dict = {}			# { chromosome: { strand: set(5'ss coordinate) } }
		for fp in fivep_sites.split(','):
			chrom, start, end, strand = fp.split(';')
			if chrom not in fivep_dict:
				fivep_dict[chrom] = {s:set() for s in ['+', '-']}
			if strand == '+':
				fivep_dict[chrom][strand].add(int(start))
			else:
				fivep_dict[chrom][strand].add(int(end)-1)

		# Loop through top alignment(s) to 3'ss's AND
		for align_info in top_alignments:
			threep_chrom, threep_strand, threep_coord, alignment_start, alignment_end, threep_is_reverse, trimmed_read_seq, _ = align_info
			threep_coord, alignment_end = int(threep_coord), int(alignment_end)
			
			# If 3'ss chromosome and orientation doesn't match at least 1 5'ss, skip alignment
			if threep_chrom not in fivep_dict:
				continue
			if len(fivep_dict[threep_chrom][threep_strand]) == 0:
				continue
			
			# Collect the 5'ss's that share a gene annotation with the 3'ss
			threep_genes = [g.data['gene_name'] for g in gene_ranges[threep_chrom][threep_strand].overlap(threep_coord, threep_coord+1)]
			same_gene_fivep = []
			for fp_coord in fivep_dict[threep_chrom][threep_strand]:
				fivep_genes = [g.data['gene_name'] for g in gene_ranges[threep_chrom][threep_strand].overlap(fp_coord, fp_coord+1)]
				gene_matches = sum(1 for g in fivep_genes if g in threep_genes) > 0
				if gene_matches:
					same_gene_fivep.append(fp_coord)

			# If not exactly one 5'ss shares a gene with the 3'ss, skip alignment
			if len(same_gene_fivep) != 1:
				continue
			
			fivep_coord = same_gene_fivep[0]

			# If the 5'ss is at the same coordinate or downstream of the 3'ss, skip alignment
			if (threep_strand == '+' and fivep_coord >= threep_coord) or (threep_strand == '-' and fivep_coord <= threep_coord):
				continue

			# If the read alignment orientation doesn't match the 3'ss's orientation, skip alignment
			if fivep_is_reverse != threep_is_reverse:
				continue
			
			# Infer branchpoint coordinate
			if threep_strand == '+':
				bp_coord = alignment_end + threep_coord-threep_lengths[(threep_chrom, threep_strand, threep_coord)]-1
			else:
				bp_coord = (threep_lengths[(threep_chrom, threep_strand, threep_coord)]-alignment_end) + threep_coord

			# If the 5'ss is at the same coordinate or downstream of the branchpoint, skip alignment
			if (threep_strand == '+' and fivep_coord >= bp_coord) or (threep_strand == '-' and fivep_coord <= bp_coord):
				continue
				
			# Add alignment to dict
			read_bp_nt = trimmed_read_seq[-1]
			if rid not in potential_alignments:
				potential_alignments[rid] = []
			potential_alignments[rid].append([read_seq, threep_chrom, threep_strand, fivep_coord, fivep_is_reverse, fivep_start, fivep_end, threep_coord, bp_coord, read_bp_nt])
							
	# Loop through reads with potential alignments, checking whether or not the branchpoint sequence in the read matches the genomic sequence, and output 1 mapping for each read
	out_rids = [rid[:-4] for rid in potential_alignments]
	temp_bp_bed, temp_bp_seq = output_base+'temp_bp_seqs.bed', output_base+'temp_bp_seqs.txt'
	with open(output_base+'_final_info_table.txt', 'w') as out_file:
		out_file.write('read_id\tread_seq\tchrom\tstrand\tfivep_site\tread_is_reverse\tfivep_read_start\tfivep_read_end\t')
		out_file.write('threep_site\tbp_site\tread_bp_nt\tgenomic_bp_nt\tgenomic_bp_window\n')

		for rid in out_rids:

			align_mismatch = {'_for':{True:[], False:[]}, '_rev':{True:[], False:[]}}
			for fivep_dir in ['_for', '_rev']:
				if rid+fivep_dir not in potential_alignments:
					continue

				# Loop through alignments, checking branchpoint read sequence vs genomic sequence
				for align_info in potential_alignments[rid+fivep_dir]:
					read_seq, chrom, strand, fivep_coord, read_is_reverse, fivep_start, fivep_end, _, bp_coord, read_bp_nt = align_info
					temp_file = open(temp_bp_bed, 'w')
					if strand == '+':
						bp_start, bp_end = bp_coord-4, bp_coord+6
					else:
						bp_start, bp_end = bp_coord-5, bp_coord+5
					temp_file.write(f'{chrom}\t{bp_start}\t{bp_end}\t{chrom};{bp_coord};{strand}\t0\t{strand}\n')
					temp_file.close()
					run(f'bedtools getfasta -fi {genome_fasta} -bed {temp_bp_bed} -fo {temp_bp_seq} -nameOnly -s -tab'.split(' '))
					temp_file = open(temp_bp_seq)
					name, genomic_bp_window = temp_file.readline().strip().split()
					temp_file.close()
					genomic_bp_window = genomic_bp_window.upper()
					genomic_bp_nt = genomic_bp_window[4]
					align_mismatch[fivep_dir][genomic_bp_nt!=read_bp_nt].append(align_info+[genomic_bp_nt, genomic_bp_window])

			# Output 1 alignment for the read, prioritizing based on mistmatch and orientation, choosing a random alignment if multiple in same category
			if len(align_mismatch['_for'][True]) > 0:
				output = [rid] + sample(align_mismatch['_for'][True], 1)[0]
			elif len(align_mismatch['_rev'][True]) > 0:
				output = [rid] + sample(align_mismatch['_rev'][True], 1)[0]
			elif len(align_mismatch['_for'][False]) > 0:
				output = [rid] + sample(align_mismatch['_for'][False], 1)[0]
			else:
				output = [rid] + sample(align_mismatch['_rev'][False], 1)[0]

			out_file.write('\t'.join([str(e) for e in output]) + '\n')
	
	run(f'rm {temp_bp_bed} {temp_bp_seq}'.split(' '))


if __name__ == '__main__':
	trimmed_reads_to_threep, threep_lengths, fivep_info_table, gtf_file, genome_fasta, output_base = sys.argv[1:]

	# Load data
	gene_ranges = parse_gene_ranges(gtf_file)
	threep_lengths = parse_threep_lengths(threep_lengths)
	read_info = parse_read_info(fivep_info_table)

	# Filter trimmed read alignments to 3'ss's and write to final_info_table.txt
	filter_threep_reads(trimmed_reads_to_threep, threep_lengths, read_info, gene_ranges, genome_fasta, output_base)

