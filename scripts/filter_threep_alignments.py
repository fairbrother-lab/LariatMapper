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
	Output format = { read id: [read sequence, 5'ss sequence, coordinates of all aligned 5'ss's that passed filtering, sequence alignment is reverse-complementary, start of alignment in read, end of alignment in read, read passed 5'ss filtering, 5'ss filter that was failed] }
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



def filter_threep_reads(trimmed_reads_to_threep, threep_lengths, read_info, gene_ranges, genome_fasta, output_base, out_file, run_data):
	'''
	
	'''
	mapped_rids = set()
	threep_info = {}	# { read id + "_for" or "_rev" : (chromosome, strand, 3'ss coordinate, alignment start in 3'ss sequence, alignment end in 3'ss sequence, alignment to 3'ss is reverse-complementary, read sequence, mapping quality) }
	with open(trimmed_reads_to_threep) as sam_file:
		# Loop through lines in SAM file
		for line in sam_file:
			alignment_info = line.strip().split('\t')
			rid, flag, threep_site, alignment_start, mapping_quality, read_cig, _, _, _, trimmed_read_seq = alignment_info[:10]
			mapped_rids.add(rid[:-4])
			
			# Parse the number of mismatches
			for alignment_tag in alignment_info[11:]:
				if alignment_tag[:2] == 'XM':
					num_mismatch = int(alignment_tag.split(':')[-1])
			
			read_cig = parse_cigar(read_cig)

			# Prep information for addition to dict
			chrom, start, end, strand = threep_site[:-3].split(';')
			threep_coord = int(end) if strand=='+' else int(start)
			alignment_start = int(alignment_start)-1
			alignment_len = sum(c[0] for c in read_cig if c[1] in ('M', 'D', 'N'))
			alignment_end = alignment_start+alignment_len
			bit_flags = bin(int(flag))
			threep_is_reverse = True if len(bit_flags)>=7 and bit_flags[-5]=='1' else False
			mapping_quality = int(mapping_quality)

			# Add alignment to dict
			if rid not in threep_info:
				threep_info[rid] = []
			threep_info[rid].append((chrom, strand, threep_coord, alignment_start, alignment_end, threep_is_reverse, trimmed_read_seq, mapping_quality, num_mismatch, read_cig))
	
	with open(run_data, 'a') as a:
		a.write(f'threep_mapped_reads\t{len(mapped_rids)}\n')
	
	# Filter alignments, grouping by read
	filtered_alignments = {}			# { read id: [read sequence, chromosome, strand, 5'ss coordinate, alignment to genome is reverse-complementary, 5'ss start in read, 5'ss end in read, 3'ss coordinate, branchpoint coordinate, branchpoint base] }
	failed_alignments = []
	for rid in threep_info:
		read_seq, fivep_seq, fivep_sites, read_is_reverse, fivep_start, fivep_end = read_info[rid]
		read_is_reverse = True if read_is_reverse == 'True' else False

		# Get best mapping quality score
		max_score = max(s[-3] for s in threep_info[rid])

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

		# Assess top alignment(s) to 3'ss's
		for threep_chrom, threep_strand, threep_coord, alignment_start, alignment_end, threep_is_reverse, trimmed_read_seq, mapping_quality, num_mismatch, read_cig in threep_info[rid]:
			passed_filtering = True
			fail_reason = None

			# Check if mapping quality is not max
			if mapping_quality < max_score:
				passed_filtering = False
				fail_reason = 'suboptimal_map_quality' if fail_reason is None else fail_reason

			# Check if too many mismatches
			if num_mismatch > MAX_MISMATCHES or num_mismatch/len(trimmed_read_seq) > MAX_MISMATCH_PERCENT:
				passed_filtering = False
				fail_reason = 'mismatches' if fail_reason is None else fail_reason

			# Check if if more than one insertion/deletion OR one insertion/deletion that's longer than 3bp
			indels = [c[0] for c in read_cig if c[1] in ('D', 'I')]
			if len(indels) > 1 or (len(indels) == 1 and indels[0] > 3):
				passed_filtering = False
				fail_reason = 'indel' if fail_reason is None else fail_reason

			# Check if the read alignment orientation doesn't match the 3'ss's orientation
			if read_is_reverse != threep_is_reverse:
				passed_filtering = False
				fail_reason = 'same_orientation' if fail_reason is None else fail_reason

			# Check if 3'ss chromosome and orientation doesn't match at least 1 5'ss
			if threep_chrom not in fivep_dict or len(fivep_dict[threep_chrom][threep_strand]) == 0:
				passed_filtering = False
				fail_reason = 'matching_chrom_strand' if fail_reason is None else fail_reason
			
			# Collect the 5'ss's that share a gene annotation with the 3'ss
			threep_genes = [g.data['gene_name'] for g in gene_ranges[threep_chrom][threep_strand].overlap(threep_coord, threep_coord+1)]
			same_gene_fivep = []
			for fp_coord in fivep_dict.get(threep_chrom, {'+':[], '-':[]})[threep_strand]:
				fivep_genes = [g.data['gene_name'] for g in gene_ranges[threep_chrom][threep_strand].overlap(fp_coord, fp_coord+1)]
				gene_matches = sum(1 for g in fivep_genes if g in threep_genes) > 0
				if gene_matches:
					same_gene_fivep.append(fp_coord)

			# Check if not exactly one 5'ss shares a gene with the 3'ss
			# If not, we gotta skip the rest because they depend on a 5'ss-3'ss pair
			if len(same_gene_fivep) != 1:
				passed_filtering = False
				fail_reason = 'one_gene' if fail_reason is None else fail_reason
				failed_alignments.append((rid, read_seq, threep_chrom, threep_strand, read_is_reverse, fivep_sites, threep_coord, fail_reason))
				continue
			fivep_coord = same_gene_fivep[0]

			# Infer branchpoint coordinate
			if threep_strand == '+':
				bp_coord = alignment_end + threep_coord-threep_lengths[(threep_chrom, threep_strand, threep_coord)]-1
			else:
				bp_coord = (threep_lengths[(threep_chrom, threep_strand, threep_coord)]-alignment_end) + threep_coord

			# Check if the 5'ss is at the same coordinate or downstream of the 3'ss
			if (threep_strand == '+' and fivep_coord >= threep_coord) or (threep_strand == '-' and fivep_coord <= threep_coord):
				passed_filtering = False
				fail_reason = '5p_3p_order' if fail_reason is None else fail_reason

			# Check if the 5'ss is at the same coordinate or downstream of the branchpoint
			if (threep_strand == '+' and fivep_coord >= bp_coord) or (threep_strand == '-' and fivep_coord <= bp_coord):
				passed_filtering = False
				fail_reason = '5p_bp_order' if fail_reason is None else fail_reason

			# Compare branchpoint base in read to genome
			read_bp_nt = trimmed_read_seq[-1]

			temp_bp_bed, temp_bp_seq = output_base+'temp_bp_seqs.bed', output_base+'temp_bp_seqs.txt'
			temp_file = open(temp_bp_bed, 'w')
			if strand == '+':
				bp_start, bp_end = bp_coord-4, bp_coord+6
			else:
				bp_start, bp_end = bp_coord-5, bp_coord+5
			temp_file.write(f'{threep_chrom}\t{bp_start}\t{bp_end}\t{threep_chrom};{bp_coord};{threep_strand}\t0\t{threep_strand}\n')
			temp_file.close()
			run(f'bedtools getfasta -fi {genome_fasta} -bed {temp_bp_bed} -fo {temp_bp_seq} -nameOnly -s -tab'.split(' '))

			temp_file = open(temp_bp_seq)
			name, genomic_bp_window = temp_file.readline().strip().split()
			temp_file.close()

			genomic_bp_window = genomic_bp_window.upper()
			genomic_bp_nt = genomic_bp_window[4]

			# Add alignment to dict
			fail_reason = 'n/a' if fail_reason is None else fail_reason
			if passed_filtering is True: 
				if rid not in filtered_alignments:
					filtered_alignments[rid] = []
				filtered_alignments[rid].append([read_seq, threep_chrom, threep_strand, fivep_coord, read_is_reverse, fivep_start, fivep_end, threep_coord, bp_coord, read_bp_nt, genomic_bp_nt, genomic_bp_window])
			else:
				failed_alignments.append((rid, read_seq, threep_chrom, threep_strand, threep_coord, read_is_reverse, fivep_sites, fail_reason))
			
	run(f'rm {temp_bp_bed} {temp_bp_seq}'.split(' '))

	with open(output_base + '_failed_threep_alignments.tsv', 'w') as w:
		w.write('read_id\tread_seq\t\tthreep_chrom\tthreep_strand\tthreep_coord\tread_is_reverse\tfivep_sites\tfail_reason\n')
		for alignment in failed_alignments:
			w.write('\t'.join([str(x) for x in alignment]) + '\n')
	
	filtered_reads = set( [rid[:-4] for rid in filtered_alignments] )
	with open(run_data, 'a') as a:
		a.write(f'threep_filtered_reads\t{len(filtered_reads)}')
	
	# Loop through reads with potential alignments, outputting 1 mapping for each read
	base_rids = [rid[:-4] for rid in filtered_alignments]
	with open(out_file, 'w') as out:
		out.write('read_id\tread_seq\tchrom\tstrand\tfivep_site\tread_is_reverse\tfivep_read_start\tfivep_read_end\t')
		out.write('threep_site\tbp_site\tread_bp_nt\tgenomic_bp_nt\tgenomic_bp_window\n')
		
		for base_rid in base_rids:
			align_mismatch = {'_for':{True:[], False:[]}, '_rev':{True:[], False:[]}}
			for orientation in ('_for', '_rev'):
				for align_info in filtered_alignments.get(base_rid+orientation, []):
					read_seq, threep_chrom, threep_strand, fivep_coord, read_is_reverse, fivep_start, fivep_end, threep_coord, bp_coord, read_bp_nt, genomic_bp_nt, genomic_bp_window = align_info

					bp_mismatch = read_bp_nt != genomic_bp_nt
					align_mismatch[orientation][bp_mismatch].append(align_info)

				# # Loop through alignments, checking branchpoint read sequence vs genomic sequence
				# for align_info in alignments[rid+fivep_dir]:
				# 	read_seq, chrom, strand, fivep_coord, read_is_reverse, fivep_start, fivep_end, _, bp_coord, read_bp_nt, passed_filtering, fail_reason = align_info
				# 	temp_file = open(temp_bp_bed, 'w')
				# 	if strand == '+':
				# 		bp_start, bp_end = bp_coord-4, bp_coord+6
				# 	else:
				# 		bp_start, bp_end = bp_coord-5, bp_coord+5
				# 	temp_file.write(f'{chrom}\t{bp_start}\t{bp_end}\t{chrom};{bp_coord};{strand}\t0\t{strand}\n')
				# 	temp_file.close()
				# 	run(f'bedtools getfasta -fi {genome_fasta} -bed {temp_bp_bed} -fo {temp_bp_seq} -nameOnly -s -tab'.split(' '))
				# 	temp_file = open(temp_bp_seq)
				# 	name, genomic_bp_window = temp_file.readline().strip().split()
				# 	temp_file.close()
				# 	genomic_bp_window = genomic_bp_window.upper()
				# 	genomic_bp_nt = genomic_bp_window[4]
				# 	align_mismatch[fivep_dir][genomic_bp_nt!=read_bp_nt].append(align_info+[genomic_bp_nt, genomic_bp_window])

			# Output 1 alignment for the read, prioritizing based on mistmatch and orientation, choosing a random alignment if multiple in same category
			if len(align_mismatch['_for'][True]) > 0:
				output = [base_rid] + sample(align_mismatch['_for'][True], 1)[0]
			elif len(align_mismatch['_rev'][True]) > 0:
				output = [base_rid] + sample(align_mismatch['_rev'][True], 1)[0]
			elif len(align_mismatch['_for'][False]) > 0:
				output = [base_rid] + sample(align_mismatch['_for'][False], 1)[0]
			else:
				output = [base_rid] + sample(align_mismatch['_rev'][False], 1)[0]

			out.write('\t'.join([str(e) for e in output]) + '\n')
	


if __name__ == '__main__':
	trimmed_reads_to_threep, threep_lengths, fivep_info_table, gtf_file, genome_fasta, output_base, run_data = sys.argv[1:]

	# Load data
	gene_ranges = parse_gene_ranges(gtf_file)
	threep_lengths = parse_threep_lengths(threep_lengths)
	read_info = parse_read_info(fivep_info_table)

	# Filter trimmed read alignments to 3'ss's and write to final_info_table.txt
	out_file = output_base + '_threep_info_table.tsv'
	filter_threep_reads(trimmed_reads_to_threep, threep_lengths, read_info, gene_ranges, genome_fasta, output_base, out_file, run_data)

