
import sys
import gzip
from intervaltree import Interval, IntervalTree
from collections import Counter
from subprocess import run
from random import sample
from os.path import join

def parse_cigar(cig_str):
	cig_tuples, curr_str = [], ''
	for c in cig_str:
		if not c.isalpha():
			curr_str = curr_str + c
		else:
			cig_tuples.append((int(curr_str), c))
			curr_str = ''
	return cig_tuples

def parse_read_info(read_info_table):

	read_info = {}
	with open(read_info_table) as in_file:
		for line in in_file:
			items = line.strip().split('\t')
			if items[0] != 'read_id':
				read_info[items[0]] = items[1:]

	return read_info

def parse_gene_ranges(gtf_file):

	if gtf_file[-2:] == 'gz':
		in_file = gzip.open(gtf_file, 'rt')
	else:
		in_file = open(gtf_file)

	gene_ranges = {}
	for line in in_file:
		if line[0] != '#':
			chrom, _, feat, start, end, _, strand, _, info = line.strip().split('\t')
			if feat == 'gene':
				info = {e.split(' ')[0]:e.split(' ')[1].replace('\"', '') for e in info[:-1].split('; ')}
				if chrom not in gene_ranges:
					gene_ranges[chrom] = {s:IntervalTree() for s in ('-', '+')}
				gene_ranges[chrom][strand].add(Interval(int(start)-1, int(end), {'gene_name':info['gene_name']}))
	in_file.close()

	return gene_ranges

def parse_threep_lengths(threep_lengths):

	lengths = {}
	with open(threep_lengths) as in_file:
		for line in in_file:
			chrom, threep_site, strand, length = line.strip().split('\t')
			lengths[(chrom, strand, int(threep_site))] = int(length)

	return lengths

def filter_threep_reads(reads_to_threep, threep_lengths, read_info, gene_ranges, genome_file, output_base):

	threep_info = {}
	with open(reads_to_threep) as read_file:
		for line in read_file:
			alignment_info = line.strip().split('\t')
			rid, flag, threep_site, reference_start, mapping_quality, read_cig, _, _, _, read_seq = alignment_info[:10]
			for alignment_tag in alignment_info[11:]:
				if alignment_tag[:2] == 'XM':
					num_mismatch = int(alignment_tag.split(':')[-1])
			read_cig = parse_cigar(read_cig)
			read_len = float(sum(c[0] for c in read_cig if c[1] in ('S', 'I', 'M')))
			if num_mismatch <= 5 and num_mismatch/read_len <= 0.1:
				indels = [c[0] for c in read_cig if c[1] in ('D', 'I')]
				if len(indels) == 0 or (len(indels) == 1 and indels[0] <= 3):
					chrom, start, end, strand = threep_site[:-3].split(';')
					threep_coord = end if strand=='+' else start
					reference_start = int(reference_start)-1
					alignment_len = sum(c[0] for c in read_cig if c[1] in ('M', 'D', 'N'))
					bit_flags = bin(int(flag))
					is_reverse = True if len(bit_flags)>=7 and bit_flags[-5]=='1' else False
					if rid not in threep_info:
						threep_info[rid] = []
					threep_info[rid].append((chrom, strand, threep_coord, reference_start, reference_start+alignment_len, is_reverse, read_seq, int(mapping_quality)))
				
	potential_alignments = {}
	for rid in threep_info:
		max_score = max(s[-1] for s in threep_info[rid])
		top_alignments = [s for s in threep_info[rid] if s[-1]==max_score]
		read_seq, fivep_seq, fivep_sites, read_is_reverse, fivep_start, fivep_end = read_info[rid]
		read_is_reverse = True if read_is_reverse=='True' else False

		fivep_dict = {}
		for fp in fivep_sites.split(','):
			chrom, start, end, strand = fp.split(';')
			if chrom not in fivep_dict:
				fivep_dict[chrom] = {s:set() for s in ['+', '-']}
			if strand == '+':
				fivep_dict[chrom][strand].add(int(start))
			else:
				fivep_dict[chrom][strand].add(int(end)-1)
		
		top_alignments_filtered = []
		for align_info in top_alignments:
			threep_chrom, threep_strand, threep_site, threep_start, threep_end, threep_is_reverse, threep_read, _ = align_info
			threep_site, threep_end = int(threep_site), int(threep_end)
			if threep_chrom in fivep_dict and len(fivep_dict[threep_chrom][threep_strand]) > 0:
				threep_genes = [g.data['gene_name'] for g in gene_ranges[threep_chrom][threep_strand].overlap(threep_site, threep_site+1)]
				same_gene_fivep = []
				for fp in fivep_dict[threep_chrom][threep_strand]:
					fivep_genes = [g.data['gene_name'] for g in gene_ranges[threep_chrom][threep_strand].overlap(fp, fp+1)]
					gene_matches = sum(1 for g in fivep_genes if g in threep_genes) > 0
					if gene_matches:
						same_gene_fivep.append(fp)

				if len(same_gene_fivep) == 1:
					fivep_site = same_gene_fivep[0]
					if (strand == '+' and fivep_site < threep_site) or (strand == '-' and fivep_site > threep_site):
						if read_is_reverse == threep_is_reverse:
							if threep_strand == '+':
								bp_site = threep_end + threep_site-threep_lengths[(threep_chrom, threep_strand, threep_site)]-1
							else:
								bp_site = (threep_lengths[(threep_chrom, threep_strand, threep_site)]-threep_end) + threep_site
							if (strand == '+' and fivep_site < bp_site) or (strand == '-' and fivep_site > bp_site):
								read_bp_nt = threep_read[-1]
								if rid not in potential_alignments:
									potential_alignments[rid] = []
								potential_alignments[rid].append([read_seq, threep_chrom, threep_strand, fivep_site, read_is_reverse, fivep_start, fivep_end, threep_site, bp_site, read_bp_nt])
							
	out_rids = [rid[:-4] for rid in potential_alignments]
	temp_bp_bed, temp_bp_seq = output_base+'temp_bp_seqs.bed', output_base+'temp_bp_seqs.txt'
	with open(output_base+'_final_info_table.txt', 'w') as out_file:
		out_file.write('read_id\tread_seq\tchrom\tstrand\tfivep_site\tread_is_reverse\tfivep_read_start\tfivep_read_end\t')
		out_file.write('threep_site\tbp_site\tread_bp_nt\tgenomic_bp_nt\tgenomic_bp_window\n')
		for rid in out_rids:
			align_mismatch = {'_for':{True:[], False:[]}, '_rev':{True:[], False:[]}}
			output_rev = False
			for fivep_dir in ['_for', '_rev']:
				if rid+fivep_dir in potential_alignments:
					for align_info in potential_alignments[rid+fivep_dir]:
						read_seq, chrom, strand, fivep_site, read_is_reverse, fivep_start, fivep_end, threep_site, bp_site, read_bp_nt = align_info
						temp_file = open(temp_bp_bed, 'w')
						if strand == '+':
							bp_start, bp_end = bp_site-4, bp_site+6
						else:
							bp_start, bp_end = bp_site-5, bp_site+5
						temp_file.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, bp_start, bp_end, '{};{};{}'.format(chrom, bp_site, strand), 0, strand))
						temp_file.close()
						run('bedtools getfasta -fi {} -bed {} -fo {} -nameOnly -s -tab'.format(genome_file, temp_bp_bed, temp_bp_seq).split(' '))
						temp_file = open(temp_bp_seq)
						name, genomic_bp_window = temp_file.readline().strip().split()
						temp_file.close()
						genomic_bp_window = genomic_bp_window.upper()
						genomic_bp_nt = genomic_bp_window[4]
						align_mismatch[fivep_dir][genomic_bp_nt!=read_bp_nt].append(align_info+[genomic_bp_nt, genomic_bp_window])
			if len(align_mismatch['_for'][True]) > 0:
				output = [rid] + sample(align_mismatch['_for'][True], 1)[0]
			elif len(align_mismatch['_rev'][True]) > 0:
				output = [rid] + sample(align_mismatch['_rev'][True], 1)[0]
			elif len(align_mismatch['_for'][False]) > 0:
				output = [rid] + sample(align_mismatch['_for'][False], 1)[0]
			else:
				output = [rid] + sample(align_mismatch['_rev'][False], 1)[0]

			out_file.write('\t'.join([str(e) for e in output]) + '\n')
			
	run('rm {} {}'.format(temp_bp_bed, temp_bp_seq).split(' '))


if __name__ == '__main__':

	reads_to_threep, threep_lengths, fivep_info_table, gtf_file, genome_file, output_base = sys.argv[1:]

	gene_ranges = parse_gene_ranges(gtf_file)
	threep_lengths = parse_threep_lengths(threep_lengths)
	read_info = parse_read_info(fivep_info_table)
	filter_threep_reads(reads_to_threep, threep_lengths, read_info, gene_ranges, genome_file, output_base)

