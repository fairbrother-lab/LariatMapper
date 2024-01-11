import sys
from pyfaidx import Fasta
from collections import Counter



comp_nts = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
def reverse_complement(seq):
	return ''.join([comp_nts[seq[i]] for i in range(len(seq)-1,-1,-1)])



def filter_fivep_reads(unmapped_fasta, fivep_to_reads, fivep_upstream, fivep_trimmed_reads_out, fivep_info_table_out, run_data):
	'''
	Filter and trim the reads to which 5'ss sequences were mapped
	Write trimmed read sequences to [NAME]_fivep_mapped_reads_trimmed.fa
	Write trimmed read information and their aligned 5' splice site(s) to [NAME]_fivep_info_table_out.txt
	fivep info table is in TSV format, values are (read id, read sequence, 5'ss sequence, coordinates of all aligned 5'ss's that passed filtering, sequence alignment is reverse-complementary, start of alignment in read, end of alignment in read)
	'''
	# Load the collection of 5bp upstream sequences
	fivep_upstream_seqs = {}
	with open(fivep_upstream) as in_file:
		for line in in_file:
			fivep_site, seq = line.strip().split('\t')
			fivep_upstream_seqs[fivep_site[:-3]] = seq.upper()

	# # Extract reads with perfect alignments from the 5'ss mapping 
	# Load alignment data
	site_coords = {}			# { read id: {first 20bp of intron sequence: (alignment start position in read, alignment end position in read, is reverse-complementary)} }
	mapped_rids = set()
	with open(fivep_to_reads) as fivep_file:
		# Loop through alignments
		for line in fivep_file:
			alignment_info = line.strip().split('\t')
			fivep_site, flag, rid, reference_start, _, read_cig = alignment_info[:6]

			mapped_rids.add(rid)

			# # Get mismatch count
			# for alignment_tag in alignment_info[11:]:
			# 	if alignment_tag[:2] == 'XM':
			# 		num_mismatch = int(alignment_tag.split(':')[-1])

			# # If it's a perfect alignment, add it to read_sites and sites_coords
			# if num_mismatch == 0 and read_cig == '20M':
			reference_start = int(reference_start)-1
			bit_flags = bin(int(flag))
			read_is_reverse = True if len(bit_flags)>=7 and bit_flags[-5]=='1' else False
			fivep_site = fivep_site[:-3]
			if rid not in site_coords:
				site_coords[rid] = {}
			site_coords[rid][fivep_site] = (reference_start, reference_start+20, read_is_reverse)

	with open(run_data, 'a') as a:
		a.write(f'fivep_mapped_reads\t{len(mapped_rids)}\n')

	# Filter, trim, and write reads
	read_fasta = Fasta(unmapped_fasta, as_raw=True)
	with open(fivep_trimmed_reads_out, 'w') as trimmed_out, open(fivep_info_table_out, 'w') as info_out:
		info_out.write('read_id\tread_seq\tfivep_seq\tfivep_sites\tread_is_reverse\tread_fivep_start\tread_fivep_end\n')
		
		# Loop through reads with extracted alignments
		passed_rids = set()
		for rid in site_coords:
			read_seq = read_fasta[rid][:]
			fivep_pass = {True:[], False:[]}	# { is reverse: [(first 20bp of intron sequence, (alignment start position in read, alignment end position in read, is reverse-complementary)...], is not reverse: [...] }
			
			# Check if the 5bp upstream of the alignment in the read matches the 5bp upstream of the 5'ss in the genome. 
			# If it does NOT, add the read alignment to fivep_pass
			for site in site_coords[rid]:
				fivep_start, fivep_end, read_is_reverse = site_coords[rid][site]
				if read_is_reverse:
					read_upstream = read_seq[fivep_end:fivep_end+5].upper()
					upstream_mismatch = read_upstream != reverse_complement(fivep_upstream_seqs[site])
				else:
					read_upstream = read_seq[fivep_start-5:fivep_start].upper()
					upstream_mismatch = read_upstream != fivep_upstream_seqs[site]
				if upstream_mismatch:
					fivep_pass[read_is_reverse].append((site, site_coords[rid][site]))

			# For each orientation, trim off the upstream-most 5'ss and everything upstream of it, then write the trimmed sequence + alignments to file
			for read_is_reverse in fivep_pass:
				out_rid = rid + '_rev' if read_is_reverse else rid + '_for'

				# Check if there are no alignments for the read in the given orientation
				if len(fivep_pass[read_is_reverse]) == 0:
					continue

				if read_is_reverse:
					# Get the start and end of the rightmost alignment in the read 
					fivep_start, fivep_end, _ = max(fivep_pass[read_is_reverse], key=lambda fp:fp[1][0])[1]
					# Keep the subset of 5'ss alignments that start at the same rightmost position
					fivep_pass_sub = [fp for fp in fivep_pass[read_is_reverse] if fp[1][0]==fivep_start]
					# Trim off the rightmost alignment and everything to the left of it
					trim_seq = read_seq[fivep_end:]
					# Get sequence of rightmost alignment
					fivep_seq = reverse_complement(read_seq[fivep_start:fivep_end])
				else:
					# Get the start and end of the leftmost alignment in the read 
					fivep_start, fivep_end, _ = min(fivep_pass[read_is_reverse], key=lambda fp:fp[1][0])[1]
					# Keep the subset of 5'ss alignments that start at the same leftmost position
					fivep_pass_sub = [fp for fp in fivep_pass[read_is_reverse] if fp[1][0]==fivep_start]
					# Trim off the leftmost alignment and everything to the right of it
					trim_seq = read_seq[:fivep_start]
					# Get sequence of leftmost alignment
					fivep_seq = read_seq[fivep_start:fivep_end]

				# Check if less than 20bp is left in the read
				if len(trim_seq) < 20:
					continue
				
				# If passed filtering, add trimmed seq to fasta
				passed_rids.add(rid)
				trimmed_out.write(f'>{out_rid}\n{trim_seq}\n')
				
				fivep_sites = sorted([fp[0] for fp in fivep_pass_sub])
				fivep_sites = ','.join(fivep_sites)
				info_out.write(f'{out_rid}\t{read_seq}\t{fivep_seq}\t{fivep_sites}\t{read_is_reverse}\t{fivep_start}\t{fivep_end}\n')
	
	with open(run_data, 'a') as a:
		a.write(f'fivep_filtered_reads\t{len(passed_rids)}\n')


if __name__ == '__main__' :
	
	unmapped_fasta, fivep_to_reads, fivep_upstream, fivep_trimmed_reads_out, fivep_info_table_out, run_data = sys.argv[1:]
	
	filter_fivep_reads(unmapped_fasta, fivep_to_reads, fivep_upstream, fivep_trimmed_reads_out, fivep_info_table_out, run_data)



			

