import sys
from pyfaidx import Fasta
from collections import Counter



comp_nts = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
def reverse_complement(seq):
	return ''.join([comp_nts[seq[i]] for i in range(len(seq)-1,-1,-1)])



def filter_fivep_reads(unmapped_fasta, fivep_to_reads, fivep_upstream, fivep_trimmed_reads_out, fivep_info_table_out):
	'''
	Filter and trim the reads that 5'ss sequences mapped to
	Write trimmed read sequences to [NAME]_fivep_mapped_reads_trimmed.fa
	Write trimmed read information and their aligned 5' splice site(s) to [NAME]_fivep_info_table_out.txt
	fivep info table is TSV format, values are (read id, read sequence, 5'ss sequence, coordinates of all aligned 5'ss's that passed filtering, sequence alignment is reverse-complementary, start of alignment in read, end of alignment in read)
	'''
	# Load the collection of 5bp upstream sequences
	fivep_upstream_seqs = {}
	with open(fivep_upstream) as in_file:
		for line in in_file:
			fivep_site, seq = line.strip().split('\t')
			fivep_upstream_seqs[fivep_site[:-3]] = seq.upper()

	# Extract reads with perfect alignments from the 5'ss mapping 
	read_sites = {}				# { read id: set(first 20bp of intron sequence) }
	site_coords = {}			# { read id: {first 20bp of intron sequence: (alignment start position in read, alignment end position in read, is reverse-complementary)} }
	with open(fivep_to_reads) as fivep_file:
		# Loop through alignments
		for line in fivep_file:
			alignment_info = line.strip().split('\t')
			fivep_site, flag, rid, reference_start, _, read_cig = alignment_info[:6]

			# Get mismatch count
			for alignment_tag in alignment_info[11:]:
				if alignment_tag[:2] == 'XM':
					num_mismatch = int(alignment_tag.split(':')[-1])

			# If it's a perfect alignment, add it to read_sites and sites_coords
			if num_mismatch == 0 and read_cig == '20M':
				reference_start = int(reference_start)-1
				bit_flags = bin(int(flag))
				is_reverse = True if len(bit_flags)>=7 and bit_flags[-5]=='1' else False
				fivep_site = fivep_site[:-3]
				if rid not in read_sites:
					read_sites[rid] = set()
					site_coords[rid] = {}
				read_sites[rid].add(fivep_site)
				site_coords[rid][fivep_site] = (reference_start, reference_start+20, is_reverse)

	# Filter, trim, and write reads
	read_fasta = Fasta(unmapped_fasta, as_raw=True)
	with open(fivep_trimmed_reads_out, 'w') as trimmed_out, open(fivep_info_table_out, 'w') as info_out:
		info_out.write('read_id\tread_seq\tfivep_seq\tfivep_sites\tfivep_first\tread_fivep_start\tread_fivep_end\n')
		
		# Loop through reads with extracted alignments
		for rid in read_sites:
			read_seq = read_fasta[rid][:]
			fivep_pass = {True:[], False:[]}	# { is reverse: [(first 20bp of intron sequence, (alignment start position in read, alignment end position in read, is reverse-complementary)...], is not reverse: [...] }
			
			# Check if the 5bp upstream of the alignment in the read matches the 5bp upstream of the 5'ss in the genome. 
			# If it does NOT, add the read alignment to fivep_pass
			for fp in site_coords[rid]:
				fivep_start, fivep_end, is_reverse = site_coords[rid][fp]
				if is_reverse:
					read_upstream = read_seq[fivep_end:fivep_end+5].upper()
					upstream_mismatch = read_upstream != reverse_complement(fivep_upstream_seqs[fp])
				else:
					read_upstream = read_seq[fivep_start-5:fivep_start].upper()
					upstream_mismatch = read_upstream != fivep_upstream_seqs[fp]
				if upstream_mismatch:
					fivep_pass[is_reverse].append((fp, site_coords[rid][fp]))

			# For each orientation, trim off the upstream-most 5'ss and everything upstream of it, then write the trimmed sequence + alignments to file
			for is_reverse in fivep_pass:
				# if there are no alignments for the read in the given orientation, skip it
				if len(fivep_pass[is_reverse]) == 0:
					continue

				if is_reverse:
					# Get the start and end of the rightmost alignment in the read 
					fivep_start, fivep_end, _ = max(fivep_pass[is_reverse], key=lambda fp:fp[1][0])[1]
					# Keep the subset of 5'ss alignments that start at the same rightmost position
					fivep_pass_sub = [fp for fp in fivep_pass[is_reverse] if fp[1][0]==fivep_start]
					# Trim off the rightmost alignment and everything to the left of it
					trim_seq = read_seq[fivep_end:]
					# Get sequence of rightmost alignment
					fivep_seq = reverse_complement(read_seq[fivep_start:fivep_end])
				else:
					# Get the start and end of the leftmost alignment in the read 
					fivep_start, fivep_end, _ = min(fivep_pass[is_reverse], key=lambda fp:fp[1][0])[1]
					# Keep the subset of 5'ss alignments that start at the same leftmost position
					fivep_pass_sub = [fp for fp in fivep_pass[is_reverse] if fp[1][0]==fivep_start]
					# Trim off the leftmost alignment and everything to the right of it
					trim_seq = read_seq[:fivep_start]
					# Get sequence of leftmost alignment
					fivep_seq = read_seq[fivep_start:fivep_end]

				# If less than 20bp is left in the read, skip it
				if len(trim_seq) < 20:
					continue

				out_rid = rid + '_rev' if is_reverse else rid + '_for'
				trimmed_out.write('>{}\n{}\n'.format(out_rid, trim_seq))

				fivep_sites = ','.join([fp[0] for fp in fivep_pass_sub])
				info_out.write(f'{out_rid}\t{read_seq}\t{fivep_seq}\t{fivep_sites}\t{is_reverse}\t{fivep_start}\t{fivep_end}\n')

if __name__ == '__main__' :
	
	unmapped_fasta, fivep_to_reads, fivep_upstream, fivep_trimmed_reads_out, fivep_info_table_out = sys.argv[1:]
	
	filter_fivep_reads(unmapped_fasta, fivep_to_reads, fivep_upstream, fivep_trimmed_reads_out, fivep_info_table_out)



			

