import sys

from pyfaidx import Fasta
import pandas as pd
import time



# =============================================================================#
#                                  Constants                                   #
# =============================================================================#
COMP_NTS = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def reverse_complement(seq):
	return ''.join([COMP_NTS[seq[i]] for i in range(len(seq)-1,-1,-1)])


def load_alignments(fivep_to_reads:str) -> dict:
	'''
	Load 5'ss alignments to dict
	Returns { read id: {first 20bp of intron sequence: (alignment start position in read, alignment end position in read, is reverse-complementary)} }
	'''
	alignments = {}	
	with open(fivep_to_reads) as fivep_file:
		# Loop through alignments
		for line in fivep_file:
			alignment_info = line.strip().split('\t')
			fivep_site, flag, rid, read_fivep_start, _, read_cig = alignment_info[:6]
			
			read_fivep_start = int(read_fivep_start)-1
			bit_flags = bin(int(flag))
			read_is_reverse = True if len(bit_flags)>=7 and bit_flags[-5]=='1' else False
			# fivep_site = fivep_site[:-3]

			if rid not in alignments:
				alignments[rid] = {}
			alignments[rid][fivep_site] = (read_fivep_start, read_fivep_start+20, read_is_reverse)

	return alignments


# def load_fivep_upstream(fivep_upstream:str) -> dict:
# 	'''
# 	Load the collection of 5bp upstream sequences
# 	Returns { 5'ss site : 5bp upstream sequence }, e.g. { "chr1;201283904;201283924;+": "TCGAG" }
# 	'''
# 	fivep_upstream_seqs = {}
# 	with open(fivep_upstream) as in_file:
# 		for line in in_file:
# 			fivep_site, seq = line.strip().split('\t')
# 			fivep_upstream_seqs[fivep_site[:-3]] = seq.upper()

# 	return fivep_upstream_seqs


def filter_fivep_reads(unmapped_fasta:str, alignments:dict, fivep_upstream_seqs:dict):
	'''
	Filter and trim the reads to which 5'ss sequences were mapped
	Write trimmed read sequences to [NAME]_fivep_mapped_reads_trimmed.fa
	Write trimmed read information and their aligned 5' splice site(s) to [NAME]_fivep_info_table_out.txt
	fivep info table is in TSV format, values are (read id, read sequence, 5'ss sequence, coordinates of all aligned 5'ss's that passed filtering, alignment is reverse-complementary, start of alignment in read, end of alignment in read)
	'''
	failed_alignments = []		# [ (read id, read_seq, 5'ss site, alignment start in read, alignment end in read, alignment is reverse-complementary, filter that it failed or None)...] }
	out_reads = []
	read_fasta = Fasta(unmapped_fasta, as_raw=True)

	# Loop through reads with extracted alignments
	for rid in alignments:
		read_seq = read_fasta[rid][:]
		
		# Check if the 5bp upstream of the alignment in the read matches the 5bp upstream of the 5'ss in the genome. 
		# If it does NOT, add the read alignment to fivep_pass
		fivep_pass = {True:[], False:[]}	# { is reverse: [(5'ss site), (alignment start position in read, alignment end position in read, is reverse-complementary)...], is not reverse: [...] }
		for site in alignments[rid]:
			read_fivep_start, read_fivep_end, read_is_reverse = alignments[rid][site]
			if read_is_reverse:
				# read_upstream = read_seq[read_fivep_end-1:read_fivep_end+4].upper()
				read_upstream = read_seq[read_fivep_end+1:read_fivep_end+6].upper()
				upstream_mismatch = read_upstream != reverse_complement(fivep_upstream_seqs[site])
			else:
				read_upstream = read_seq[read_fivep_start-5:read_fivep_start].upper()
				upstream_mismatch = read_upstream != fivep_upstream_seqs[site]
				
			if upstream_mismatch:
				fivep_pass[read_is_reverse].append((site, alignments[rid][site]))
			else:
				failed_alignments.append((rid, read_seq, site, read_fivep_start, read_fivep_end, read_is_reverse, '5bp_up_match'))

		# For each orientation, trim off the upstream-most 5'ss and everything upstream of it, then write the trimmed sequence + alignments to file
		for read_is_reverse in fivep_pass:
			# Check if there are no alignments for the read in the given orientation
			if len(fivep_pass[read_is_reverse]) == 0:
				continue

			if read_is_reverse:
				# Get the start and end of the rightmost alignment in the read 
				read_fivep_start, read_fivep_end, _ = max(fivep_pass[read_is_reverse], key=lambda fp:fp[1][0])[1]
				# Trim off the rightmost alignment and everything to the left of it
				# trim_seq = read_seq[read_fivep_end-1:]
				trim_seq = read_seq[read_fivep_end:]
				# Get sequence of rightmost alignment
				fivep_seq = reverse_complement(read_seq[read_fivep_start:read_fivep_end])
			else:
				# Get the start and end of the leftmost alignment in the read 
				read_fivep_start, read_fivep_end, _ = min(fivep_pass[read_is_reverse], key=lambda fp:fp[1][0])[1]
				# Trim off the leftmost alignment and everything to the right of it
				trim_seq = read_seq[:read_fivep_start]
				# Get sequence of leftmost alignment
				fivep_seq = read_seq[read_fivep_start:read_fivep_end]

			# Keep the subset of 5'ss alignments that start at the upstream-most position and fail the rest
			fivep_pass_sub = []
			for fp in fivep_pass[read_is_reverse]:
				if fp[1][0] == read_fivep_start:
					fivep_pass_sub.append(fp)
				else:
					failed_alignments.append((rid, read_seq, fp[0], fp[1][0], fp[1][1], fp[1][2], 'furthest_upstream'))

			# Check if less than 20bp is left in the read
			if len(trim_seq) < 20:
				for fp in fivep_pass_sub:
					failed_alignments.append((rid, read_seq, fp[0], fp[1][0], fp[1][1], fp[1][2], 'enough_trim_seq'))
				continue
			
			# Add reads + alignment(s) that passed filtering to out_reads 
			out_rid = rid + '_rev' if read_is_reverse else rid + '_for'
			fivep_sites = sorted([fp[0] for fp in fivep_pass_sub])
			fivep_sites = ','.join(fivep_sites)
			out_reads.append((trim_seq, out_rid, read_seq, fivep_seq, fivep_sites, read_is_reverse, read_fivep_start, read_fivep_end))
		
	return out_reads, failed_alignments


def write_out_reads(out_reads:list, fivep_trimmed_reads_out:str, fivep_info_table_out:str) -> None:
	'''
	Write 
	'''
	with open(fivep_trimmed_reads_out, 'w') as trimmed_out, open(fivep_info_table_out, 'w') as info_out:
		info_out.write('read_id\tread_seq\tfivep_seq\tfivep_sites\tread_is_reverse\tread_fivep_start\tread_fivep_end\n')

		for info in out_reads:
			trim_seq = info[0]
			out_rid = info[1]
			trimmed_out.write(f'>{out_rid}\n{trim_seq}\n')

			row = '\t'.join([str(x) for x in info[1:]])
			info_out.write(row + '\n')
	

# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__' :
	print(time.strftime('%m/%d/%y - %H:%M:%S') + f' | Arguments recieved: {sys.argv[1:]}')
	unmapped_fasta, fivep_to_reads, fivep_upstream, fivep_trimmed_reads_out, fivep_info_table_out, output_base = sys.argv[1:]

	alignments = load_alignments(fivep_to_reads)
	# fivep_upstream_seqs = load_fivep_upstream(fivep_upstream)
	fivep_upstream_seqs = pd.read_csv(fivep_upstream, sep='\t', index_col='fivep_site').upstream_sequence.to_dict()

	out_reads, failed_alignments = filter_fivep_reads(unmapped_fasta, alignments, fivep_upstream_seqs)

	write_out_reads(out_reads, fivep_trimmed_reads_out, fivep_info_table_out)

	with open(f'{output_base}failed_fivep_alignments.tsv', 'w') as w:
		w.write('read_id\tread_seq\tfivep_site\tread_fivep_start\tread_fivep_end\tread_is_reverse\tfail_reason\n')
		for info in failed_alignments:
			row = '\t'.join([str(x) for x in info])
			w.write(row + '\n')

	out_rids = set([x[1][:-4] for x in out_reads])
	with open(f'{output_base}run_data.tsv', 'a') as a:
		a.write(f'fivep_mapped_reads\t{len(alignments.keys())}\n')
		a.write(f'fivep_filtered_reads\t{len(out_rids)}\n')



			

