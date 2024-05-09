
import sys
from pyfaidx import Fasta
import subprocess
import threading

# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
COMP_NTS = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}


# =============================================================================#
#                                   Classes                                    #
# =============================================================================#
class filter_thread (threading.Thread):
	def __init__(self, unmapped_fasta, alignments, genome_fasta):
		threading.Thread.__init__(self)
		self.unmapped_fasta, self.alignments, self.genome_fasta = unmapped_fasta, alignments, genome_fasta
	def run(self):
		self.out_reads, self.failed_alignments = filter_fivep_reads(self.unmapped_fasta, self.alignments, self.genome_fasta)


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
		for line in fivep_file:
			alignment_info = line.strip().split('\t')
			fivep_site, flag, rid, read_fivep_start, _, read_cig = alignment_info[:6]
			
			read_fivep_start = int(read_fivep_start)-1
			bit_flags = bin(int(flag))
			read_is_reverse = True if len(bit_flags)>=7 and bit_flags[-5]=='1' else False

			if rid not in alignments:
				alignments[rid] = {}
			alignments[rid][fivep_site] = (read_fivep_start, read_fivep_start+20, read_is_reverse)

	return alignments


def filter_fivep_reads(unmapped_fasta:str, alignments:dict, genome_fasta:str):
	'''
	Filter and trim the reads to which 5'ss sequences were mapped
	Write trimmed read sequences to [NAME]_fivep_mapped_reads_trimmed.fa
	Write trimmed read information and their aligned 5' splice site(s) to [NAME]_fivep_info_table_out.txt
	fivep info table is in TSV format, values are (read id, read sequence, 5'ss sequence, coordinates of all aligned 5'ss's that passed filtering, alignment is reverse-complementary, start of alignment in read, end of alignment in read)
	'''
	failed_alignments = []		# [ (read id, read_seq, 5'ss site, alignment start in read, alignment end in read, alignment is reverse-complementary, filter that it failed or None)...] }
	out_reads = []
	read_fasta = Fasta(unmapped_fasta, as_raw=True)

	# Check if the 5bp upstream of the alignment in the read matches the 5bp upstream of the 5'ss in the genome. 
	# If it does NOT, add the read alignment to fivep_pass
	# { is reverse: [(5'ss site), (alignment start position in read, alignment end position in read, is reverse-complementary)...], is not reverse: [...] }
	bedtools_input = ''
	for rid in alignments:
		for site in alignments[rid]:
			chrom, fivep_pos, strand = site.split(';')
			fivep_pos = int(fivep_pos)
			bedtools_input += f'{chrom}\t{fivep_pos-5}\t{fivep_pos}\t{site}\t0\t{strand}\n' if strand == '+' else f'{chrom}\t{fivep_pos+1}\t{fivep_pos+6}\t{site}\t0\t{strand}\n'

	bedtools_call = f'bedtools getfasta -s -tab -nameOnly -fi {genome_fasta} -bed -'
	bedtools_output = subprocess.run(bedtools_call.split(' '), input=bedtools_input, check=True, capture_output=True, text=True).stdout.strip().split('\n')
	fivep_upstream_seqs = dict([(l.split('\t')[0][:-3], l.split('\t')[1].upper()) for l in bedtools_output])
	
	for rid in alignments:
		read_seq = read_fasta[rid][:]
		fivep_pass = {True:[], False:[]}
		for site in alignments[rid]:
			read_fivep_start, read_fivep_end, read_is_reverse = alignments[rid][site]
			if read_is_reverse:
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
				read_bp_pos = read_fivep_end
				# Trim off the rightmost alignment and everything to the left of it
				trim_seq = read_seq[read_fivep_end:]
				# Get sequence of rightmost alignment
				fivep_seq = reverse_complement(read_seq[read_fivep_start:read_fivep_end])
			else:
				# Get the start and end of the leftmost alignment in the read 
				read_fivep_start, read_fivep_end, _ = min(fivep_pass[read_is_reverse], key=lambda fp:fp[1][0])[1]
				read_bp_pos = read_fivep_start - 1
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
			out_reads.append((trim_seq, out_rid, read_seq, fivep_seq, fivep_sites, read_is_reverse, read_fivep_start, read_fivep_end, read_bp_pos))
		
	return out_reads, failed_alignments


def write_out_reads(out_reads:list, fivep_trimmed_reads_out:str, fivep_info_table_out:str) -> None:
	'''
	Write 
	'''
	with open(fivep_trimmed_reads_out, 'w') as trimmed_out, open(fivep_info_table_out, 'w') as info_out:
		info_out.write('read_id\tread_seq\tfivep_seq\tfivep_sites\tread_is_reverse\tread_fivep_start\tread_fivep_end\tread_bp_pos\n')

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

	threads, unmapped_fasta, fivep_to_reads, genome_fasta, fivep_trimmed_reads_out, fivep_info_table_out, output_base = sys.argv[1:]

	alignments = load_alignments(fivep_to_reads)
	read_ids, threads = list(alignments.keys()), int(threads)

	thread_list = []
	for thread_num in range(threads):
		thread_alignments = {rid:alignments[rid] for rid in read_ids[thread_num::threads]}
		new_thread = filter_thread(unmapped_fasta, thread_alignments, genome_fasta)
		new_thread.start()
		thread_list.append(new_thread)

	for t in thread_list:
		t.join()

	out_reads, failed_alignments = [], []
	for t in thread_list:
		out_reads, failed_alignments = out_reads+t.out_reads, failed_alignments+t.failed_alignments

	write_out_reads(out_reads, fivep_trimmed_reads_out, fivep_info_table_out)

	with open(f'{output_base}failed_fivep_alignments.tsv', 'w') as failed_out:
		failed_out.write('read_id\tread_seq\tfivep_site\tread_fivep_start\tread_fivep_end\tread_is_reverse\tfail_reason\n')
		for info in failed_alignments:
			failed_out.write('\t'.join([str(x) for x in info]) + '\n')

	out_rids = set([x[1][:-4] for x in out_reads])
	with open(f'{output_base}run_data.tsv', 'a') as stats_out:
		stats_out.write(f'fivep_mapped_reads\t{len(alignments.keys())}\n')
		stats_out.write(f'fivep_filtered_reads\t{len(out_rids)}\n')



			

