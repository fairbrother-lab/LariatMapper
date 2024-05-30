import sys
import subprocess
import itertools as it
import multiprocessing
import time

from pyfaidx import Fasta
import pandas as pd



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
COMP_NTS = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
FIVEP_INFO_TABLE_COLS = ['read_id',
						'read_seq',
						'fivep_seq',
						'fivep_sites',
						'read_is_reverse',
						'read_fivep_start',
						'read_fivep_end',
						'read_bp_pos',
						]
FAILED_ALIGNMENTS_COLS = ['read_id',
						'read_seq',
						'fivep_site',
						'read_fivep_start',
						'read_fivep_end',
						'read_is_reverse',
						'fail_reason',
						]

out_lock = multiprocessing.Lock()
failed_out_lock = multiprocessing.Lock()



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def reverse_complement(seq):
	return ''.join([COMP_NTS[seq[i]] for i in range(len(seq)-1,-1,-1)])


def parse_alignments(fivep_to_reads:str) -> dict:
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


def get_read_seqs(unmapped_fasta:str) -> dict:
	read_fasta = Fasta(unmapped_fasta, as_raw=True)
	read_seqs = {}
	for rid in read_fasta.keys():
		read_seqs[rid] = read_fasta[rid][:]
	
	return read_seqs


def get_fivep_upstream_seqs(alignments:dict, genome_fasta:str) -> dict:
	# Prepare 5' sites for input into bedtools getfasta
	# We're retrieving the 5bp upstream of the 5'ss, which is the last 5bp in the upstream exon
	bedtools_input = ''
	for rid in alignments:
		for site in alignments[rid]:
			chrom, fivep_pos, strand = site.split(';')
			fivep_pos = int(fivep_pos)
			bedtools_input += f'{chrom}\t{fivep_pos-5}\t{fivep_pos}\t{site}\t0\t{strand}\n' if strand == '+' else f'{chrom}\t{fivep_pos+1}\t{fivep_pos+6}\t{site}\t0\t{strand}\n'

	# Call bedtools getfasta 
	bedtools_call = f'bedtools getfasta -s -tab -nameOnly -fi {genome_fasta} -bed -'
	bedtools_output = subprocess.run(bedtools_call.split(' '), input=bedtools_input, check=True, capture_output=True, text=True).stdout.strip().split('\n')

	# Parse output
	fivep_upstream_seqs = dict([(l.split('\t')[0][:-3], l.split('\t')[1].upper()) for l in bedtools_output])

	return fivep_upstream_seqs


def filter_fivep_reads(alignments:dict, read_seqs:dict, fivep_upstream_seqs:dict, output_base:str) -> None:
	'''
	Filter and trim the reads to which 5'ss sequences were mapped
	Write trimmed read information and their aligned 5' splice site(s) to fivep_info_table_out.txt
	Write trimmed read sequences to fivep_mapped_reads_trimmed.fa
	'''
	failed_alignments = []		# [ (read id, read_seq, 5'ss site, alignment start in read, alignment end in read, alignment is reverse-complementary, filter that it failed or None)...] }
	out_reads = []

	# Check if the 5bp upstream of the alignment in the read matches the 5bp upstream of the 5'ss in the genome. 
	# If it does NOT, add the read alignment to fivep_pass
	# { is reverse: [(5'ss site), (alignment start position in read, alignment end position in read, is reverse-complementary)...], is not reverse: [...] }
	for rid in alignments:
		read_seq = read_seqs[rid]
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

	# Write the filtered alignments and trimmed sequences to file
	with out_lock:
		with open(f'{output_base}fivep_info_table.tsv', 'a') as info_out, open(f'{output_base}fivep_mapped_reads_trimmed.fa', 'a') as trimmed_out:
			for row in out_reads:
				row = [str(item) for item in row]
				info_out.write('\t'.join(row[1:]) + '\n')

				trim_seq = row[0]
				out_rid = row[1]
				trimmed_out.write(f'>{out_rid}\n{trim_seq}\n')
	
	# Write the failed alignments to file
	with failed_out_lock:
		failed_alignments = pd.DataFrame(failed_alignments, columns=FAILED_ALIGNMENTS_COLS)
		failed_alignments.to_csv(f'{output_base}failed_fivep_alignments.tsv', mode='a', sep='\t', index=False, header=False)



# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__' :
	threads, genome_fasta, output_base = sys.argv[1:]
	threads = int(threads)

	alignments = parse_alignments(f'{output_base}fivep_to_reads.sam')
	read_seqs = get_read_seqs(f'{output_base}unmapped_reads.fa')
	fivep_upstream_seqs = get_fivep_upstream_seqs(alignments, genome_fasta)

	# Record read count
	rids = set(rid.split('/')[0] for rid in alignments.keys())
	with open(f'{output_base}read_counts.tsv', 'a') as a:
		a.write(f'fivep_mapped_reads\t{len(rids)}\n')	

	if len(alignments)==0:
		print(time.strftime('%m/%d/%y - %H:%M:%S') + '| No reads remaining')
		exit()

	# Prep out files with a header row
	with open(f'{output_base}fivep_info_table.tsv', 'w') as w:
		w.write('\t'.join(FIVEP_INFO_TABLE_COLS) + '\n')
	with open(f'{output_base}fivep_mapped_reads_trimmed.fa', 'w') as w:
		pass
	with open(f'{output_base}failed_fivep_alignments.tsv', 'w') as w:
		w.write('\t'.join(FAILED_ALIGNMENTS_COLS) + '\n')

	rid_list = list(alignments.keys())
	processes = []
	for i in range(threads):
		rid_subset = rid_list[i::threads]
		alignments_subset = {rid:alignments[rid] for rid in rid_subset}
		if len(alignments_subset)==0:
			continue
		subset_process = multiprocessing.Process(target=filter_fivep_reads, args=(alignments_subset, read_seqs, fivep_upstream_seqs, output_base,))
		subset_process.start()
		processes.append(subset_process)
		
	for p in processes:
		p.join()