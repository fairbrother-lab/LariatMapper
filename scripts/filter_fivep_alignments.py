from dataclasses import dataclass, field
import sys
import subprocess
import collections
import multiprocessing
import time
import math
import logging

from pyfaidx import Fasta
import pandas as pd
import pysam
import numpy as np



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
COMP_NTS = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
# READ_CHUNKSIZE = 1_000

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



	

	




# def parse_alignments(fivep_to_reads:str) -> dict:
# 	'''
# 	Load 5'ss alignments to dict
# 	Returns { read id: {first 20bp of intron sequence: (alignment start position in read, alignment end position in read, is reverse-complementary)} }
# 	'''
# 	alignments = {}	
# 	with open(fivep_to_reads) as fivep_file:
# 		for line in fivep_file:
# 			alignment_info = line.strip().split('\t')
# 			fivep_site, flag, read_id, read_fivep_start, _, read_cig = alignment_info[:6]
			
# 			read_fivep_start = int(read_fivep_start)-1
# 			bit_flags = bin(int(flag))
# 			read_is_reverse = True if len(bit_flags)>=7 and bit_flags[-5]=='1' else False

# 			if read_id not in alignments:
# 				alignments[read_id] = {}
# 			alignments[read_id][fivep_site] = (read_fivep_start, read_fivep_start+20, read_is_reverse)

# 	return alignments


# def parse_alignments(fivep_to_reads:str) -> dict:
# 	'''
# 	Load 5'ss alignments to dict
# 	Returns { read id: {first 20bp of intron sequence: (alignment start position in read, alignment end position in read, is reverse-complementary)} }
# 	'''
# 	align_iterator = pysam.AlignmentFile(fivep_to_reads, 'r')
# 	alignments = collections.defaultdict(dict)
# 	i = 0
# 	for align in align_iterator:
# 		alignments[align.query_name][align.reference_name] = (align.query_alignment_start, align.query_alignment_end, align.is_reverse)

# 		i += 1
# 		if i == ALIGN_CHUNKSIZE:
# 			# Make sure you're getting all the alignments for the current read id
# 			yield alignments
# 			# Reset
# 			i = 0
# 			alignments = collections.defaultdict(dict)
			

# 	return alignments


# def parse_alignments(reads_to_fivep:str) -> dict:
# 	'''
# 	Load 5'ss alignments to dict
# 	Returns { read id: {fivep site: (alignment start position in read, alignment end position in read, is reverse-complementary)} }
# 	'''
# 	alignments = collections.defaultdict(dict)
# 	for align in pysam.AlignmentFile(reads_to_fivep, 'r'):
# 		alignments[align.query_name][align.reference_name] = (align.query_alignment_start, align.query_alignment_end, align.is_reverse)

# 	return alignments


def get_read_seqs(unmapped_fasta:str) -> dict:
	read_fasta = Fasta(unmapped_fasta, as_raw=True)
	read_seqs = {}
	for read_id in read_fasta.keys():
		read_seqs[read_id] = read_fasta[read_id][:]
	
	return read_seqs


# def get_fivep_upstream_seqs(alignments:dict, genome_fasta:str) -> dict:
# 	# Prepare 5' sites for input into bedtools getfasta
# 	# We're retrieving the 5bp upstream of the 5'ss, which is the last 5bp in the upstream exon
# 	bedtools_input = ''
# 	for read_id in alignments:
# 		for site in alignments[read_id]:
# 			chrom, fivep_pos, strand = site.split(';')
# 			fivep_pos = int(fivep_pos)
# 			bedtools_input += f'{chrom}\t{fivep_pos-5}\t{fivep_pos}\t{site}\t0\t{strand}\n' if strand == '+' else f'{chrom}\t{fivep_pos+1}\t{fivep_pos+6}\t{site}\t0\t{strand}\n'

# 	# Call bedtools getfasta 
# 	bedtools_call = f'bedtools getfasta -s -tab -nameOnly -fi {genome_fasta} -bed -'
# 	bedtools_output = subprocess.run(bedtools_call.split(' '), input=bedtools_input, check=True, capture_output=True, text=True).stdout.strip().split('\n')

# 	# Parse output
# 	fivep_upstream_seqs = dict([(l.split('\t')[0][:-3], l.split('\t')[1].upper()) for l in bedtools_output])

# 	return fivep_upstream_seqs


def get_fivep_upstream_seqs(fivep_fasta:str, genome_fasta:str) -> dict:
	# Prepare 5' sites for input into bedtools getfasta
	# We're retrieving the 5bp upstream of the 5'ss, which is the last 5bp in the upstream exon
	bedtools_input = ''
	fivep_sites = Fasta(fivep_fasta).keys()
	for site in fivep_sites:
		chrom, fivep_pos, strand = site.split(';')
		fivep_pos = int(fivep_pos)
		bedtools_input += f'{chrom}\t{fivep_pos-5}\t{fivep_pos}\t{site}\t0\t{strand}\n' if strand == '+' else f'{chrom}\t{fivep_pos+1}\t{fivep_pos+6}\t{site}\t0\t{strand}\n'

	# Call bedtools getfasta 
	bedtools_call = f'bedtools getfasta -s -tab -nameOnly -fi {genome_fasta} -bed -'
	bedtools_output = subprocess.run(bedtools_call.split(' '), input=bedtools_input, check=True, capture_output=True, text=True).stdout.strip().split('\n')

	# Parse output
	fivep_upstream_seqs = dict([(l.split('\t')[0][:-3], l.split('\t')[1].upper()) for l in bedtools_output])

	return fivep_upstream_seqs


# def filter_fivep_reads(alignments:dict, read_seqs:dict, fivep_upstream_seqs:dict, output_base:str) -> None:
# 	'''
# 	Filter and trim the reads to which 5'ss sequences were mapped
# 	Write trimmed read information and their aligned 5' splice site(s) to fivep_info_table_out.txt
# 	Write trimmed read sequences to fivep_mapped_reads_trimmed.fa
# 	'''
# 	failed_alignments = []		# [ (read id, read_seq, 5'ss site, alignment start in read, alignment end in read, alignment is reverse-complementary, filter that it failed or None)...] }
# 	out_reads = []
# 	# read_fasta = Fasta(f'{output_base}unmapped_reads.fa', as_raw=True)

# 	# Check if the 5bp upstream of the alignment in the read matches the 5bp upstream of the 5'ss in the genome. 
# 	# If it does NOT, add the read alignment to fivep_pass
# 	# { is reverse: [(5'ss site), (alignment start position in read, alignment end position in read, is reverse-complementary)...], is not reverse: [...] }
# 	for read_id in alignments:
# 		read_seq = read_seqs[read_id]
# 		fivep_pass = {True:[], False:[]}
# 		for site in alignments[read_id]:
# 			read_fivep_start, read_fivep_end, read_is_reverse = alignments[read_id][site]
# 			if read_is_reverse:
# 				read_upstream = read_seq[read_fivep_end+1:read_fivep_end+6].upper()
# 				upstream_mismatch = read_upstream != reverse_complement(fivep_upstream_seqs[site])
# 			else:
# 				read_upstream = read_seq[read_fivep_start-5:read_fivep_start].upper()
# 				upstream_mismatch = read_upstream != fivep_upstream_seqs[site]
				
# 			if upstream_mismatch:
# 				fivep_pass[read_is_reverse].append((site, alignments[read_id][site]))
# 			else:
# 				failed_alignments.append((read_id, read_seq, site, read_fivep_start, read_fivep_end, read_is_reverse, '5bp_up_match'))

# 		# For each orientation, trim off the upstream-most 5'ss and everything upstream of it, then write the trimmed sequence + alignments to file
# 		for read_is_reverse in fivep_pass:
# 			# Check if there are no alignments for the read in the given orientation
# 			if len(fivep_pass[read_is_reverse]) == 0:
# 				continue

# 			if read_is_reverse:
# 				# Get the start and end of the rightmost alignment in the read 
# 				read_fivep_start, read_fivep_end, _ = max(fivep_pass[read_is_reverse], key=lambda fp:fp[1][0])[1]
# 				read_bp_pos = read_fivep_end
# 				# Trim off the rightmost alignment and everything to the left of it
# 				trim_seq = read_seq[read_fivep_end:]
# 				# Get sequence of rightmost alignment
# 				fivep_seq = reverse_complement(read_seq[read_fivep_start:read_fivep_end])
# 			else:
# 				# Get the start and end of the leftmost alignment in the read 
# 				read_fivep_start, read_fivep_end, _ = min(fivep_pass[read_is_reverse], key=lambda fp:fp[1][0])[1]
# 				read_bp_pos = read_fivep_start - 1
# 				# Trim off the leftmost alignment and everything to the right of it
# 				trim_seq = read_seq[:read_fivep_start]
# 				# Get sequence of leftmost alignment
# 				fivep_seq = read_seq[read_fivep_start:read_fivep_end]

# 			# Keep the subset of 5'ss alignments that start at the upstream-most position and fail the rest
# 			fivep_pass_sub = []
# 			for fp in fivep_pass[read_is_reverse]:
# 				if fp[1][0] == read_fivep_start:
# 					fivep_pass_sub.append(fp)
# 				else:
# 					failed_alignments.append((read_id, read_seq, fp[0], fp[1][0], fp[1][1], fp[1][2], 'furthest_upstream'))

# 			# Check if less than 20bp is left in the read
# 			if len(trim_seq) < 20:
# 				for fp in fivep_pass_sub:
# 					failed_alignments.append((read_id, read_seq, fp[0], fp[1][0], fp[1][1], fp[1][2], 'enough_trim_seq'))
# 				continue
			
# 			# Add reads + alignment(s) that passed filtering to out_reads 
# 			out_rid = read_id + '_rev' if read_is_reverse else read_id + '_for'
# 			fivep_sites = sorted([fp[0] for fp in fivep_pass_sub])
# 			fivep_sites = ','.join(fivep_sites)
# 			out_reads.append((trim_seq, out_rid, read_seq, fivep_seq, fivep_sites, read_is_reverse, read_fivep_start, read_fivep_end, read_bp_pos))

# 	# Write the filtered alignments and trimmed sequences to file
# 	with out_lock:
# 		with open(f'{output_base}fivep_info_table.tsv', 'a') as info_out, open(f'{output_base}fivep_mapped_reads_trimmed.fa', 'a') as trimmed_out:
# 			for row in out_reads:
# 				row = [str(item) for item in row]
# 				info_out.write('\t'.join(row[1:]) + '\n')

# 				trim_seq = row[0]
# 				out_rid = row[1]
# 				trimmed_out.write(f'>{out_rid}\n{trim_seq}\n')
	
# 	# Write the failed alignments to file
# 	with failed_out_lock:
# 		failed_alignments = pd.DataFrame(failed_alignments, columns=FAILED_ALIGNMENTS_COLS)
# 		failed_alignments.to_csv(f'{output_base}failed_fivep_alignments.tsv', mode='a', sep='\t', index=False, header=False)

def decide_chunk_ranges(n_aligns:int, threads:int):
	# If there are only a handful of alignments just go through them all in one thread
	if n_aligns < 2*threads:
		return [[0, n_aligns]]
	
	# Determine the range of lines to assign to each thread
	chunk_size = int(n_aligns / threads)
	chunk_ranges = [[t*chunk_size, (t+1)*chunk_size - 1]  for t in range(threads)]
	chunk_ranges[-1][-1] = n_aligns

	return chunk_ranges


def parse_line(sam_line:str):
	fivep_site, flag, read_id, read_fivep_start = sam_line.split('\t')[:4]

	read_fivep_start = int(read_fivep_start)-1
	read_fivep_end = read_fivep_start+20
	bit_flags = bin(int(flag))
	read_is_reverse = True if len(bit_flags)>=7 and bit_flags[-5]=='1' else False

	return read_id, fivep_site, read_fivep_start, read_fivep_end, read_is_reverse


def yield_read_aligns(fivep_to_reads:str, chunk_start:int, chunk_end:int, n_aligns:int):
	with open(fivep_to_reads) as align_file:
		# Run to the starting line
		for _ in range(chunk_start-2):
			next(align_file)
		
		# Check the line before the starting line to see if the chunk starts in the middle of a read's collection of alignments  
		# If it is, move the starting line up until it reaches the next read's alignments
		previous_line_rid = align_file.readline().split('\t')[2]
		start_line = align_file.readline()
		while start_line.split('\t')[2] == previous_line_rid:
			start_line = align_file.readline()
			chunk_start += 1

		# Add the relevant fivep site info to fivep_sites
		read_id, fivep_site, read_fivep_start, read_fivep_end, read_is_reverse = parse_line(start_line)
		# This dict will hold all of the read's 5'ss alignments, seperated into reverse(True) and forward(False) alignments
		fivep_sites = {True: [], False: []}
		fivep_sites[read_is_reverse].append((fivep_site, read_fivep_start, read_fivep_end))

		line_num = chunk_start
		while line_num < n_aligns:
			line_num += 1
			# Get the next line's alignment info
			next_line_rid, fivep_site, read_fivep_start, read_fivep_end, read_is_reverse = parse_line(align_file.readline())

			# If we're still in the same read's collection of alignments, add the info to fivep_sites
			if next_line_rid == read_id:
				fivep_sites[read_is_reverse].append((fivep_site, read_fivep_start, read_fivep_end))

			# If we've reached the first alignment for a new read...
			else:
				# Yield the reverse-aligning 5'ss for filtering
				yield read_id, True, fivep_sites[True]
				# Then yield the forward-aligning 5'ss for filtering
				yield read_id, False, fivep_sites[False]

				# If we're at or have passed the end of the assigned chunk, we're done
				# We don't do this check until we know we got all of the last read's alignments
				if line_num >= chunk_end:
					break

				# Set to processing next read's alignments
				read_id = next_line_rid
				fivep_sites = {True:[], False:[]}
				fivep_sites[read_is_reverse].append((fivep_site, read_fivep_start, read_fivep_end))

		# If this is the last chunk and it reached the end of the file, make sure to yield the last read's alignments
		if line_num == n_aligns:
				yield read_id, True, fivep_sites[True]
				yield read_id, False, fivep_sites[False]



def filter_reads_chunk(fivep_to_reads:str, chunk_start:int, chunk_end:int, n_aligns:int, read_seqs:dict, fivep_upstream_seqs:dict, output_base:str) -> None:
	'''
	Filter and trim the reads to which 5'ss sequences were mapped
	Write trimmed read information and their aligned 5' splice site(s) to fivep_info_table_out.txt
	Write trimmed read sequences to fivep_mapped_reads_trimmed.fa
	'''
	failed_alignments = []		
	out_reads = []
	
	# Go through the assigned chunk of lines, reading and processing all 5'ss alignments to each read together
	# We will split each read into a forward read (with all the 5'ss aligning to the forward sequence) and a reverse read (with all the 5'ss alinging to the reverse sequence)
	for read_id, read_is_reverse, fivep_sites in yield_read_aligns(fivep_to_reads, chunk_start, chunk_end, n_aligns):
		if len(fivep_sites) == 0:
			continue

		read_seq = read_seqs[read_id]
		fivep_pass = []
		for fivep_site, read_fivep_start, read_fivep_end in fivep_sites:
			# Check if the 5bp upstream of the alignment in the read matches the 5bp upstream of the 5'ss in the genome. 
			# If it does NOT, add the read alignment to fivep_pass
			if read_is_reverse:
				read_upstream = read_seq[read_fivep_end+1:read_fivep_end+6].upper()
				upstream_mismatch = read_upstream != reverse_complement(fivep_upstream_seqs[fivep_site])
			else:
				read_upstream = read_seq[read_fivep_start-5:read_fivep_start].upper()
				upstream_mismatch = read_upstream != fivep_upstream_seqs[fivep_site]
				
			if upstream_mismatch:
				fivep_pass.append((fivep_site, read_fivep_start, read_fivep_end))
			else:
				failed_alignments.append((read_id, read_seq, fivep_site, read_fivep_start, read_fivep_end, read_is_reverse, '5bp_up_match'))

		# Check if there are no alignments for the read in the given orientation
		if len(fivep_pass) == 0:
			continue
			
		# For each orientation, trim off the upstream-most 5'ss and everything upstream of it, then write the trimmed sequence + alignments to file
		if read_is_reverse:
			# Get the start and end of the rightmost alignment in the read 
			_, furthest_fivep_start, furthest_fivep_end = max(fivep_pass, key=lambda fp:fp[1])
			read_bp_pos = furthest_fivep_end
			# Trim off the rightmost alignment and everything to the left of it
			trim_seq = read_seq[furthest_fivep_end:]
			# Get sequence of rightmost alignment
			fivep_seq = reverse_complement(read_seq[furthest_fivep_start:furthest_fivep_end])
		else:
			# Get the start and end of the leftmost alignment in the read 
			_, furthest_fivep_start, furthest_fivep_end = min(fivep_pass, key=lambda fp:fp[1])
			read_bp_pos = furthest_fivep_start - 1
			# Trim off the leftmost alignment and everything to the right of it
			trim_seq = read_seq[:furthest_fivep_start]
			# Get sequence of leftmost alignment
			fivep_seq = read_seq[furthest_fivep_start:furthest_fivep_end]

		# Keep the subset of 5'ss alignments that start at the upstream-most position and fail the rest
		fivep_pass_sub = []
		for fp in fivep_pass:
			if fp[1] == furthest_fivep_start:
				fivep_pass_sub.append(fp)
			else:
				failed_alignments.append((read_id, read_seq, *fp, read_is_reverse, 'furthest_upstream'))

		# Check if less than 20bp is left in the read
		if len(trim_seq) < 20:
			for fp in fivep_pass_sub:
				failed_alignments.append((read_id, read_seq, *fp, read_is_reverse, 'enough_trim_seq'))
			continue
		
		# Add reads + alignment(s) that passed filtering to out_reads 
		out_rid = read_id + '_rev' if read_is_reverse else read_id + '_for'
		fivep_sites = sorted([fp[0] for fp in fivep_pass_sub])
		fivep_sites = ','.join(fivep_sites)
		out_reads.append((trim_seq, out_rid, read_seq, fivep_seq, fivep_sites, read_is_reverse, furthest_fivep_start, furthest_fivep_end, read_bp_pos))

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
	# Get logger
	log = logging.getLogger() 
	log.setLevel('DEBUG')
	handler = logging.StreamHandler(sys.stdout)
	handler.setLevel('DEBUG')
	log.addHandler(handler)

	threads, genome_fasta, fivep_fasta, output_base = sys.argv[1:]
	log.debug(f'Args recieved: {sys.argv[1:]}')
	threads = int(threads)
	fivep_to_reads = f'{output_base}fivep_to_reads.sam'

	# read_ids = pysam.AlignmentFile(f'{output_base}fivep_to_reads.bam', 'r').references
	# Break the read ids tuple into chunks for processing
	# n_chunks = math.ceil(len(read_ids)/READ_CHUNKSIZE)
	# read_ids = np.array_split(read_ids, n_chunks)

	# Load reference info
	# alignments = parse_alignments(f'{output_base}fivep_to_reads.sam')
	# alignments = parse_alignments(f'{output_base}reads_to_fivep.sam')
	read_seqs = get_read_seqs(f'{output_base}unmapped_reads.fa')
	# fivep_upstream_seqs = get_fivep_upstream_seqs(alignments, genome_fasta)
	fivep_upstream_seqs = get_fivep_upstream_seqs(fivep_fasta, genome_fasta)

	# Get the total number of alignments
	with open(fivep_to_reads) as sam:
		n_aligns = sum(1 for _ in sam)
	log.debug(f'{n_aligns:,} read-fivep alignments')

	chunk_ranges = decide_chunk_ranges(n_aligns, threads)
	log.debug(f'chunk ranges: {chunk_ranges}')

	# # Record read count
	# # rids = set(read_id.split('/')[0] for read_id in alignments.keys())
	# with open(f'{output_base}read_counts.tsv', 'a') as a:
	# 	a.write(f'fivep_mapped_reads\t{len(read_ids)}\n')

	# if len(read_ids)==0:
	# 	print(time.strftime('%m/%d/%y - %H:%M:%S') + '| No reads remaining')
	# 	exit()

	# Prep out files with a header row
	with open(f'{output_base}fivep_info_table.tsv', 'w') as w:
		w.write('\t'.join(FIVEP_INFO_TABLE_COLS) + '\n')
	with open(f'{output_base}fivep_mapped_reads_trimmed.fa', 'w') as w:
		pass
	with open(f'{output_base}failed_fivep_alignments.tsv', 'w') as w:
		w.write('\t'.join(FAILED_ALIGNMENTS_COLS) + '\n')

	# rid_list = list(alignments.keys())
	# processes = []
	# for i in range(threads):
	# 	rid_subset = rid_list[i::threads]
	# 	alignments_subset = {read_id:alignments[read_id] for read_id in rid_subset}
	# 	if len(alignments_subset)==0:
	# 		continue 
	# 	subset_process = multiprocessing.Process(target=filter_fivep_reads, args=(alignments_subset, read_seqs, fivep_upstream_seqs, output_base,))
	# 	subset_process.start()
	# 	processes.append(subset_process)
		
	# for p in processes:
	# 	p.join()

	log.debug(f'Processing')
	processes = []
	for chunk_start, chunk_end in chunk_ranges:
		process = multiprocessing.Process(target=filter_reads_chunk, args=(fivep_to_reads, chunk_start, chunk_end, n_aligns, read_seqs, fivep_upstream_seqs, output_base,))
		process.start()
		processes.append(process)

	for process in processes:
		process.join()
		if process.exitcode != 0:
			raise RuntimeError()
		
	log.debug('End of script')