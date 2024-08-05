from dataclasses import dataclass, field
import sys
import subprocess
import multiprocessing
import os

from pyfaidx import Fasta
import pandas as pd
import numpy as np

import functions



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
TAILS_COLS = ['read_id',
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
						'filter_failed',
						]

out_lock = multiprocessing.Lock()
failed_out_lock = multiprocessing.Lock()



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def get_read_seqs(unmapped_fasta:str) -> dict:
	read_fasta = Fasta(unmapped_fasta, as_raw=True)
	read_seqs = {}
	for read_id in read_fasta.keys():
		read_seqs[read_id] = read_fasta[read_id][:]
	
	return read_seqs


def get_fivep_upstream_seqs(fivep_fasta:str, genome_fasta:str, log) -> dict:
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
	bedtools_output = functions.run_command(bedtools_call, input=bedtools_input, log=log)
	bedtools_output = bedtools_output.split('\n')

	# Parse output
	fivep_upstream_seqs = dict([(l.split('\t')[0][:-3], l.split('\t')[1].upper()) for l in bedtools_output])

	return fivep_upstream_seqs


def decide_chunk_ranges(n_aligns:int, threads:int):
	'''
	Returns [[chunk_1_start, chunk_1_end],... [chunk_t_start, n_aligns]] where t = threads
	Start and end positions are 1-based inclusive
	'''
	# If there are only a handful of alignments just go through them all in one thread
	if n_aligns <= 2*threads:
		return [[1, n_aligns]]
	
	# Determine the range of lines to assign to each thread
	chunk_size = int(n_aligns / threads)
	chunk_ranges = [[t*chunk_size+1, (t+1)*chunk_size]  for t in range(threads)]
	chunk_ranges[-1][-1] = n_aligns

	return chunk_ranges


def parse_line(sam_line:str):
	fivep_site, flag, read_id, read_fivep_start = sam_line.split('\t')[:4]

	read_fivep_start = int(read_fivep_start)-1
	read_fivep_end = read_fivep_start+20
	read_is_reverse = functions.align_is_reverse(flag)

	return read_id, fivep_site, read_fivep_start, read_fivep_end, read_is_reverse


def yield_read_aligns(fivep_to_reads:str, chunk_start:int, chunk_end:int, n_aligns:int):
	with open(fivep_to_reads) as align_file:
		# Run to the line directly before the starting line
		for _ in range(chunk_start-2):
			next(align_file)
		
		if chunk_start == 1:
			start_line = align_file.readline()
		# Check the line before the starting line to see if the chunk starts in the middle of a read's collection of alignments  
		# If it is, move the starting line up until it reaches the next read's alignments
		else:
			previous_line_rid = align_file.readline().split('\t')[2]
			start_line = align_file.readline()
			while start_line.split('\t')[2] == previous_line_rid:
				start_line = align_file.readline()
				chunk_start += 1

		# Add the relevant fivep site info to fivep_sites
		current_read_id, fivep_site, read_fivep_start, read_fivep_end, read_is_reverse = parse_line(start_line)
		# This dict will hold all of the read's 5'ss alignments, seperated into reverse(True) and forward(False) alignments
		fivep_sites = {True: [], False: []}
		fivep_sites[read_is_reverse].append((fivep_site, read_fivep_start, read_fivep_end))

		line_num = chunk_start
		while line_num < n_aligns:
			line_num += 1
			# Get the next line's alignment info
			line_rid, fivep_site, read_fivep_start, read_fivep_end, read_is_reverse = parse_line(align_file.readline())

			# If we're still in the same read's collection of alignments, add the info to fivep_sites
			if line_rid == current_read_id:
				fivep_sites[read_is_reverse].append((fivep_site, read_fivep_start, read_fivep_end))

			# If we've reached the first alignment for a new read...
			else:
				# Yield the reverse-aligning 5'ss for filtering
				yield current_read_id, True, fivep_sites[True]
				# Then yield the forward-aligning 5'ss for filtering
				yield current_read_id, False, fivep_sites[False]

				# Set to processing next read's alignments
				current_read_id = line_rid
				fivep_sites = {True:[], False:[]}
				fivep_sites[read_is_reverse].append((fivep_site, read_fivep_start, read_fivep_end))

				# If we're at or have passed the end of the assigned chunk, we're done
				# We don't do this check until we know we got all of the last read's alignments, 
				# so we might process a few lines after the assigned chunk_end 
				if line_num >= chunk_end:
					break

		# Yield the last read's alignments
		yield current_read_id, True, fivep_sites[True]
		yield current_read_id, False, fivep_sites[False]


def filter_reads_chunk(fivep_to_reads:str, chunk_start:int, chunk_end:int, n_aligns:int, read_seqs:dict, fivep_upstream_seqs:dict, output_base:str, log_level) -> None:
	'''
	Filter and trim the reads to which 5'ss sequences were mapped, the trimmed section being the read "head" and the trimmed-off section being the read "tail" 
	Write read information and their aligned 5' splice site(s) to tails.tsv
	Write head sequences to heads.fa
	'''
	# We have to set the log level in each process because the children don't inherit the log level from their parent,
	# even if you pass the log object itself
	log = functions.get_logger(log_level)
	log.debug(f'Process {os.getpid()}: Born and assigned lines {chunk_start:,}-{chunk_end:,}')

	# Go through the assigned chunk of lines, reading and processing all 5'ss alignments to each read together
	# We will split each read into a forward read (with all the 5'ss aligning to the forward sequence) and a reverse read (with all the 5'ss alinging to the reverse sequence)
	failed_alignments = []		
	out_reads = []
	for read_id, read_is_reverse, fivep_sites in yield_read_aligns(fivep_to_reads, chunk_start, chunk_end, n_aligns):
		if len(fivep_sites) == 0:
			continue

		read_seq = read_seqs[read_id]
		fivep_pass = []
		for fivep_site, read_fivep_start, read_fivep_end in fivep_sites:
			# Check if the 5bp upstream of the alignment in the read matches the 5bp upstream of the 5'ss in the genome. 
			# If it does NOT, add the read alignment to fivep_pass
			if read_is_reverse:
				read_upstream = read_seq[read_fivep_end:read_fivep_end+5].upper()
				upstream_mismatch = read_upstream != functions.reverse_complement(fivep_upstream_seqs[fivep_site])
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
			
		# For each orientation, trim off the upstream-most 5'ss and everything upstream of it, then write the sequence + alignments to file
		if read_is_reverse:
			# Get the start and end of the rightmost alignment in the read 
			_, furthest_fivep_start, furthest_fivep_end = max(fivep_pass, key=lambda fp:fp[1])
			read_bp_pos = furthest_fivep_end
			# Trim off the rightmost alignment and everything to the left of it
			head_seq = read_seq[furthest_fivep_end:]
			# Get sequence of rightmost alignment
			fivep_seq = functions.reverse_complement(read_seq[furthest_fivep_start:furthest_fivep_end])
		else:
			# Get the start and end of the leftmost alignment in the read 
			_, furthest_fivep_start, furthest_fivep_end = min(fivep_pass, key=lambda fp:fp[1])
			read_bp_pos = furthest_fivep_start - 1
			# Trim off the leftmost alignment and everything to the right of it
			head_seq = read_seq[:furthest_fivep_start]
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
		if len(head_seq) < 20:
			for fp in fivep_pass_sub:
				failed_alignments.append((read_id, read_seq, *fp, read_is_reverse, 'enough_head_seq'))
			continue
		
		# Add reads + alignment(s) that passed filtering to out_reads 
		out_rid = read_id + '_rev' if read_is_reverse else read_id + '_for'
		fivep_sites = sorted([fp[0] for fp in fivep_pass_sub])
		fivep_sites = ','.join(fivep_sites)
		out_reads.append((head_seq, out_rid, read_seq, fivep_seq, fivep_sites, read_is_reverse, furthest_fivep_start, furthest_fivep_end, read_bp_pos))

	# Write the filtered alignments and head sequences to file
	with out_lock:
		with open(f'{output_base}tails.tsv', 'a') as info_out, open(f'{output_base}heads.fa', 'a') as seq_out:
			for row in out_reads:
				row = [str(item) for item in row]
				info_out.write('\t'.join(row[1:]) + '\n')

				head_seq = row[0]
				out_rid = row[1]
				seq_out.write(f'>{out_rid}\n{head_seq}\n')
	
	# Write the failed alignments to file
	with failed_out_lock:
		failed_alignments = pd.DataFrame(failed_alignments, columns=FAILED_ALIGNMENTS_COLS)
		failed_alignments.to_csv(f'{output_base}failed_fivep_alignments.tsv', mode='a', sep='\t', index=False, header=False)

	log.debug(f'Process {os.getpid()}: Chunk finished')



# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__' :
	# Get args
	threads, genome_fasta, fivep_fasta, output_base, log_level = sys.argv[1:]

	# Get logger
	log = functions.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}')
	
	threads = int(threads)
	fivep_to_reads = f'{output_base}fivep_to_reads.sam'

	# Load reference info
	read_seqs = get_read_seqs(f'{output_base}unmapped_reads.fa')
	fivep_upstream_seqs = get_fivep_upstream_seqs(fivep_fasta, genome_fasta, log)

	# Get the total number of alignments
	with open(fivep_to_reads) as sam:
		n_aligns = sum(1 for _ in sam)
	log.debug(f'{n_aligns:,} read-fivep alignments')

	if n_aligns == 0:
		log.info('No reads remaining')
		exit()

	chunk_ranges = decide_chunk_ranges(n_aligns, threads)
	log.debug(f'chunk ranges: {chunk_ranges}')

	# Prep out files with a header row
	with open(f'{output_base}tails.tsv', 'w') as w:
		w.write('\t'.join(TAILS_COLS) + '\n')
	with open(f'{output_base}heads.fa', 'w') as w:
		pass
	with open(f'{output_base}failed_fivep_alignments.tsv', 'w') as w:
		w.write('\t'.join(FAILED_ALIGNMENTS_COLS) + '\n')

	log.debug(f'Parallel processing {len(chunk_ranges):,} chunks...')
	processes = []
	for chunk_start, chunk_end in chunk_ranges:
		process = multiprocessing.Process(target=filter_reads_chunk, args=(fivep_to_reads, chunk_start, chunk_end, n_aligns, read_seqs, fivep_upstream_seqs, output_base, log_level, ))
		process.start()
		processes.append(process)

	# Check if any processes hit an error
	for process in processes:
		process.join()
		if process.exitcode != 0:
			raise RuntimeError()
		
	log.debug('End of script')