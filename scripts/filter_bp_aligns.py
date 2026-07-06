import sys
import os
import multiprocessing as mp

from pyfaidx import Fasta
import pandas as pd

import utils





# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
MIN_TAIL_LEN = 20

# In files
UNMAPPED_READS_FILE = "{}unmapped_reads.fa"
BP_TO_READS_FILE = "{}bp_to_reads.sam"
# Out files
FAILED_BPS_FILE = "{}failed_bp_alignments.tsv"
HEADS_FILE = "{}heads.tsv"
TAILS_SEQ_FILE = "{}tails.fa"

HEADS_COLS = ['read_id',
			  'read_seq',
			  'bp_seq',
			  'bp_sites',
			  'read_orient_to_gene',
			  'read_bp_start',
			  'read_bp_end',
			  'read_bp_pos',
			  ]
FAILED_BPS_COLS = ['read_id',
				   'read_seq',
				   'bp_site',
				   'read_bp_start',
				   'read_bp_end',
				   'alignment_score',
				   'read_orient_to_gene',
				   'filter_failed',
				   ]


out_lock = mp.Lock()
failed_out_lock = mp.Lock()



# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def get_read_seqs(unmapped_fasta:str, read_ids_to_parse: set) -> dict:

	read_seqs = {}
	with open(unmapped_fasta) as in_file:
		read_id, read_seq = in_file.readline().strip(), in_file.readline().strip()
		while read_id:
			read_id = read_id[1:]
			if read_id in read_ids_to_parse:
				read_seqs[read_id] = read_seq.upper()
			read_id, read_seq = in_file.readline().strip(), in_file.readline().strip()
	
	return read_seqs


def get_bp_downstream_seqs(bp_sites_to_parse:set, genome_fasta:str, log) -> dict:
	# Prepare bp sites for input into bedtools getfasta
	# We're retrieving the 5bp downstream of the branchpoint
	bedtools_input = ''
	chrom_lengths = utils.get_chrom_lengths(f'{genome_fasta}.fai')
	for site in bp_sites_to_parse:
		chrom, bp_pos, strand = site.split(';')
		bp_pos = int(bp_pos)
		if strand == '+':
			window_start_pos, window_end_pos = bp_pos-4, min(bp_pos+6, chrom_lengths[chrom])
		else:
			window_start_pos, window_end_pos = max(bp_pos-5, 0), bp_pos+5
		bedtools_input += f'{chrom}\t{window_start_pos}\t{window_end_pos}\t{site}\t0\t{strand}\n'

	# Call bedtools getfasta 
	bp_downstream_seqs = utils.getfasta(genome_fasta, bedtools_input, log)
	bp_downstream_seqs = bp_downstream_seqs.set_index('name', drop=True)
	# Make it a dict
	bp_downstream_seqs['downstream_seq'] = bp_downstream_seqs['seq'].str.slice(start=5)
	bp_downstream_seqs['bp_seq'] = bp_downstream_seqs['seq'].str.slice(stop=5)
	bp_downstream_dict = bp_downstream_seqs['downstream_seq'].to_dict()
	bp_seq_dict = bp_downstream_seqs['bp_seq'].to_dict()

	return bp_seq_dict, bp_downstream_dict


def parse_line(sam_line:str):
	alignment_info = sam_line.split('\t')
	bp_site, flag, read_id, read_bp_start, _, cigar_str = alignment_info[:6]
	alignment_score = int(alignment_info[11].removeprefix('AS:i:'))

	cigar_tuples = utils.get_cigar_tuples(cigar_str)
	read_bp_start = int(read_bp_start)-1
	read_bp_end = utils.get_align_end(read_bp_start, cigar_tuples)
	left_clip = cigar_tuples[0][0] if cigar_tuples[0][1]=='S' else 0
	right_clip = cigar_tuples[-1][0] if cigar_tuples[-1][1]=='S' else 0
	read_orient_to_gene = utils.align_orient(flag)

	return read_id, bp_site, read_bp_start, read_bp_end, read_orient_to_gene, left_clip, right_clip, alignment_score

def yield_read_aligns(bp_to_reads:str, chunk_start:int, chunk_end:int, n_aligns:int):
	with open(bp_to_reads) as align_file:
		# Run to the line directly before the starting line
		for _ in range(chunk_start-2):
			next(align_file)
		
		if chunk_start == 1:
			start_align = align_file.readline()
		# Check the line before the starting line to see if the chunk starts in the middle of a read's collection of alignments  
		# If it is, move the starting line up until it reaches the next read's alignments
		else:
			previous_align_rid = align_file.readline().split('\t')[2]
			start_align = align_file.readline()
			while start_align.split('\t')[2] == previous_align_rid:
				start_align = align_file.readline()
				chunk_start += 1

		# Add the relevant BP site info to bp_sites
		current_read_id, bp_site, read_bp_start, read_bp_end, read_orient_to_gene, left_clip, right_clip, alignment_score = parse_line(start_align)
		# This dict will hold all of the read's BP alignments, seperated into reverse and forward alignments
		bp_sites = {"Forward": [], "Reverse": []}
		bp_sites[read_orient_to_gene].append((bp_site, read_bp_start, read_bp_end, left_clip, right_clip, alignment_score))

		align_num = chunk_start
		while align_num < n_aligns:
			align_num += 1
			# Get the next line's alignment info
			align_rid, bp_site, read_bp_start, read_bp_end, read_orient_to_gene, left_clip, right_clip, alignment_score = parse_line(align_file.readline())

			# If we're still in the same read's collection of alignments, add the info to bp_sites
			if align_rid == current_read_id:
				bp_sites[read_orient_to_gene].append((bp_site, read_bp_start, read_bp_end, left_clip, right_clip, alignment_score))

			# If we've reached the first alignment for a new read...
			else:
				# Yield the reverse-aligning BPs for filtering
				yield current_read_id, "Reverse", bp_sites["Reverse"]
				# Then yield the forward-aligning BPs for filtering
				yield current_read_id, "Forward", bp_sites["Forward"]

				# Set to processing next read's alignments
				current_read_id = align_rid
				bp_sites = {"Forward":[], "Reverse":[]}
				bp_sites[read_orient_to_gene].append((bp_site, read_bp_start, read_bp_end, left_clip, right_clip, alignment_score))

				# If we're at or have passed the end of the assigned chunk, we're done
				# We don't do this check until we know we got all of the last read's alignments, 
				# so we might process a few lines after the assigned chunk_end 
				if align_num > chunk_end:
					break

		# Yield the last read's alignments
		if align_num == n_aligns:
			yield current_read_id, "Reverse", bp_sites["Reverse"]
			yield current_read_id, "Forward", bp_sites["Forward"]

def adjust_align_coords(read_bp_start, read_bp_end, read_adj_pos, read_orient_to_gene):
	if read_orient_to_gene == 'Forward':
		return read_bp_start, read_adj_pos
	else:
		return read_adj_pos, read_bp_end
	
DNTS = set(['GT'])
def find_fivep_dnt(bp_window_seq, read_seq, read_bp_start, read_bp_end, read_orient_to_gene):

	if read_orient_to_gene == 'Forward':
		read_pos_range, read_window_seq_for = range(read_bp_end, read_bp_end+5), read_seq[read_bp_end:read_bp_end+7]
		adj_offset_inc, adj_offset_exc = 0, 1
	else:
		read_pos_range, read_window_seq_for = range(read_bp_start-1, read_bp_start-6, -1), utils.reverse_complement(read_seq[read_bp_start-7:read_bp_start])
		adj_offset_inc, adj_offset_exc = 1, 0
	adj_fail, out_bp_rel_pos = False, None
	for read_pos, window_index, read_nt, bp_window_nt, bp_rel_pos in zip(read_pos_range, range(5),
																		 read_window_seq_for[:5],
																		 bp_window_seq, range(-4, 1)):
		nt_match = read_nt == bp_window_nt
		if nt_match:
			if bp_rel_pos == 0:
				if read_window_seq_for[window_index+1:window_index+3] in DNTS:
					read_bp_start, read_bp_end = adjust_align_coords(read_bp_start, read_bp_end, read_pos+adj_offset_exc, read_orient_to_gene)
					out_bp_rel_pos = 0
				else:
					adj_fail = True
		else:
			if bp_rel_pos == 0:
				if read_window_seq_for[window_index+1:window_index+3] in DNTS:
					read_bp_start, read_bp_end = adjust_align_coords(read_bp_start, read_bp_end, read_pos+adj_offset_exc, read_orient_to_gene)
					out_bp_rel_pos = 0
				elif read_window_seq_for[window_index:window_index+2] in DNTS:
					read_bp_start, read_bp_end = adjust_align_coords(read_bp_start, read_bp_end, read_pos+adj_offset_inc, read_orient_to_gene)
					out_bp_rel_pos = -1
				else:
					adj_fail = True
			else:
				if read_window_seq_for[window_index:window_index+2] in DNTS:
					read_bp_start, read_bp_end = adjust_align_coords(read_bp_start, read_bp_end, read_pos+adj_offset_inc, read_orient_to_gene)
					out_bp_rel_pos = bp_rel_pos-1
				elif read_window_seq_for[window_index+1:window_index+3] in DNTS:
					read_bp_start, read_bp_end = adjust_align_coords(read_bp_start, read_bp_end, read_pos+adj_offset_exc, read_orient_to_gene)
					out_bp_rel_pos = bp_rel_pos
				else:
					adj_fail = True
				break
	
	return read_bp_start, read_bp_end, adj_fail, out_bp_rel_pos

def filter_reads_chunk(chunk_start:int, chunk_end:int, n_aligns:int, read_seqs:dict, bp_window_seqs:dict, bp_downstream_seqs:dict, strand:str, output_base:str, log_level:str) -> None:
	'''
	Filter and trim the reads to which BP sequences were mapped, the trimmed section being the read "head" and the trimmed-off section being the read "tail" 
	Write read information and their aligned branchpoint site(s) to heads.tsv
	Write head sequences to heads.fa
	'''
	# We have to set the log level in each process because the children don't inherit the log level from their parent,
	# even if you pass the log object itself
	log = utils.get_logger(log_level)
	log.debug(f'Process {os.getpid()}: Born and assigned lines {chunk_start:,}-{chunk_end:,}')

	# Go through the assigned chunk of lines, reading and processing all BP alignments to each read together
	# We will split each read into a forward read (with all the BP aligning to the forward sequence) and a reverse read (with all the BP aligning to the reverse sequence)
	failed_alignments = []		
	out_reads = []
	for read_id, read_orient_to_gene, bp_sites in yield_read_aligns(BP_TO_READS_FILE.format(output_base), chunk_start, chunk_end, n_aligns):
		if len(bp_sites) == 0:
			continue

		read_num = read_id[-1]
		read_seq = read_seqs[read_id]
		bp_pass = []
		for bp_site, read_bp_start, read_bp_end, left_clip, right_clip, alignment_score in bp_sites:
			# If the sequencing data is strand-specific, we can rule out alignments on 
			# the wrong strand for the read mate
			if strand == 'First':
				if (read_num=='1' and read_orient_to_gene == 'Reverse') or (read_num=='2' and read_orient_to_gene == 'Forward'):
					failed_alignments.append((read_id, read_seq, bp_site, read_bp_start, read_bp_end, alignment_score, read_orient_to_gene, 'wrong_strand'))
			elif strand == 'Second':
				if (read_num=='1' and read_orient_to_gene == 'Forward') or (read_num=='2' and read_orient_to_gene == 'Reverse'):
					failed_alignments.append((read_id, read_seq, bp_site, read_bp_start, read_bp_end, alignment_score, read_orient_to_gene, 'wrong_strand'))

			if read_orient_to_gene == "Forward" and len(read_seq)-read_bp_end<MIN_TAIL_LEN:
				failed_alignments.append((read_id, read_seq, bp_site, read_bp_start, read_bp_end, alignment_score, read_orient_to_gene, 'not_enough_tail_seq'))
				continue
			elif read_orient_to_gene == "Reverse" and read_bp_start<MIN_TAIL_LEN:
				failed_alignments.append((read_id, read_seq, bp_site, read_bp_start, read_bp_end, alignment_score, read_orient_to_gene, 'not_enough_tail_seq'))
				continue

			# Check if the 5bp downstream of the alignment in the read matches the 5bp downstream of the BP in the genome. 
			# If it does NOT, add the read alignment to bp_pass
			# If it does, add the read alignment to failed_alignments
			# In almost all cases we're comparing the 5bp downstream in the genome to the 5bp 
			# downstream of in the read, but if the BP is at the very edge of its chromosome there 
			# may be less than 5bp of sequence downstream of the BP, so we use the length of 
			# bp_downstream_seqs[bp_site] to determine how many bp to reach downstream in the read. 
			# In almost all cases it'll be 5bp

			read_bp_start, read_bp_end, adj_fail, bp_rel_pos = find_fivep_dnt(bp_window_seqs[bp_site], read_seq, read_bp_start, read_bp_end, read_orient_to_gene)
			if adj_fail:
				failed_alignments.append((read_id, read_seq, bp_site, read_bp_start, read_bp_end, alignment_score, read_orient_to_gene, 'BP_dNT_adj_fail'))
				continue

			if read_orient_to_gene == "Forward":
				segment_end = read_bp_end + len(bp_downstream_seqs[bp_site])
				read_downstream = read_seq[read_bp_end:segment_end]
				downstream_mismatch = read_downstream != bp_downstream_seqs[bp_site]
			elif read_orient_to_gene == "Reverse":
				segment_start = read_bp_start - len(bp_downstream_seqs[bp_site])
				read_downstream = read_seq[segment_start:read_bp_start]
				downstream_mismatch = read_downstream != utils.reverse_complement(bp_downstream_seqs[bp_site])
				
			if downstream_mismatch:
				#clip_length = right_clip if read_orient_to_gene=='Forward' else left_clip
				bp_site += f';{abs(bp_rel_pos)}'
				bp_pass.append((bp_site, read_bp_start, read_bp_end, alignment_score))
			else:
				failed_alignments.append((read_id, read_seq, bp_site, read_bp_start, read_bp_end, alignment_score, read_orient_to_gene, 'BP_downstream_match'))

		# Check if there are no alignments for the read in the given orientation
		if len(bp_pass) == 0:
			continue
			
		# For each orientation, trim off the most downstream BP and everything upstream of it, then write the sequence + alignments to file
		if read_orient_to_gene == "Forward":
			# Get the start and end of the rightmost alignment in the read 
			furthest_bp_site, furthest_bp_start, furthest_bp_end, furthest_alignment_score = max(bp_pass, key=lambda bp:(bp[2], bp[3]))
			read_bp_pos = furthest_bp_end-1
			# Trim off the leftmost alignment and everything to the right of it
			tail_seq = read_seq[furthest_bp_end:]
			# Get sequence of leftmost alignment
			bp_seq = read_seq[furthest_bp_start:furthest_bp_end]
		elif read_orient_to_gene == "Reverse":
			# Get the start and end of the leftmost alignment in the read 
			furthest_bp_site, furthest_bp_start, furthest_bp_end, furthest_alignment_score = min(bp_pass, key=lambda bp:(bp[1], -bp[3]))
			read_bp_pos = furthest_bp_start
			# Trim off the rightmost alignment and everything to the left of it
			tail_seq = read_seq[:furthest_bp_start]
			# Get sequence of rightmost alignment
			bp_seq = utils.reverse_complement(read_seq[furthest_bp_start:furthest_bp_end])

		# Keep the subset of BP alignments that end at the downstream-most position and fail the rest
		bp_pass_sub = []
		for bp in bp_pass:
			if (read_orient_to_gene == 'Forward' and bp[2] == furthest_bp_end) or (read_orient_to_gene == 'Reverse' and bp[1] == furthest_bp_start):
				#if bp[-1] == furthest_alignment_score:
				bp_pass_sub.append(bp)
			else:
				failed_alignments.append((read_id, read_seq, *bp, read_orient_to_gene, 'not_furthest_downstream'))

		
		# Add reads + alignment(s) that passed filtering to out_reads 
		out_rid = read_id + '_' + read_orient_to_gene[:3].lower()
		bp_sites = sorted([bp[0] for bp in bp_pass_sub])
		bp_sites = ','.join(bp_sites)
		out_reads.append((tail_seq, out_rid, read_seq, bp_seq, bp_sites, read_orient_to_gene, furthest_bp_start, furthest_bp_end, read_bp_pos))

	# Write the filtered alignments and head sequences to file
	with out_lock:
		with open(HEADS_FILE.format(output_base), 'a') as info_out, open(TAILS_SEQ_FILE.format(output_base), 'a') as seq_out:
			for row in out_reads:
				row = [str(item) for item in row]
				info_out.write('\t'.join(row[1:]) + '\n')

				tail_seq = row[0]
				out_rid = row[1]
				seq_out.write(f'>{out_rid}\n{tail_seq}\n')
	
	# Write the failed alignments to file
	with failed_out_lock:
		failed_alignments = pd.DataFrame(failed_alignments, columns=FAILED_BPS_COLS)
		failed_alignments.to_csv(FAILED_BPS_FILE.format(output_base), mode='a', sep='\t', index=False, header=False)

	log.debug(f'Process {os.getpid()}: Chunk finished')





# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__' :
	# Get args
	output_base, log_level, genome_fasta, strand, threads, = sys.argv[1:]

	# Get logger
	log = utils.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}')
	
	threads = int(threads)

	# Get the total number of alignments and which reads to pull sequences for
	read_ids_to_parse, bp_sites_to_parse, n_aligns = set(), set(), 0
	with open(BP_TO_READS_FILE.format(output_base)) as sam:
		for line in sam:
			bp_site, _, read_id = line.split('\t')[:3]
			read_ids_to_parse.add(read_id)
			bp_sites_to_parse.add(bp_site)
			n_aligns += 1
	log.debug(f'{n_aligns:,} read-bp alignments')

	# If there are no alignments, end the run early
	if n_aligns == 0:
		sys.exit(4)

	# Load reference info
	read_seqs = get_read_seqs(UNMAPPED_READS_FILE.format(output_base), read_ids_to_parse)
	del read_ids_to_parse
	bp_window_seqs, bp_downstream_seqs = get_bp_downstream_seqs(bp_sites_to_parse, genome_fasta, log)
	del bp_sites_to_parse

	# Decide how to divide the collection of alignments into chunks
	# for parallel processing, so each thread gets a roughly equal number
	chunk_ranges = utils.decide_chunk_ranges(n_aligns, threads)
	log.debug(f'chunk ranges: {chunk_ranges}')

	# Prep out files with a header row
	with open(HEADS_FILE.format(output_base), 'w') as w:
		w.write('\t'.join(HEADS_COLS) + '\n')
	with open(TAILS_SEQ_FILE.format(output_base), 'w') as w:
		pass
	with open(FAILED_BPS_FILE.format(output_base), 'w') as w:
		w.write('\t'.join(FAILED_BPS_COLS) + '\n')

	# Start parallel processes, leaving the first chunk for the main process
	log.debug(f'Parallel processing {len(chunk_ranges):,} chunks...')
	processes = []
	for chunk_start, chunk_end in chunk_ranges[1:]:
		process = mp.Process(target=filter_reads_chunk, 
							args=(chunk_start, chunk_end, n_aligns, read_seqs, 
								  bp_window_seqs, bp_downstream_seqs, strand, output_base, log_level,))
		process.start()
		processes.append(process)
		
	# Process the first chunk in the main process
	filter_reads_chunk(chunk_ranges[0][0], chunk_ranges[0][1], n_aligns, read_seqs,
					   bp_window_seqs, bp_downstream_seqs, strand, output_base, log_level) 

	# Check if any processes hit an error
	for process in processes:
		process.join()
		if process.exitcode != 0:
			raise RuntimeError(f'Error in process {process.pid}')

		
	log.debug('End of script')