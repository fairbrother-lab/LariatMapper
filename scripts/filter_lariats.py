from sys import argv
from os import rmdir, mkdir
from os.path import join, isfile, exists
from intervaltree import Interval, IntervalTree
from collections import Counter
from subprocess import run, DEVNULL
from gzip import open as gzipopen
from time import strftime



# =============================================================================#
#                                  Constants                                   #
# =============================================================================#
BASE_RESULTS_COLUMNS = ('gene',
                        'gene_ensembl_id',
                        'gene_type',
                        'read_id',
                        'read_seq',
                        'chrom',
                        'strand',
                        'fivep_pos',
                        'threep_pos',
                        'bp_pos',
                        'read_bp_nt',
                        'genomic_bp_nt',
                        'genomic_bp_context',
                        'bp_dist_to_threep',
                        'total_mapped_reads')



# =============================================================================#
#                                  Functions                                  #
# =============================================================================#
def load_intron_info(ref_introns) -> tuple:
    '''
    Returns a dict formatted as follows:
    {Chromosome: {Strand(+ or -): Intervaltree(StartPosition(int), EndPosition(int))}}
    for 3' splice sites (+/- 2 bases), 5' splice sites (+/- 2 bases), and the introns they come from (start to end)
    '''
    threep_sites, fivep_sites, introns = {}, {}, {}

    if ref_introns[-2:] == 'gz':
        intron_file = gzipopen(ref_introns, 'rt')
    else:
        intron_file = open(ref_introns)

    introns_done = set()
    for line in intron_file:
        chrom, start, end, _, _, strand = line.strip().split()
        if chrom not in threep_sites:
            threep_sites[chrom] = {s: IntervalTree() for s in ['+', '-']}
            fivep_sites[chrom] = {s: IntervalTree() for s in ['+', '-']}
            introns[chrom] = {s: IntervalTree() for s in ['+', '-']}
        intron_id = '{}_{}_{}_{}'.format(chrom, strand, start, end)
        if intron_id not in introns_done:
            start, end = int(start), int(end)
            if strand == '+':
                threep_sites[chrom][strand].add(Interval(end-2, end+2))
                fivep_sites[chrom][strand].add(Interval(start-2, start+2))
            else:
                threep_sites[chrom][strand].add(Interval(start-2, start+2))
                fivep_sites[chrom][strand].add(Interval(end-2, end+2))
            introns[chrom][strand].add(Interval(start, end))
            introns_done.add(intron_id)

    intron_file.close()

    return threep_sites, fivep_sites, introns


# def load_gene_info(ref_gtf) -> dict:
#     '''
#     Get whole-genome gene info from annotation file
#     '''
#     if ref_gtf[-2:] == 'gz':
#         gtf_file = gzipopen(ref_gtf, 'rt')
#     else:
#         gtf_file = open(ref_gtf)

#     gene_info = {}
#     for line in gtf_file:
#         if line[0] != '#':
#             chrom, _, feat, start, end, _, strand, _, annotations = line.strip().split('\t')
#             if feat == 'gene':
#                 annotations = {a.split(' ')[0]: a.split(' ')[1].replace(
#                     '\"', '') for a in annotations[:-1].split('; ')}
#                 if chrom not in gene_info:
#                     gene_info[chrom] = {s: IntervalTree() for s in ('-', '+')}
#                 gene_int_data = {
#                     'gene_name': annotations['gene_name'], 'gene_type': annotations['gene_type'], 'ensembl_id': annotations['gene_id']}
#                 gene_info[chrom][strand].add(
#                     Interval(int(start)-1, int(end), gene_int_data))
#     gtf_file.close()

#     return gene_info


def load_lariat_table(final_info_table: str) -> dict:
	'''
	For a given lariat-mapping of a fastq file, retrieve all the lariat reads from the XXX_final_info_table.tsv and put them in a dict, which
	can then be added to the experiment-wide superset dict
	'''
	lariats = {}
	with open(final_info_table, 'r') as lariat_file:
		next(lariat_file)  # Skip the header line
		for line in lariat_file:
			rid, read_seq, chrom, strand, gene_name, gene_id, gene_type, fivep_site, read_is_reverse, fivep_read_start, fivep_read_end, bp_site, read_bp_nt, genomic_bp_nt, genomic_bp_window, threep_site, bp_dist_to_threep = line.strip().split('\t')			
			fivep_site = int(fivep_site)
			threep_site = int(threep_site)
			bp_site = int(bp_site)
			bp_dist_to_threep = int(bp_dist_to_threep)

			lariats.append([rid, read_seq, chrom, strand, gene_name, gene_id, gene_type, fivep_site, read_is_reverse, fivep_read_start, fivep_read_end, bp_site, read_bp_nt, genomic_bp_nt, genomic_bp_window, threep_site, bp_dist_to_threep])
	
	return lariats


def get_mapped_reads(output_dir, output_base_name) -> int:
	'''
	Get the number of reads that mapped to the reference genome from the *_run_data.tsv file 
	'''
	read_count_path = join(output_dir, f'{output_base_name}_run_data.tsv')
	if not isfile(read_count_path):
		raise FileNotFoundError(f'No mapped read count found for "{output_base_name}"')
	with open(read_count_path, 'r') as file:
		line = file.readline()
		sample_read_count = int(line.split('\t')[1])

	return sample_read_count



def filter_lariats(lariats: dict, threep_sites: dict, fivep_sites: dict, output_dir:str, num_cpus:str, ref_b2index:str, ref_repeatmasker: str) -> dict:
	'''
	Filter the candidate lariat reads in order to exclude any that meet the following criteria:
			- BP is within 2bp of a splice site (likely an intron circle, not a lariat)
			- Read maps to UBB or UBC (likely false positive due to the repetitive nature of the genes)
			- There is a valid aligment for the 3' segment upstream of the 5' segment (likely a pre-mRNA read)
			- Both the 5'SS and the BP overlap with repetitive regions from RepeatMasker (likely false positive due to sequence repetition)
	'''
	# create temporary directory in the output dir
	temp_dir = output_dir + '/temp'
	if not exists(temp_dir):
		mkdir(temp_dir)

	# Check if both the 5'SS and the BP overlap with a repetitive region
	# Write the 5'ss and BP coordinates to BED files
	fivep_tmp_bed, bp_tmp_bed = f'{temp_dir}/fivep_tmp.bed', f'{temp_dir}/bp_tmp.bed'
	with open(fivep_tmp_bed, 'w') as fivep_out, open(bp_tmp_bed, 'w') as bp_out:
		for rid, read_seq, chrom, strand, gene_name, gene_id, gene_type, fivep_site, read_is_reverse, fivep_read_start, fivep_read_end, bp_site, read_bp_nt, genomic_bp_nt, genomic_bp_window, threep_site, bp_dist_to_threep in lariat_reads
			fivep_out.write('{}\t{}\t{}\t{}\n'.format(
				chrom, fivep_site-1, fivep_site+1, rid))
			bp_out.write('{}\t{}\t{}\t{}\n'.format(
				chrom, bp_site-1, bp_site+1, rid))
				
	# Identify 5'ss's that overlap a repeat region
	fivep_overlap_bed, bp_overlap_bed = f'{temp_dir}/fivep_repeat_overlaps.bed', f'{temp_dir}/bp_repeat_overlaps.bed'
	with open(fivep_overlap_bed, 'w') as out_file:
		run(f'bedtools intersect -u -a {fivep_tmp_bed} -b {ref_repeatmasker}'.split(' '),
			stdout=out_file)
	# Identify BPs that overlap a repeat region
	with open(bp_overlap_bed, 'w') as out_file:
		run(f'bedtools intersect -u -a {bp_tmp_bed} -b {ref_repeatmasker}'.split(' '),
			stdout=out_file)

	# Add reads where both sites overlapped to repeat_rids for removal
	fivep_repeat_rids, bp_repeat_rids = set(), set()
	with open(fivep_overlap_bed) as in_file:
		for line in in_file:
			_, _, _, rid = line.strip().split()
			fivep_repeat_rids.add(rid)
	with open(bp_overlap_bed) as in_file:
		for line in in_file:
			_, _, _, rid = line.strip().split()
			bp_repeat_rids.add(rid)
	repeat_rids = fivep_repeat_rids.intersection(bp_repeat_rids)

	# Now loop through reads and check them against filters
	failed_reads = {}
	filtered_lariats = {}
	for rid, read_seq, chrom, strand, gene_name, gene_id, gene_type, fivep_site, read_is_reverse, fivep_read_start, fivep_read_end, bp_site, read_bp_nt, genomic_bp_nt, genomic_bp_window, threep_site, bp_dist_to_threep in lariat_reads:
		passed_filtering = True
		fail_reason = None

		# Check if BP is within 2bp of an annotated splice site
		if fivep_sites[chrom][strand].overlaps(bp_site) or threep_sites[chrom][strand].overlaps(bp_site):
			passed_filtering = False
			fail_reason = 'near_ss' if fail_reason is None else fail_reason

		# Check if read mapped to a ubiquitin gene
		if gene_name in ('UBC', 'UBB'):
			passed_filtering = False
			fail_reason = 'ubiquitin_gene' if fail_reason is None else fail_reason
		
		if rid in repeat_rids:
			passed_filtering = False
			fail_reason = 'in_repeat' if fail_reason is None else fail_reason

		# If passed
		if passed_filtering:
			filtered_lariats.append([gene_name, gene_id, gene_type, rid, read_seq, chrom, strand, fivep_site, threep_site, bp_site, read_bp_nt, genomic_bp_nt, genomic_bp_window, bp_dist_to_threep])
		else:
			failed_reads.append([gene_name, gene_id, gene_type, rid, read_seq, chrom, strand, fivep_site, threep_site, bp_site, read_bp_nt, genomic_bp_nt, genomic_bp_window, bp_dist_to_threep, fail_reason])

	print(strftime('%m/%d/%y - %H:%M:%S') + f' | Pre-filter read count = {len(lariat_reads)}')
	print(strftime('%m/%d/%y - %H:%M:%S') + f' | Post-filter read count = {len(filtered_lariats)}')

	# Delete all temporary files
	# run(f'rm {seq_tmp_fa} {seq_tmp_sam} {fivep_tmp_bed} {bp_tmp_bed} {fivep_overlap_bed} {bp_overlap_bed}'.split(' '))
	# rmdir(temp_dir)

	return filtered_lariats, failed_reads


# 	with open(out_file, 'w') as out:
# 		out.write('read_id\tread_seq\tchrom\tstrand\tfivep_site\tread_is_reverse\tread_fivep_start\tread_fivep_end\tthreep_site\tbp_site\tread_bp_nt\tgenomic_bp_nt\tgenomic_bp_window\n')
		
# 		for base_rid in base_rids:
# 			align_mismatch = {'_for':{True:[], False:[]}, '_rev':{True:[], False:[]}}
# 			for orientation in ('_for', '_rev'):
# 				for align_info in filtered_alignments.get(base_rid+orientation, []):
# 					read_seq, threep_chrom, threep_strand, fivep_site, read_is_reverse, fivep_start, fivep_end, threep_site, bp_site, read_bp_nt, genomic_bp_nt, genomic_bp_window = align_info

# 					bp_mismatch = read_bp_nt != genomic_bp_nt
# 					align_mismatch[orientation][bp_mismatch].append(align_info)

# 			# Output 1 alignment for the read, prioritizing based on mistmatch and orientation, choosing a random alignment if multiple in same category
# 			if len(align_mismatch['_for'][True]) > 0:
# 				output = [base_rid] + sample(align_mismatch['_for'][True], 1)[0]
# 			elif len(align_mismatch['_rev'][True]) > 0:
# 				output = [base_rid] + sample(align_mismatch['_rev'][True], 1)[0]
# 			elif len(align_mismatch['_for'][False]) > 0:
# 				output = [base_rid] + sample(align_mismatch['_for'][False], 1)[0]
# 			else:
# 				output = [base_rid] + sample(align_mismatch['_rev'][False], 1)[0]

# 			out.write('\t'.join([str(e) for e in output]) + '\n')



# =============================================================================#
#                                    Main                                      #
# =============================================================================#

if __name__ == '__main__':
	output_dir, output_base_name, num_cpus, ref_b2index, ref_gtf, ref_introns, ref_repeatmasker = argv[1:]
	keep_failed_reads = True

	run_data = join(output_dir, f'{output_base_name}_run_data.tsv')

	# Load intron, splice site, and gene coordinates
	print(strftime('%m/%d/%y - %H:%M:%S | Parsing intron info...'))
	threep_sites, fivep_sites, introns = load_intron_info(ref_introns)
	# print(strftime('%m/%d/%y - %H:%M:%S | Parsing gene info...'))
	# gene_info = load_gene_info(ref_gtf)

	print(strftime('%m/%d/%y - %H:%M:%S | Parsing lariat reads...'))
	final_info_table = join(output_dir, f'{output_base_name}_threep_info_table.tsv')
	lariat_reads = load_lariat_table(final_info_table)

	# Parse counts of linearly aligned reads 
	print(strftime('%m/%d/%y - %H:%M:%S | Retrieving total mapped reads...'))
	sample_read_count = get_mapped_reads(output_dir, output_base_name)

	# Filter lariat reads
	print(strftime('%m/%d/%y - %H:%M:%S | Filtering lariat reads...'))
	filtered_lariats, failed_reads = filter_lariats(lariat_reads, threep_sites, fivep_sites, output_dir, num_cpus, ref_b2index, ref_repeatmasker)
	
	with open(run_data, 'a') as a:
		a.write(f'filtered_lariat_reads\t{len(filtered_lariats)}\n')

	# Now write it all to file
	print(strftime('%m/%d/%y - %H:%M:%S | Writing results to output file...'))
	with open(join(output_dir, f'{output_base_name}_lariat_reads.tsv'), 'w') as results_file:
		# make and write the header row
		header = '\t'.join(BASE_RESULTS_COLUMNS) + '\n'
		results_file.write(header)
			
		for read_info in filtered_lariats:
			read_output = '\t'.join([str(e) for e in read_info]) + f'\t{sample_read_count}'
			results_file.write(read_output + '\n')

	with open(join(output_dir, f'{output_base_name}_failed_lariat_reads.tsv'), 'w') as failed_file:
		header = '\t'.join(BASE_RESULTS_COLUMNS) + '\tfail_reason' + '\n'
		failed_file.write(header)

		for read_info in failed_reads:
			read_output = '\t'.join([str(e) for e in read_info]) + f'\t{sample_read_count}'
			failed_file.write(read_output + '\n')
























	# # Check if the 3' segment has a valid alignment upstream of the 5' segment
	# # Write trimmed sequences to fasta
	# trim_seqs = {rid: vals[0] for rid, vals in lariat_reads.items()}
	# seq_tmp_fa, seq_tmp_sam = f'{temp_dir}/trim_seq_tmp.fa', f'{temp_dir}/trim_seq_tmp.sam'
	# with open(seq_tmp_fa, 'w') as out_file:
	# 	for rid in trim_seqs:
	# 		out_file.write(f'>{rid}\n{trim_seqs[rid]}\n')
	#
	# # Map trimmed sequences to genome
	# map_call = f'bowtie2 --end-to-end --no-unal --threads {num_cpus} -k 30 -f -x {ref_b2index} -U {seq_tmp_fa} -S {seq_tmp_sam}'
	# run(map_call.split(' '), stdout=DEVNULL, stderr=DEVNULL)
	#
	# # Check if the trimmed sequence (which mapped to the 3'ss sequence) maps upstream of the 5'ss
	# # If it does, it's probably a intron-exon junction read from a pre-mRNA instead of a lariat
	# upstream_rids = set()
	# with open(seq_tmp_sam) as read_file:
	# 	for line in read_file:
	# 		if line[0] == '@':
	# 			continue
	#
	# 		# Load alignment info and the read's info
	# 		rid, _, trim_chrom, reference_start = line.strip().split('\t')[:4]
	# 		reference_start = int(reference_start)-1
	# 		_, _, chrom, strand, fivep_site, _, _, _, _, bp_site = lariat_reads[rid][:10]
	# 		gene_data = gene_info[chrom][strand].at(bp_site)
	# 		gene_names = [g.data['gene_name'] for g in gene_data]
	#
	# 		if trim_chrom != chrom:
	# 			continue
	#
	# 		# Check if where the trimmed sequence mapped is upstream of the 5'ss
	# 		trim_genes = [g.data['gene_name']
	# 						for g in gene_info[chrom][strand].at(reference_start)]
	# 		genes_overlap = sum(
	# 			1 for g in trim_genes if g in gene_names) > 0
	# 		if not genes_overlap:
	# 			continue
	#
	# 		# If it is, add it to upstream_rids for removal
	# 		if strand == '+' and reference_start <= fivep_site:
	# 			upstream_rids.add(rid)
	# 		elif strand == '-' and reference_start >= fivep_site:
	# 			upstream_rids.add(rid)
			# if rid in upstream_rids:
		# 	passed_filtering = False
		# 	fail_reason = 'upstream_found' if fail_reason is None else fail_reason

