from sys import argv
from os import rmdir, mkdir
from os.path import join, isfile, exists
from intervaltree import Interval, IntervalTree
from collections import Counter
from subprocess import run, DEVNULL
from gzip import open as gzipopen
from time import strftime

# =============================================================================#
#                                  Constants                                  #
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
                        'total_mapped_reads', 
						'passed_filtering',
						'fail_reason')



# =============================================================================#
#                                  Functions                                  #
# =============================================================================#
def get_intron_info(ref_introns) -> tuple:
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


def get_gene_info(ref_gtf) -> dict:
    '''
    Get whole-genome gene info from annotation file
    '''
    if ref_gtf[-2:] == 'gz':
        gtf_file = gzipopen(ref_gtf, 'rt')
    else:
        gtf_file = open(ref_gtf)

    gene_info = {}
    for line in gtf_file:
        if line[0] != '#':
            chrom, _, feat, start, end, _, strand, _, annotations = line.strip().split('\t')
            if feat == 'gene':
                annotations = {a.split(' ')[0]: a.split(' ')[1].replace(
                    '\"', '') for a in annotations[:-1].split('; ')}
                if chrom not in gene_info:
                    gene_info[chrom] = {s: IntervalTree() for s in ('-', '+')}
                gene_int_data = {
                    'gene_name': annotations['gene_name'], 'gene_type': annotations['gene_type'], 'ensembl_id': annotations['gene_id']}
                gene_info[chrom][strand].add(
                    Interval(int(start)-1, int(end), gene_int_data))
    gtf_file.close()

    return gene_info


def parse_lariat_table(threep_info_table: str) -> dict:
	'''
	For a given lariat-mapping of a fastq file, retrieve all the lariat reads from the XXX_threep_info_table.tsv and put them in a dict, which
	can then be added to the experiment-wide superset dict
	'''
	reads = {}
	with open(threep_info_table, 'r') as lariat_file:
		next(lariat_file)  # Skip the header line
		for line in lariat_file:
			read_id, read_seq, chrom, strand, fivep_site, read_is_reverse, fivep_read_start, fivep_read_end, threep_site, bp_site, read_bp_nt, genomic_bp_nt, genomic_bp_window = line.strip().split('\t')			

			trim_seq = read_seq[int(fivep_read_end):] if read_is_reverse == 'True' else read_seq[:int(fivep_read_start)]
			fivep_site = int(fivep_site)
			threep_site = int(threep_site)
			bp_site = int(bp_site)

			reads[read_id] = [trim_seq, read_seq, chrom, strand, fivep_site, read_is_reverse, fivep_read_start, fivep_read_end]
			reads[read_id] += [threep_site, bp_site, read_bp_nt, genomic_bp_nt, genomic_bp_window]
	return reads


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


def merge_opposite_reads(lariat_reads: dict, output_base_name) -> dict:
    '''
    Merge the lariats reads of each pair of mates into a single list of unique lariats reads from the given sample
    '''
    merged_reads = {}
    read_one_output_name = f'{output_base_name}_R1'
    read_two_output_name = f'{output_base_name}_R2'
    read_one_rids, read_two_rids, repeat_rids = set(), set(), set()
    for rid in lariat_reads[read_one_output_name]:
        merged_reads[rid] = lariat_reads[read_one_output_name][rid]
        read_one_rids.add(rid)
    for rid in lariat_reads[read_two_output_name]:
        if rid not in read_one_rids:
            merged_reads[rid] = lariat_reads[read_two_output_name][rid]
            read_two_rids.add(rid)
        else:
            repeat_rids.add(rid)

    return merged_reads



def filter_lariat_reads(lariat_reads: dict, threep_sites: dict, fivep_sites: dict, introns: dict, gene_info: dict, output_dir:str, num_cpus:str, ref_b2index:str, ref_repeatmasker: str) -> dict:
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

	trim_seqs = {}
	filtered_reads = {}
	for rid in lariat_reads:
		trim_seq, read_seq, chrom, strand, fivep_site, read_is_reverse, fivep_read_start, fivep_read_end, threep_site, bp_site, read_bp_nt, genomic_bp_nt, genomic_bp_window = lariat_reads[rid]

		gene_data = gene_info[chrom][strand].at(bp_site)
		gene_names = [g.data['gene_name'] for g in gene_data]

		passed_filtering = True
		fail_reason = None
		# Check if BP is within 2bp of an annotated splice site
		if fivep_sites[chrom][strand].overlaps(bp_site) or threep_sites[chrom][strand].overlaps(bp_site):
			passed_filtering = False
			fail_reason = 'near_ss' if fail_reason is None else fail_reason
		# Check if read mapped to a ubiquitin gene
		if 'UBB' in gene_names or 'UBC' in gene_names:
			passed_filtering = False
			fail_reason = 'ubiquitin_gene' if fail_reason is None else fail_reason
		
		# ???
		overlap_introns = list(introns[chrom][strand].overlap(bp_site, bp_site+1))
		# if len(overlap_introns) > 0:
		assert overlap_introns != [], f'{lariat_reads[rid]}'
		if strand == '+':
			threep_site = min(overlap_introns, key=lambda s: s.end-bp_site).end
		else:
			threep_site = min(overlap_introns, key=lambda s: bp_site-s.begin).begin
		dist_to_threep = bp_site-threep_site if strand == '+' else threep_site-bp_site

		trim_seqs[rid] = trim_seq
		filtered_reads[rid] = [gene_data, read_seq, chrom, strand, fivep_site,
										threep_site, bp_site, read_bp_nt, genomic_bp_nt, genomic_bp_window, dist_to_threep, passed_filtering, fail_reason]

	# Check if the 3' segment has a valid alignment upstream of the 5' segment
	# Write trimmed sequences to fasta
	seq_tmp_fa, seq_tmp_sam = f'{temp_dir}/trim_seq_tmp.fa', f'{temp_dir}/trim_seq_tmp.sam'
	with open(seq_tmp_fa, 'w') as out_file:
		for rid in trim_seqs:
			out_file.write(f'>{rid}\n{trim_seqs[rid]}\n')

	# Map trimmed sequences to genome
	map_call = f'bowtie2 --end-to-end --no-unal --threads {num_cpus} -k 30 -f -x {ref_b2index} -U {seq_tmp_fa} -S {seq_tmp_sam}'
	run(map_call.split(' '), stdout=DEVNULL, stderr=DEVNULL)

	# Check if the trimmed sequence (which mapped to the 3'ss sequence) maps upstream of the 5'ss
	# If it does, it's probably a intron-exon junction read from a pre-mRNA instead of a lariat
	upstream_rids = set()
	with open(seq_tmp_sam) as read_file:
		for line in read_file:
			if line[0] == '@':
				continue

			# Load alignment info and the read's info
			rid, _, trim_chrom, reference_start = line.strip().split('\t')[:4]
			reference_start = int(reference_start)-1
			gene_data, _, chrom, strand, fivep, threep = filtered_reads[rid][:6]
			gene_names = [g.data['gene_name'] for g in gene_data]

			if trim_chrom != chrom:
				continue

			# Check if where the trimmed sequence mapped is upstream of the 5'ss
			trim_genes = [g.data['gene_name']
							for g in gene_info[chrom][strand].at(reference_start)]
			genes_overlap = sum(
				1 for g in trim_genes if g in gene_names) > 0
			if not genes_overlap:
				continue

			# If it is, add it to upstream_rids for removal
			if strand == '+' and reference_start <= fivep:
				upstream_rids.add(rid)
			elif strand == '-' and reference_start >= fivep:
				upstream_rids.add(rid)

	# Mark reads as failed
	for rid in upstream_rids:
		filtered_reads[rid][-2] = False
		filtered_reads[rid][-1] = 'upstream_found'

	# Check if both the 5'SS and the BP overlap with a repetitive region
	# Write the 5'ss and BP coordinates to BED files
	fivep_tmp_bed, bp_tmp_bed = f'{temp_dir}/fivep_tmp.bed', f'{temp_dir}/bp_tmp.bed'
	with open(fivep_tmp_bed, 'w') as fivep_out, open(bp_tmp_bed, 'w') as bp_out:
		for rid in filtered_reads:
			if filtered_reads[rid][-2] is True:
				_, _, chrom, _, fivep_site, _, bp_site = filtered_reads[rid][:7]
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

	# Add reads where both sites overlapped to repeat_rids for remoal
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

	# Mark reads as failed
	for rid in repeat_rids:
		filtered_reads[rid][-2] = False
		filtered_reads[rid][-1] = 'in_repeat'

	# Add gene info to reads
	for rid in filtered_reads:
		gene_data = filtered_reads[rid][0].pop()
		gene_data = [gene_data.data['gene_name'], gene_data.data['ensembl_id'], gene_data.data['gene_type'], rid]
		filtered_reads[rid] = gene_data + filtered_reads[rid][1:]

	passed_reads = len([val for val in filtered_reads.values() if val[-2] is True])
	print(strftime('%m/%d/%y - %H:%M:%S') + f' | Pre-filter read count = {len(lariat_reads)}')
	print(strftime('%m/%d/%y - %H:%M:%S') + f' | Post-filter read count = {passed_reads}')

	# Delete all temporary files
	# run(f'rm {seq_tmp_fa} {seq_tmp_sam} {fivep_tmp_bed} {bp_tmp_bed} {fivep_overlap_bed} {bp_overlap_bed}'.split(' '))
	# rmdir(temp_dir)

	return filtered_reads



# =============================================================================#
#                                    Main                                     #
# =============================================================================#

if __name__ == '__main__':
	output_dir, output_base_name, num_cpus, ref_b2index, ref_gtf, ref_introns, ref_repeatmasker = argv[1:]
	keep_failed_reads = True

	run_data = join(output_dir, f'{output_base_name}_run_data.tsv')

	# Load intron, splice site, and gene coordinates
	print(strftime('%m/%d/%y - %H:%M:%S | Parsing intron info...'))
	threep_sites, fivep_sites, introns = get_intron_info(ref_introns)
	print(strftime('%m/%d/%y - %H:%M:%S | Parsing gene info...'))
	gene_info = get_gene_info(ref_gtf)

	print(strftime('%m/%d/%y - %H:%M:%S | Parsing lariat reads...'))
	lariat_reads = {}
	threep_info_table = join(output_dir, f'{output_base_name}_threep_info_table.tsv')
	lariat_reads = parse_lariat_table(threep_info_table)

	# Parse counts of linearly aligned reads 
	print(strftime('%m/%d/%y - %H:%M:%S | Retrieving total mapped reads...'))
	sample_read_count = get_mapped_reads(output_dir, output_base_name)

	# Filter lariat reads
	print(strftime('%m/%d/%y - %H:%M:%S | Filtering lariat reads...'))
	filtered_lariats = filter_lariat_reads(lariat_reads, threep_sites, fivep_sites, introns, gene_info, output_dir, num_cpus, ref_b2index, ref_repeatmasker)
	
	passed_all = sum(1 for val in filtered_lariats.values() if val[-2] is True)
	with open(run_data, 'a') as a:
		a.write(f'passed_all_filters_reads\t{passed_all}\n')
	if not keep_failed_reads:
		filtered_lariats = {rid: values for rid, values in filtered_lariats.items() if values[-2] is True}

	# Now write it all to file
	print(strftime('%m/%d/%y - %H:%M:%S | Writing results to output file...'))
	with open(join(output_dir, f'{output_base_name}_lariat_reads.tsv'), 'w') as results_file:
		# make and write the header row
		header = '\t'.join(BASE_RESULTS_COLUMNS) + '\n'
		results_file.write(header)
			
		for read_info in filtered_lariats.values():
			read_output = read_info[:-2] + [sample_read_count] + read_info[-2:]
			results_file.write('\t'.join([str(e) for e in read_output]) + '\n')