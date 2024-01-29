from sys import argv
from os import rmdir, mkdir
from os.path import join, isfile, exists
from subprocess import run, DEVNULL
import gzip
from time import strftime

from intervaltree import Interval, IntervalTree
import pandas as pd


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
        intron_file = gzip.open(ref_introns, 'rt')
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


def load_lariat_table(out_dir: str) -> pd.DataFrame:
	'''
	For a given lariat-mapping of a fastq file, retrieve all the lariat reads from the XXX_lariat_info_table.tsv and put them in a dict, which
	can then be added to the experiment-wide superset dict
	'''
	lariat_reads = pd.read_csv(f'{out_dir}/lariat_info_table.tsv', sep='\t')

	return lariat_reads


def get_mapped_reads(out_dir:str) -> int:
	'''
	Get the number of reads that mapped to the reference genome from the *_run_data.tsv file 
	'''
	with open(f'{out_dir}/run_data.tsv', 'r') as file:
		line = file.readline()
		sample_read_count = int(line.split('\t')[1])

	return sample_read_count


def filt(row:pd.Series, fivep_sites:dict, threep_sites:dict, repeat_rids:set):
	# Check if BP is within 2bp of an annotated splice site
	if fivep_sites[row['chrom']][row['strand']].overlaps(row['bp_site']) or threep_sites[row['chrom']][row['strand']].overlaps(row['bp_site']):
		return 'near_ss'

	# Check if read mapped to a ubiquitin gene
	if row['gene_name'] in ('UBC', 'UBB'):
		return 'ubiquitin_gene'
	
	# Check if read mapped to a repetitive region
	if row['rid'] in repeat_rids:
		return 'in_repeat'

	return pd.NA


def filter_lariats(lariat_reads: pd.DataFrame, threep_sites: dict, fivep_sites: dict, output_dir:str, num_cpus:str, ref_b2index:str, ref_repeatmasker: str) -> dict:
	'''
	Filter the candidate lariat reads in order to exclude any that meet the following criteria:
			- BP is within 2bp of a splice site (likely an intron circle, not a lariat)
			- Read maps to UBB or UBC (likely false positive due to the repetitive nature of the genes)
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
		for i, row in lariat_reads.iterrows():
			fivep_out.write(f"{row['chrom']}\t{row['fivep_site']-1}\t{row['fivep_site']+1}\t{row['read_id']}\n")
			bp_out.write(f"{row['chrom']}\t{row['bp_site']-1}\t{row['bp_site']+1}\t{row['read_id']}\n")
				
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
	lariat_reads['filter_failed'] = lariat_reads.apply(filt, repeat_rids=repeat_rids, fivep_sites=fivep_sites, threep_sites=threep_sites, axis=1)
	
	choose 1 lariat per read
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

	failed_reads = lariat_reads[lariat_reads.filter_failed.notna()].copy()
	filtered_lariats = lariat_reads[lariat_reads.filter_failed.isna()].drop(columns='filter_failed')

	print(strftime('%m/%d/%y - %H:%M:%S') + f' | Pre-filter read count = {len(lariat_reads)}')
	print(strftime('%m/%d/%y - %H:%M:%S') + f' | Post-filter read count = {len(filtered_lariats)}')

	# Delete all temporary files
	# run(f'rm {seq_tmp_fa} {seq_tmp_sam} {fivep_tmp_bed} {bp_tmp_bed} {fivep_overlap_bed} {bp_overlap_bed}'.split(' '))
	# rmdir(temp_dir)

	return filtered_lariats, failed_reads





# =============================================================================#
#                                    Main                                      #
# =============================================================================#

if __name__ == '__main__':
	num_cpus, ref_b2index, ref_gtf, ref_introns, ref_repeatmasker, out_dir = argv[1:]
	keep_failed_reads = True

	# Load intron, splice site, and gene coordinates
	print(strftime('%m/%d/%y - %H:%M:%S | Parsing intron info...'))
	threep_sites, fivep_sites, introns = load_intron_info(ref_introns)
	# print(strftime('%m/%d/%y - %H:%M:%S | Parsing gene info...'))
	# gene_info = load_gene_info(ref_gtf)

	print(strftime('%m/%d/%y - %H:%M:%S | Parsing lariat reads...'))
	lariat_reads = load_lariat_table(out_dir)

	# Parse counts of linearly aligned reads 
	print(strftime('%m/%d/%y - %H:%M:%S | Retrieving total mapped reads...'))
	sample_read_count = get_mapped_reads(out_dir)

	# Filter lariat reads
	print(strftime('%m/%d/%y - %H:%M:%S | Filtering lariat reads...'))
	filtered_lariats, failed_reads = filter_lariats(lariat_reads, threep_sites, fivep_sites, out_dir, num_cpus, ref_b2index, ref_repeatmasker)

	with open(f'{out_dir}/run_data.tsv', 'a') as a:
		a.write(f'filtered_lariat_reads\t{len(filtered_lariats)}\n')

	# Now write it all to file
	print(strftime('%m/%d/%y - %H:%M:%S | Writing results to output file...'))
	filtered_lariats.to_csv(f'{out_dir}/lariat_reads.tsv', sep='\t', index=False)
	failed_reads.to_csv(f'{out_dir}/failed_lariat_reads.tsv', sep='\t', index=False)
	






















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

