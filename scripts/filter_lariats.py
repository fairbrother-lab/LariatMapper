
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
                        'total_mapped_reads')

RMAP = {'one': 'R1', 'two': 'R2'}


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


def parse_lariat_table(final_info_table: str) -> dict:
    '''
    For a given lariat-mapping of a fastq file, retrieve all the lariat reads from the XXX_final_info_table.txt and put them in a dict, which
    can then be added to the experiment-wide superset dict
    '''
    reads = {}
    with open(final_info_table, 'r') as lariat_file:
        next(lariat_file)  # Skip the header line
        for line in lariat_file:
            read_id, read_seq, chrom, strand, fivep_site, read_is_reverse, fivep_read_start, fivep_read_end, threep_site, bp_site, read_bp_nt, genomic_bp_nt, genomic_bp_window = line.strip().split()
            if read_is_reverse == 'False':
                trim_seq = read_seq[:int(fivep_read_start)]
            else:
                trim_seq = read_seq[int(fivep_read_end):]
            reads[read_id] = [trim_seq, read_seq, chrom, strand, int(fivep_site), read_is_reverse, fivep_read_start, fivep_read_end]
            reads[read_id] += [int(threep_site), int(bp_site), read_bp_nt, genomic_bp_nt, genomic_bp_window]
    return reads


def get_aligned_reads(output_dir, output_base_name) -> int:
    '''
    Get total aligned reads for each sample from the *_total_reads.txt file in its results directory
    '''
    sample_read_count = 0
    for read_num in ['one', 'two']:
        read_output_name = f'{output_base_name}_{RMAP[read_num]}'
        read_count_path = join(output_dir, read_output_name, f'{read_output_name}_total_reads.txt')
        if not isfile(read_count_path):
            raise FileNotFoundError(f'No total reads found for sample "{read_output_name}"')
        with open(read_count_path, 'r') as count_file:
            # line format: "total_reads=[total aligned reads]"
            aligned_reads = int(count_file.readline().split('=')[1])
            sample_read_count = max(sample_read_count, aligned_reads)

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


def filter_lariat_reads(merged_lariats: dict, threep_sites: dict, fivep_sites: dict, introns: dict, gene_info: dict, output_dir:str, num_cpus:str, ref_b2index:str, ref_repeatmasker: str) -> dict:
    '''
    Filter the candidate lariat reads in order to exclude any that meet the following criteria:
            - BP is within 2bp of a splice site (likely an intron circle, not a lariat)
            - 5'SS and 3'SS are not in the correct order
            - Read maps to a Ubiquitin gene (likely false positive due to repetitive nature of gene)
            - There is a valid aligment for the 3' segment upstream of the 5' segment
            - Both the 5'SS and the BP overlap with repetitive regions from RepeatMasker (likely false positive)
    '''

    # create temporary directory in the output dir
    temp_dir = output_dir + '/temp'
    if not exists(temp_dir):
        mkdir(temp_dir)

    filtered_reads = {}
    trim_reads = {}
    fail_reason = Counter()
    for rid in merged_lariats:
        trim_seq, read_seq, chrom, strand, fivep_site, read_is_reverse, fivep_read_start, fivep_read_end, threep_site, bp_site, read_bp_nt, genomic_bp_nt, genomic_bp_window = merged_lariats[rid]
        gene_data = gene_info[chrom][strand].at(bp_site)
        gene_names = [g.data['gene_name'] for g in gene_data]

        # Primary exclusion criteria:
        # Branchpoint is not within 2bp of an annotated splice site (filters intron circles)
        outside_ss = not fivep_sites[chrom][strand].overlaps(bp_site) and not threep_sites[chrom][strand].overlaps(bp_site)
        # The branchpoint is downstream of the 5' splice site
        correct_order = (strand == '+' and fivep_site < bp_site) or (strand == '-' and fivep_site > bp_site)
        # Read did not map to a ubiquitin gene
        outside_ubiquitin = 'UBB' not in gene_names and 'UBC' not in gene_names

        if outside_ss and correct_order and outside_ubiquitin:
            overlap_introns = list(introns[chrom][strand].overlap(bp_site, bp_site+1))
            if len(overlap_introns) > 0:
                if strand == '+':
                    threep_site = min(overlap_introns, key=lambda s: s.end-bp_site).end
                else:
                    threep_site = min(overlap_introns, key=lambda s: bp_site-s.begin).begin
                dist_to_threep = bp_site-threep_site if strand == '+' else threep_site-bp_site
                trim_reads[rid] = trim_seq
                filtered_reads[rid] = [gene_data, read_seq, chrom, strand, fivep_site,
                                               threep_site, bp_site, read_bp_nt, genomic_bp_nt, genomic_bp_window, dist_to_threep]
        elif not outside_ss:
            fail_reason['near_ss'] += 1
        elif not correct_order:
            fail_reason['fivep_bp_disordered'] += 1
        elif not outside_ubiquitin:
            fail_reason['ubiquitin_gene'] += 1

    # Filter reads were the 3' segment has a valid alignment upstream of the 5' segment
    seq_tmp_fa, seq_tmp_sam = f'{temp_dir}/trim_seq_tmp.fa', f'{temp_dir}/trim_seq_tmp.sam'
    with open(seq_tmp_fa, 'w') as out_file:
        for rid in trim_reads:
            out_file.write('>{}\n{}\n'.format(rid, trim_reads[rid]))

    map_call = f'bowtie2 --end-to-end --no-unal --threads {num_cpus} -k 30 -f -x {ref_b2index} -U {seq_tmp_fa} -S {seq_tmp_sam}'
    run(map_call.split(' '), stdout=DEVNULL, stderr=DEVNULL)

    upstream_rids = set()
    with open(seq_tmp_sam) as read_file:
        for line in read_file:
            if line[0] != '@':
                rid, _, trim_chrom, reference_start = line.strip().split('\t')[
                    :4]
                gene_data, _, chrom, strand, fivep, threep = filtered_reads[rid][:6]
                gene_names = [g.data['gene_name'] for g in gene_data]
                if trim_chrom == chrom:
                    reference_start = int(reference_start)-1
                    trim_genes = [g.data['gene_name']
                                  for g in gene_info[chrom][strand].at(reference_start)]
                    genes_overlap = sum(
                        1 for g in trim_genes if g in gene_names) > 0
                    if genes_overlap:
                        if strand == '+' and reference_start <= fivep:
                            upstream_rids.add(rid)
                        elif strand == '-' and reference_start >= fivep:
                            upstream_rids.add(rid)
    fail_reason['upstream_found'] += len(upstream_rids)

    # Filter reads where both the 5'SS and the BP overlap with repetitive regions
    fivep_tmp_bed, bp_tmp_bed = f'{temp_dir}/fivep_tmp.bed', f'{temp_dir}/bp_tmp.bed'
    with open(fivep_tmp_bed, 'w') as fivep_out:
        with open(bp_tmp_bed, 'w') as bp_out:
            for rid in filtered_reads:
                if rid not in upstream_rids:
                    _, _, chrom, _, fivep_site, _, bp_site = filtered_reads[rid][:7]
                    fivep_out.write('{}\t{}\t{}\t{}\n'.format(
                        chrom, fivep_site-1, fivep_site+1, rid))
                    bp_out.write('{}\t{}\t{}\t{}\n'.format(
                        chrom, bp_site-1, bp_site+1, rid))

    fivep_overlap_bed, bp_overlap_bed = f'{temp_dir}/fivep_repeat_overlaps.bed', f'{temp_dir}/bp_repeat_overlaps.bed'
    with open(fivep_overlap_bed, 'w') as out_file:
        run(f'bedtools intersect -u -a {fivep_tmp_bed} -b {ref_repeatmasker}'.split(' '),
            stdout=out_file)
    with open(bp_overlap_bed, 'w') as out_file:
        run(f'bedtools intersect -u -a {bp_tmp_bed} -b {ref_repeatmasker}'.split(' '),
            stdout=out_file)

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
    fail_reason['in_repeat'] += len(repeat_rids)

    discard_rids = upstream_rids.union(repeat_rids)

    # Discard filtered reads and output final set of lariat reads
    filtered_reads = {rid:filtered_reads[rid] for rid in filtered_reads if rid not in discard_rids}
    for rid in filtered_reads:
        gene_data = filtered_reads[rid][0].pop()
        gene_data = [gene_data.data['gene_name'], gene_data.data['ensembl_id'], gene_data.data['gene_type'], rid]
        filtered_reads[rid] = gene_data + filtered_reads[rid][1:]

    fail_reason_out = ', '.join([f'{r} ({fail_reason[r]})' for r in fail_reason])
    print(strftime('%m/%d/%y - %H:%M:%S') + f' | Pre-filter read count = {len(merged_lariats)}')
    print(strftime('%m/%d/%y - %H:%M:%S') + f' | Post-filter read count = {len(filtered_reads)}')
    print(strftime('%m/%d/%y - %H:%M:%S') + f' | Discarded read counts = {fail_reason_out}')

    # Delete all temporary files
    run(f'rm {seq_tmp_fa} {seq_tmp_sam} {fivep_tmp_bed} {bp_tmp_bed} {fivep_overlap_bed} {bp_overlap_bed}'.split(' '))
    rmdir(temp_dir)

    return filtered_reads


# =============================================================================#
#                                    Main                                     #
# =============================================================================#
if __name__ == '__main__':

    fastq_dir, read_one_file, read_two_file, output_dir, output_base_name, num_cpus, ref_b2index, \
    ref_fasta, ref_gtf, ref_5p_fasta, ref_3p_b2index, ref_3p_lengths, ref_introns, ref_repeatmasker = argv[1:]

    # Load intron, splice site, and gene coordinates
    print(strftime('%m/%d/%y - %H:%M:%S | Parsing intron info...'))
    threep_sites, fivep_sites, introns = get_intron_info(ref_introns)
    print(strftime('%m/%d/%y - %H:%M:%S | Parsing gene info...'))
    gene_info = get_gene_info(ref_gtf)

    print(strftime('%m/%d/%y - %H:%M:%S | Parsing lariat reads...'))
    lariat_reads = {}
    for read_num in ['one', 'two']:
        read_output_name = f'{output_base_name}_{RMAP[read_num]}'
        final_info_table = join(output_dir, read_output_name, f'{read_output_name}_final_info_table.txt')
        lariat_reads[read_output_name] = parse_lariat_table(final_info_table)

    # Parse counts of linearly aligned reads for each sample
    print(strftime('%m/%d/%y - %H:%M:%S | Parsing total mapped reads...'))
    sample_read_count = get_aligned_reads(output_dir, output_base_name)

    # Merge results from read one (R1) and read two (R2) files
    print(strftime('%m/%d/%y - %H:%M:%S | Merging R1 and R2 lariat reads...'))
    merged_lariats = merge_opposite_reads(lariat_reads, output_base_name)

    # Filter merged lariat reads
    print(strftime('%m/%d/%y - %H:%M:%S | Filtering lariat reads...'))
    merged_filtered_lariats = filter_lariat_reads(merged_lariats, threep_sites, fivep_sites, introns, gene_info, output_dir, num_cpus, ref_b2index, ref_repeatmasker)

    # Now write it all to file
    print(strftime('%m/%d/%y - %H:%M:%S | Writing results to output file...'))
    with open(join(output_dir, f'{output_base_name}_lariat_reads.txt'), 'w') as results_file:
        # make and write the header row
        header = '\t'.join(BASE_RESULTS_COLUMNS) + '\n'
        results_file.write(header)
        for read_info in merged_filtered_lariats.values():
            read_output = read_info + [sample_read_count]
            results_file.write('\t'.join([str(e) for e in read_output]) + '\n')
