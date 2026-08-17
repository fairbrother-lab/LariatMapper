
import os
import sys
import pandas
import subprocess
from collections import defaultdict, Counter
from intervaltree import Interval, IntervalTree

import utils

# Input files
OUTPUT_SJ_FILE = "{}SJ.out.tab"
PUTATITVE_LARIATS_FILE = "{}putative_lariats.tsv"
# Output file
FILTERED_LARIATS_FILE = "{}reverse_mapped_lariats_filtered.tsv"

def parse_genes(ref_genes:str, log) -> dict[str, dict[str, IntervalTree]]:
	genes_df, genes = pandas.read_csv(ref_genes, sep='\t'), {}
	for _, gene_info in genes_df.iterrows():
		gene_id, chrom, strand, start, end = [gene_info[c] for c in ['gene_id', 'chrom', 'strand', 'start', 'end']]
		if chrom not in genes:
			genes[chrom] = {'+': IntervalTree(), '-': IntervalTree()}
		genes[chrom][strand].add(Interval(start, end, {'gene_id':gene_id}))
	
	return genes

def parse_introns(ref_introns:str, log) -> dict[str, dict[str, IntervalTree]]:
											
	introns_df = pandas.read_csv(ref_introns, sep='\t', na_filter=False)
	introns_df['gene_id'] = introns_df.gene_id.transform(lambda gene_id: set(gene_id.split(',')))

	chromosomes = set(introns_df.chrom.unique())
	introns = {chrom:{'+': IntervalTree(), '-': IntervalTree()} for chrom in chromosomes}
	for _, intron_info in introns_df.iterrows():
		chrom, strand, start, end = [intron_info[c] for c in ['chrom', 'strand', 'start', 'end']]
		intron_id = f'{chrom}_{strand}_{start}_{end}'
		introns[chrom][strand].add(Interval(start, end, {'gene_id': intron_info['gene_id'], 'intron_id':intron_id}))
	
	return introns

def parse_splice_junction_count(sj_file:str, log) -> int:
	'''
	Get the total count of splice junctions from the SJ.out.tab file output by STAR
	'''
	total_sj_count = 0
	with open(sj_file) as in_file:
		for line in in_file:
			uniq_read_count = line.strip().split('\t')[6]
			total_sj_count += int(uniq_read_count)

	return total_sj_count

def parse_line(sam_line:str):
	alignment_info = sam_line.split('\t')
	read_id, flag, chrom, start_pos, _, cigar_str = alignment_info[:6]
	tags = [tag_info.split(':') for tag_info in alignment_info[11:]]
	tags = {tag_info[0]:tag_info[2] for tag_info in tags}
	alignment_score, mismatch_count = int(tags['AS']), int(tags['NM'])

	cigar_tuples = utils.get_cigar_tuples(cigar_str)
	start_pos = int(start_pos)-1
	end_pos = utils.get_align_end(start_pos, cigar_tuples)
	gaps = [length for length, op in cigar_tuples if op in 'ID']
	query_length = sum([length for length, op in cigar_tuples if op in 'MIS=X'])
	mismatch_p = mismatch_count/float(query_length)
	read_orient_to_chrom = utils.align_orient(flag)

	return read_id, chrom, start_pos, end_pos, read_orient_to_chrom, cigar_str, mismatch_count, mismatch_p, gaps

def verify_bp_alignment(lariat_df:pandas.DataFrame, output_base:str, hisat2_index:str, genes:dict, threads:str, log):

	tmp_fasta = f"{output_base}head_full_seq.tmp.fa"
	read_info_dict = {}
	with open(tmp_fasta, 'w') as out_file:
		for _, read_info in lariat_df.iterrows():
			bp_id = f"{read_info['chrom']};{read_info['bp_pos']};{read_info['strand']};{read_info['bp_align_pos']}"
			tmp_read_id = f"{read_info['read_id']}|{bp_id}"
			tmp_read_seq = read_info['read_seq_forward'][:(read_info['read_bp_pos']+1)]
			out_file.write(f'>{tmp_read_id}\n{tmp_read_seq}\n')
			read_info_dict[tmp_read_id] = read_info
		
	hisat2_cmd = f"hisat2 --no-softclip --no-spliced-alignment --no-unal --sensitive -k 5 --threads {threads} -f -x {hisat2_index} -U {tmp_fasta}"
	hisat2_output = subprocess.run(hisat2_cmd, shell=True, capture_output=True, text=True).stdout.strip().split('\n')
	passed_alignments = []
	read_failures = defaultdict(list)
	for sam_line in hisat2_output:
		if sam_line[0] != '@':
			tmp_read_id, chrom, start_pos, end_pos, read_orient_to_chrom, cigar_str, mismatch_count, mismatch_p, gaps = parse_line(sam_line)
			read_id, bp_id = tmp_read_id.split('|')
			read_mate_id = read_id[:-4]
			bp_chrom, bp_pos, bp_strand, bp_align_pos = bp_id.split(';')

			if mismatch_p > 0.1:
				read_failures[read_mate_id].append('mismatch_rate')
				continue

			if len(gaps) > 1:
				read_failures[read_mate_id].append('multiple_gaps')
				continue

			if sum(gaps) > 3:
				read_failures[read_mate_id].append('long_gap')
				continue

			if bp_chrom != chrom or (bp_strand == '+' and read_orient_to_chrom == 'Reverse') or (bp_strand == '-' and read_orient_to_chrom == 'Forward'):
				read_failures[read_mate_id].append('chrom_strand_mismatch')
				continue
			
			align_genes = set([g.data['gene_id'] for g in genes[chrom][bp_strand].overlap(start_pos, end_pos)])
			orig_genes = set(read_info_dict[tmp_read_id]['gene_id'].split(','))
			overlap_genes = align_genes.intersection(orig_genes)
			if len(overlap_genes) == 0:
				read_failures[read_mate_id].append('gene_mismatch')
				continue

			full_align_bp_pos = end_pos-1 if bp_strand == '+' else start_pos
			if full_align_bp_pos != int(bp_align_pos):
				read_failures[read_mate_id].append('bp_align_pos_mismatch')
				continue

			gaps = '0' if len(gaps)==0 else str(gaps[0])

			passed_alignments.append([read_id, chrom, bp_strand, int(bp_pos), int(bp_align_pos), start_pos, end_pos, cigar_str, mismatch_count, gaps])
	
	failure_score = {'mismatch_rate':0, 'multiple_gaps':1, 'long_gap':2,
				     'chrom_strand_mismatch':3, 'gene_mismatch':4, 'bp_align_pos_mismatch':5}
	read_failure_counts = Counter([max(read_failures[rid], key=lambda f:failure_score[f]) for rid in read_failures])
	log.info(f"Putative lariat branchpoint alignment filter failures: {read_failure_counts}")
	
	passed_alignments = pandas.DataFrame(passed_alignments, columns=['read_id', 'chrom', 'strand', 'bp_pos', 'bp_align_pos', 'head_align_start_pos', 'head_align_end_pos',
																     'head_cigar_str', 'head_mismatches', 'head_gaps'])
	
	lariat_pass_df = pandas.merge(lariat_df, passed_alignments, on=['read_id', 'chrom', 'strand', 'bp_pos', 'bp_align_pos'], how='right')

	os.remove(tmp_fasta)
	
	return lariat_pass_df

def filter_nested_reads(lariat_df:pandas.DataFrame, log):

	align_pass = []
	for _, read_info in lariat_df.iterrows():
		head_align_start_pos, head_align_end_pos, tail_align_start_pos, tail_align_end_pos = [read_info[c] for c in ['head_align_start_pos', 'head_align_end_pos', 'tail_align_start_pos', 'tail_align_end_pos']]

		if head_align_start_pos >= tail_align_start_pos and head_align_end_pos <= tail_align_end_pos:
			align_pass.append(False)
		elif tail_align_start_pos >= head_align_start_pos and tail_align_end_pos <= head_align_end_pos:
			align_pass.append(False)
		else:
			align_pass.append(True)

	log.info(f'Putative lariat nested read filter failures: {len(align_pass)-sum(align_pass)}')
	
	lariat_df = lariat_df[align_pass]

	return lariat_df

def filter_multi_alignments(lariat_df:pandas.DataFrame, introns:dict, log):
	
	align_pass, align_tail_introns, align_head_introns = [], [], []
	for _, read_info in lariat_df.iterrows():
		chrom, strand, head_align_start_pos, head_align_end_pos, tail_align_start_pos, tail_align_end_pos = [read_info[c] for c in ['chrom', 'strand', 'head_align_start_pos', 'head_align_end_pos', 'tail_align_start_pos', 'tail_align_end_pos']]
		# Confirm tail and head alignments exist entirely within introns
		tail_introns = [i for i in introns[chrom][strand].overlap(tail_align_start_pos, tail_align_end_pos-1) if i.begin<=tail_align_start_pos and i.end>=tail_align_end_pos]
		head_introns = [i for i in introns[chrom][strand].overlap(head_align_start_pos, head_align_end_pos-1) if i.begin<=head_align_start_pos and i.end>=head_align_end_pos]
		pass_align = False
		if len(tail_introns) > 0 and len(head_introns) > 0:
			# Limit to introns within the same genes
			tail_genes, head_genes = set.union(*[i.data['gene_id'] for i in tail_introns]), set.union(*[i.data['gene_id'] for i in head_introns])
			overlap_genes = tail_genes.intersection(head_genes)
			tail_introns = [i for i in tail_introns if len(i.data['gene_id'].intersection(overlap_genes))>0]
			head_introns = [i for i in head_introns if len(i.data['gene_id'].intersection(overlap_genes))>0]
			if len(tail_introns) > 0 and len(head_introns) > 0:
				pass_align = True
				
		align_pass.append(pass_align)
		if pass_align:
			align_tail_introns.append(';'.join([i.data['intron_id'] for i in tail_introns]))
			align_head_introns.append(';'.join([i.data['intron_id'] for i in head_introns]))
		else:
			align_tail_introns.append('NA')
			align_head_introns.append('NA')
	
	lariat_df['tail_introns'] = align_tail_introns
	lariat_df['head_introns'] = align_head_introns
	lariat_df = lariat_df[align_pass]

	lariat_df.read_id = lariat_df.read_id.str.slice(0,-4)
	lariat_df[['read_id', 'read_num']] = lariat_df.read_id.str.split('/', expand=True)
	lariat_df.read_num = lariat_df.read_num.astype(int)
	lariat_df.insert(1, 'read_num', lariat_df.pop('read_num'))

	align_counts = Counter(lariat_df.groupby('read_id').size())

	log.info(f"Putative lariat intron filter failures: {align_counts}")
	log.info(f"Remaining alignment counts after intron filter: {align_counts}")

	lariat_df = lariat_df.groupby('read_id').sample(n=1, random_state=1)

	return lariat_df


def adjust_bp_fivep(lariat_df:pandas.DataFrame, log):

	def move_values(ldf:pandas.DataFrame, index, s, dist):
		if s=='+':
			ldf.at[index, 'fivep_pos'] -= dist
			ldf.at[index, 'bp_align_pos'] -= dist
			ldf.at[index, 'fivep_dist_to_anno'] -= dist
		else:
			ldf.at[index, 'fivep_pos'] += dist
			ldf.at[index, 'bp_align_pos'] += dist
			ldf.at[index, 'fivep_dist_to_anno'] += dist
		ldf.at[index, 'read_bp_pos'] -= dist
		ldf.at[index, 'read_fivep_pos'] -= dist
		return ldf
		
	for i, lar_info in lariat_df.iterrows():

		strand, fivep_pos, fivep_nearest_anno_pos = [lar_info[c] for c in ['strand', 'fivep_pos', 'fivep_nearest_anno_pos']]
		
		fivep_rel_pos = fivep_pos-fivep_nearest_anno_pos if strand=='+' else fivep_nearest_anno_pos-fivep_pos
		if fivep_rel_pos <= 3:
			lariat_df = move_values(lariat_df, i, strand, abs(fivep_rel_pos))
			continue
	
	return lariat_df
			

def overlap_snaptron(lariat_pass_df:pandas.DataFrame, snaptron_file:str, log):

	snaptron_columns = defaultdict(list)
	juncs_to_process = sum((lariat_pass_df['fivep_pos']!=lariat_pass_df['fivep_nearest_anno_pos']))
	log.info(f'Pulling Snaptron junction info for {juncs_to_process} alignments')
	aligns_done = 0
	for _, read_info in lariat_pass_df.iterrows():
		chrom, strand, fivep_pos, fivep_nearest_anno_pos = [read_info[c] for c in ['chrom', 'strand', 'fivep_pos', 'fivep_nearest_anno_pos']]
		if fivep_pos != fivep_nearest_anno_pos:
			tabix_cmd = f'tabix {snaptron_file} {chrom}:{fivep_pos-1}-{fivep_pos+1}'
			junction_lines = subprocess.run(tabix_cmd, shell=True, capture_output=True, text=True).stdout.strip().split('\n')
			found_fivep, fivep_cov_sum, recursive_threep = False, 0, False
			fivep_junc_threep, recursive_junc_fivep = [], []
			for jx_line in junction_lines:
				if jx_line != '':
					_, _, jx_start, jx_end, _, jx_strand, _, _, _, _, _, samples, _, coverage_sum, _, _, _ = jx_line.split('\t')
					if jx_strand == strand:
						if strand == '+':
							fivep, threep = int(jx_start)-1, int(jx_end)-1
						else:
							fivep, threep = int(jx_end)-1, int(jx_start)-1
						if fivep == fivep_pos:
							found_fivep = True
							fivep_cov_sum += int(coverage_sum)
							fivep_junc_threep.append(f'{threep}:{coverage_sum}')
						if threep == fivep_pos-1:
							recursive_threep = True
							recursive_junc_fivep.append(f'{fivep}:{coverage_sum}')

			fivep_junc_threep = ';'.join(fivep_junc_threep) if len(fivep_junc_threep)>0 else None
			recursive_junc_fivep = ';'.join(recursive_junc_fivep) if len(recursive_junc_fivep)>0 else None
			snaptron_columns['snaptron_fivep'].append(found_fivep)
			snaptron_columns['snaptron_fivep_cov'].append(fivep_cov_sum)
			snaptron_columns['snaptron_recursive_threep'].append(recursive_threep)
			snaptron_columns['snaptron_downstream_threep'].append(fivep_junc_threep)
			snaptron_columns['snaptron_upstream_recursive_fivep'].append(recursive_junc_fivep)
			
			aligns_done += 1
			if aligns_done % 25 == 0:
				log.info(f'Done pulling junctions for {aligns_done} alignments')
		else:
			snaptron_columns['snaptron_fivep'].append(None)
			snaptron_columns['snaptron_fivep_cov'].append(None)
			snaptron_columns['snaptron_recursive_threep'].append(None)
			snaptron_columns['snaptron_downstream_threep'].append(None)
			snaptron_columns['snaptron_upstream_recursive_fivep'].append(None)
	
	lariat_pass_df.loc[:, 'snaptron_fivep'] = snaptron_columns['snaptron_fivep']
	lariat_pass_df.loc[:, 'snaptron_recursive_threep'] = snaptron_columns['snaptron_recursive_threep']
	lariat_pass_df.loc[:, 'snaptron_downstream_threep'] = snaptron_columns['snaptron_downstream_threep']
	lariat_pass_df.loc[:, 'snaptron_upstream_recursive_fivep'] = snaptron_columns['snaptron_upstream_recursive_fivep']
	
	return lariat_pass_df



if __name__ == '__main__':

	threads, output_base, hisat2_index, ref_genes, ref_introns, log_level = sys.argv[1:]

	log = utils.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}')

	lariat_df = pandas.read_csv(PUTATITVE_LARIATS_FILE.format(output_base), sep='\t')

	read_mate_ids = lariat_df['read_id'].str.slice(stop=-4)
	read_mate_id_count = read_mate_ids.unique().size
	read_id_count = read_mate_ids.str.slice(stop=-2).unique().size
	log.info(f'Processing putative lariat alignments for {read_id_count} reads ({read_mate_id_count} mates)')

	genes = parse_genes(ref_genes, log)
	introns = parse_introns(ref_introns, log)

	lariat_pass_df = adjust_bp_fivep(lariat_df, log)

	lariat_pass_df = verify_bp_alignment(lariat_df, output_base, hisat2_index, genes, threads, log)

	lariat_pass_df = filter_nested_reads(lariat_pass_df, log)

	lariat_pass_df = filter_multi_alignments(lariat_pass_df, introns, log)

	# Get splice junction count
	lariat_df['sj_count_total'] = parse_splice_junction_count(OUTPUT_SJ_FILE.format(output_base), log)

	snaptron_path = '/users/jbuerer/scratch/genomes/snaptron/srav3_hg38_junctions.bgz'
	lariat_pass_df = overlap_snaptron(lariat_pass_df, snaptron_path, log)

	lariat_pass_df.to_csv(FILTERED_LARIATS_FILE.format(output_base), sep='\t', index=False)








