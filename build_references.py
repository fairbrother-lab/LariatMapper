#!/usr/bin/env python3

import os
import subprocess
import time
import gzip
import argparse

#=============================================================================#
#                                  Functions                                  #
#=============================================================================#
def parse_attributes(attribute_string:str, file_type:str) -> dict:
	if file_type == 'gtf':
		attributes = attribute_string.rstrip('";').split('; ')
		attributes = [attr.split(' ') for attr in attributes]
		tags = [attr_val.strip('"') for attr_name, attr_val in attributes if attr_name=='tag']
		attributes = {attr_name: attr_val.strip('"') for attr_name, attr_val in attributes if attr_name!='tag'}
		attributes['tags'] = tags
	else:
		attributes = [attr.split('=') for attr in attribute_string.split(';')]
		attributes = [(attr[0].lstrip(), attr[1]) for attr in attributes]
		attributes = dict(attributes)
		if 'tag' in attributes:
			attributes['tags'] = attributes['tag'].split(',')

	return attributes

def parse_transcript_info(ref_anno:str, anno_type:str, gunzip:bool) -> dict:
	if gunzip:
		in_file = gzip.open(ref_anno, 'rt')
	else:
		in_file = open(ref_anno)
	
	transcripts = {}
	for line in in_file:
		if line[0] != '#':
			chrom, _, feature, start, end, _, strand, _, attributes = line.strip().split('\t')
			attributes = parse_attributes(attributes, anno_type)
			if feature == 'transcript':
				transcripts[attributes['transcript_id']] = {'attributes':attributes, 'coords':(chrom, strand, int(start), int(end)), 'exons':[]}
			elif feature == 'exon':
				transcripts[attributes['transcript_id']]['exons'].append((int(start), int(end)))
	in_file.close()

	return transcripts


def build_exons_introns(transcripts:dict, out_exons_bed:str, out_introns_bed:str) -> dict:
	exon_genes, intron_genes = {}, {}
	for transcript_id in transcripts:
		chrom, strand, _, _ = transcripts[transcript_id]['coords']
		reverse_orientation = strand=='-'
		gene_orientation_exons = sorted(transcripts[transcript_id]['exons'], key=lambda e:e[0], reverse=reverse_orientation)
		for i in range(len(gene_orientation_exons)):
			exon_start, exon_end = gene_orientation_exons[i]
			exon = (chrom, strand, exon_start, exon_end)
			if exon not in exon_genes:
				exon_genes[exon] = []
			exon_genes[exon].append(f'{transcripts[transcript_id]["attributes"]["gene_id"]}-{i+1}')
			if i < len(gene_orientation_exons)-1:
				if strand == '+':
					intron_start, intron_end = gene_orientation_exons[i][1]+1, gene_orientation_exons[i+1][0]-1
				else:
					intron_start, intron_end = gene_orientation_exons[i+1][1]+1, gene_orientation_exons[i][0]-1
				intron = (chrom, strand, intron_start, intron_end)
				if intron not in intron_genes:
					intron_genes[intron] = []
				intron_genes[intron].append(f'{transcripts[transcript_id]["attributes"]["gene_id"]}-{i+1}')
	
	for region, region_genes, out_bed in zip(('exon', 'intron'), (exon_genes, intron_genes), (out_exons_bed, out_introns_bed)):
		with open(out_bed, 'w') as out_file:
			for region_id in region_genes:
				chrom, strand, start, end = region_id
				start -= 1
				region_gene_str = '|'.join(region_genes[region_id])
				bed_info = f'{region};{chrom};{start};{end};{strand};{region_gene_str}'
				out_file.write(f'{chrom}\t{start}\t{end}\t{bed_info}\t0\t{strand}\n')
	with open(out_introns_bed, 'w') as out_file:
		for intron in intron_genes:
			chrom, strand, start, end = intron
			start -= 1
			gene_str = '|'.join(intron_genes[intron])
			bed_info = f'{chrom};{start};{end};{strand};{gene_str}'
			out_file.write(f'{chrom}\t{start}\t{end}\t{bed_info}\t0\t{strand}\n')

	return intron_genes


def build_fivep(intron_genes:dict, ref_fasta:str, out_fivep_fasta:str) -> None:
	fivep_sites = set()
	for intron in intron_genes:
		chrom, strand, intron_start, intron_end = intron
		if intron_end-intron_start+1 > 20:
			fivep_pos = intron_start-1 if strand=='+' else intron_end-1
			fivep_sites.add((chrom, strand, fivep_pos))
	
	bedtools_input = ''
	for fivep in fivep_sites:
		chrom, strand, fivep_pos = fivep
		if strand == '+':
			window_start, window_end = fivep_pos-5, fivep_pos+20
		else:
			window_start, window_end = fivep_pos-19, fivep_pos+6
		fivep_site = f'{chrom};{fivep_pos};{strand}'
		fivep_line = f'{chrom}\t{window_start}\t{window_end}\t{fivep_site}\t0\t{strand}\n'
		bedtools_input += fivep_line
	
	bedtools_call = f'bedtools getfasta -s -tab -nameOnly -fi {ref_fasta} -bed -'
	bedtools_output = subprocess.run(bedtools_call.split(' '), input=bedtools_input, check=True, capture_output=True, text=True).stdout.strip()
	with open(out_fivep_fasta, 'w') as fasta_out:
		for fivep_line in bedtools_output.split('\n'):
			fivep_site, fivep_seq = fivep_line.split('\t')
			fivep_site, fivep_seq = fivep_site[:-3], fivep_seq.upper()
			fasta_out.write(f'>{fivep_site}\n{fivep_seq[5:]}\n')


#=============================================================================#
#                                    Main                                     #
#=============================================================================#
def main():
	# Parse arguments
	parser = argparse.ArgumentParser(prog='build_references',
								  	description='Build custom reference files and create reference directory for lariat mapping')
	
	parser.add_argument('-f', '--ref_fasta', required=True, help='Path to reference genome fasta file')
	parser.add_argument('-a', '--ref_anno', required=True, help='Path to reference gene annotation file in GTF or GFF format\n(may be gzipped with .gz extension)')
	parser.add_argument('-r', '--ref_repeatmasker', required=True, help='Path to BED file with RepeatMasker annotation of reference genome')
	parser.add_argument('-i', '--hisat2_index', required=True, help='Path to base name of hisat2 index of reference genome (i.e. everything before .1.ht2 extension)')
	parser.add_argument('-o', '--out_dir', required=True, help='Path to directory where reference files will be output (will be created if it does not exist)')

	args = parser.parse_args()
	ref_fasta, ref_anno, ref_repeatmasker, hisat2_index, out_dir = args.ref_fasta, args.ref_anno, args.ref_repeatmasker, args.hisat2_index, args.out_dir
	
	if os.path.isfile(f'{hisat2_index}.1.ht2'):
		hisat2_file_extensions = [f'.{i}.ht2' for i in range(1,9)]
	elif os.path.isfile(f'{hisat2_index}.1.ht2l'):
		hisat2_file_extensions = [f'.{i}.ht2l' for i in range(1,9)]
	else:
		raise FileNotFoundError(f'hisat2 index does not exist at {hisat2_index}')
	
	# Check input files
	ref_names = ['Genome fasta', 'Reference annotation', 'RepeatMasker BED'] + ['hisat2 index']*8
	ref_files = [ref_fasta, ref_anno, ref_repeatmasker] + [f'{hisat2_index}{ext}' for ext in hisat2_file_extensions]
	for rn, rf in zip(ref_names, ref_files):
		if not os.path.isfile(rf):
			raise FileNotFoundError(f'{rn} file does not exist at {rf}')
	
	# Make dir
	if not os.path.isdir(out_dir):
		os.mkdir(out_dir)

	# Make symbolic links to input reference files
	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Creating symbolic links to input reference files...')
	
	prev_ext, last_ext = ref_anno.split('.')[-2:]
	if last_ext == 'gz':
		gunzip, anno_type = True, prev_ext
		anno_file = f'annotation.{anno_type}.gz'
	else:
		gunzip, anno_type = False, last_ext
		anno_file = f'annotation.{anno_type}'
	
	ref_paths = [ref_fasta, ref_anno, ref_repeatmasker] + [f'{hisat2_index}{ext}' for ext in hisat2_file_extensions]
	symlink_paths = [os.path.join(out_dir, 'genome.fa'), os.path.join(out_dir, anno_file), os.path.join(out_dir, 'repeatmasker.bed')]
	symlink_paths += [os.path.join(out_dir, f'hisat2_index{ext}') for ext in hisat2_file_extensions]
	for rp, sp in zip(ref_paths, symlink_paths):
		if os.path.exists(sp):
			os.remove(sp)
		os.symlink(rp, sp)
	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Linking done.')

	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Parsing transcripts from annotation file...')
	transcripts = parse_transcript_info(ref_anno, anno_type, gunzip)
	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Done parsing transcripts.')
	
	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Processing exons and introns...')
	out_exons_bed = os.path.join(out_dir, 'exons.bed') 
	out_introns_bed = os.path.join(out_dir, 'introns.bed') 
	intron_genes = build_exons_introns(transcripts, out_exons_bed, out_introns_bed)
	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Done building exon and intron files.')

	print(time.strftime('%m/%d/%y - %H:%M:%S') + " | Processing 5' splice sites...")
	out_fivep_fasta = os.path.join(out_dir, 'fivep_sites.fa')
	build_fivep(intron_genes, ref_fasta, out_fivep_fasta)
	print(time.strftime('%m/%d/%y - %H:%M:%S') + " | Done building 5' splice site files.")

	print(time.strftime('%m/%d/%y - %H:%M:%S') + ' | Reference building complete.')

if __name__ == '__main__':
	main()
