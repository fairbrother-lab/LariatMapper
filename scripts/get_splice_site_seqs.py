
import sys
import gzip
from subprocess import run

chroms = set(['chr'+str(i) for i in range(1,23)] + ['chrY', 'chrX', 'chrM'])

def output_splice_site_seqs(intron_file, genome_file, out_base):

	if intron_file[-2:] == 'gz':
		in_file = gzip.open(intron_file, 'rt')
	else:
		in_file = open(intron_file)

	fivep_done, threep_done = set(), set()
	with open(out_base+'.threep_sites.bed', 'w') as threep_out, open(out_base+'.fivep_sites.bed', 'w') as fivep_out:
		with open(out_base+'.threep_seq_lens.txt', 'w') as len_out, open(out_base+'.fivep_upstream.bed', 'w') as upstream_out:
			for line in in_file:
				chrom, start, end, _, _, strand = line.strip().split()
				start, end = int(start), int(end)
				if chrom in chroms and end-start >= 20:
					if strand == '+':
						fivep_start, fivep_end, fivep_site = start, start+20, start
						upstream_start, upstream_end = start-5, start
					else:
						fivep_start, fivep_end, fivep_site = end-20, end, end
						upstream_start, upstream_end = end, end+5
					if (chrom, fivep_site, strand) not in fivep_done:
						fivep_out.write('{0}\t{1}\t{2}\t{0};{1};{2};{3}\t0\t{3}\n'.format(chrom, fivep_start, fivep_end, strand))
						upstream_out.write('{0}\t{1}\t{2}\t{0};{3};{4};{5}\t0\t{5}\n'.format(chrom, upstream_start, upstream_end, fivep_start, fivep_end, strand))
					if end-start+1 > 250:
						if strand == '+':
							threep_start, threep_end, threep_site = end-250, end, end
						else:
							threep_start, threep_end, threep_site = start, start+250, start
					else:
						threep_start, threep_end = start, end
						threep_site = end if strand == '+' else start
					if (chrom, threep_site, strand) not in threep_done:
						threep_out.write('{0}\t{1}\t{2}\t{0};{1};{2};{3}\t0\t{3}\n'.format(chrom, threep_start, threep_end, strand))
						len_out.write('{}\t{}\t{}\t{}\n'.format(chrom, threep_site, strand, threep_end-threep_start))
					fivep_done.add((chrom, fivep_site, strand))
					threep_done.add((chrom, threep_site, strand))

	fivep_cmd = f'bedtools getfasta -nameOnly -s -fo {out_base}.fivep_sites.fa -fi {genome_file} -bed {out_base}.fivep_sites.bed'
	#run(fivep_cmd.split(' '))

	threep_cmd = f'bedtools getfasta -nameOnly -s -fo {out_base}.threep_sites.fa -fi {genome_file} -bed {out_base}.threep_sites.bed'
	#run(threep_cmd.split(' '))

	upstream_cmd = f'bedtools getfasta -tab -nameOnly -s -fo {out_base}.fivep_upstream.txt -fi {genome_file} -bed {out_base}.fivep_upstream.bed'
	run(upstream_cmd.split(' '))

	run(f'rm {out_base}.fivep_sites.bed {out_base}.threep_sites.bed {out_base}.fivep_upstream.bed'.split(' '))

if __name__ == '__main__' :

	intron_file, genome_file, out_base = sys.argv[1:]
	output_splice_site_seqs(intron_file, genome_file, out_base)




