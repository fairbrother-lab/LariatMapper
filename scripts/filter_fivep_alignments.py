
import sys
from pyfaidx import Fasta
from collections import Counter

comp_nts = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
def reverse_complement(seq):
	return ''.join([comp_nts[seq[i]] for i in range(len(seq)-1,-1,-1)])

def filter_fivep_reads(fivep_to_reads, fivep_upstream, read_file, trimmed_out_path, info_out_path):

	fivep_upstream_seqs = {}
	with open(fivep_upstream) as in_file:
		for line in in_file:
			fivep_site, seq = line.strip().split('\t')
			fivep_upstream_seqs[fivep_site[:-3]] = seq.upper()

	read_sites, site_coords = {}, {}
	with open(fivep_to_reads) as fivep_file:
		for line in fivep_file:
			alignment_info = line.strip().split('\t')
			fivep_site, flag, rid, reference_start, _, read_cig = alignment_info[:6]
			for alignment_tag in alignment_info[11:]:
				if alignment_tag[:2] == 'XM':
					num_mismatch = int(alignment_tag.split(':')[-1])
			if num_mismatch == 0 and read_cig == '20M':
				reference_start = int(reference_start)-1
				bit_flags = bin(int(flag))
				is_reverse = True if len(bit_flags)>=7 and bit_flags[-5]=='1' else False
				fivep_site = fivep_site[:-3]
				if rid not in read_sites:
					read_sites[rid] = set()
					site_coords[rid] = {}
				read_sites[rid].add(fivep_site)
				site_coords[rid][fivep_site] = (reference_start, reference_start+20, is_reverse)

	read_fasta = Fasta(read_file, as_raw=True)
	with open(trimmed_out_path, 'w') as trimmed_out_file, open(info_out_path, 'w') as info_out_file:
		info_out_file.write('read_id\tread_seq\tfivep_seq\tfivep_sites\tfivep_first\tread_fivep_start\tread_fivep_end\n')
		for rid in read_sites:
			read_seq = read_fasta[rid][:]
			fivep_pass = {True:[], False:[]}
			for fp in site_coords[rid]:
				fivep_start, fivep_end, is_reverse = site_coords[rid][fp]
				if is_reverse:
					read_upstream = read_seq[fivep_end:fivep_end+5].upper()
					upstream_mismatch = read_upstream != reverse_complement(fivep_upstream_seqs[fp])
				else:
					read_upstream = read_seq[fivep_start-5:fivep_start].upper()
					upstream_mismatch = read_upstream != fivep_upstream_seqs[fp]
				if upstream_mismatch:
					fivep_pass[is_reverse].append((fp, site_coords[rid][fp]))

			for is_reverse in fivep_pass:
				if len(fivep_pass[is_reverse]) > 0:
					if is_reverse:
						fivep_start, fivep_end, _ = max(fivep_pass[is_reverse], key=lambda fp:fp[1][0])[1]
						fivep_pass_sub = [fp for fp in fivep_pass[is_reverse] if fp[1][0]==fivep_start]
						trim_seq, fivep_seq = read_seq[fivep_end:], reverse_complement(read_seq[fivep_start:fivep_end])
					else:
						fivep_start, fivep_end, _ = min(fivep_pass[is_reverse], key=lambda fp:fp[1][0])[1]
						fivep_pass_sub = [fp for fp in fivep_pass[is_reverse] if fp[1][0]==fivep_start]
						trim_seq, fivep_seq = read_seq[:fivep_start], read_seq[fivep_start:fivep_end]

					if len(trim_seq) >= 20:
						out_rid = rid + {True:'_rev', False:'_for'}[is_reverse]
						trimmed_out_file.write('>{}\n{}\n'.format(out_rid, trim_seq))
						fivep_sites = ','.join([fp[0] for fp in fivep_pass_sub])
						info_out_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(out_rid, read_seq, fivep_seq, fivep_sites, is_reverse, fivep_start, fivep_end))

if __name__ == '__main__' :
	
	read_file, fivep_to_reads, fivep_upstream, trimmed_out_path, info_out_path = sys.argv[1:]
	
	filter_fivep_reads(fivep_to_reads, fivep_upstream, read_file, trimmed_out_path, info_out_path)



			

