import os
import sys
import json

import pandas as pd
import pyfaidx
import pysam

import functions



#=============================================================================#
#                                   Globals                                   #
#=============================================================================#
# In files
ARGS_FILE = "{}settings.json"
OUTPUT_BAM_FILE = "{}output.bam"
LINEAR_COUNTS_FILE = "{}output.bam_count.tsv"
READ_CLASSES_FILE = "{}read_classes.tsv.gz"
# Out files
SUMMARY_FILE = "{}summary.txt"
READ_COUNTS_FILE = "{}read_counts.tsv"

SUMMARY_TEMPLATE = (
					"----------------------------------------\n"
					"                Metadata                \n"
					"----------------------------------------\n"
					"Version:\t{version}\n"
					#TODO: Time and resources
					"\n"
					"----------------------------------------\n"
					"                Settings                \n"
					"----------------------------------------\n"
					"Input reads:\t{input_reads}\n"
					"Input type:\t{seq_type}\n"
					"Input strandedness:\t{strand}\n"
					"Reference directory:\t{ref_dir}\n"
					"Reference HISAT2 index:\t{ref_h2index}\n"
					"Reference genome FASTA:\t{ref_fasta}\n"
					"Reference 5'ss FASTA:\t{ref_5p_fasta}\n"
					"Reference introns:\t{ref_introns}\n"
					"Reference RepeatMasker:\t{ref_repeatmasker}\n"
					"Output path:\t{output_base}\n"
					"Threads:\t{threads}\n"
					"Make UCSC track:\t{ucsc_track}\n"
					"Keep read classes file:\t{keep_classes}\n"
					"Keep temporary files:\t{keep_temp}\n"
					"Skip version check:\t{skip_version_check}\n"
					"Log level:\t{log_level}\n"
					"\n"
					"----------------------------------------\n"
					"              Read classes              \n"
					"----------------------------------------\n"
					"mRNA:\t{exon_exon_junc}\n"
					"pre-mRNA:\t{exon_intron_junc}\n"
					"Exonic:\t{exon_only}\n"
					"Intronic:\t{intron_only}\n"
					# "Genic, ambiguous:\t{ambig}\n"
					# "Intergenic:\t{intergenic}\n"
					"Unmapped:\t{Unmapped}\n"
					"Unmapped with 5'ss alignment:\t{Unmapped_with_5ss_alignment}\n"
					"Template-switching:\t{Template_switching}\n"
					"Circularized intron:\t{Circularized_intron}\n"
					"In repetitive region:\t{In_repetitive_region}\n"
					"Lariat:\t{Lariat}\n"
					"\n"
					"----------------------------------------\n"
					"      Read count after each stage       \n"
					"----------------------------------------\n"
					"Input:\t{input_count}\n"
					"Linear mapping:\t{Linear_mapping}\n"
					"5'ss mapping:\t{5ss_mapping}\n"
					"5'ss alignment filtering:\t{5ss_alignment_filtering}\n"
					"Head mapping:\t{Head_mapping}\n"
					"Head alignment filtering:\t{Head_alignment_filtering}\n"
					"Lariat filtering:\t{Lariat_filtering}\n"
					"\n"
					"----------------------------------------\n"
					"          Additional statistics         \n"
					"----------------------------------------\n"
					"pre-mRNA/mRNA = {pre_ratio:.1%}\n"
					"Lariat RPM = {lariat_rpm:.3g}\n"
					"Circularized intron RPM = {circ_rpm:.3g}\n"
					"Lariat reads, genomic_bp_nt = A:\t{A:.1%}\n"
					"Lariat reads, genomic_bp_nt = C:\t{C:.1%}\n"
					"Lariat reads, genomic_bp_nt = G:\t{G:.1%}\n"
					"Lariat reads, genomic_bp_nt = T:\t{T:.1%}\n"
					"Lariat reads, genomic_bp_nt = N:\t{N:.1%}\n"
					"Lariat reads, genomic_bp_nt ≠ read_bp_nt:\t{bp_mismatch:.1%}\n"
					"Lariat reads, |bp_dist_to_threep| ≤ 70:\t{within_70:.1%}\n"
)
READ_COUNTS_TEMPLATE = (
						"Category\tSubcategory\tReads\n"
						"Input\tTotal\t{input_count}\n"
						"Linearly mapped\tTotal\t{Linear}\n"
						"Linearly mapped\tmRNA\t{exon_exon_junc}\n"
						"Linearly mapped\tpre-mRNA\t{exon_intron_junc}\n"
						"Linearly mapped\tExonic\t{exon_only}\n"
						"Linearly mapped\tIntronic\t{intron_only}\n"
						# "Linearly mapped\tGenic, ambiguous\t{ambig}\n"
						# "Linearly mapped\tIntergenic\t{intergenic}\n"
						"Not linearly mapped\tTotal\t{not_linear}\n"
						"Not linearly mapped\tUnmapped\t{Unmapped}\n"
						"Not linearly mapped\tUnmapped with 5'ss alignment\t{Unmapped_with_5ss_alignment}\n"
						"Not linearly mapped\tTemplate-switching\t{Template_switching}\n"
						"Not linearly mapped\tCircularized intron\t{Circularized_intron}\n"
						"Not linearly mapped\tIn repetitive region\t{In_repetitive_region}\n"
						"Not linearly mapped\tLariat\t{Lariat}\n"
						"Only one mate linearly mapped\tTotal\t{mixed_pairs}\n"
						"Read count after stage\tLinear mapping:\t{Linear_mapping}\n"
						"Read count after stage\t5'ss mapping:\t{5ss_mapping}\n"
						"Read count after stage\t5'ss alignment filtering:\t{5ss_alignment_filtering}\n"
						"Read count after stage\tHead mapping:\t{Head_mapping}\n"
						"Read count after stage\tHead alignment filtering:\t{Head_alignment_filtering}\n"
						"Read count after stage\tLariat filtering:\t{Lariat_filtering}\n"
)

NOLINEAR_READ_CLASSES = ("Unmapped", "Unmapped_with_5ss_alignment", 'In_repetitive_region', 
						'Template_switching', 'Circularized_intron', 'Lariat')
#=============================================================================#
#                               Functions                                     #
#=============================================================================#
def add_mapped_reads(output_base:str, seq_type:str, log) -> int:
	'''
	Get the number of reads that mapped to the reference genome
	'''
	command = f'samtools view --count --exclude-flags 12 {OUTPUT_BAM_FILE.format(output_base)}'
	count = int(functions.run_command(command, log=log))
	if seq_type == 'paired':
		count = count//2

	return count



#=============================================================================#
#                                    Main                                     #
#=============================================================================#
if __name__ == '__main__':
	# Get args
	output_base, log_level, seq_type = sys.argv[1:]

	# Get logger
	log = functions.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}')

	# Initialise stats dict
	stats = {}

	# Run information
	stats['version'] = functions.version()
	with open(ARGS_FILE.format(output_base), 'r') as json_file:
		settings = json.load(json_file)
	stats.update(settings)

	# Add input read count
	r1 = pysam.FastxFile(settings['input_reads'].split(',')[0])
	stats['input_count'] = sum(1 for read in r1)

	# Add linear mapping count
	stats['Linear'] = add_mapped_reads(output_base, seq_type, log)

	# Add linear read class counts
	linear_counts = pd.read_csv(LINEAR_COUNTS_FILE.format(output_base), sep='\t', index_col=0)
	# linear_counts['rowsum'] = linear_counts.sum(axis=1) - linear_counts.gene
	# linear_counts['ambig'] = linear_counts.gene - linear_counts.rowsum
	stats.update(linear_counts.apply(sum).to_dict())
	# stats['intergenic'] = stats['Linear'] - stats['gene']

	# Add nonlinear read class counts
	read_classes = pd.read_csv(READ_CLASSES_FILE.format(output_base), sep='\t', na_filter=False)
	read_classes.read_class = (read_classes.read_class
								.str.replace(' ', '_')
								.str.replace('-', '_')
								.str.replace("'", ''))
	classes = (pd.Categorical(read_classes.read_class, categories=NOLINEAR_READ_CLASSES, ordered=True)
				.value_counts()
				.to_dict())
	log.debug(f'Read classes: {classes}')
	stats.update(classes)
	stats['not_linear'] = stats['input_count'] - stats['Linear']

	# Add read counts after each stage
	unmapped = set()
	for rid in pyfaidx.Fasta(f'{output_base}unmapped_reads.fa', as_raw=True):
		unmapped.add(rid.name[:-2])
	stats['Linear_mapping'] = len(unmapped)

	fivep_maps = set()
	with open(f'{output_base}fivep_to_reads.sam', 'r') as r:
		for line in r:
			rid = line.split('\t')[2][:-2]
			fivep_maps.add(rid)
	stats['5ss_mapping'] = len(fivep_maps)

	tails = pd.read_csv(f'{output_base}tails.tsv', sep='\t', usecols=['read_id']).read_id
	stats['5ss_alignment_filtering'] = tails.str.slice(0,-6).nunique()

	head_maps = set()
	with open(f'{output_base}heads_to_genome.sam', 'r') as r:
		for line in r:
			rid = line.split('\t')[0][:-6]
			head_maps.add(rid)
	stats['Head_mapping'] = len(head_maps)

	putative_lariats = pd.read_csv(f'{output_base}putative_lariats.tsv', sep='\t', usecols=['read_id']).read_id
	stats['Head_alignment_filtering'] = putative_lariats.str.slice(0,-6).nunique()

	filtered_lariats = pd.read_csv(f'{output_base}lariat_reads.tsv', sep='\t', usecols=['read_id']).read_id
	stats['Lariat_filtering'] = filtered_lariats.nunique()

	# For paired-end data, add count of reads where one mate mapped linearly in the 
	# initial mapping and the other didn't
	if seq_type == 'single':
		stats['mixed_pairs'] = 'N/A'
	elif seq_type == 'paired':
		mp = 0
		for align in pysam.AlignmentFile(OUTPUT_BAM_FILE.format(output_base), 'rb'):
			if align.is_mapped and align.mate_is_unmapped:
				mp += 1
		stats['mixed_pairs'] = mp

	# Add additional info
	stats['pre_ratio'] = stats['exon_intron_junc'] / stats['exon_exon_junc']
	stats['lariat_rpm'] = stats['Lariat'] / stats['Linear'] * 1e6
	stats['circ_rpm'] = stats['Circularized_intron'] / stats['Linear'] * 1e6
	lariat_reads = pd.read_csv(f'{output_base}lariat_reads.tsv', sep='\t')
	lariat_reads.genomic_bp_nt = pd.Categorical(lariat_reads.genomic_bp_nt, categories=['A', 'C', 'G', 'T', 'N'])
	stats.update(lariat_reads.genomic_bp_nt.value_counts(normalize=True).to_dict())
	stats['within_70'] = lariat_reads.bp_dist_to_threep.abs().le(70).sum()/lariat_reads.shape[0]
	stats['bp_mismatch'] = lariat_reads.genomic_bp_nt.ne(lariat_reads.read_bp_nt).sum()/lariat_reads.shape[0]
	
	
	# Write summary info to file
	log.debug(f'Summary stats: {stats}')
	with open(SUMMARY_FILE.format(output_base), 'w') as w:
		w.write(SUMMARY_TEMPLATE.format(**stats))

	# Write read counts to file
	read_count_stats = {key: stats[key] for key in stats if '{'+key+'}' in READ_COUNTS_TEMPLATE}
	log.debug(f'Read count stats: {read_count_stats}')
	with open(READ_COUNTS_FILE.format(output_base), 'w') as w:
		w.write(READ_COUNTS_TEMPLATE.format(**read_count_stats))

	log.debug('End of script')

