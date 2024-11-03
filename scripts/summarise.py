import os
import sys
import json

import pandas as pd
import pyfaidx
import pysam

import functions
from filter_head_aligns import CIRCULARS_COLS, TEMP_SWITCH_COLS
from filter_lariats import FINAL_RESULTS_COLS



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
					"Other linear (placeholder):\t{linear_other}\n"
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
					"mRNA/pre-mRNA = {pre_ratio}\n"
					"Lariat RPM = {lariat_rpm}\n"
					"Circularized intron RPM = {circ_rpm}\n"
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
						"Linearly mapped\tOther (placeholder)\t{linear_other}\n"
						"Not linearly mapped\tTotal\t{not_linear}\n"
						"Not linearly mapped\tUnmapped\t{Unmapped}\n"
						"Not linearly mapped\tUnmapped with 5'ss alignment\t{Unmapped_with_5ss_alignment}\n"
						"Not linearly mapped\tTemplate-switching\t{Template_switching}\n"
						"Not linearly mapped\tCircularized intron\t{Circularized_intron}\n"
						"Not linearly mapped\tIn repetitive region\t{In_repetitive_region}\n"
						"Not linearly mapped\tLariat\t{Lariat}\n"
						"Other\tOne mate linearly mapped\t{mixed_pairs}\n"
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

	# Make neccesary output files if they are absent due to the run ending early
	if not os.path.isfile(f'{output_base}lariat_reads.tsv'):
		with open(f'{output_base}lariat_reads.tsv', 'w') as w:
			w.write('\t'.join(FINAL_RESULTS_COLS))
	if not os.path.isfile(f'{output_base}circularized_intron_reads.tsv'):
		with open(f'{output_base}circularized_intron_reads.tsv', 'w') as w:
			w.write('\t'.join(CIRCULARS_COLS))
	if not os.path.isfile(f'{output_base}template_switching_reads.tsv'):
		with open(f'{output_base}template_switching_reads.tsv', 'w') as w:
			w.write('\t'.join(TEMP_SWITCH_COLS))


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
	stats['linear_other'] = stats['Linear'] - stats['exon_exon_junc'] - stats['exon_intron_junc'] - stats['exon_only'] - stats['intron_only']

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
	if os.path.isfile(f'{output_base}unmapped_reads.fa'):
		for rid in pyfaidx.Fasta(f'{output_base}unmapped_reads.fa', as_raw=True):
			unmapped.add(rid.name[:-2])
	stats['Linear_mapping'] = len(unmapped)

	fivep_maps = set()
	if os.path.isfile(f'{output_base}fivep_to_reads.sam'):
		with open(f'{output_base}fivep_to_reads.sam', 'r') as r:
			for line in r:
				rid = line.split('\t')[2][:-2]
				fivep_maps.add(rid)
	stats['5ss_mapping'] = len(fivep_maps)

	if os.path.isfile(f'{output_base}tails.tsv'):
		tails = pd.read_csv(f'{output_base}tails.tsv', sep='\t', usecols=['read_id']).read_id
		stats['5ss_alignment_filtering'] = tails.str.slice(0,-6).nunique()
	else:
		stats['5ss_alignment_filtering'] = 0

	head_maps = set()
	if os.path.isfile(f'{output_base}heads_to_genome.sam'):
		with open(f'{output_base}heads_to_genome.sam', 'r') as r:
			for line in r:
				rid = line.split('\t')[0][:-6]
				head_maps.add(rid)
	stats['Head_mapping'] = len(head_maps)

	if os.path.isfile(f'{output_base}putative_lariats.tsv'):
		putative_lariats = pd.read_csv(f'{output_base}putative_lariats.tsv', sep='\t', usecols=['read_id']).read_id
		stats['Head_alignment_filtering'] = putative_lariats.str.slice(0,-6).nunique()
	else:
		stats['Head_alignment_filtering'] = 0

	lariat_reads = pd.read_csv(f'{output_base}lariat_reads.tsv', sep='\t')
	stats['Lariat_filtering'] = lariat_reads.read_id.nunique()

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
	if stats['exon_intron_junc'] > 0:
		stats['pre_ratio'] = f"{stats['exon_exon_junc'] / stats['exon_intron_junc']:.4g}"
	else:
		stats['pre_ratio'] = 'N/A'

	if stats['Linear'] > 0:
		stats['lariat_rpm'] = f"{stats['Lariat'] / stats['Linear'] * 1e6:.4g}"
		stats['circ_rpm'] = f"{stats['Circularized_intron'] / stats['Linear'] * 1e6:.4g}"
	else:
		stats['lariat_rpm'] = 'N/A'
		stats['circ_rpm'] = 'N/A'

	if lariat_reads.shape[0]==0:
		stats.update({bp: 0 for bp in ['A', 'C', 'G', 'T', 'N']})
		stats['within_70'] = 0
		stats['bp_mismatch'] = 0
	else:
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


	# Check read counts output for consistency, warn if inconsistent
	counts = pd.read_csv(READ_COUNTS_FILE.format(output_base), sep='\t', index_col=[0,1])
	lin_tot = counts.at[('Linearly mapped', 'Total'), 'Reads']
	lin_subs = counts.loc['Linearly mapped'].Reads.sum() - lin_tot
	nonlin_tot = counts.at[('Not linearly mapped', 'Total'), 'Reads']
	nonlin_subs = counts.loc['Not linearly mapped'].Reads.sum() - nonlin_tot
	inp = counts.at[('Input', 'Total'), 'Reads'] 
	if lin_tot != lin_subs:
		raise RuntimeError(f'Linearly mapped read subcategories do not sum to total ({lin_subs:,} ≠ {lin_tot:,})')
	if nonlin_tot != nonlin_subs:
		raise RuntimeError(f'Not linearly mapped read subcategories do not sum to total ({nonlin_subs:,} ≠ {nonlin_tot:,})')
	if inp != lin_tot + nonlin_tot:
		raise RuntimeError(f'Input reads do not sum to linearly and non-linearly mapped reads ({inp:,} ≠ {lin_tot + nonlin_tot:,})')
	survival_counts = counts.loc[counts.index.get_level_values('Category')=='Read count after stage'].Reads.values
	if not all(survival_counts[i] >= survival_counts[i+1] for i in range(len(survival_counts)-1)):
		raise RuntimeError(f'Read counts after each stage are not in decreasing order')


	log.debug('End of script')

