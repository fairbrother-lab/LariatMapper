import os
import sys
import json

import pandas as pd
import pyfaidx
import pysam

import utils
from filter_head_aligns import TEMP_SWITCH_COLS, TRANS_SPLICING_COLS, CIRCULARS_COLS
from filter_lariats import FINAL_RESULTS_COLS





#=============================================================================#
#                                   Globals                                   #
#=============================================================================#
# In files
ARGS_FILE = "{}settings.json"
OUTPUT_BAM_FILE = "{}output.bam"
LINEAR_COUNTS_FILE = "{}output.bam_summary_count.tsv"
READ_CLASSES_FILE = "{}read_classes.tsv.gz"
# Out files
SUMMARY_FILE = "{}summary.txt"
READ_COUNTS_FILE = "{}read_counts.tsv"

SUMMARY_TEMPLATE = (
					"----------------------------------------\n"
					"                Metadata                \n"
					"----------------------------------------\n"
					"Version: {version}\n"
					#TODO: Time and resources
					"\n"
					"----------------------------------------\n"
					"                Settings                \n"
					"----------------------------------------\n"
					"Input reads: {input_reads}\n"
					"Input type: {seq_type}\n"
					"Input strandedness: {strand}\n"
					"Reference directory: {ref_dir}\n"
					"Reference HISAT2 index: {ref_h2index}\n"
					"Reference genome FASTA: {ref_fasta}\n"
					"Reference 5'ss FASTA: {ref_5p_fasta}\n"
					"Reference introns: {ref_introns}\n"
					"Strandedness: {strand}\n"
					"Reference RepeatMasker: {ref_repeatmasker}\n"
					"BP position correction: {bp_correction}\n"
					"BP position correction files: {bp_correction_files}\n"
					"Output path: {output_base}\n"
					"Threads: {threads}\n"
					"Make UCSC track: {ucsc_track}\n"
					"Keep read classes file: {keep_classes}\n"
					"Keep temporary files: {keep_temp}\n"
					"Log level: {log_level}\n"
					"\n"
					"----------------------------------------\n"
					"              Read classes              \n"
					"----------------------------------------\n"
					"Linear, exon-exon junction: {exon_exon_junc}\n"
					"Linear, exon-intron junction: {exon_intron_junc}\n"
					"Linear, exon only: {exon_only}\n"
					"Linear, intron only: {intron_only}\n"
					"Linear, intergenic or ambiguous: {intergenic_ambiguous}\n"
					"No alignment: {No_alignment}\n"
					"Fivep alignment: {Fivep_alignment}\n"
					"Template-switching: {Template_switching}\n"
					"Trans-splicing: {Trans_splicing}\n"
					"Circularized intron: {Circularized_intron}\n"
					"In repetitive region: {In_repetitive_region}\n"
					"Lariat: {Lariat}\n"
					"\n"
					"----------------------------------------\n"
					"      Read count after each stage       \n"
					"----------------------------------------\n"
					"Input: {input_count}\n"
					"Linear mapping: {Linear_mapping}\n"
					"Fivep mapping: {5ss_mapping}\n"
					"Fivep alignment filtering: {5ss_alignment_filtering}\n"
					"Head mapping: {Head_mapping}\n"
					"Head alignment filtering: {Head_alignment_filtering}\n"
					"Lariat filtering: {Lariat_filtering}\n"
					"\n"
					"----------------------------------------\n"
					"          Additional statistics         \n"
					"----------------------------------------\n"
					"Exon-exon/exon-intron ratio: {pre_ratio}\n"
					"Lariat RPM: {lariat_rpm}\n"
					"Circularized intron RPM: {circ_rpm}\n"
					"Lariat reads, genomic_bp_nt = A: {A:.1%}\n"
					"Lariat reads, genomic_bp_nt = C: {C:.1%}\n"
					"Lariat reads, genomic_bp_nt = G: {G:.1%}\n"
					"Lariat reads, genomic_bp_nt = T: {T:.1%}\n"
					"Lariat reads, genomic_bp_nt = N: {N:.1%}\n"
					"Lariat reads, genomic_bp_nt ≠ read_bp_nt: {bp_mismatch:.1%}\n"
					"Lariat reads, |bp_dist_to_threep| ≤ 70: {within_70:.1%}\n"
					"Lariat reads, branchpoint position corrected: {bp_corrected_prop}\n"
)
READ_COUNTS_TEMPLATE = (
						"Category\tSubcategory\tReads\n"
						"Input\tTotal\t{input_count}\n"
						"Linearly_mapped\tTotal\t{Linear}\n"
						"Linearly_mapped\tExon_exon_junction\t{exon_exon_junc}\n"
						"Linearly_mapped\tExon_intron_junction\t{exon_intron_junc}\n"
						"Linearly_mapped\tExon_only\t{exon_only}\n"
						"Linearly_mapped\tIntron_only\t{intron_only}\n"
						"Linearly_mapped\tIntergenic_or_ambiguous\t{intergenic_ambiguous}\n"
						"Not_linearly_mapped\tTotal\t{not_linear}\n"
						"Not_linearly_mapped\tNo_alignment\t{No_alignment}\n"
						"Not_linearly_mapped\tFivep_alignment\t{Fivep_alignment}\n"
						"Not_linearly_mapped\tTemplate_switching\t{Template_switching}\n"
						"Not_linearly_mapped\tTrans_splicing\t{Trans_splicing}\n"
						"Not_linearly_mapped\tCircularized_intron\t{Circularized_intron}\n"
						"Not_linearly_mapped\tIn_repetitive_region\t{In_repetitive_region}\n"
						"Not_linearly_mapped\tLariat\t{Lariat}\n"
						"Other\tOne_mate_linearly_mapped\t{mixed_pairs}\n"
						"Read_count_after_stage\tLinear_mapping\t{Linear_mapping}\n"
						"Read_count_after_stage\tFivep_mapping\t{5ss_mapping}\n"
						"Read_count_after_stage\tFivep_alignment_filtering\t{5ss_alignment_filtering}\n"
						"Read_count_after_stage\tHead_mapping\t{Head_mapping}\n"
						"Read_count_after_stage\tHead_alignment_filtering\t{Head_alignment_filtering}\n"
						"Read_count_after_stage\tLariat_filtering\t{Lariat_filtering}\n"
)

NONLINEAR_READ_CLASSES = ("No_alignment", "Fivep_alignment", 'In_repetitive_region', 
						'Template_switching', 'Trans_splicing', 'Circularized_intron', 'Lariat')





#=============================================================================#
#                                    Main                                     #
#=============================================================================#
if __name__ == '__main__':
	# Get args
	output_base, log_level, seq_type = sys.argv[1:]

	# Get logger
	log = utils.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}')

	# Make neccesary output files if they are absent due to the run ending early
	if not os.path.isfile(f'{output_base}lariat_reads.tsv'):
		with open(f'{output_base}lariat_reads.tsv', 'w') as w:
			w.write('\t'.join(FINAL_RESULTS_COLS) + '\n')
	if not os.path.isfile(f'{output_base}trans_splicing_reads.tsv'):
		with open(f'{output_base}trans_splicing_reads.tsv', 'w') as w:
			w.write('\t'.join(TRANS_SPLICING_COLS) + '\n')
	if not os.path.isfile(f'{output_base}template_switching_reads.tsv'):
		with open(f'{output_base}template_switching_reads.tsv', 'w') as w:
			w.write('\t'.join(TEMP_SWITCH_COLS) + '\n')
	if not os.path.isfile(f'{output_base}circularized_intron_reads.tsv'):
		with open(f'{output_base}circularized_intron_reads.tsv', 'w') as w:
			w.write('\t'.join(CIRCULARS_COLS) + '\n')


	# Initialise stats dict
	stats = {}

	# Run information
	stats['version'] = utils.version()
	with open(ARGS_FILE.format(output_base), 'r') as json_file:
		settings = json.load(json_file)
	stats.update(settings)
	if settings['pwm_correction'] != "":
		stats['bp_correction'] = "Position weight matrix"
		stats['bp_correction_files'] = settings['pwm_correction']
	elif settings['model_correction'] != "":
		stats['bp_correction'] = "Machine-learning model"
		stats['bp_correction_files'] = settings['model_correction']
	else:
		stats['bp_correction'] = "False"
		stats['bp_correction_files'] = ""

	# Add input read count
	r1 = pysam.FastxFile(settings['input_reads'].split(',')[0])
	stats['input_count'] = sum(1 for read in r1)

	# Add linear read class counts
	linear_counts = pd.read_csv(LINEAR_COUNTS_FILE.format(output_base), sep='\t')
	stats.update(linear_counts.iloc[0].to_dict())
	stats['Linear'] = stats['total_reads']
	
	# Add nonlinear read class counts
	read_classes = pd.read_csv(READ_CLASSES_FILE.format(output_base), sep='\t', na_filter=False)
	read_classes.read_class = (read_classes.read_class
								.str.replace(' ', '_')
								.str.replace('-', '_')
								.str.replace("'", ''))
	classes = (pd.Categorical(read_classes.read_class, categories=NONLINEAR_READ_CLASSES, ordered=True)
				.value_counts()
				.to_dict())
	log.debug(f'Read classes: {classes}')
	stats.update(classes)
	stats['not_linear'] = stats['input_count'] - stats['Linear']

	# Add read counts after each stage
	unmapped = set()
	if os.path.isfile(f'{output_base}unmapped_reads.fa'):
		for rid in pyfaidx.Fasta(f'{output_base}unmapped_reads.fa', as_raw=True, rebuild=False):
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
		stats['pre_ratio'] = f"{(2*stats['exon_exon_junc']) / stats['exon_intron_junc']:.4g}"
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
		stats['bp_corrected_prop'] = 'N/A'
	else:
		lariat_reads.genomic_bp_nt = pd.Categorical(lariat_reads.genomic_bp_nt, categories=['A', 'C', 'G', 'T', 'N'])
		stats.update(lariat_reads.genomic_bp_nt.value_counts(normalize=True).to_dict())
		stats['within_70'] = lariat_reads.bp_dist_to_threep.abs().le(70).sum()/lariat_reads.shape[0]
		stats['bp_mismatch'] = lariat_reads.genomic_bp_nt.ne(lariat_reads.read_bp_nt).sum()/lariat_reads.shape[0]
		if 'corrected' in lariat_reads.columns:
			stats['bp_corrected_prop'] = f"{lariat_reads.corrected.sum()/lariat_reads.shape[0]:.1%}"
		else:
			stats['bp_corrected_prop'] = 'N/A'
		

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

