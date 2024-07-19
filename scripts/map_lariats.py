import sys
import os
import subprocess

from scripts import functions



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
OUTPUT_BAM='{output_base}output.bam'
MAPPED_BAM='{output_base}mapped_reads.bam'
UNMAPPED_BAM='{output_base}unmapped_reads.bam'
RUN_DATA='{output_base}read_counts.tsv'
UNMAPPED_FASTA='{output_base}unmapped_reads.fa'
FIVEP_TO_READS='{output_base}fivep_to_reads.sam'
HEADS_FASTA='{output_base}heads.fa'
HEADS_TO_GENOME='{output_base}heads_to_genome.sam'
TAILS='{output_base}tails.tsv'
PUTATIVE_LARIATS='{output_base}putative_lariats.tsv'
FAILED_FIVEP='{output_base}failed_fivep_alignments.tsv'
FAILED_HEAD='{output_base}failed_head_alignments.tsv'
FAILED_LARIAT='{output_base}failed_lariat_alignments.tsv'

TEMP_FILES = (
				OUTPUT_BAM, MAPPED_BAM, UNMAPPED_BAM, UNMAPPED_FASTA,
				'{output_base}unmapped_reads.fa', '{output_base}.1.bt2l', '{output_base}.2.bt2l', '{output_base}.3.bt2l', '{output_base}.4.bt2l', '{output_base}.rev.1.bt2l', '{output_base}.rev.2.bt2l',
				FIVEP_TO_READS, HEADS_FASTA, HEADS_TO_GENOME, TAILS, PUTATIVE_LARIATS, 
				FAILED_FIVEP, FAILED_HEAD, FAILED_LARIAT, 
)




# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def delete_temp(log):
	log.debug('Deleting temporary files...')

	temp_files = [path_skeleton.format(output_base) for path_skeleton in TEMP_FILES]
	for file in temp_files:
		if os.path.isfile(file):
			os.remove(file)


def try_run_command(command, log, input=None):
	try:
		functions.run_command(command, input=input, log=log)
	except subprocess.CalledProcessError as error:
		if keep_temp is False:
			delete_temp(output_base, log)
		raise error


def linear_map(threads, genome_index, seq_type, read_files: list, log):
	log.info('Mapping reads to genome...')

	command= 'hisat2 --no-softclip -k 1 --max-seeds 20 --pen-noncansplice 0 --n-ceil L,0,0.05 --score-min L,0,-0.24 --bowtie2-dp 1' \
		f' --threads {threads} -x {genome_index}'
	if seq_type == 'single':
		command += f' -U {read_files[0]} | samtools view --bam --with-header --add-flags PAIRED,READ1'
	elif seq_type == 'paired':
		command += f' -1 {read_files[0]} -2 {read_files[1]} | samtools view --bam --with-header'
	command += ' > $output_bam'
	
	response = try_run_command(command, log)
	log.info(response)


def build_unmapped_reads_index(threads, log):
	log.info('Extracting unmapped reads...')

	# Make seperate BAM files
	command = f'samtools view --bam --with-header --exclude-flags 4 {OUTPUT_BAM.format(output_base)}' \
		 	f' > {MAPPED_BAM.format(output_base)}'
	try_run_command(command, log)
	command= f'samtools view --bam --with-header --require-flags 4 {OUTPUT_BAM.format(output_base)}' \
		 	f'> {UNMAPPED_BAM.format(output_base)}'
	try_run_command(command, log)

	# Check if 0 reads were left unmapped
	# If so, end early
	command = f'samtools view --count {UNMAPPED_BAM.format(output_base)}'
	unmapped_count = try_run_command(command, log)
	if unmapped_count == '0':
		log.warning('All reads mapped linearly, no reads remaining')
		if keep_temp is False:
			delete_temp(output_base, log)
		exit()

	# Create fasta file of unmapped reads 
	command = f'samtools fasta -N -o {UNMAPPED_FASTA.format(output_base)} {UNMAPPED_BAM.format(output_base)}'
	try_run_command(command, log)
	command = f'samtools faidx {UNMAPPED_FASTA.format(output_base)} {UNMAPPED_BAM.format(output_base)}'
	try_run_command(command, log)

	# Build a bowtie2 index 
	log.info('Building bowtie2 index of unmapped fasta...')
	command = f'bowtie2-build --quiet --large --threads {threads} {UNMAPPED_FASTA.format(output_base)} {UNMAPPED_FASTA.format(output_base)}'
	try_run_command(command, log)
	

def map_fivep_to_reads(threads, fivep_fasta, log):
	log.info("Mapping 5' splice sites to reads...")

	command = f'bowtie2 --end-to-end --sensitive --no-unal -f -k 10000 --score-min C,0,0' \
			f' --threads {threads} -x {UNMAPPED_FASTA.format(output_base)} -U {fivep_fasta}' \
			f' | samtools sort --verbosity 0 --output-fmt SAM -M  --threads {threads} ' \
			  '| samtools view' \
			f'> {FIVEP_TO_READS.format(output_base)}'
	response = try_run_command(command, log)
	log.info(response)


def map_heads_to_genome(threads, genome_index, log):
	log.info('Mapping heads to genome...')

	command = 'hisat2 --no-softclip --no-spliced-alignment --very-sensitive -f -k 100 --no-unal' \
			f' --threads {threads} -x {genome_index} -U {HEADS_FASTA.format(output_base)}' \
			f' | samtools sort --verbosity 0 --output-fmt SAM -n --threads {threads}' \
			 ' | samtools view' \
			f' > {HEADS_TO_GENOME.format(output_base)}'
	response = try_run_command(command, log)
	log.info(response)



# =============================================================================#
#                                    Main                                      #
# =============================================================================#
if __name__ == '__main__':
	# Get args and make keep_temp and output_base global variables
	# for reading ONLY so we don't need to pass them into basically every function
	global output_base
	global keep_temp
	output_base, output_prefix, threads, genome_index, genome_fasta, fivep_fasta, exons_tsv, introns_tsv, repeats_bed, keep_temp, ucsc_track, pipeline_dir, log_level, seq_type, read_files = sys.argv[1:]

	# Get logger and make it global for logging ONLY 
	# so we don't need to pass them into basically every function
	# global log
	log = functions.get_logger(log_level)
	log.debug(f'Args recieved: {sys.argv[1:]}')

	# Process args
	keep_temp = True if keep_temp=='True' else False
	ucsc_track = True if ucsc_track=='True' else False
	read_files = read_files.split(',')

	# Map reads to genome and keep unmapped reads. Lariat reads crossing the brachpoint will not be able to map to the gene they're from
	linear_map(threads, genome_index, seq_type, read_files, log)

	# Extract unmapped reads from linear mapping and build a bowtie2 index of them
	# If there are 0 unmapped reads, warn and end early
	build_unmapped_reads_index(threads, log)

	# Align fasta file of all 5' splice sites (first 20nts of introns) to unmapped reads index
	# We need to order the output SAM by reference (the read id, in this case) for the following filtering process
	map_fivep_to_reads(threads, fivep_fasta, log)

	# Extract reads with a mapped 5' splice site and trim it off
	log.info("Finding 5' read alignments and trimming reads...")
	command = f'python -u {pipeline_dir}/scripts/filter_fivep_aligns.py {threads} {genome_fasta} {fivep_fasta} {output_base} {log_level}'
	try_run_command(command, log)

	# Map read heads to genome
	map_heads_to_genome(threads, genome_index, log)

	# Filter head alignments
	log.info('Analyzing head alignments and outputting lariat table...')
	command = f'python -u {pipeline_dir}/scripts/filter_head_aligns.py {threads} {introns_tsv} {genome_fasta} {output_base} {log_level}'
	try_run_command(command, log)

	# Filter lariat mappings and choose 1 for each read
	log.info('Filtering putative lariat alignments...')
	command = f'python -u {pipeline_dir}/scripts/filter_lariats.py {genome_fasta} {repeats_bed} {output_base} {log_level}'
	try_run_command(command, log)

	# Make a custom track BED file of identified lariats 
	if ucsc_track is True:
		log.info('Making UCSC Genome Browser track...')
		command = f'python -u {pipeline_dir}/scripts/make_track.py {output_base} {log_level}'
		try_run_command(command, log)

	# Classify reads
	command = f'python -u {pipeline_dir}/scripts/classify_linear.py {output_base} {exons_tsv} {introns_tsv} {seq_type} {log_level}'
	try_run_command(command, log)
	command = f'python -u {pipeline_dir}/scripts/classify_nonlinear.py {output_base} {seq_type} {log_level}'
	try_run_command(command, log)

	# Delete the temporary files 
	if keep_temp is False:
		delete_temp(output_base, log)
	
	# End
	if output_prefix == 'None':
		log.info('Lariat mapping complete.')
	else:
		log.info(f'Lariat mapping complete for {output_prefix}.')