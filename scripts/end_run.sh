#!/bin/bash

#=============================================================================#
#                                  Arguments                                  #
#=============================================================================#
# Reference directory
REF_DIR="${1}"
# Output directory 
OUTPUT_BASE="${2}"
# Keep read_classes.tsv.gz. true or false, default false
KEEP_CLASSES="${3}"
# Keep the temp files created during the run. true or false, default false
KEEP_TEMP="${4}"
# Type of sequencing data. "single" or "paired"
SEQ_TYPE="${5}"
# File for writing logs
LOG_FILE="${6}"
# Level for python logging. "DEBUG", "INFO", "WARNING", or "ERROR"
LOG_LEVEL="${7}"
# Directory containing lariat mapping pipeline files
PIPELINE_DIR="${8}"

#=============================================================================#
#                                  Variables                                  #
#=============================================================================#
output_bam="$OUTPUT_BASE"output.bam
unmapped_fasta="$OUTPUT_BASE"unmapped_reads.fa
fivep_to_reads="$OUTPUT_BASE"fivep_to_reads.sam
heads_fasta="$OUTPUT_BASE"heads.fa
heads_to_genome="$OUTPUT_BASE"heads_to_genome.sam
tails="$OUTPUT_BASE"tails.tsv
putative_lariats="$OUTPUT_BASE"putative_lariats.tsv
failed_fiveps="$OUTPUT_BASE"failed_fivep_alignments.tsv
failed_heads="$OUTPUT_BASE"failed_head_alignments.tsv
failed_lariat="$OUTPUT_BASE"failed_lariat_alignments.tsv
settings_file="$OUTPUT_BASE"settings.json

# Have to use * because bowtie2 index could be small (X.bt2) or large (X.bt2l)
temp_files=(
	$output_bam $output_bam.bai $unmapped_fasta $unmapped_fasta.fai
	$unmapped_fasta.1.bt* $unmapped_fasta.2.bt* $unmapped_fasta.3.bt* $unmapped_fasta.4.bt* 
	$unmapped_fasta.rev.1.bt* $unmapped_fasta.rev.2.bt* $fivep_to_reads $heads_fasta 
	$heads_to_genome $tails $putative_lariats $failed_fiveps $failed_heads $failed_lariat 
	$settings_file
)

#=============================================================================#
#                                  Functions                                  #
#=============================================================================#
end_run () {
	### Classify reads
	Rscript $PIPELINE_DIR/scripts/linear_mapping_wrapper.R -i $output_bam -f $PIPELINE_DIR/scripts/linear_mapping.R -r $REF_DIR -l $SEQ_TYPE -o $OUTPUT_BASE
	
	python -u $PIPELINE_DIR/scripts/classify_nonlinear.py $OUTPUT_BASE $SEQ_TYPE $LOG_FILE $LOG_LEVEL 

	### Summarise results
	### Also create lariat, circ, and temp-switch output files if they don't exist
	python -u $PIPELINE_DIR/scripts/summarise.py $OUTPUT_BASE $LOG_FILE $LOG_LEVEL $SEQ_TYPE

	### Delete the temporary files 
	if ! $KEEP_TEMP; then
		printf "$(date +'%d/%b/%Y %H:%M:%S') | Deleting temporary files...\n" >> $LOG_FILE
		for file in "${temp_files[@]}"; do
			rm -f $file
		done
		if ! $KEEP_CLASSES; then
			rm -f $OUTPUT_BASE"read_classes.tsv.gz"
		fi
	fi 

	### Print completion message
	if [ "${OUTPUT_BASE: -1}" == "/" ];then
		printf "$(date +'%d/%b/%Y %H:%M:%S') | Lariat mapping complete.\n" >> $LOG_FILE
	else
		printf "$(date +'%d/%b/%Y %H:%M:%S') | Lariat mapping complete for "$(echo ${OUTPUT_BASE:0:-1} | sed "s:.*/::")".\n" >> $LOG_FILE
	fi

	### Exit
	exit 0
}

end_run