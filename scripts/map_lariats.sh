#!/bin/bash



#=============================================================================#
#                                  Arguments                                  #
#=============================================================================#
# RNA-seq fastq read file(s)
INPUT_FILES="${1}"
# Reference directory
REF_DIR="${2}"
# hisat2 index of reference genome
GENOME_INDEX="${3}"
# Reference genome FASTA file
GENOME_FASTA="${4}"
# FASTA file of 5' splice sites (first 20nts of all introns)
FIVEP_FASTA="${5}"
# Annotated exons
EXONS_TSV="${6}"
# Annotated introns
INTRONS_TSV="${7}"
# Strand-specificity of the RNA-seq data. "Unstranded, "Forward", or "Reverse"
STRAND="${8}"
# Parameters of the template-switching filter
TEMP_SWITCH_FILTER="${9}"
# Annotated repeat regions
REPEATS_BED="${10}"
# A comma-seperated list of PWM files or a single model file
PWM_FILES="${11}"
MODEL_FILE="${12}"
# Run make_track.py after filter_lariats.py. true or false, default false
UCSC_TRACK="${13}"
# Keep output.bam
KEEP_BAM="${14}"
# Keep read_classes.tsv.gz. true or false, default false
KEEP_CLASSES="${15}"
# Keep the temp files created during the run. true or false, default false
KEEP_TEMP="${16}"
# Number of threads to use
THREADS="${17}"
# Type of sequencing data. "single" or "paired"
SEQ_TYPE="${18}"
# Output directory 
OUTPUT_BASE="${19}"
# Level for python logging. "DEBUG", "INFO", "WARNING", or "ERROR"
LOG_LEVEL="${20}"
# Directory containing lariat mapping pipeline files
PIPELINE_DIR="${21}"



#=============================================================================#
#                                  Variables                                  #
#=============================================================================#
output_sam="$OUTPUT_BASE"output.sam
output_bam="$OUTPUT_BASE"output.bam
unmapped_fasta="$OUTPUT_BASE"unmapped_reads.fa
fivep_to_reads="$OUTPUT_BASE"fivep_to_reads.sam
heads_fasta="$OUTPUT_BASE"heads.fa
heads_to_genome="$OUTPUT_BASE"heads_to_genome.sam

if [ "$STRAND" == "Unstranded" ]; then
	hisat2_strand_arg=""
elif [ "$STRAND" == "First" ] && [ "$SEQ_TYPE" == "single" ]; then
	hisat2_strand_arg="--rna-strandness F"
elif [ "$STRAND" == "First" ] && [ "$SEQ_TYPE" == "paired" ]; then
	hisat2_strand_arg="--rna-strandness FR"
elif [ "$STRAND" == "Second" ] && [ "$SEQ_TYPE" == "single" ]; then
	hisat2_strand_arg="--rna-strandness R"
elif [ "$STRAND" == "Second" ] && [ "$SEQ_TYPE" == "paired" ]; then
	hisat2_strand_arg="--rna-strandness RF"
fi

# List of temporary files to delete at the end of the run if KEEP_TEMP is false
temp_files=(
	$output_sam $output_bam $output_bam.bai $unmapped_fasta $unmapped_fasta.fai
	# Have to use * because bowtie2 index could be small (X.bt2) or large (X.bt2l)
	$unmapped_fasta.1.bt* $unmapped_fasta.2.bt* $unmapped_fasta.3.bt* $unmapped_fasta.4.bt* 
	$unmapped_fasta.rev.1.bt* $unmapped_fasta.rev.2.bt* $fivep_to_reads $fivep_to_reads.tmp 
	$heads_fasta $heads_to_genome $heads_to_genome.tmp "$OUTPUT_BASE"tails.tsv "$OUTPUT_BASE"putative_lariats.tsv
	$OUTPUT_BASE"settings.json" $OUTPUT_BASE"read_classes.tsv.gz" $OUTPUT_BASE"output.bam_summary_count.tsv"
	"$OUTPUT_BASE"failed_fivep_alignments.tsv "$OUTPUT_BASE"failed_head_alignments.tsv "$OUTPUT_BASE"failed_lariat_alignments.tsv
)



#=============================================================================#
#                                  Functions                                  #
#=============================================================================#
print_message(){
	text="$1"
	out_msg="$(date +'%d/%b/%Y %H:%M:%S') | $text \n"
	printf "$out_msg"
}

delete_temp_files(){
	print_message "Deleting temporary files..."
	for file in "${temp_files[@]}"; do
		if [ $file == $OUTPUT_BASE"read_classes.tsv.gz" ] & $KEEP_CLASSES; then
			continue
		elif [ $file == $OUTPUT_BASE"output.bam" ] & $KEEP_BAM; then
			continue
		else
			rm $file
		fi
	done
}

end_run() {
	### Classify reads
	print_message "Classifying linearly-aligned reads..."
	Rscript $PIPELINE_DIR/scripts/linear_mapping_wrapper.R -i $output_bam -f $PIPELINE_DIR/scripts/linear_mapping.R -r $REF_DIR -l $SEQ_TYPE -s $STRAND -o $OUTPUT_BASE
	check_exitcode
	print_message "Classifying non-linearly-aligned reads..."
	python -u $PIPELINE_DIR/scripts/classify_nonlinear.py $OUTPUT_BASE $SEQ_TYPE $LOG_LEVEL 
	check_exitcode

	### Summarise results
	### Also create lariat, circ, and temp-switch output files if they don't exist
	print_message "Summarising results..."
	python -u $PIPELINE_DIR/scripts/summarise.py $OUTPUT_BASE $LOG_LEVEL $SEQ_TYPE
	check_exitcode

	### Create a subdir named "plots" and make summary plots in it
	print_message "Building plots..."
	Rscript $PIPELINE_DIR/scripts/summary_plots.R -o $OUTPUT_BASE
	check_exitcode

	### Delete the temporary files 
	if ! $KEEP_TEMP; then
		delete_temp_files
	fi
	check_exitcode

	### Print completion message
	if [ "${OUTPUT_BASE: -1}" == "/" ];then
		print_message "Lariat mapping complete."
	else
		print_message "Lariat mapping complete for "$(echo ${OUTPUT_BASE:0:-1} | sed "s:.*/::")"."
	fi

	### Exit
	exit 0
}

check_exitcode() {
	exit_code=$?
	### Exit code 0 = success
	if [ $exit_code -eq 0 ];then
		return 0
	### Exit code 4 = No reads or alignments remain for processing part-way through the pipeline
	### End the run early
	elif [ $exit_code -eq 4 ];then
		print_message "No reads remaining."
		end_run
	else
		exit $exit_code
	fi
}



#=============================================================================#
#                                    Calls                                    #
#=============================================================================#
### Map filtered reads to genome and keep unmapped reads. Lariat reads crossing the brachpoint will not be able to map to the gene they're from
print_message "Mapping reads to genome..."
if [ "$SEQ_TYPE" == "single" ]; then
	# Map
	hisat2 --no-softclip -k 1 --max-seeds 20 --pen-noncansplice 0 --n-ceil L,0,0.05 --score-min L,0,-0.24 --bowtie2-dp 1 \
	       $hisat2_strand_arg --threads $THREADS -x $GENOME_INDEX -U $INPUT_FILES \
		> $output_sam
	check_exitcode
	samtools view --bam --with-header --add-flags PAIRED,READ1 $output_sam \
		| samtools sort --threads $THREADS --verbosity 0 \
		> $output_bam 
	check_exitcode

elif [ "$SEQ_TYPE" == "paired" ]; then
	# Get the two read files from the comma-separated list
	IFS=',' read -r read_one read_two <<< "$INPUT_FILES"
	# Map
	hisat2 --no-softclip -k 1 --max-seeds 20 --pen-noncansplice 0 --n-ceil L,0,0.05 --score-min L,0,-0.24 --bowtie2-dp 1 \
			$hisat2_strand_arg --threads $THREADS -x $GENOME_INDEX -1 $read_one -2 $read_two \
		> $output_sam
	check_exitcode
	samtools view --bam --with-header $output_sam \
		| samtools sort --threads $THREADS --verbosity 0 \
		> $output_bam
	check_exitcode

fi

samtools index $output_bam
check_exitcode

unmapped_read_count=$(samtools view --count --require-flags 4 $output_bam)
if [ $unmapped_read_count == 0 ];then
	print_message "All reads mapped linearly to genome."
	end_run
fi


### Create fasta file of unmapped reads 
print_message "Creating fasta file of unmapped reads..."
samtools fasta -N --require-flags 4 -o $unmapped_fasta $output_bam >/dev/null 2>&1
check_exitcode
print_message "Indexing unmapped reads fasta file..."
samtools faidx $unmapped_fasta 
check_exitcode

### Build a bowtie2 index of the unmapped reads
print_message "Building bowtie2 index of unmapped fasta..."
bowtie2-build --quiet --threads $THREADS $unmapped_fasta $unmapped_fasta
check_exitcode 


## Align fasta file of all 5' splice sites (first 20nts of introns) to unmapped reads index
# We need to order the output SAM by read_id (which is the reference in this case)
# so we can process all the alignments to each read together in the following filtering 
print_message "Mapping 5' splice sites to reads..."
bowtie2 --end-to-end --sensitive --no-unal -f -k 10000 --score-min C,0,0 --threads $THREADS -x $unmapped_fasta -U $FIVEP_FASTA \
	> $fivep_to_reads.tmp
check_exitcode
samtools sort --threads $THREADS --verbosity 0 --output-fmt SAM $fivep_to_reads.tmp \
	| samtools view \
	> $fivep_to_reads
check_exitcode


## Extract reads with a mapped 5' splice site and trim it off
print_message "Finding 5' read alignments and trimming reads..."
# scalene --html --outfile "$OUTPUT_BASE"filter_fivep_aligns.html \
# 	$PIPELINE_DIR/scripts/filter_fivep_aligns.py $OUTPUT_BASE $LOG_LEVEL $GENOME_FASTA $FIVEP_FASTA $STRAND $THREADS
python -u $PIPELINE_DIR/scripts/filter_fivep_aligns.py $OUTPUT_BASE $LOG_LEVEL $GENOME_FASTA $FIVEP_FASTA $STRAND $THREADS
check_exitcode


### Map read heads to genome
print_message "Mapping heads to genome..."
hisat2 --no-softclip --no-spliced-alignment --very-sensitive -k 100 \
	   --no-unal --threads $THREADS -f -x $GENOME_INDEX -U $heads_fasta \
	> $heads_to_genome.tmp
check_exitcode
# We need to order the output SAM by read_id (which is the reference in this case)
# so we can process all the alignments to each read together in the following filtering 
samtools sort --threads $THREADS --verbosity 0 --output-fmt SAM -n $heads_to_genome.tmp \
	> $heads_to_genome
check_exitcode


### Filter head alignments
print_message "Analyzing head alignments and outputting lariat table..."
# scalene --html --outfile "$OUTPUT_BASE"filter_head_aligns.html \
# 	$PIPELINE_DIR/scripts/filter_head_aligns.py $THREADS $EXONS_TSV $INTRONS_TSV $GENOME_FASTA $TEMP_SWITCH_FILTER $OUTPUT_BASE $LOG_LEVEL 
python -u $PIPELINE_DIR/scripts/filter_head_aligns.py $THREADS $EXONS_TSV $INTRONS_TSV $GENOME_FASTA $TEMP_SWITCH_FILTER $OUTPUT_BASE $LOG_LEVEL 
check_exitcode


### Filter lariat mappings and choose 1 for each read
print_message "Filtering putative lariat alignments..."
python -u $PIPELINE_DIR/scripts/filter_lariats.py $OUTPUT_BASE $LOG_LEVEL $SEQ_TYPE $GENOME_FASTA $REPEATS_BED
check_exitcode


### Correct branchpoint positions, if extra files are provided
if ! [ "$PWM_FILES" == "" ]; then
	print_message "Correcting branchpoints with PWMs..."
	Rscript $PIPELINE_DIR/scripts/bp_correction_wrapper.R \
		--input "$OUTPUT_BASE"lariat_reads.tsv \
		--ref_fasta $GENOME_FASTA \
		--file $PIPELINE_DIR/scripts/bp_correction.R \
		--method PWM \
		--PWM_path $PWM_FILES \
		--log_level $LOG_LEVEL \
		--output_base $OUTPUT_BASE 
	check_exitcode
elif ! [ "$MODEL_FILE" == "" ]; then
	print_message "Correcting branchpoints with DeepEnsemble..."
	Rscript $PIPELINE_DIR/scripts/bp_correction_wrapper.R \
		--input "$OUTPUT_BASE"lariat_reads.tsv \
		--ref_fasta $GENOME_FASTA \
		--file $PIPELINE_DIR/scripts/bp_correction.R \
		--method Model-based \
		--model_path $MODEL_FILE \
		--log_level $LOG_LEVEL \
		--output_base $OUTPUT_BASE
	check_exitcode
fi


### Make a custom track BED file of identified lariats 
if $UCSC_TRACK; then
	print_message "Making UCSC Genome Browser track..."
	python -u $PIPELINE_DIR/scripts/make_track.py $OUTPUT_BASE $LOG_LEVEL 
	check_exitcode
fi



wait
end_run