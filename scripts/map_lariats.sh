#!/bin/bash



#=============================================================================#
#                                  Arguments                                  #
#=============================================================================#
# RNA-seq fastq read file(s)
INPUT_FILES="${1}"
# Type of sequencing data. "single" or "paired"
SEQ_TYPE="${2}"
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
# Output directory 
OUTPUT_BASE="${8}"
# Strand-specificity of the RNA-seq data. "Unstranded, "Forward", or "Reverse"
STRAND="${9}"
# Annotated repeat regions
REPEATS_BED="${10}"
# Run make_track.py after filter_lariats.py. true or false, default false
UCSC_TRACK="${11}"
# Keep read_classes.tsv.gz. true or false, default false
KEEP_CLASSES="${12}"
# Keep the temp files created during the run. true or false, default false
KEEP_TEMP="${13}"
# Number of threads to use
THREADS="${14}"
# Level for python logging. "DEBUG", "INFO", "WARNING", or "ERROR"
LOG_LEVEL="${15}"
# Directory containing lariat mapping pipeline files
PIPELINE_DIR="${16}"



#=============================================================================#
#                                  Variables                                  #
#=============================================================================#
output_bam="$OUTPUT_BASE"output.bam
unmapped_fasta="$OUTPUT_BASE"unmapped_reads.fa
fivep_to_reads="$OUTPUT_BASE"fivep_to_reads.sam
heads_fasta="$OUTPUT_BASE"heads.fa
heads_to_genome="$OUTPUT_BASE"heads_to_genome.sam
# reads_to_fivep="$OUTPUT_BASE"reads_to_fivep.sam
tails="$OUTPUT_BASE"tails.tsv
putative_lariats="$OUTPUT_BASE"putative_lariats.tsv
failed_fiveps="$OUTPUT_BASE"failed_fivep_alignments.tsv
failed_heads="$OUTPUT_BASE"failed_head_alignments.tsv
failed_lariat="$OUTPUT_BASE"failed_lariat_alignments.tsv

if [ "$STRAND" == "Unstranded" ]; then
	hisat2_strand_arg=""
elif [ "$STRAND" == "Forward" ] & [ "$SEQ_TYPE" == "single" ]; then
	hisat2_strand_arg="--rna-strandness F"
elif [ "$STRAND" == "Forward" ] & [ "$SEQ_TYPE" == "paired" ]; then
	hisat2_strand_arg="--rna-strandness FR"
elif [ "$STRAND" == "Reverse" ] & [ "$SEQ_TYPE" == "single" ]; then
	hisat2_strand_arg="--rna-strandness R"
elif [ "$STRAND" == "Reverse" ] & [ "$SEQ_TYPE" == "paired" ]; then
	hisat2_strand_arg="--rna-strandness RF"
fi



#=============================================================================#
#                                    Calls                                    #
#=============================================================================#
### Map filtered reads to genome and keep unmapped reads. Lariat reads crossing the brachpoint will not be able to map to the gene they're from
printf "$(date +'%d/%b/%y %H:%M:%S') | Mapping reads to genome...\n"
if [ "$SEQ_TYPE" == "single" ]; then
	hisat2 --no-softclip -k 1 --max-seeds 20 --pen-noncansplice 0 --n-ceil L,0,0.05 --score-min L,0,-0.24 --bowtie2-dp 1 \
	       $hisat2_strand_arg --threads $THREADS -x $GENOME_INDEX -U $INPUT_FILES \
		| samtools view --bam --with-header --add-flags PAIRED,READ1 \
		| samtools sort --threads $THREADS --verbosity 0 \
		> $output_bam \
		|| exit 1
elif [ "$SEQ_TYPE" == "paired" ]; then
	# Get the two read files from the comma-separated list
	IFS=',' read -r read_one read_two <<< "$INPUT_FILES"
	# Map
	hisat2 --no-softclip -k 1 --max-seeds 20 --pen-noncansplice 0 --n-ceil L,0,0.05 --score-min L,0,-0.24 --bowtie2-dp 1 \
		   $hisat2_strand_arg --threads $THREADS -x $GENOME_INDEX -1 $read_one -2 $read_two \
		| samtools view --bam --with-header \
		| samtools sort --threads $THREADS --verbosity 0 \
		> $output_bam \
		|| exit 1
fi
samtools index $output_bam

unmapped_read_count=$(samtools view --count --require-flags 4 $output_bam)
if [ $unmapped_read_count == 0 ];then
	printf "$(date +'%d/%b/%y %H:%M:%S') | No reads remaining\n"
	python -u $PIPELINE_DIR/scripts/classify_nonlinear.py $OUTPUT_BASE $SEQ_TYPE $LOG_LEVEL \
		|| exit 1
	python -u $PIPELINE_DIR/scripts/summarise.py $OUTPUT_BASE $LOG_LEVEL $SEQ_TYPE \
		|| exit 1
fi

### Create fasta file of unmapped reads 
printf "$(date +'%d/%b/%y %H:%M:%S') | Creating fasta file of unmapped reads...\n"
samtools fasta -N --require-flags 4 -o $unmapped_fasta $output_bam >/dev/null 2>&1 || exit 1
printf "$(date +'%d/%b/%y %H:%M:%S') | Indexing unmapped reads fasta file...\n"
samtools faidx $unmapped_fasta || exit 1

### Build a bowtie2 index of the unmapped reads
printf "$(date +'%d/%b/%y %H:%M:%S') | Building bowtie2 index of unmapped fasta...\n"
bowtie2-build --quiet --threads $THREADS $unmapped_fasta $unmapped_fasta || exit 1 

## Align fasta file of all 5' splice sites (first 20nts of introns) to unmapped reads index
# We need to order the output SAM by reference (the read id, in this case) for the following filtering process
printf "$(date +'%d/%b/%y %H:%M:%S') | Mapping 5' splice sites to reads...\n"
bowtie2 --end-to-end --sensitive --no-unal -f -k 10000 --score-min C,0,0 --threads $THREADS -x $unmapped_fasta -U $FIVEP_FASTA \
	| samtools sort --threads $THREADS --verbosity 0 --output-fmt SAM -M \
	| samtools view \
	> $fivep_to_reads \
	|| exit 1 

# # ## Align unmapped reads to index of all 5' splice sites (first 20nts of introns)
# printf "$(date +'%d/%b/%y %H:%M:%S') | Mapping 5' splice sites to reads...\n"
# bowtie2 --local -k 1000 -L 20 -i C,1,0 --ma 1 --mp 1,1 --np 1 --rdg 1,1 --rfg 1,1 --score-min C,20,0 \
# 		--no-unal --no-head --threads $THREADS -x $FIVEP_INDEX -f -U $unmapped_fasta \
# 	> $reads_to_fivep 

## Extract reads with a mapped 5' splice site and trim it off
printf "$(date +'%d/%b/%y %H:%M:%S') | Finding 5' read alignments and trimming reads...\n"
python -u $PIPELINE_DIR/scripts/filter_fivep_aligns.py $OUTPUT_BASE $LOG_LEVEL $GENOME_FASTA $FIVEP_FASTA $STRAND $THREADS \
	|| exit 1 

### Map read heads to genome
printf "$(date +'%d/%b/%y %H:%M:%S') | Mapping heads to genome...\n"
hisat2 --no-softclip --no-spliced-alignment --very-sensitive -k 100 \
	   --no-unal --threads $THREADS -f -x $GENOME_INDEX -U $heads_fasta \
	| samtools sort --threads $THREADS --verbosity 0 --output-fmt SAM -n \
	| samtools view \
	> $heads_to_genome \
	|| exit 1

### Filter head alignments
printf "$(date +'%d/%b/%y %H:%M:%S') | Analyzing head alignments and outputting lariat table...\n"
python -u $PIPELINE_DIR/scripts/filter_head_aligns.py $THREADS $INTRONS_TSV $GENOME_FASTA $OUTPUT_BASE $LOG_LEVEL \
	|| exit 1 

### Filter lariat mappings and choose 1 for each read
printf "$(date +'%d/%b/%y %H:%M:%S') | Filtering putative lariat alignments...\n"
python -u $PIPELINE_DIR/scripts/filter_lariats.py $OUTPUT_BASE $LOG_LEVEL $SEQ_TYPE $GENOME_FASTA $REPEATS_BED \
	|| exit 1 

### Classify reads
python -u $PIPELINE_DIR/scripts/classify_nonlinear.py $OUTPUT_BASE $SEQ_TYPE $LOG_LEVEL \
	|| exit 1
python -u $PIPELINE_DIR/scripts/summarise.py $OUTPUT_BASE $LOG_LEVEL $SEQ_TYPE \
	|| exit 1

### Make a custom track BED file of identified lariats 
if $UCSC_TRACK; then
	printf "$(date +'%d/%b/%y %H:%M:%S') | Making UCSC Genome Browser track...\n"
	python -u $PIPELINE_DIR/scripts/make_track.py $OUTPUT_BASE $LOG_LEVEL \
		|| exit 1
fi



wait
### Delete the temporary files 
if ! $KEEP_TEMP; then
	printf "$(date +'%d/%b/%y %H:%M:%S') | Deleting temporary files...\n"
	rm $output_bam
	rm $output_bam.bai
	rm $unmapped_fasta 
	rm $unmapped_fasta.fai
	# Have to use * because bowtie2 index could be small (X.bt2) or large (X.bt2l)
	rm $unmapped_fasta.1.bt*
	rm $unmapped_fasta.2.bt*
	rm $unmapped_fasta.3.bt*
	rm $unmapped_fasta.4.bt*
	rm $unmapped_fasta.rev.1.bt*
	rm $unmapped_fasta.rev.2.bt*
	rm $fivep_to_reads
	rm $heads_fasta
	rm $heads_to_genome
	rm $tails
	rm $putative_lariats
	rm $failed_fiveps
	rm $failed_heads
	rm $failed_lariat
	rm $OUTPUT_BASE"settings.json"
	if ! $KEEP_CLASSES; then
		rm $OUTPUT_BASE"read_classes.tsv.gz"
	fi
fi
fi 

if [ "${OUTPUT_BASE: -1}" == "/" ];then
	printf "$(date +'%d/%b/%y %H:%M:%S') | Lariat mapping complete.\n"
else
	printf "$(date +'%d/%b/%y %H:%M:%S') | Lariat mapping complete for "$(echo ${OUTPUT_BASE:0:-1} | sed "s:.*/::")".\n"
fi