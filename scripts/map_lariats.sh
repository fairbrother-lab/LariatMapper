#!/bin/bash



#=============================================================================#
#                                  Arguments                                  #
#=============================================================================#
# Output directory 
OUTPUT_BASE=$1
# Number of threads to use
THREADS=$2
# hisat2 index of reference genome
GENOME_INDEX=$3
# Reference genome FASTA file
GENOME_FASTA=$4
# Gene annotation of the reference genome in GTF or GFF format
# ANNO_FILE=$5
# # FASTA file of 5' splice sites (first 20nts of all introns)
FIVEP_FASTA=$5
# Annotated introns
INTRONS_BED=$6
# Annotated repeat regions
REPEATS_BED=$7
# Keep the intermediate files created during the run ("True" or "False", default "False")
KEEP_INTERMEDIATES=$8
# Run make_lariat_track.py after filter_lariats.py
UCSC_TRACK="${9}"
# Directory containing lariat mapping pipeline files
PIPELINE_DIR="${10}"
# Directory containing lariat mapping pipeline files
LOG_LEVEL="${11}"
# RNA-seq fastq read file(s)
if [ $# == 12 ]; then
	READ_FILE="${12}"
	single_end=true
else
	READ_ONE="${12}"
	READ_TWO="${13}"
	single_end=false
fi

#=============================================================================#
#                                  Variables                                  #
#=============================================================================#
output_bam="$OUTPUT_BASE"output.bam
mapped_bam="$OUTPUT_BASE"mapped_reads.bam
unmapped_bam="$OUTPUT_BASE"unmapped_reads.bam
run_data="$OUTPUT_BASE"read_counts.tsv
unmapped_fasta="$OUTPUT_BASE"unmapped_reads.fa
fivep_to_reads="$OUTPUT_BASE"fivep_to_reads.sam
fivep_trimmed_reads="$OUTPUT_BASE"fivep_mapped_reads_trimmed.fa
trimmed_reads_to_genome="$OUTPUT_BASE"trimmed_reads_to_genome.sam
reads_to_fivep="$OUTPUT_BASE"reads_to_fivep.sam
FIVEP_INDEX=/datasets2/lariat_mapping/testing/output/test_ref/fivep_sites

fivep_info_table="$OUTPUT_BASE"fivep_info_table.tsv
trimmed_info_table="$OUTPUT_BASE"trimmed_info_table.tsv
failed_fivep="$OUTPUT_BASE"failed_fivep_alignments.tsv
failed_trimmed="$OUTPUT_BASE"failed_trimmed_alignments.tsv
failed_lariat="$OUTPUT_BASE"failed_lariat_alignments.tsv





#=============================================================================#
#                                    Calls                                    #
#=============================================================================#
### Map filtered reads to genome and keep unmapped reads. Lariat reads crossing the brachpoint will not be able to map to the gene they're from
printf "$(date +'%d/%b/%y %H:%M:%S') | Mapping reads to genome...\n"
if $single_end; then
	hisat2 --no-softclip -k 1 --max-seeds 20 --pen-noncansplice 0 --n-ceil L,0,0.05 --score-min L,0,-0.24 --bowtie2-dp 1 \
	       --threads $THREADS -x $GENOME_INDEX -U $READ_FILE \
		| samtools view --bam --with-header --add-flags PAIRED,READ1 \
		> $output_bam \
		|| exit 1
else
	hisat2 --no-softclip -k 1 --max-seeds 20 --pen-noncansplice 0 --n-ceil L,0,0.05 --score-min L,0,-0.24 --bowtie2-dp 1 \
		   --threads $THREADS -x $GENOME_INDEX -1 $READ_ONE -2 $READ_TWO \
		| samtools view --bam --with-header \
		> $output_bam \
		|| exit 1
fi
printf "$(date +'%d/%b/%y %H:%M:%S') | Extracting unmapped reads...\n"
samtools view --bam --with-header --exclude-flags 4 $output_bam > $mapped_bam
samtools view --bam --with-header --require-flags 4 $output_bam > $unmapped_bam
mapped_read_count=$(samtools view --count $mapped_bam)
unmapped_read_count=$(samtools view --count $unmapped_bam)
if ! $single_end; then
	mapped_read_count=$((mapped_read_count/2))
	unmapped_read_count=$((unmapped_read_count/2))
fi
echo -e "linear_mapped\t$mapped_read_count" > $run_data
echo -e "linear_unmapped\t$unmapped_read_count" >> $run_data

### Check if 0 reads were left unmapped
if [ $unmapped_read_count == 0 ];then
	printf "$(date +'%d/%b/%y %H:%M:%S') | No reads remaining"
	exit
fi

### Create fasta file of unmapped reads 
printf "$(date +'%d/%b/%y %H:%M:%S') | Creating fasta file of unmapped reads...\n"
samtools fasta -N -o $unmapped_fasta $unmapped_bam >/dev/null 2>&1 || exit 1
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
# scalene --html --outfile "$OUTPUT_BASE"filter_fivep_profile.html $PIPELINE_DIR/scripts/filter_fivep_alignments.py $THREADS $GENOME_FASTA $FIVEP_FASTA $OUTPUT_BASE \
python -u $PIPELINE_DIR/scripts/filter_fivep_alignments.py $THREADS $GENOME_FASTA $FIVEP_FASTA $OUTPUT_BASE $LOG_LEVEL \
	|| exit 1 





### Map 5' trimmed reads to genome
printf "$(date +'%d/%b/%y %H:%M:%S') | Mapping 5' trimmed reads to genome...\n"
hisat2 --no-softclip --no-spliced-alignment --very-sensitive -k 100 \
	   --no-unal --threads $THREADS -f -x $GENOME_INDEX -U $fivep_trimmed_reads \
	| samtools sort --threads $THREADS --verbosity 0 --output-fmt SAM -n \
	| samtools view \
	> $trimmed_reads_to_genome \
	|| exit 1

### Filter trimmed alignments
printf "$(date +'%d/%b/%y %H:%M:%S') | Analyzing trimmed alignments and outputting lariat table...\n"
# scalene --html --outfile "$OUTPUT_BASE"filter_trim_profile.html $PIPELINE_DIR/scripts/filter_trimmed_alignments.py $THREADS $INTRONS_BED $GENOME_FASTA $OUTPUT_BASE \
python -u $PIPELINE_DIR/scripts/filter_trimmed_alignments.py $THREADS $INTRONS_BED $GENOME_FASTA $OUTPUT_BASE $LOG_LEVEL \
	|| exit 1 

### Filter lariat mappings and choose 1 for each read
printf "$(date +'%d/%b/%y %H:%M:%S') | Filtering putative lariat alignments...\n"
python -u $PIPELINE_DIR/scripts/filter_lariats.py $GENOME_FASTA $REPEATS_BED $OUTPUT_BASE $LOG_LEVEL \
	|| exit 1 

### Make a custom track BED file of identified lariats 
if $UCSC_TRACK; then
	printf "$(date +'%d/%b/%y %H:%M:%S') | Making UCSC Genome Browser track...\n"
	python -u $PIPELINE_DIR/scripts/make_lariat_track.py $OUTPUT_BASE $LOG_LEVEL \
		|| exit 1
fi

### Classify reads
python -u $PIPELINE_DIR/scripts/classify_linear_reads.py $EXONS_TSV $INTRONS_TSV $OUTPUT_BASE $LOG_LEVEL \
	|| exit 1
python -u $PIPELINE_DIR/scripts/classify_nonlinear_reads.py $OUTPUT_BASE $LOG_LEVEL \
	|| exit 1



wait
### Delete the intermediate files 
if ! $KEEP_INTERMEDIATES; then
	printf "$(date +'%d/%b/%y %H:%M:%S') | Deleting intermediate files...\n"
	rm $output_bam
	rm $mapped_bam
	rm $unmapped_bam
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
	rm $fivep_trimmed_reads
	rm $trimmed_reads_to_genome
	rm $fivep_info_table
	rm $trimmed_info_table
	rm $failed_fivep
	rm $failed_trimmed
	rm $failed_lariat
fi 

if [ "${OUTPUT_BASE: -1}" == "/" ];then
	printf "$(date +'%d/%b/%y %H:%M:%S') | Lariat mapping complete.\n"
else
	printf "$(date +'%d/%b/%y %H:%M:%S') | Lariat mapping complete for "$(echo ${OUTPUT_BASE:0:-1} | sed "s:.*/::")".\n"
fi