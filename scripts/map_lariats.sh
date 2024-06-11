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
ANNO_FILE=$5
# FASTA file of 5' splice sites (first 20nts of all introns)
FIVEP_FASTA=$6
# Annotated introns
INTRONS_BED=$7
# Annotated repeat regions
REPEATS_BED=$8
# Keep the intermediate files created during the run ("True" or "False", default "False")
KEEP_INTERMEDIATES=$9
# Run make_lariat_track.py after filter_lariats.py
UCSC_TRACK="${10}"
# Directory containing lariat mapping pipeline files
PIPELINE_DIR="${11}"
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
output_bam="$OUTPUT_BASE"mapped_reads.bam
unmapped_bam="$OUTPUT_BASE"unmapped_reads.bam
run_data="$OUTPUT_BASE"read_counts.tsv
unmapped_fasta="$OUTPUT_BASE"unmapped_reads.fa
fivep_to_reads="$OUTPUT_BASE"fivep_to_reads.sam
fivep_trimmed_reads="$OUTPUT_BASE"fivep_mapped_reads_trimmed.fa
trimmed_reads_to_genome="$OUTPUT_BASE"trimmed_reads_to_genome.sam



#=============================================================================#
#                                    Calls                                    #
#=============================================================================#
### Map filtered reads to genome and keep unmapped reads. Lariat reads crossing the brachpoint will not be able to map to the gene they're from
printf "$(date +'%m/%d/%y - %H:%M:%S') | Mapping reads and extracting unmapped reads...\n"
if $single_end; then
	hisat2 --no-softclip --max-seeds 20 --bowtie2-dp 1 --pen-noncansplice 0 -k 1 --n-ceil L,0,0.05 --score-min L,0,-0.24 \
	       --threads $THREADS -x $GENOME_INDEX -U $READ_FILE \
		| samtools view --bam --with-header --add-flags PAIRED,READ1 \
		> $output_bam \
		|| exit 1 
else
	hisat2 --no-softclip --max-seeds 20 --bowtie2-dp 1 --pen-noncansplice 0 -k 1 --n-ceil L,0,0.05 --score-min L,0,-0.24 \
		   --threads $THREADS -x $GENOME_INDEX -1 $READ_ONE -2 $READ_TWO \
		| samtools view --bam --with-header \
		> $output_bam \
		|| exit 1 
fi
samtools view --bam --with-header --require-flags 4 $output_bam > $unmapped_bam
mapped_read_count=$(samtools view --count --exclude-flags 4 $output_bam)
unmapped_read_count=$(samtools view --count $unmapped_bam)
if ! $single_end; then
	mapped_read_count=$((mapped_read_count/2))
	unmapped_read_count=$((unmapped_read_count/2))
fi
run_data="$OUTPUT_BASE"read_counts.tsv
echo -e "linear_mapped\t$mapped_read_count" > $run_data
echo -e "linear_unmapped\t$unmapped_read_count" >> $run_data

# Save linear-mapped read alignments for later classification
printf "$(date +'%m/%d/%y - %H:%M:%S') | Writing mapped reads to BED file...\n"
bedtools bamtobed -split -i $output_bam | gzip > "$OUTPUT_BASE"mapped_reads.bed.gz 

### Create fasta file of unmapped reads 
printf "$(date +'%m/%d/%y - %H:%M:%S') | Creating fasta file of unmapped reads...\n"
samtools fasta -N -o $unmapped_fasta $unmapped_bam >/dev/null || exit 1 
printf "$(date +'%m/%d/%y - %H:%M:%S') | Indexing unmapped reads fasta file...\n"
samtools faidx $unmapped_fasta

### Build a bowtie2 index of the unmapped reads
printf "$(date +'%m/%d/%y - %H:%M:%S') | Building bowtie2 index of unmapped fasta...\n"
bowtie2-build --threads $THREADS $unmapped_fasta $unmapped_fasta >/dev/null \
	|| exit 1 

## Align unmapped reads to fasta file of all 5' splice sites (first 20nts of introns)
printf "$(date +'%m/%d/%y - %H:%M:%S') | Mapping 5' splice sites to reads...\n"
bowtie2 --end-to-end --sensitive --no-unal -f -k 10000 --score-min C,0,0 --threads $THREADS -x $unmapped_fasta -U $FIVEP_FASTA \
	| samtools view \
	> $fivep_to_reads \
	|| exit 1

## Extract reads with a mapped 5' splice site and trim it off
printf "$(date +'%m/%d/%y - %H:%M:%S') | Finding 5' read alignments and trimming reads...\n"
python $PIPELINE_DIR/scripts/filter_fivep_alignments.py $THREADS $GENOME_FASTA $OUTPUT_BASE \
	|| exit 1

### Map 5' trimmed reads to genome
printf "$(date +'%m/%d/%y - %H:%M:%S') | Mapping 5' trimmed reads to genome...\n"
hisat2 --no-softclip --no-spliced-alignment --very-sensitive -k 10 --no-unal --threads $THREADS -f -x $GENOME_INDEX -U $fivep_trimmed_reads \
	> $trimmed_reads_to_genome \
	|| exit 1

### Filter trimmed alignments
printf "$(date +'%m/%d/%y - %H:%M:%S') | Analyzing trimmed alignments and outputting lariat table...\n"
# scalene --html --outfile "$OUTPUT_BASE"filter_trim_profile.html $PIPELINE_DIR/scripts/filter_trimmed_alignments.py $THREADS $ANNO_FILE $INTRONS_BED $GENOME_FASTA $OUTPUT_BASE \
python -u $PIPELINE_DIR/scripts/filter_trimmed_alignments.py $THREADS $ANNO_FILE $INTRONS_BED $GENOME_FASTA $OUTPUT_BASE \
	|| exit 1

### Filter lariat mappings and choose 1 for each read
printf "$(date +'%m/%d/%y - %H:%M:%S') | Filtering putative lariat alignments...\n"
python -u $PIPELINE_DIR/scripts/filter_lariats.py $GENOME_FASTA $REPEATS_BED $OUTPUT_BASE \
	|| exit 1

### Make a custom track BED file of identified lariats 
if $UCSC_TRACK; then
	printf "$(date +'%m/%d/%y - %H:%M:%S') | Making UCSC Genome Browser track...\n"
	python -u $PIPELINE_DIR/scripts/make_lariat_track.py $OUTPUT_BASE \
		|| exit 1
fi


wait
### Delete the intermediate files 
if ! $KEEP_INTERMEDIATES; then
	printf "$(date +'%m/%d/%y - %H:%M:%S') | Deleting intermediate files...\n"
	rm $output_bam
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
	rm $fivep_info_table
	rm $trimmed_reads_to_genome
fi 

if [ "${OUTPUT_BASE: -1}" == "/" ];then
	printf "$(date +'%m/%d/%y - %H:%M:%S') | Lariat mapping complete.\n"
else
	printf "$(date +'%m/%d/%y - %H:%M:%S') | Lariat mapping complete for "$(echo ${OUTPUT_BASE:0:-1} | sed "s:.*/::")".\n"
fi
