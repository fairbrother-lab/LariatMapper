#!/bin/bash

#=============================================================================#
#                                  Arguments                                  #
#=============================================================================#

# RNA-seq fastq read file
READ_FILE=$1
# Output directory with sample name prefix
OUTPUT_BASE=$2
# Number of CPUs to use
CPUS=$3
# Genome Bowtie2 index base name
GENOME_BOWTIE2_INDEX=$4
# Reference genome FASTA file
GENOME_FASTA=$5
# GTF file containing gene annotatin for mapping genome
GTF_FILE=$6
# FASTA file of 5' splice sites (first 20nts of all introns)
FIVEP_FASTA=$7
# Bowtie2 index of 3' splice sites genome (last 250nts of all introns)
THREEP_BOWTIE2_INDEX=$8
# 3' splice site lengths
THREEP_LENGTHS=$9
# Python scripts to filter 5' and 3' alignments
FIVEP_SCRIPT="${10}"
THREEP_SCRIPT="${11}"
# Bash script with the timestamp function, which wraps a print statement with a timestamp for logging
TIMESTAMP="${12}"


#=============================================================================#
#                                    Setup                                    #
#=============================================================================#

source $TIMESTAMP

printf "Run parameters:\nREAD_FILE: $READ_FILE \nOUTPUT_BASE: $OUTPUT_BASE \nCPUS: $CPUS \nGENOME_BOWTIE2_INDEX: $GENOME_BOWTIE2_INDEX \nGENOME_FASTA: $GENOME_FASTA \nGTF_FILE: $GTF_FILE \nFIVEP_FASTA: $FIVEP_FASTA \nTHREEP_BOWTIE2_INDEX: $THREEP_BOWTIE2_INDEX \nTHREEP_LENGTHS: $THREEP_LENGTHS\n"


#=============================================================================#
#                                    Calls                                    #
#=============================================================================#
### Map filtered reads to genome and keep unmapped reads. Lariat reads crossing the brachpoint will not be able to map to the gene they're from
timestamp "(1/7) Mapping reads and extracting unmapped reads" -1
output_bam=$OUTPUT_BASE"_mapped_reads.bam"
unmapped_bam=$OUTPUT_BASE"_unmapped_reads.bam"
bowtie2 --end-to-end --sensitive --score-min L,0,-0.24 -k 1 --n-ceil L,0,0.05 --threads $CPUS -x $GENOME_BOWTIE2_INDEX -U $READ_FILE \
	| samtools view --bam --with-header > $output_bam
samtools view --bam --with-header --require-flags 4 $output_bam > $unmapped_bam
mapped_read_count=$(samtools view --count --exclude-flags 4 $output_bam)
echo "total_reads=$mapped_read_count" > $OUTPUT_BASE"_total_reads.txt"
# rm $output_bam

### Create fasta file of unmapped reads 
timestamp "(2/7) Creating fasta file of unmapped reads" -1
unmapped_fasta=$OUTPUT_BASE"_unmapped_reads.fa"
samtools fasta $unmapped_bam > $unmapped_fasta
samtools faidx $unmapped_fasta
# rm $unmapped_bam

### Build a bowtie index of the unmapped reads
timestamp "(3/7) Building bowtie index of unmapped fasta" -1
bowtie2-build --large-index --threads $CPUS $unmapped_fasta $unmapped_fasta > /dev/null

### Align unmapped reads to fasta file of all 5' splice sites (first 20nts of introns)
timestamp "(4/7) Mapping 5' splice sites to reads" -1
fivep_to_reads=$OUTPUT_BASE"_fivep_to_reads.sam"
bowtie2 --end-to-end --sensitive --no-unal -f --k 10000 --threads $CPUS -x $unmapped_fasta -U $FIVEP_FASTA \
	| samtools view > $fivep_to_reads

### Extract reads with a mapped 5' splice site and trim it off
timestamp "(5/7) Finding 5' read alignments and trimming reads" -1
fivep_info_table=$OUTPUT_BASE"_fivep_info_table.txt"
fivep_trimmed_reads=$OUTPUT_BASE"_fivep_mapped_reads_trimmed.fa"
python $FIVEP_SCRIPT $unmapped_fasta $fivep_to_reads $fivep_trimmed_reads $fivep_info_table

### Map 5' trimmed reads to 3' sites (last 250nts of introns)
timestamp "(6/7) Mapping 5' trimmed reads to 3' sites" -1
trimmed_reads_to_threep=$OUTPUT_BASE"_fivep_reads_trimmed_mapped_to_threep.sam"
bowtie2 --end-to-end --sensitive -k 10 --no-unal --threads $CPUS -f -x $THREEP_BOWTIE2_INDEX -U $fivep_trimmed_reads \
	| samtools view > $trimmed_reads_to_threep

### Filter 3' splice site alignments and output info table, including the branchpoint site
timestamp "(7/7) Filtering 3' alignments and outputting final table" -1
python $THREEP_SCRIPT $trimmed_reads_to_threep $THREEP_LENGTHS $fivep_info_table $GTF_FILE $GENOME_FASTA $OUTPUT_BASE

### Delete all intermediate/uneeded files that were created throughout this process
wait
rm $output_bam
rm $unmapped_bam
rm $unmapped_fasta* $fivep_to_reads* $fivep_trimmed_reads $trimmed_reads_to_threep*  
