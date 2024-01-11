#!/bin/bash

#=============================================================================#
#                                  Arguments                                  #
#=============================================================================#

# RNA-seq fastq read file
READ_FILE=$1
# Output directory 
OUTPUT_DIR=$2
# Output base name 
NAME=$3
# Number of CPUs to use
CPUS=$4
# Genome Bowtie2 index base name
GENOME_BOWTIE2_INDEX=$5
# Reference genome FASTA file
GENOME_FASTA=$6
# GTF file containing gene annotatin for mapping genome
GTF_FILE=$7
# FASTA file of 5' splice sites (first 20nts of all introns)
FIVEP_FASTA=$8
# Custom file of sequences in 5nt window upstream of 5'ss
FIVEP_UPSTREAM=$9
# Bowtie2 index of 3' splice sites genome (last 250nts of all introns)
THREEP_BOWTIE2_INDEX="${10}"
# TSV file with 3' splice site coordinates and lengths (max 250)
THREEP_LENGTHS="${11}"

#=============================================================================#
#                                    Calls                                    #
#=============================================================================#
### Map filtered reads to genome and keep unmapped reads. Lariat reads crossing the brachpoint will not be able to map to the gene they're from
echo ""
printf "$(date +'%m/%d/%y - %H:%M:%S') | Mapping reads and extracting unmapped reads...\n"
output_bam=$OUTPUT_DIR/$NAME"_mapped_reads.bam"
unmapped_bam=$OUTPUT_DIR/$NAME"_unmapped_reads.bam"
bowtie2 --end-to-end --sensitive --score-min L,0,-0.24 -k 1 --n-ceil L,0,0.05 --threads $CPUS -x $GENOME_BOWTIE2_INDEX -U $READ_FILE \
	| samtools view --bam --with-header > $output_bam
samtools view --bam --with-header --require-flags 4 $output_bam > $unmapped_bam
mapped_read_count=$(samtools view --count --exclude-flags 4 $output_bam)
unmapped_read_count=$(samtools view --count $unmapped_bam)
run_data=$OUTPUT_DIR/$NAME"_run_data.tsv"
echo -e "ref_mapped_reads\t$mapped_read_count" > $run_data
echo -e "ref_unmapped_reads\t$unmapped_read_count" >> $run_data

### Create fasta file of unmapped reads 
echo ""
printf "$(date +'%m/%d/%y - %H:%M:%S') | Creating fasta file of unmapped reads...\n"
unmapped_fasta=$OUTPUT_DIR/$NAME"_unmapped_reads.fa"
samtools fasta $unmapped_bam > $unmapped_fasta
samtools faidx $unmapped_fasta

### Build a bowtie index of the unmapped reads
echo ""
printf "$(date +'%m/%d/%y - %H:%M:%S') | Building bowtie index of unmapped fasta...\n"
bowtie2-build --large-index --threads $CPUS $unmapped_fasta $unmapped_fasta > /dev/null

### Align unmapped reads to fasta file of all 5' splice sites (first 20nts of introns)
echo ""
printf "$(date +'%m/%d/%y - %H:%M:%S') | Mapping 5' splice sites to reads...\n"
fivep_to_reads=$OUTPUT_DIR/$NAME"_fivep_to_reads.sam"
bowtie2 --end-to-end --sensitive --no-unal -f -k 10000 --score-min C,0,0 --threads $CPUS -x $unmapped_fasta -U $FIVEP_FASTA \
	| samtools view > $fivep_to_reads

### Extract reads with a mapped 5' splice site and trim it off
echo ""
printf "$(date +'%m/%d/%y - %H:%M:%S') | Finding 5' read alignments and trimming reads...\n"
fivep_trimmed_reads=$OUTPUT_DIR/$NAME"_fivep_mapped_reads_trimmed.fa"
fivep_info_table=$OUTPUT_DIR/$NAME"_fivep_info_table.tsv"
python scripts/filter_fivep_alignments.py $unmapped_fasta $fivep_to_reads $FIVEP_UPSTREAM $fivep_trimmed_reads $fivep_info_table $run_data

### Map 5' trimmed reads to 3' sites (last 250nts of introns)
echo ""
printf "$(date +'%m/%d/%y - %H:%M:%S') | Mapping 5' trimmed reads to 3' sites...\n"
trimmed_reads_to_threep=$OUTPUT_DIR/$NAME"_fivep_reads_trimmed_mapped_to_threep.sam"
bowtie2 --end-to-end --sensitive -k 10 --no-unal --threads $CPUS -f -x $THREEP_BOWTIE2_INDEX -U $fivep_trimmed_reads \
	| samtools view > $trimmed_reads_to_threep

### Filter 3' splice site alignments and output info table, including the branchpoint site
echo ""
printf "$(date +'%m/%d/%y - %H:%M:%S') | Analyzing 3' alignments and outputting lariat table...\n"
python scripts/filter_threep_alignments.py $trimmed_reads_to_threep $THREEP_LENGTHS $fivep_info_table $GTF_FILE $GENOME_FASTA $OUTPUT_DIR/$NAME $run_data


### Delete all intermediate/uneeded files that were created throughout this process
wait
rm $output_bam
rm $unmapped_bam
# rm $unmapped_fasta* $fivep_to_reads* $fivep_trimmed_reads $trimmed_reads_to_threep*  
