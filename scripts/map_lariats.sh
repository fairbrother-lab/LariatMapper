#!/bin/bash

#=============================================================================#
#                                  Arguments                                  #
#=============================================================================#
# RNA-seq fastq read file
READ_FILE=$1
# Output directory 
OUTPUT_BASE=$2
# Number of CPUs to use
CPUS=$3
# Genome Bowtie2 index base name
GENOME_INDEX=$4
# Reference genome FASTA file
GENOME_FASTA=$5
# GTF file containing gene annotatin for mapping genome
GTF_FILE=$6
# FASTA file of 5' splice sites (first 20nts of all introns)
FIVEP_FASTA=$7
# Custom file of sequences in 5nt window upstream of 5'ss
FIVEP_UPSTREAM=$8
# Annotated transcripts
TRANSCRIPTS_BED="${9}"
# Annotated introns
INTRONS_BED="${10}"
# Annotated exons
EXONS_BED="${11}"
# Annotated repeat regions
REPEATS_BED="${12}"
# Keep the intermediate files created during the run ("True" or "False", default "False")
KEEP_INTERMEDIATES="${13}"
# REFERENCE_DIR=



#=============================================================================#
#                                    Variables                                #
#=============================================================================#
# genome_index=$REFERENCE_DIR/



#=============================================================================#
#                                    Calls                                    #
#=============================================================================#
### Map filtered reads to genome and keep unmapped reads. Lariat reads crossing the brachpoint will not be able to map to the gene they're from
# echo ""
# printf "$(date +'%m/%d/%y - %H:%M:%S') | Mapping reads and extracting unmapped reads...\n"
# output_bam="$OUTPUT_BASE"mapped_reads.bam
# unmapped_bam="$OUTPUT_BASE"unmapped_reads.bam
# bowtie2 --end-to-end --sensitive --score-min L,0,-0.24 -k 1 --n-ceil L,0,0.05 --threads $CPUS -x $GENOME_INDEX -U $READ_FILE \
# 	| samtools view --bam --with-header \
# 	> $output_bam
# samtools view --bam --with-header --require-flags 4 $output_bam > $unmapped_bam
# mapped_read_count=$(samtools view --count --exclude-flags 4 $output_bam)
# unmapped_read_count=$(samtools view --count $unmapped_bam)
# run_data="$OUTPUT_BASE"run_data.tsv
# echo -e "ref_mapped_reads\t$mapped_read_count" > $run_data
# echo -e "ref_unmapped_reads\t$unmapped_read_count" >> $run_data

# ### Create fasta file of unmapped reads 
# echo ""
# printf "$(date +'%m/%d/%y - %H:%M:%S') | Creating fasta file of unmapped reads...\n"
# unmapped_fasta="$OUTPUT_BASE"unmapped_reads.fa
# samtools fasta $unmapped_bam > $unmapped_fasta
# samtools faidx $unmapped_fasta

# ### Build a bowtie index of the unmapped reads
# echo ""
# printf "$(date +'%m/%d/%y - %H:%M:%S') | Building bowtie index of unmapped fasta...\n"
# bowtie2-build --large-index --threads $CPUS $unmapped_fasta $unmapped_fasta > /dev/null

### Align unmapped reads to fasta file of all 5' splice sites (first 20nts of introns)
# echo ""
# printf "$(date +'%m/%d/%y - %H:%M:%S') | Mapping 5' splice sites to reads...\n"
# fivep_to_reads="$OUTPUT_BASE"fivep_to_reads.sam
# bowtie2 --end-to-end --sensitive --no-unal -f -k 10000 --score-min C,0,0 --threads $CPUS -x $unmapped_fasta -U $FIVEP_FASTA \
# 	| samtools view \
# 	> $fivep_to_reads

### Extract reads with a mapped 5' splice site and trim it off
# echo ""
# printf "$(date +'%m/%d/%y - %H:%M:%S') | Finding 5' read alignments and trimming reads...\n"
fivep_trimmed_reads="$OUTPUT_BASE"fivep_mapped_reads_trimmed.fa
# fivep_info_table="$OUTPUT_BASE"fivep_info_table.tsv
# python scripts/filter_fivep_alignments.py $unmapped_fasta $fivep_to_reads $FIVEP_UPSTREAM $fivep_trimmed_reads $fivep_info_table $OUTPUT_BASE

### Map 5' trimmed reads to genome
printf "$(date +'%m/%d/%y - %H:%M:%S') | Mapping 5' trimmed reads to genome...\n"
trimmed_reads_to_genome="$OUTPUT_BASE"trimmed_reads_to_genome.bam
bowtie2 --end-to-end --very-sensitive -k 10 --no-unal --threads $CPUS -f -x $GENOME_INDEX -U $fivep_trimmed_reads \
	| samtools view --bam \
	| bedtools tag -names -f 1 -i - -files $INTRONS_BED \
	> $trimmed_reads_to_genome

### Filter trimmed alignments
echo ""
printf "$(date +'%m/%d/%y - %H:%M:%S') | Analyzing trimmed alignments and outputting lariat table...\n"
# python -u scripts/filter_trimmed_alignments.py $GTF_FILE $INTRONS_BED $GENOME_FASTA $OUTPUT_BASE $KEEP_INTERMEDIATES
scalene --html --outfile "$OUTPUT_BASE"filter_trim_profile.html scripts/filter_trimmed_alignments.py $GTF_FILE $INTRONS_BED $GENOME_FASTA $OUTPUT_BASE $KEEP_INTERMEDIATES

### Filter lariat mappings and choose 1 for each read
printf "$(date +'%m/%d/%y - %H:%M:%S') | Filtering putative lariat mappings...\n"
python -u scripts/filter_lariats.py $GTF_FILE $INTRONS_BED $REPEATS_BED $OUTPUT_BASE $KEEP_INTERMEDIATES

wait
### Delete the intermediate files 
if [ "$KEEP_INTERMEDIATES" = "False" ]; then
	echo "Deleting intermediate files..."
	rm $output_bam
	rm $unmapped_bam
	rm $fivep_to_reads 
	rm $fivep_trimmed_reads 
	rm $trimmed_reads_to_genome
	rm $unmapped_fasta 
	rm $unmapped_fasta.fai
	rm $unmapped_fasta.1.bt2l 
	rm $unmapped_fasta.2.bt2l 
	rm $unmapped_fasta.3.bt2l 
	rm $unmapped_fasta.4.bt2l 
	rm $unmapped_fasta.rev.1.bt2l 
	rm $unmapped_fasta.rev.2.bt2l 
fi 

printf "$(date +'%m/%d/%y - %H:%M:%S') | Lariat mapping complete.\n"