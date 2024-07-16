import argparse
import collections
import logging

from intervaltree import Interval, IntervalTree
import pandas as pd
import numpy as np
# import psutil

import functions



# =============================================================================#
#                                  Globals                                     #
# =============================================================================#
INTERMEDIATE_OUTFILES = {
						'output_bam': '${output_base}output.bam',
						'mapped_bam': '${output_base}mapped_reads.bam',
						'unmapped_bam': '${output_base}unmapped_reads.bam',
						'run_data': '${output_base}read_counts.tsv',
						'unmapped_fasta': '${output_base}unmapped_reads.fa',
						'unmapped_fasta': '${output_base}unmapped_reads.fai',
						'unmapped_fasta': '${output_base}unmapped_reads.fa.1.bt2',
						'unmapped_fasta': '${output_base}unmapped_reads.fa.2.bt2',
						'unmapped_fasta': '${output_base}unmapped_reads.fa.3.bt2',
						'unmapped_fasta': '${output_base}unmapped_reads.fa.4.bt2',
						'unmapped_fasta': '${output_base}unmapped_reads.fa.rev.1.bt2',
						'unmapped_fasta': '${output_base}unmapped_reads.fa.rev.2.bt2',
						'fivep_to_reads': '${output_base}fivep_to_reads.sam',
						'fivep_trimmed_reads': '${output_base}fivep_mapped_reads_trimmed.fa',
						'trimmed_reads_to_genome': '${output_base}trimmed_reads_to_genome.sam',
						'reads_to_fivep': '${output_base}reads_to_fivep.sam',
						'fivep_info_table': '${output_base}fivep_info_table.tsv',
						'trimmed_info_table': '${output_base}trimmed_info_table.tsv',
						'failed_fivep': '${output_base}failed_fivep_alignments.tsv',
						'failed_trimmed': '${output_base}failed_trimmed_alignments.tsv',
						'failed_lariat': '${output_base}failed_lariat_alignments.tsv',
						}
# 'FIVEP_INDEX': '/datasets2/lariat_mapping/testing/output/test_ref/fivep_sites'
RESULTS_OUTFILES = {
	
}


# =============================================================================#
#                                  Functions                                   #
# =============================================================================#
def linear_map_to_genome(seq_type, threads, genome_index, read_file=None, read_one=None, read_two=None):
	command = 'hisat2 --no-softclip -k 1 --max-seeds 20 --pen-noncansplice 0 --n-ceil L,0,0.05 --score-min L,0,-0.24 --bowtie2-dp 1'
	command += f'--threads {threads}  -x ${genome_index}'
	if seq_type == 'single-end':
		command += '-U $READ_FILE' 
		command += ' \ \n| samtools view --bam --with-header' '--add-flags PAIRED,READ1'
			> $output_bam \
			|| exit 1
	elif seq_type == 'paired-end':
		command += '-1 ${read_one} -2 ${read_two}' \
			| samtools view --bam --with-header \
			> $output_bam \
			|| exit 1
	fi

	printf "$(date +'%d/%b/%y %H:%M:%S') | Extracting unmapped reads...\n"
samtools view --bam --with-header --exclude-flags 4 $output_bam > $mapped_bam
samtools view --bam --with-header --require-flags 4 $output_bam > $unmapped_bam


def classify_reads():



def end_early():


def end_error():
	
# =============================================================================#
#                                    Main                                      #
# =============================================================================#
def main(**kwargs):

	linear_map_to_genome()
	### Map filtered reads to genome and keep unmapped reads. Lariat reads crossing the brachpoint will not be able to map to the gene they're from
printf "$(date +'%d/%b/%y %H:%M:%S') | Mapping reads to genome...\n"
if $single_end; then
	# hisat2 --no-softclip -k 1 --max-seeds 20 --pen-noncansplice 0 --n-ceil L,0,0.05 --score-min L,0,-0.24 --bowtie2-dp 1 \
	hisat2 --no-softclip -k 5 --max-seeds 20 --pen-noncansplice 0 --n-ceil L,0,0.05 --score-min L,0,-0.24 --bowtie2-dp 1 \
	       --threads $THREADS -x $GENOME_INDEX -U $READ_FILE \
		| samtools view --bam --with-header --add-flags PAIRED,READ1 \
		> $output_bam \
		|| exit 1
else
	# hisat2 --no-softclip -k 1 --max-seeds 20 --pen-noncansplice 0 --n-ceil L,0,0.05 --score-min L,0,-0.24 --bowtie2-dp 1 \
	hisat2 --no-softclip -k 5 --max-seeds 20 --pen-noncansplice 0 --n-ceil L,0,0.05 --score-min L,0,-0.24 --bowtie2-dp 1 \
		   --threads $THREADS -x $GENOME_INDEX -1 $READ_ONE -2 $READ_TWO \
		| samtools view --bam --with-header \
		> $output_bam \
		|| exit 1
fi


	process_linear_map()
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


	build_read_index()
### Build a bowtie2 index of the unmapped reads
printf "$(date +'%d/%b/%y %H:%M:%S') | Building bowtie2 index of unmapped fasta...\n"
bowtie2-build --quiet --threads $THREADS $unmapped_fasta $unmapped_fasta || exit 1 

	map_fivep_to_reads()
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

	filter_fivep_aligns()
## Extract reads with a mapped 5' splice site and trim it off
printf "$(date +'%d/%b/%y %H:%M:%S') | Finding 5' read alignments and trimming reads...\n"
python -u $PIPELINE_DIR/scripts/filter_fivep_aligns.py $THREADS $GENOME_FASTA $FIVEP_FASTA $OUTPUT_BASE $LOG_LEVEL \
	|| exit 1 

	map_trimmed_to_genome()
### Map 5' trimmed reads to genome
printf "$(date +'%d/%b/%y %H:%M:%S') | Mapping 5' trimmed reads to genome...\n"
hisat2 --no-softclip --no-spliced-alignment --very-sensitive -k 100 \
	   --no-unal --threads $THREADS -f -x $GENOME_INDEX -U $fivep_trimmed_reads \
	| samtools sort --threads $THREADS --verbosity 0 --output-fmt SAM -n \
	| samtools view \
	> $trimmed_reads_to_genome \
	|| exit 1

	
	filter_trim_aligns()
### Filter trimmed alignments
printf "$(date +'%d/%b/%y %H:%M:%S') | Analyzing trimmed alignments and outputting lariat table...\n"
python -u $PIPELINE_DIR/scripts/filter_trim_aligns.py $THREADS $INTRONS_TSV $GENOME_FASTA $OUTPUT_BASE $LOG_LEVEL \
	|| exit 1 

	filter_trim_aligns()
### Filter lariat mappings and choose 1 for each read
printf "$(date +'%d/%b/%y %H:%M:%S') | Filtering putative lariat alignments...\n"
python -u $PIPELINE_DIR/scripts/filter_lariats.py $GENOME_FASTA $REPEATS_BED $OUTPUT_BASE $LOG_LEVEL \
	|| exit 1 

	make_track()
### Make a custom track BED file of identified lariats 
if $UCSC_TRACK; then
	printf "$(date +'%d/%b/%y %H:%M:%S') | Making UCSC Genome Browser track...\n"
	python -u $PIPELINE_DIR/scripts/make_track.py $OUTPUT_BASE $LOG_LEVEL \
		|| exit 1
fi

	classify_reads()
### Classify reads
python -u $PIPELINE_DIR/scripts/classify_linear.py $EXONS_TSV $INTRONS_TSV $OUTPUT_BASE $LOG_LEVEL \
	|| exit 1
python -u $PIPELINE_DIR/scripts/classify_nonlinear.py $OUTPUT_BASE $LOG_LEVEL \
	|| exit 1



wait
	end()
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


if __name__ == '__main__':
	main()


