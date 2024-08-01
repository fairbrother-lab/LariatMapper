# LariatMapper | The Fairbrother Lab

## Overview

A pipeline for extracting lariats and their branchpoints from RNA-sequencing data. 

Input: RNA-sequencing reads and a reference genome

Output: Lariats, circularized introns, and template-switching events

## Setup

### Dependencies
This pipeline has the following dependencies:
- bowtie2
- hisat2
- samtools
- bedtools
- python3
- numpy
- pandas
- pyarrow
- [pyfaidx](https://pypi.org/project/pyfaidx/)
- [intervaltree](https://pypi.org/project/intervaltree/)

These dependencies are included in the file `environment.yaml` which can be used to make a conda environment for the pipeline by running

	conda env create -f environment.yaml

The environment can then be activated by running 

	conda activate larmap_env
 
before running scripts in the pipeline (RECOMMENDED).

For M1 mac users: please install packages `bowtie2`, `bedtools`, and `samtools` using the command `arch -arm64 brew install [package]` before running `conda`, if any of the above pacakges has not previously been installed.

### Reference files
Create a directory to store the neccessary reference files with `build_references.py`. Alternatively, you can skip this step and input the reference files at runtime.

You will need:
- A FASTA file of reference genome sequences
- A GTF or GFF file of annotations for the reference genome
- A hisat2 index of the reference genome

Run `build_references.py` with the following arguments:

	python build_references.py -f <REF_FASTA> -a <REF_ANNO> -i <HISAT2_INDEX> -o <OUT_DIR>

If you have a BED file of repetitive regions from RepeatMasker, you can include the argument `-r <REF_REPEATMASKER>` to copy it to the reference directory. You can find such a file for several reference genomes on the UCSC Genome Browser (https://genome.ucsc.edu/cgi-bin/hgTables) in group "Repeats", track "RepeatMasker".

You can then use `OUT_DIR` as the reference files directory when running the pipeline (argument `-r, --ref_dir`)

### Input
LariatMapper accepts FASTQ-format RNA-sequencing data. 

The data should be preprocessed to remove low-quality reads, adapter sequences, and unique molecular identifiers (UMIs). De-duplication is not required but *is* recommended.


## Running the Pipeline
Run `python larmap.py` with the following arguments:


## Output
All output will be written in the directory specified with `-o, --output_dir`. This includes:
	
- `lariat_reads.tsv`: A table of lariat reads and their alignments
- `circularized_intron_reads.tsv`: A table of circularized intron reads and their alignments
- `template_switching_reads.tsv`: A table of 
- `summary_statistics.tsv`: A collection of read counts for different categories of RNA, read counts for each step in the pipeline, and performance measures

`lariat_reads.tsv` and `circularized_intron_reads.tsv` columns:

- read_id: The read's ID (unique)
- gene_name*: The name of the gene that produced the lariat 
- gene_id*: The Ensembl ID of the gene that produced the lariat
- gene_type*: The type of the gene that produced the lariat 
- chrom: The chromosome of the gene that produced the lariat
- strand: The strand of the gene that produced the lariat. "+" for the forward strand and "-" for the reverse strand
- fivep_pos: The genomic position of the lariat's 5' splice site 
- threep_pos: The genomic position of the closest 3' splice site that is downstream of the branchpoint
- bp_pos: The genomic position of the lariat's branchpoint 
- read_bp_nt: The nucleotide of the branchpoint according to the read's sequence. Reverse-complemented if strand is "-"
- bp_dist_to_threep: The distance of the branchpoint to the 3' splice site in nucleotides
- read_alignment: "forward" if the RNA-seq read's sequence matches the source intron's sequence, "reverse" if it matches the reverse-complement
- read_bp_pos: The position of the branchpoint in the read
- read_seq: The read's DNA sequence 
- read_bp_nt: The nucleotide of the branchpoint according to the read. Reverse-complemented if read_alignment is "reverse"
- genomic_bp_nt: The nucleotide of the branchpoint according to the reference genome. Reverse-complemented if strand is "-"
- genomic_bp_context: The genomic sequence from positions -4 to +5 of the branchpoint. Reverse-complemented if strand is "-"
- total_mapped_reads: The number of input reads that mapped linearly to the reference genome (identical across all rows)

*may be multiple comma-delimited values

`template_switching_reads.tsv` columns:

- read_id: The read's ID (unique)
- fivep_sites*: The 5' splice sites that mapped to the read. Format is \<chromosome>;\<strand>;\<position>,...
- temp_switch_sites*: The location where the reverse transcriptase transfered to. Format is \<chromosome>;\<position>,...
- read_seq*: The read's DNA sequence
- fivep_seq*: The 5' splice sites' DNA sequence
- genomic_bp_context*: The genomic sequence from positions -4 to +5 of the branchpoint. Reverse-complemented if strand is "-"
- read_bp_pos*: The position of the branchpoint in the read 

*may be multiple comma-delimited values

All position values (`fivep_pos`, `threep_pos`, `bp_pos`) are 0-based and inclusive. If a lariat's 5' splice site and branchpoint could be attributed to multiple gene annotations, the gene values will appear like so:

	gene_name	gene_id	gene_type
	PCBP2,ENSG00000257379	ENSG00000197111.16,ENSG00000257379.1	lncRNA,protein_coding