# WARNING: WORK IN PROGRESS, NOT CURRENTLY RELIABLE


# Lariat Mapping | The Fairbrother Lab

## Overview

A pipeline for mapping lariat-derived reads present in RNA-seq data through the identification of reads with gapped, inverted alignments to introns.

## Setup

This pipeline has the following dependencies:
- python3
- bowtie2
- samtools
- bedtools
- numpy
- [pyfaidx](https://pypi.org/project/pyfaidx/) (tested with v0.7.2.1)
- [intervaltree](https://pypi.org/project/intervaltree/) (tested with v3.1.0)

These dependencies are included in the file `environment.yaml` which can be used to make a conda environment for the pipeline by running `conda env create -f environment.yaml`. The environment can then be activated by running `conda activate larmap_env` before running scripts in the pipeline (RECOMMENDED).
For M1 mac users: please install packages `bowtie2`, `bedtools`, and `samtools` using the command `arch -arm64 brew install [package]` before running `conda`, if any of the above pacakges has not previously been installed.

The pipeline requires the following standard reference files: 
- FASTA file of the reference genome
- bowtie2 index of the reference genome
- GTF file containing gene annotations for the reference genome
- BED file containing intron annotations for the reference genome (available on the [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables) in group "Genes and Gene Predictions" by selecting output format)
- BED file containing the RepeatMasker annotation for the mapping genome (available on the [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables) in group "Repeats", track "RepeatMasker")

To produce the required custom reference files, run `python get_splice_site_seqs.py [introns BED file] [genome fasta file] [output file prefix]`, then make a bowtie2 index of the threep_sites.fa file by running `bowtie2-build [output file prefix].threep_sites.fa [output file prefix].threep_sites]`. 

## Running the Pipeline
Run `bash larmap_run.sh` with the following arguments:

      -r, --read_file           FASTQ file
      -o, --output_dir          Directory for output files
      -e, --output_base_name    Prefix to add to output files
      -c, --num_cpus            Number of CPUs available
      -i, --ref_b2index         Bowtie2 index of the full reference genome
      -f, --ref_fasta           FASTA file of the full reference genome
      -g, --ref_gtf             GTF file with gene annotation of the reference genome
      -5, --ref_5p_fasta        FASTA file with sequences of first 20nt from reference 5' splice sites (first 20nt of introns)
      -u, --ref_5p_upstream     Custom file of sequences in 5nt window upstream of 5' splice sites
      -3, --ref_3p_b2index      Bowtie2 index file of last 250nt from reference 3' splice sites (last 250nt of introns)
      -l, --ref_3p_lengths      Custom file with the lengths of the sequences in ref_3p_b2index (some introns are <250nt)
      -n, --ref_introns         BED file of all introns in the reference genome
      -m, --ref_repeatmasker    BED file of repetitive elements from RepeatMasker

A directory named `[output_base_name]_lariat_mapping` will be created in `[output_dir]`. Upon completion of the pipeline, this directory will contain a results file named `[output_base_name]_lariat_reads.tsv` .

## Output
`[output_base_name]_lariat_reads.tsv` contains a table in tab-separated values format. Each row is an RNA-seq read that has a valid lariat mapping.

The columns are as follows:

	gene_name            The name of the gene that produced the lariat
	gene_id              The Ensembl ID of the gene that produced the lariat
	gene_type            The type of the gene that produced the lariat 
	read_id              The RNA-seq read's ID
	read_seq             The RNA-seq read's DNA sequence
	chrom                The chromosome of the gene that produced the lariat
	strand               The strand of the gene that produced the lariat. "+" for the forward strand and "-" for the reverse strand
	fivep_pos            The genomic position of the lariat's 5' splice site
	threep_pos           The genomic position of the closest 3' splice site that is downstream of the branchpoint
	bp_pos               The genomic position of the lariat's branchpoint 
	read_bp_nt           The nucleotide of the branchpoint according to the RNA-seq read's sequence. Reverse-complemented if strand is "-"
	genomic_bp_nt        The nucleotide of the branchpoint according to the reference genome. Reverse-complemented if strand is "-"
	genomic_bp_context   The genomic sequence from positions -4 to +5 of the branchpoint. Reverse-complemented if strand is "-"
	bp_dist_to_threep    The distance of the branchpoint to the 3' splice site in nucleotides
	total_mapped_reads   The number of input reads that mapped linearly to the reference genome. Identical across all rows



## Pipeline Workflow

1. `larmap_run.sh` calls `map_lariats.sh` on the FASTQ file. This will produce three files in the output subdirectory for the read file:
    - `[output_base_name]_fivep_info_table.tsv` (intermediate file containing info on the mapping of the 5'SS sequences to the unmapped reads)
    - `[output_base_name]_final_info_table.tsv` (results file containing candidate lariat reads obtained after mapping the 5'SS trimmed reads to the 3'SS region sequences)
    -`[output_base_name]_run_data.tsv` (miscellaneous read counts)

    The mapping script `map_lariats.sh` will:
    - Align reads to the reference genome with bowtie2; save mapped read count and proceed with unmapped reads
    - Convert the unmapped reads bam file to FASTA format with samtools
    - Build a bowtie2 index of the unmapped reads FASTA file
    - Align a FASTA file of 5'SS to the unmapped reads index
    - Trim reads with 5'SS alignments and write trimmed reads to FASTA file
    - Align the trimmed reads to a Bowtie2 index of 3'SS regions
    - Take the mapped trimmed reads from and create an output file containing candidate lariat reads

3. The `filter_lariats.py` script loads intron and gene information from provided annotation files and performs post-mapping filtering before outputting the final lariat mapping results. 

    The candidate lariat reads are filtered based on the following criteria:
   - BP is within 2bp of a splice site (likely from an intron circle, not a lariat)
   - 5'SS and 3'SS are not in the correct order
   - Read maps to a Ubiquitin gene (likely false positive due to repetitive nature of gene)
   - There is a valid aligment for the 3' segment upstream of the 5' segment
   - Both the 5'SS and the BP overlap with repetitive regions from RepeatMasker (likely false positive)
