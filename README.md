# Lariat-Mapping | The Fairbrother Lab

## Overview

A splicing branchpoint reporter pipeline.

## Dependencies

This pipeline has the following dependencies:
- python3 (tested with v3.10.6)
- bowtie2 (tested with v2.4.5)
- samtools (tested with v1.15.1)
- bedtools (tested with v2.30.0)
- numpy (tested with v1.23.2)
- [pyfaidx](https://pypi.org/project/pyfaidx/) (tested with v0.7.2.1)
- [intervaltree](https://pypi.org/project/intervaltree/) (tested with v3.1.0)

These dependencies are included in the file `environment.yaml` which can be used to make a conda environment for the pipeline by running `conda env create -f environment.yaml`. Then, activate the environment with `conda activate larmap_env` 

For M1 mac users: please install packages `bowtie2`, `bedtools`, and `samtools` using the command `arch -arm64 brew install [package]` before running `conda`, if any of the above pacakges has not previously been installed.

## Additional Reference File

Prior to running the pipeline, please use `python get_splice_site_seqs.py [introns].bed [genome].fa [output_name_prefix]` to produce the reference file of sequences in 5nt window upstream of 5' splice sites; the output is to be used as the input of the `-u`, or `--ref_5p_upstream` option in the pipeline.

## Running the Pipeline

To run the larmap pipeline, use `./larmap_run.sh` with the following arguments:

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
      -m, --ref_repeatmasker    BED file of repetitive elements from repeatmasker

A directory named `[output_base_name]_lariat_mapping` will be created in `output_dir`. Upon completion of the pipeline, this directory will contain a tab-separated results file with lariat read info called `[output_base_name]_lariat_reads.txt`.

## Pipeline Workflow

1. `larmap_run.sh` calls `map_lariats.sh` sequentially for the read one (R1) and read two (R2) FASTQ files. This will produce three files in the output subdirectory for the read file:
    -`[output_base_name]_total_reads.txt` (one line file containing count of linearly-aligned reads from read file)
    - `[output_base_name]_fivep_info_table.txt` (intermediate file containing info on the mapping of the 5'SS sequences to the unmapped reads)
    - `[output_base_name]_final_info_table.txt` (results file containing candidate lariat reads obtained after mapping the 5'SS trimmed reads to the 3'SS region sequences)

    The mapping script `map_lariats.sh` will:
    - Align reads to the reference genome with bowtie2; save mapped read count and proceed with unmapped reads
    - Convert the unmapped reads bam file to FASTA format with samtools
    - Build a bowtie2 index of the unmapped reads FASTA file
    - Align a FASTA file of 5'SS to the unmapped reads index
    - Trim reads with 5'SS alignments and write trimmed reads to FASTA file
    - Align the trimmed reads to a Bowtie2 index of 3'SS regions
    - Take the mapped trimmed reads from and create an output file containing candidate lariat reads

3. The `filter_lariats.py` script loads intron and gene information from provided annotation files and performs post-mapping filtering before outputting the final lariat mapping results. 

    The candidate lariat reads are filtered based on the following criteria: - BP is within 2bp of a splice site (likely from an intron circle, not a lariat) - 5'SS and 3'SS are not in the correct order - Read maps to a Ubiquitin gene (likely false positive due to repetitive nature of gene) - There is a valid aligment for the 3' segment upstream of the 5' segment - Both the 5'SS and the BP overlap with repetitive regions from RepeatMasker (likely false positive) - NEEDS TO BE IMPLEMENTED: Filter out template switching reads (eg. 5'SS 6bp sequence matches 6bp sequence downstream of BP) - NEEDS TO BE IMPLEMENTED: Correct BP position to account for RT skipping

    After filtering, the final set of lariat reads are output to the `results_path` file. The initial columns of this tab-delimited file are any additional custom columns which were added to the info file. After these columns, the file contains the following data for each lariat read: - Gene name - Gene Ensembl ID - Gene type (e.g. protein-coding, lncRNA) - Read ID - Read sequence - Chromosome - Strand - 5' splice site coordinate - 3' splice site coordinate - Branchpoint coordinate - Read branchpoint nucleotide - Genomic branchpoint nucleotide - Genomic branchpoint sequence context (10bp window with the BP at position 5) - Branchpoint position relative to nearest 3' splice site - Total count of linearly-aligned reads (useful for normalization when comparing lariat levels between samples)
