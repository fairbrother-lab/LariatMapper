# LariatMapper | The Fairbrother Lab

## Overview

A pipeline for mapping lariat-derived reads (AKA lariat reads) that are present in RNA-seq data through the identification of reads with gapped, inverted alignments to introns. Secondarily identifies reads derived from circularized intron RNA.

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
LariatMapper accepts FASTQ-format read data from 2nd-generation RNA-sequencing experiments, both paired-end and single-end. It does not currently support strand-specific or 3rd-generation sequencing data. 

Sequencing data should be preprocessed to remove low-quality reads, adapter sequences, and unique molecular identifiers (UMIs). De-duplication is not required but *is* recommended.


## Running the Pipeline
Run `python larmap.py` with the following arguments:

  	-r, --ref_dir				Directory with reference files for lariat mapping. Create by running build_references.py
	-o, --output_dir			Directory for output files (will be created if it does not exist)
	For paired-end sequencing data
	  -1, --read_one			Read 1 input FASTQ file when processing paired-end RNA-seq data. Can be uncompressed or gzip-compressed. 
	  -2, --read_two			Read 2 input FASTQ file when processing paired-end RNA-seq data. Can be uncompressed or gzip-compressed. 
	For single-end sequencing data
	  -f, --read_file			Input FASTQ file when processing single-end RNA-seq data. Can be uncompressed or gzip-compressed.
	
If you did not create a reference directory and instead want to input the neccesary reference files at runtime, use the following arguments instead of `-r, --ref_dir`:

	  -i, --ref_h2index			hisat2 index of the reference genome
	  -g, --ref_fasta			FASTA file of the reference genome
	  -a, --ref_anno			Gene annotation of the reference genome in GTF or GFF format (may be gzipped with .gz extension)
	  -5, --ref_5p_fasta		FASTA file with sequences of first 20nt of annotated introns
	  -n, --ref_introns			BED file of all annotated introns
	  -m, --ref_repeatmasker	BED file of repetitive regions annotated by RepeatMasker. Putative lariats that map to a repetitive region will be filtered out as false positives (Optional)

Optional arguments:

	  -q, --quiet				Don't print any status messages (work in progress)
	  -d, --debug				Print extensive any status messages (work in progress)
	  -t, --threads				Number of threads to use for parallel processing (default=1)
	  -p, --output_prefix		Add a prefix to output file names (-o OUT -p ABC -> OUT/ABC_lariat_reads.tsv)
	  -u, --ucsc_track			Add an output file named "lariat_reads.bed" which can be used as a custom track in the UCSC Genome Browser (https://www.genome.ucsc.edu/cgi-bin/hgCustom) to visualize lariat alignments
	  -k, --keep_intermediates	Don't delete the intermediate files created while running the pipeline (default=delete)
     

## Output
All output will be written in the directory specified with `-o, --output_dir`. For most users, the main output of interest will be `lariat_reads.tsv`, which is a tab-delimited table of RNA-sequencing reads identified as lariat reads. 

`lariat_reads.tsv` will contain the following columns:

	read_id           		The read's ID (unique)
	gene_name         		The name of the gene that produced the lariat (may be multiple comma-delimited values)
	gene_id           		The Ensembl ID of the gene that produced the lariat (may be multiple comma-delimited values)
	gene_type         		The type of the gene that produced the lariat (may be multiple comma-delimited values)
	chrom             		The chromosome of the gene that produced the lariat
	strand            		The strand of the gene that produced the lariat. "+" for the forward strand and "-" for the reverse strand
	fivep_pos         		The genomic position of the lariat's 5' splice site 
	threep_pos        		The genomic position of the closest 3' splice site that is downstream of the branchpoint
	bp_pos            		The genomic position of the lariat's branchpoint 
	read_bp_nt        		The nucleotide of the branchpoint according to the read's sequence. Reverse-complemented if strand is "-"
	bp_dist_to_threep 		The distance of the branchpoint to the 3' splice site in nucleotides
	read_alignment			"forward" if the RNA-seq read's sequence matches the source intron's sequence, "reverse" if it matches the reverse-complement
	read_bp_pos				The position of the branchpoint in the read
	read_seq          		The read's DNA sequence 
	read_bp_nt				The nucleotide of the branchpoint according to the read. Reverse-complemented if read_alignment is "reverse"
	genomic_bp_nt     		The nucleotide of the branchpoint according to the reference genome. Reverse-complemented if strand is "-"
	genomic_bp_context		The genomic sequence from positions -4 to +5 of the branchpoint. Reverse-complemented if strand is "-"
	total_mapped_reads		The number of input reads that mapped linearly to the reference genome. Identical across all rows

All position values (`fivep_pos`, `threep_pos`, `bp_pos`) are 0-based and inclusive. If a lariat's 5' splice site and branchpoint could be attributed to multiple gene annotations, the gene values will appear like so:

	gene_name	gene_id	gene_type
	PCBP2,ENSG00000257379	ENSG00000197111.16,ENSG00000257379.1	lncRNA,protein_coding