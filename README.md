# LariatMapper | The Fairbrother Lab

## Overview

A pipeline for extracting lariats and their branchpoints from RNA-sequencing data. 

Input: RNA-sequencing reads and a reference genome

Output: Lariats, circularized introns, and template-switching events

## Setup

### Dependencies
The software dependencies are detailed in `requirements.txt`. We recommend creating a dedicated programming environment with [mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) for LariatMapper to avoid dependency-related problems during use.

If you have mamba installed, you can create a new environment named "lariat_mapper" by running

	mamba create --name "lariat_mapper" --file requirements.txt --channel conda-forge --channel bioconda

The environment can then be activated by running

	mamba activate lariat_mapper
 
before running scripts in the pipeline.

### Reference files
LariatMapper needs a set of reference files to run. 

You will need:
- A FASTA file of the reference genome (`GENOME_FASTA`)
- A GTF or GFF file of annotations for the reference genome (`GENOME_ANNO`)
- The name of the attribute for a feature's gene ID in the annotation file (`G_ATTR`, defaults to "gene_id")
- The name of the attribute for a feature's transcript ID in the annotation file (`T_ATTR`, defaults to "transcript_id")
- A hisat2 index of the reference genome (`HISAT2_INDEX`)

Run `build_references.py` with the paths to each file and the desired output path:

	python build_references.py -f GENOME_FASTA -a <GENOME_ANNO> -i <HISAT2_INDEX> -o <OUT_DIR> -g <G_ATTR> -t <T_ATTR>

You can then use `OUT_DIR` as the reference files directory when running the pipeline (argument `-r, --ref_dir`)

If you have a BED file of repetitive regions from RepeatMasker, you can include the argument `-r REPEATMASKER_BED` to copy it to the reference directory. You can find RepeatMasker files for several reference genomes on the UCSC Genome Browser (https://genome.ucsc.edu/cgi-bin/hgTables) in group "Repeats", track "RepeatMasker". 

If a RepeatMasker file is included in a run, LariatMapper will check putative lariat alignments for false positives which arise from repetitive regions. If the 5' splice site and branchpoint are both located in a reptitive region, the alignment will be filtered out.


## Running the Pipeline
### Required arguments 
For single-end sequencing data, run

	python larmap.py -f READ_FILE -r REF_DIR -o OUTPUT_DIR

For paired-end sequencing data, run

	python larmap.py -1 READ_ONE -2 READ_TWO -r REF_DIR -o OUTPUT_DIR

LariatMapper accepts FASTQ-format files, uncompressed or gzip-compressed. The data should be preprocessed to remove low-quality reads, adapter sequences, and unique molecular identifiers (UMIs) for reliable results. 


### All arguments
```
Input read files:
  Provide either two paired-end read files or one single-end read file. Files can be uncompressed or gzip-compressed
	-1 READ_ONE, --read_one READ_ONE
                        Read 1 input FASTQ file when processing paired-end RNA-seq data. Mutually exclusive with -f
	-2 READ_TWO, --read_two READ_TWO
                        Read 2 input FASTQ file when processing paired-end RNA-seq data. Mutually exclusive with -f
	-f READ_FILE, --read_file READ_FILE
                        Input FASTQ file when processing single-end RNA-seq data. Mutually exclusive with -1 and -2

Reference files:
	-r REF_DIR, --ref_dir REF_DIR
                        Directory with reference files created by build_references.py

Output:
	-o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory for output files. Will be created if it does not exist

Optional arguments:
	-s {Unstranded,First,Second}, --strand {Unstranded,First,Second}
                        WARNING, EXPERIMENTAL FEATURE STILL IN DEVELOPMENT! Strandedness of the input reads. Choices: Unstranded = Library preparation wasn't strand-specific; First = READ_ONE/READ_FILE reads match the RNA sequence (i.e. 2nd cDNA synthesis strand); Second = READ_ONE/READ_FILE reads are reverse-complementary to the RNA sequence (i.e. 1st cDNA synthesis strand) (Default = Unstranded)
	-m REF_REPEATMASKER, --ref_repeatmasker REF_REPEATMASKER
                        BED file of repetitive regions in the genome. Putative lariats that map to a repetitive region will be filtered out as false positives (Default = REF_DIR/repeatmasker.bed if it's an existing file, otherwise skip repetitive region filtering
	-i REF_H2INDEX, --ref_h2index REF_H2INDEX
                        hisat2 index of the reference genome (Default = REF_DIR/hisat2_index)
	-g REF_FASTA, --ref_fasta REF_FASTA
                        FASTA file of the reference genome (Default = REF_DIR/genome.fa)
	-5 REF_5P_FASTA, --ref_5p_fasta REF_5P_FASTA
                        FASTA file with sequences of first 20nt of annotated introns (Default = REF_DIR/fivep_sites.fa)
	-n REF_INTRONS, --ref_introns REF_INTRONS
                        TSV file of all annotated introns (Default = REF_DIR/introns.tsv.gz)
	-p OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Add a prefix to output file names (-o OUT -p ABC -> OUT/ABC_lariat_reads.tsv)
	-u, --ucsc_track      
  						Add an output file named "lariat_reads.bed" which can be used as a custom track in the UCSC Genome Browser (https://www.genome.ucsc.edu/cgi-bin/hgCustom) to visualize lariat alignments
	-c, --keep_classes    
						Keep a file with per-read classification named "read_classes.tsv.gz" in the output (Default = delete)
	-k, --keep_temp       
						Keep all temporary files created while running the pipeline. Forces -c/--keep_classes (Default = delete)
	-t THREADS, --threads THREADS
                        Number of threads to use for parallel processing (Default = 1)
	-x, --skip_version_check
                        Don't check if LariatMapper is up-to-date with the main branch on GitHub (Default = check and warn if not up-to-date)
	-q, --quiet           
						Only print fatal error messages (sets logging level to ERROR, Default = INFO). Mutually exclusive with -w and -d
	-w, --warning         
						Print warning messages and fatal error messages (sets logging level to WARNING, Default = INFO). Mutually exclusive with -q and -d
	-d, --debug           
						Print extensive status messages (sets logging level to DEBUG, Default = INFO). Mutually exclusive with -q and -w
```

## Output
All output will be written in the directory `OUT_DIR`. This includes:
	
- `lariat_reads.tsv`: A table of lariat reads and their alignments
- `circularized_intron_reads.tsv`: A table of circularized intron reads and their alignments
- `template_switching_reads.tsv`: A table of template-switching reads and their alignments
- `summary.txt`: A collection of metadata and summary statistics for the run

`lariat_reads.tsv` and  columns:

- read_id: The read's ID (unique)
- gene_id*: The ID of the gene(s) that produced the lariat
- chrom: The chromosome of the gene that produced the lariat
- strand: The strand of the gene that produced the lariat. "+" for the forward strand and "-" for the reverse strand
- fivep_pos: The genomic position of the lariat's 5' splice site 
- bp_pos: The genomic position of the branchpoint 
- threep_pos: The genomic position of the closest 3' splice site that is downstream of the branchpoint
- bp_dist_to_threep: The distance of the branchpoint from the 3' splice site in nucleotides. Always negative
- read_is_reverse: "False" if the RNA-seq read's sequence matches the source intron's sequence, "True" if it matches the reverse-complement
- read_bp_pos: The position of the branchpoint in the read
- read_seq: The read's sequence 
- read_bp_nt: The nucleotide of the branchpoint according to the read. Reverse-complemented if read_is_reverse is "True"
- genomic_bp_nt: The nucleotide of the branchpoint according to the reference genome. Reverse-complemented if strand is "-"
- genomic_bp_context: The genomic sequence from positions -8 to +8 of the branchpoint. Reverse-complemented if strand is "-"
- total_mapped_reads: The number of input reads that mapped linearly to the reference genome (identical across all rows)
*May be multiple comma-delimited values

`circularized_intron_reads.tsv` columns:
- read_id: The read's ID (unique)
- chrom: The intron's chromosome of the intron
- strand: The intron's strand. "+" for the forward strand and "-" for the reverse strand
- fivep_pos: The genomic position of the intron's 5' splice site 
- head_end_pos: The genomic position of the read head's end
- threep_pos: The genomic position of intron's 3' splice site 
- head_dist_to_threep: The distance of the read head's end from the 3' splice site in nucleotides. Always negative or 0
- read_is_reverse: "False" if the RNA-seq read's sequence matches the source intron's sequence, "True" if it matches the reverse-complement
- read_seq: The read's sequence 
- read_head_end_pos: The position of the read head's end in the read
- read_head_end_nt: The nucleotide of the read head's end according to the read. Reverse-complemented if read_is_reverse is "True"
- head_end_nt: The nucleotide of the read head's end according to the reference genome. Reverse-complemented if strand is "-"
- head_end_context: The genomic sequence from positions -8 to +8 of the read head's end. Reverse-complemented if strand is "-"
- gene_id*: The intron's gene ID(s)
*May be multiple comma-delimited values


`template_switching_reads.tsv` columns:

- read_id: The read's ID (unique)
- fivep_sites*: The 5' splice site(s) that mapped to the read. Format is \<chromosome>;\<strand>;\<position>,...
- temp_switch_sites*: The possible location(s) to which the reverse transcriptase transfered. Format is \<chromosome>;\<position>,...
- read_seq*: The read's sequence(s)
- fivep_seq*: The 5' splice sites' sequence(s)
- genomic_bp_context*: The genomic sequence(s) from positions -8 to +8 of the branchpoint. Reverse-complemented if strand is "-"
- read_bp_pos*: The position(s) of the branchpoint(s) in the read 

*May be multiple comma-delimited values

All position values are 0-based and inclusive. 

If a lariat's 5' splice site and branchpoint could be attributed to multiple gene annotations, the gene values will appear like so:

	gene_name	gene_id	gene_type
	PCBP2,ENSG00000257379	ENSG00000197111.16,ENSG00000257379.1	lncRNA,protein_coding


## Additional information
See DEMO.md for a demonstration of a basic LariatMapper run.

See DESIGN.md for an overiview of LariatMapper's design and the theory behind it.