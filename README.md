# LariatMapper (beta) | The Fairbrother Lab
## Overview

A pipeline for extracting lariats and their branchpoints from RNA-sequencing data. 

NOTE: LariatMapper is currently in development, and may display unexpected or erroneous behavior. If you encounter any problems while using it, please let us know by [creating an issue on GitHub](https://github.com/fairbrother-lab/LariatMapper/issues/new?template=bug-report.md).

## Table of Contents
- [Setup](#setup)
	- [Dependencies](#dependencies)
	- [Reference files](#reference-files)
- [Running the Pipeline](#running-the-pipeline)
	- [Required arguments](#required-arguments)
	- [Putative branchpoint correction](#putative-branchpoint-correction)
	- [All options](#all-options)
- [Output](#output)
- [Additional information](#additional-information)
	- [Software attributions](#software-attributions)

## Setup

### Dependencies
The software dependencies are detailed in `requirements.txt`. We recommend creating a dedicated programming environment for LariatMapper with [mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) to avoid dependency-related problems during use.

If you have mamba installed, you can create a new environment named "larmap" by running

	mamba create --name "larmap" --file requirements.txt --channel conda-forge --channel bioconda

The environment can then be activated by running

	mamba activate larmap
 
before running scripts in the pipeline.

### Reference files
LariatMapper needs a set of reference files to run. 

You will need:
- A FASTA file of the reference genome (`GENOME_FASTA`)
- A GTF or GFF file of annotations for the reference genome (`GENOME_ANNO`)
- A hisat2 index of the reference genome (`HISAT2_INDEX`)

Run `build_references.py` with the paths to each file and the desired output path:

	python build_references.py -f GENOME_FASTA -a <GENOME_ANNO> -i <HISAT2_INDEX> -o <OUT_DIR>

You can then use `OUT_DIR` as the reference files directory when running LariatMapper (argument `-r, --ref_dir`)

If you have a BED file of repetitive regions from RepeatMasker, you can include the argument `-r REPEATMASKER_BED` to copy it to the reference directory. You can find RepeatMasker files for several reference genomes on the [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTables) in group "Repeats", track "RepeatMasker". 

If a RepeatMasker file is included in a run, LariatMapper will check putative lariat alignments for false positives which arise from repetitive regions. If the 5' splice site and branchpoint are both located in a reptitive region, the alignment will be filtered out.


## Running the Pipeline
### Required arguments 
For single-end sequencing data, run

	python larmap.py -f READ_FILE -r REF_DIR -o OUTPUT_DIR

For paired-end sequencing data, run

	python larmap.py -1 READ_ONE -2 READ_TWO -r REF_DIR -o OUTPUT_DIR

LariatMapper accepts FASTQ-format files, uncompressed or gzip-compressed. The data should be preprocessed to remove low-quality reads, adapter sequences, and unique molecular identifiers for reliable results. 

### Putative branchpoint correction
*To be added*

### All options
	-m REF_REPEATMASKER, --ref_repeatmasker REF_REPEATMASKER
                        BED file of repetitive regions in the genome. Putative lariats that map to a repetitive region will be filtered out as false positives. (Default = REF_DIR/repeatmasker.bed if it's an existing file, otherwise skip repetitive region filtering)
	-i REF_H2INDEX, --ref_h2index REF_H2INDEX
                        HISAT2 index of the reference genome. (Default = REF_DIR/hisat2_index)
	-g REF_FASTA, --ref_fasta REF_FASTA
                        FASTA file of the reference genome. (Default = REF_DIR/genome.fa)
	-5 REF_5P_FASTA, --ref_5p_fasta REF_5P_FASTA
                        FASTA file of 5' splice site sequences, i.e. the first 20nt of all annotated introns. (Default = REF_DIR/fivep_sites.fa)
	-n REF_INTRONS, --ref_introns REF_INTRONS
                        TSV file of all annotated introns. (Default = REF_DIR/introns.tsv.gz)
		
	-P PWM_CORRECTION, --pwm_correction PWM_CORRECTION
                        RDS file with a position weight matrix to correct apparent branchpoint positions. Multiple files can be provided in comma-seperated format. Mutually exclusive with --model_correction. See scripts/pwm_build.R to build a custom matrix (Default = no correction)
	-M PWM_CORRECTION, --model_correction MODEL_CORRECTION
                        RDS file with predictions from DeepEnsemble, a deep-learning-based branchpoint prediction model. Mutually exclusive with --pwm_correction. See <ZENODO_LINK_TO_BE_ADDED> to download predictions for specific reference genomes. (Default = no correction)
						
	-p OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Add a prefix to output file names (-o OUT -p ABC -> OUT/ABC_lariat_reads.tsv). (Default = no prefix)
						
	-u, --ucsc_track    Add an output file named "lariat_reads.bed". This can be used as a custom track in the UCSC Genome Browser to visualize lariat read alignments

	-b, --keep_bam      Keep the BAM file produced in the initial linear mapping step (Default = delete)
	-c, --keep_classes  Keep a file with per-read classification of non-linearly-aligned reads named "read_classes.tsv.gz" in the output (Default = delete)
	-k, --keep_temp     Keep all temporary files created while running the pipeline. Forces -c and -b (Default = delete)

	-t THREADS, --threads THREADS
                        Number of threads to use. (Default = 1)

	-q, --quiet         Only print fatal error messages. Mutually exclusive with -w and -d
	-w, --warning       Print warning messages and fatal error messages. Mutually exclusive with -q and -d
	-d, --debug         Print extensive status messages. Mutually exclusive with -q and -w


## Output
All output will be written in the directory `OUT_DIR`. This includes:

- `lariat_reads.tsv`: A table of the reads identified as lariats
- `circularized_intron_reads.tsv`: A table of the reads identified as circularized introns
- `template_switching_reads.tsv`: A table of reads identified as containing a reverse-transcriptase template-switching event
- `summary.txt`: A collection of metadata and summary statistics for the run
- `read_counts.tsv`: A table of counts for various read classes in machine-friendly format

<details>
	<summary><code>lariat_reads.tsv</code> columns</summary>
<i>To be added</i>
</details>



<details>
	<summary><code>circularized_intron_reads.tsv</code> columns</summary>
*To be added*
</details>

<details>
	<summary><code>template_switching_reads.tsv</code> columns</summary>
<i>To be added</i>
</details>


All position values are 0-based inclusive. 


## Additional information
See DEMO.md for a demonstration of a basic LariatMapper run.

See DESIGN.md for an overiview of LariatMapper's design and the theory behind it.

<details>
	<summary><strong>Software attributions</strong></summary>

- <strong>bedtoolsr</strong>: Patwardhan, Mayura; Wenger, Craig D.; Davis, Eric S.; Phanstiel, Douglas H.: "Bedtoolsr: An R package for genomic data analysis and manipulation" (in preparation).

- <strong>Biostrings</strong>: Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). Biostrings: Efficient manipulation of biological strings. https://bioconductor.org/packages/Biostrings.

- <strong>Bowtie2</strong>: Langmead, B., Salzberg, S. Fast gapped-read alignment with Bowtie 2. Nat Methods 9, 357–359 (2012). https://doi.org/10.1038/nmeth.1923

- <strong>GenomicAlignments</strong>, **GenomicFeatures**, **GenomicRanges**: Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan M, Carey V (2013). “Software for Computing and Annotating Genomic Ranges.” PLoS Computational Biology, 9. doi:10.1371/journal.pcbi.1003118, http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003118. 

- <strong>HISAT2</strong>: Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019). https://doi.org/10.1038/s41587-019-0201-4

- <strong>pyfaidx</strong>: Shirley MD, Ma Z, Pedersen B, Wheelan S. Efficient "pythonic" access to FASTA files using pyfaidx. PeerJ PrePrints 3:e1196. 2015. 

- <strong>rtracklayer</strong>: Lawrence M, Gentleman R, Carey V (2009). “rtracklayer: an R package for interfacing with genome browsers.” Bioinformatics, 25, 1841-1842. doi:10.1093/bioinformatics/btp328, http://bioinformatics.oxfordjournals.org/content/25/14/1841.abstract.

- <strong>tidyverse</strong>:   Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” _Journal of Open Source Software_, *4*(43), 1686. doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>.

- <strong>txdbmaker</strong>: Pagès H, Carlson M, Aboyoun P, Falcon S, Morgan M (2024). txdbmaker: Tools for making TxDb objects from genomic annotations. https://bioconductor.org/packages/txdbmaker

- <strong>samtools</strong>: Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li, Twelve years of SAMtools and BCFtools, GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008

LariatMapper was developed from an in-house analysis pipeline which was first publicized in "Large-scale mapping of branchpoints in human pre-mRNA transcripts in vivo" by Taggart et al. (2012). 
</details>