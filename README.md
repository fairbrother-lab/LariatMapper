# LariatMapper (beta) | The Fairbrother Lab
## Overview

A pipeline for extracting lariats and their branchpoints from RNA-sequencing data.

NOTE: LariatMapper is currently in development. If you encounter any problems while using it, please let us know by [creating an issue on GitHub](https://github.com/fairbrother-lab/LariatMapper/issues/new?template=bug-report.md).


## Table of Contents
- [Installation](#installation)
- [Setup](#setup)
	- [Dependencies](#dependencies)
	- [Reference files](#reference-files)
- [Running the Pipeline](#running-the-pipeline)
	- [Required arguments](#required-arguments)
	- [Putative branchpoint correction](#putative-branchpoint-correction)
	- [All options](#all-options)
- [Output](#output)
	- [Read classes](#read-classes)
	- [Output files](#output-files)
	- [Output file columns](#output-file-columns)
- [Additional information](#additional-information)
- [Software attributions](#software-attributions)


## Installation
LariatMapper works in a Linux or MacOS system. It will not work in Windows. You can download it on the command line with `git`:

```
git clone https://github.com/fairbrother-lab/LariatMapper
```


## Setup
### Dependencies
The software dependencies are detailed in `requirements.txt`. We recommend creating a dedicated programming environment for LariatMapper with [mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) to avoid dependency-related problems during use. 

If you have mamba installed, you can create a new environment named "larmap" by running
```
mamba create --name "larmap" --file requirements.txt --channel conda-forge --channel bioconda
```

The environment can then be activated by running
```
mamba activate larmap
```
before running scripts in the pipeline.


### Reference files
LariatMapper needs a set of reference files to run. 

You will need:
- A FASTA file of the reference genome (`GENOME_FASTA`)
- A GTF or GFF file of annotations for the reference genome (`GENOME_ANNO`)
- A hisat2 index of the reference genome (`HISAT2_INDEX`)

Run `build_references.py` with the paths to each file and the desired output path:

```
python build_references.py -f GENOME_FASTA -a GENOME_ANNO -i HISAT2_INDEX -o OUT_DIR
```

You can then use `OUT_DIR` as the reference files directory when running LariatMapper (argument `-r, --ref_dir`)

If you have a BED file of repetitive regions from RepeatMasker, you can include the argument `-r REPEATMASKER_BED` to copy it to the reference directory. You can find RepeatMasker files for several reference genomes on the [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTables) in group "Repeats", track "RepeatMasker". 

If a RepeatMasker file is included in a run, LariatMapper will check putative lariat alignments for false positives which arise from repetitive regions. If the 5' splice site and branchpoint are both located in a reptitive region, the alignment will be filtered out.


## Running the Pipeline
### Required arguments 
For single-end sequencing data, run
```
python larmap.py -f READ_FILE -r REF_DIR -o OUTPUT_DIR
```
For paired-end sequencing data, run
```
python larmap.py -1 READ_ONE -2 READ_TWO -r REF_DIR -o OUTPUT_DIR
```
LariatMapper accepts FASTQ-format files, uncompressed or gzip-compressed. The data should be preprocessed to remove low-quality reads, adapter sequences, and unique molecular identifiers for reliable results. 


### Branchpoint correction
LariatMapper includes an option to try correcting the apparent branchpoint positions of lariat reads to account for sequencing errors. When applied, LariatMapper will check a 3nt window downstream of the head's end position (the apparent branchpoint) to see if any of them are more likely to be the true branchpoint. 

LariatMapper can use different 2 approaches to branchpoint correction:

<details>
<summary><strong>Position Weight Matrix (PWM)-based correction.</strong></summary>

This uses a PWM of the expected branchpoint motif, including a mark to indicate the location of the branchpoint within the PWM. The PWM is matched against the genomic sequence of the apparent branchpoint and each other position within the 3nt window. If a position in the window gets a higher match score than the apparent branchpoint position and the score is 0.8 or greater, it is deemed the correct branchpoint position. Multiple PWMs can be used, in which case the highest-scoring position across all PWMs is chosen. 

See <ZENODO_LINK_TO_BE_ADDED> to download prebuilt matrices for commonly-used reference genomes. See `scripts/pwm_build.R` to build a custom matrix from a FASTA file of branchpoint sequences.

Applied with `-P, --pwm_correction`

</details>

<details>
<summary><strong>Model-based correction</summary></strong>

This uses predictions from DeepEnsemble, a deep-learning-based branchpoint prediction model. A 3nt window downstream of the apparent branchpoint is checked, and if the apparent branchpoint is *not* predicted to be a branchpoint but a position within the 3nt window *is*, the later is deemed the correct branchpoint. 

See <ZENODO_LINK_TO_BE_ADDED> to download precalculated predictions for commonly-used reference genomes. 

Applied with `-M, --model_correction`

</details>


Including branchpoint correction in a LariatMapper run will add the following columns to `lariat_reads.tsv`:

 - `corrected_bp_pos`: The genomic position of the corrected branchpoint position. Identical to `bp_pos` if `corrected` = `False`
 - `corrected_bp_dist_to_threep`: The genomic position of the closest 3' splice site that is downstream of the branchpoint, using the corrected position. Identical to `bp_dist_to_threep` if `corrected` = `False`
 - `corrected_genomic_bp_nt`: The nucleotide of the branchpoint in the genome, using the corrected position. Reverse-complemented if `strand` = `-`. Identical to `genomic_bp_nt` if `corrected` = `False`
 - `corrected_genomic_bp_context`: The genomic sequence from positions -8 to +8 of the branchpoint, using the corrected position. Reverse-complemented if `strand` = `-`. Identical to `genomic_bp_context` if `corrected` = `False`
 - `corrected_bp_mismatch `: `True` if `read_bp_nt` ≠ `corrected_genomic_bp_nt`, otherwise `False`. Identical to `bp_mismatch` if `corrected` = `False`
 - `corrected`: `True` if a more likely candidate for the branchpoint position was discovered, otherwise `False`




### All options
<details>
<summary> Expand </summary>

```
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
					RDS file with a position weight matrix to correct apparent branchpoint positions. Multiple files can be provided in comma-seperated format. Mutually exclusive with --model_correction. See <ZENODO_LINK_TO_BE_ADDED> to download prebuilt matrices. See scripts/pwm_build.R to build a custom matrix (Default = no correction)
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
```

</details>

## Output
### Read classes
LariatMapper aims to assign each input read a "read class", which denotes the type of RNA template from which the read originated. Read classes are mututally exclusive, and can be divided into 2 groups:

1. The read *does* have a valid linear alignment to the genome, and...
	- **Linear, exon-exon junction**: ... contains an exon-exon splice junction
	- **Linear, exon-intron junction**: ... contains and contains an exon-intron junction.
	- **Linear, exon only**: ... only overlaps exons.
	- **Linear, intron only**: ... only overlaps introns.
	- **Linear, intergenic or ambiguous**: ... is within either zero genes or multiple overlapping genes.
2. The read does *NOT* have a valid linear alignment to the genome, and...
	- **No alignment**: ... no 5' splice sites map to it.
	- **Fivep alignment**: ... at least one 5' splice site maps to it. 
	- **Template-switching**: ... contains a reverse transcriptase template-switching event.
	- **Circularized intron**: ... has an inverted gapped alignment with a tail that starts at a 5' splice site, and a head that ends within 2nt of a downstream 3' splice site.
	- **Lariat**: ... has an inverted gapped alignment with a tail that starts at a 5' splice site, and a head that's within a intron and downstream of the tail.
	- **In repetitive region**: ... has a lariat alignment to a repetitive region


### Output files
All output will be written in the directory `OUT_DIR`. This includes:

- `lariat_reads.tsv`: A table of the reads classified as lariats**
- `circularized_intron_reads.tsv`: A table of the reads classified as circularized introns
- `template_switching_reads.tsv`: A table of reads that contain contains a reverse transcriptase template-switching event
- `output.bam_count.tsv`: A table of linearly-mapped read counts for each gene, stratified by read class
- `summary.txt`: A collection of metadata and summary statistics for the run
- `read_counts.tsv`: A table of counts for various read classes in machine-friendly format
- `plots/Branchpoint_base_composition.png`: A plot of the distribution of branchpoint nucleotides across all lariat reads
- `plots/Branchpoint_threep_distance.png`: A plot of the distribution of distances between the branchpoint and the 3' splice site across all lariat reads


### Output file columns 
<details>
<summary><code> lariat_reads.tsv </code></summary>

- `read_id`: The read's ID (unique)
- `gene_id`<sup>*</sup>: The Ensembl ID of the gene that produced the lariat
- `chrom`: The chromosome of the gene that produced the lariat
- `strand`: The strand of the gene that produced the lariat. `+` for the forward strand and `-` for the reverse strand
- `fivep_pos`: The genomic position of the lariat's 5' splice site 
- `bp_pos`: The genomic position of the lariat's branchpoint 
- `threep_pos`: The genomic position of the closest 3' splice site that is downstream of the branchpoint
- `bp_dist_to_threep`: The distance of the branchpoint to the 3' splice site
- `read_orient_to_gene`:  `Forward` if the RNA-seq read's sequence matches the source intron's sequence, `Reverse` if it matches the reverse-complement
- `read_seq_forward`: The read's DNA sequence. Reverse-complemented if `read_orient_to_gene` = `Reverse`
- `read_bp_pos`: The position of the branchpoint in the read
- `read_bp_nt`: The nucleotide of the branchpoint in the read
- `genomic_bp_nt`: The nucleotide of the branchpoint in the genome. Reverse-complemented if `strand` = `-`
- `bp_mismatch`: `True` if `read_bp_nt` ≠ `genomic_bp_nt`, otherwise `False`
- `genomic_bp_context`: The genomic sequence from positions -8 to +8 of the branchpoint. Reverse-complemented if `strand` = `-`
- `total_mapped_reads`: The total number of input reads that mapped linearly to the reference genome (identical across all rows)

<sup>*</sup> can be multiple comma-delimited values
</details>

<details>
<summary><code> circularized_intron_reads.tsv </code></summary>

- `read_id`: The read's ID (unique)
- `gene_id`<sup>*</sup>: The Ensembl ID of the gene that produced the lariat
- `chrom`: The chromosome of the gene that produced the lariat
- `strand`: The strand of the gene that produced the lariat. `+` for the forward strand and `-` for the reverse strand
- `fivep_pos`: The genomic position of the lariat's 5' splice site 
- `head_end_pos`: The genomic position of the end of the head alignment
- `threep_pos`: The genomic position of the closest 3' splice site that is downstream of the branchpoint
- `head_end_dist_to_threep`: The distance of the end of the head alignment to the 3' splice site
- `read_orient_to_gene`:  `Forward` if the RNA-seq read's sequence matches the source intron's sequence, `Reverse` if it matches the reverse-complement
- `read_seq_forward`: The read's DNA sequence. Reverse-complemented if `read_orient_to_gene` = `Reverse`
- `read_head_end_pos`: The end position of the head in the read
- `read_head_end_nt`: The nucleotide of the end position of the head in the read
- `genomic_head_end_nt`: The nucleotide of the end of the head alignment in the genome. Reverse-complemented if `strand` = `-`
- `genomic_head_end_context`: The genomic sequence from positions -8 to +8 of the end of the head alignment. Reverse-complemented if `strand` = `-`

<sup>*</sup> can be multiple comma-delimited values
</details>

<details>
<summary><code> template_switching_reads.tsv </code></summary>

- `read_id`: The read's ID (unique)
- `read_orient_to_gene`:  `Forward` if the RNA-seq read's sequence matches the source intron's sequence, `Reverse` if it matches the reverse-complement
- `read_seq_forward`: The read's DNA sequence. Reverse-complemented if `read_orient_to_gene` = `Reverse`
- `read_bp_pos`: The position of the branchpoint in the read
- `fivep_seq`<sup>*</sup>: The 5' splice sites' DNA sequence. Reverse-complemented if `strand` = `-`
- `fivep_sites`<sup>*</sup>: The 5' splice site(s) that mapped to the read. Format is `CHROMOSOME;STRAND;POSITION`
- `genomic_bp_context`: The genomic sequence from positions -8 to +8 of the branchpoint. Reverse-complemented if `strand` = `-`
- `temp_switch_sites`<sup>*</sup>: The genomic location to which the reverse transcriptase transfered. Format is `CHROMOSOME;POSITION`

<sup>*</sup> can be multiple comma-delimited values
</details>

<details>
<summary><code> output.bam_count.tsv </code></summary>

- `gene_id`: The gene's ID (unique)
- `gene`: The total number of reads mapped to the gene. 
- `exon_intron_junc`: The number of reads mapped to the gene that contain an exon-intron junction 
- `exon_exon_junc`: The number of reads mapped to the gene that contain an exon-exon splice junction 
- `exon_only`: The number of reads mapped to the gene that only overlap exons
- `intron_only`: The number of reads mapped to the gene that only overlap introns

`gene` = `exon_intron_junc`+`exon_only`+`exon_exon_junc`+`intron_only`

</details>

All position values are 0-based and inclusive. 


## Additional information
See DEMO.md for a short demonstration of setting up and then running LariatMapper.

See DESIGN.md for an overiview of LariatMapper's design and the theory behind it.


## Contact us
You can reach us via email at 

- truman_mooney@brown.edu 
- jeremiah_buerer@brown.edu


## Software attributions

- **bedtoolsr**: Patwardhan, Mayura; Wenger, Craig D.; Davis, Eric S.; Phanstiel, Douglas H.: "Bedtoolsr: An R package for genomic data analysis and manipulation" (in preparation).

- **Biostrings**: Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). Biostrings: Efficient manipulation of biological strings. https://bioconductor.org/packages/Biostrings.

- **Bowtie2**: Langmead, B., Salzberg, S. Fast gapped-read alignment with Bowtie 2. Nat Methods 9, 357–359 (2012). https://doi.org/10.1038/nmeth.1923

- **GenomicAlignments**, **GenomicFeatures**, **GenomicRanges**: Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan M, Carey V (2013). “Software for Computing and Annotating Genomic Ranges.” PLoS Computational Biology, 9. doi:10.1371/journal.pcbi.1003118, http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003118. 

- **HISAT2**: Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019). https://doi.org/10.1038/s41587-019-0201-4

- **pyfaidx**: Shirley MD, Ma Z, Pedersen B, Wheelan S. Efficient "pythonic" access to FASTA files using pyfaidx. PeerJ PrePrints 3:e1196. 2015. 

- **rtracklayer**: Lawrence M, Gentleman R, Carey V (2009). “rtracklayer: an R package for interfacing with genome browsers.” Bioinformatics, 25, 1841-1842. doi:10.1093/bioinformatics/btp328, http://bioinformatics.oxfordjournals.org/content/25/14/1841.abstract.

- **tidyverse**:   Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” _Journal of Open Source Software_, *4*(43), 1686. doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>.

- **txdbmaker**: Pagès H, Carlson M, Aboyoun P, Falcon S, Morgan M (2024). txdbmaker: Tools for making TxDb objects from genomic annotations. https://bioconductor.org/packages/txdbmaker

- **samtools**: Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li, Twelve years of SAMtools and BCFtools, GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008


Development of the original method for alignment of lariat reads from RNA-seq data was pioneered by the Fairbrother lab ([Taggart *et al.* 2012](https://doi.org/10.1038/nsmb.2327)). The current LariatMapper algorithm was extended and refined from an in-house analysis pipeline ([Townley *et al.* 2023](https://doi.org/10.1016/j.molcel.2023.06.011)) based on the method described in [Pineda and Bradley 2018](https://doi.org/10.1101/gad.312058.118).