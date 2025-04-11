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
	- [All optional arguments](#all-optional-arguments)
- [Concepts](#concepts)
	- [Head and tail](#head-and-tail)
	- [Read classes](#read-classes)
- [Output](#output)
	- [Output files](#output-files)
	- [Output file columns](#output-file-columns)
- [Additional information](#additional-information)
- [Contact us](#contact-us)
- [Software attributions](#software-attributions)


## Installation
LariatMapper can be installed in a Linux or macOS system. It does not work in Windows. 

You can download it in a command line terminal with [git](https://git-scm.com/):

```
git clone https://github.com/fairbrother-lab/LariatMapper
```


## Setup
### Dependencies
The software dependencies are detailed in `requirements.txt`. We recommend creating a dedicated programming environment for LariatMapper with [mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) to avoid dependency-related problems during use. 

If you have mamba installed, you can create a new environment named "larmap" by running
```
mamba create --name "larmap" --file LariatMapper/requirements.txt --channel conda-forge --channel bioconda
```

The environment can then be activated by running
```
mamba activate larmap
```
before using LariatMapper.


### Reference files
LariatMapper needs a set of reference files to run.

You will need:
- A FASTA file of the reference genome (`GENOME_FASTA`)
- A GTF or GFF file of annotations for the reference genome (`GENOME_ANNO`)
- A hisat2 index of the reference genome (`HISAT2_INDEX`)

Run `build_references.py` with the paths to each file and the desired output path (`REF_DIR`):

```
python build_references.py -f GENOME_FASTA -a GENOME_ANNO -i HISAT2_INDEX -o REF_DIR
```

`build_references.py` creates symbolic links to the input files by default. To copy the input files to `REF_DIR` instead, you can use the argument `-c, --copy`.

You can then use `REF_DIR` as the reference files directory when running LariatMapper.

If you have a BED file of repetitive regions from RepeatMasker, you can use the argument `-r, --repeatmasker_bed` to copy it to the reference directory. You can find RepeatMasker files for several reference genomes on the [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTables) in group "Repeats", track "RepeatMasker". 

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
LariatMapper includes an option to try to correct the apparent branchpoint position of lariat reads to account for sequencing errors. When applied, LariatMapper will check a 3nt window downstream of the head's end position (the apparent branchpoint) to see if any of them are more likely to be the true branchpoint position. If there *is* a more likely branchpoint position, it will record that corrected position. 

Including branchpoint correction in a LariatMapper run will add the following columns to `lariat_reads.tsv`:

 - `corrected_bp_pos`: The genomic position of the corrected branchpoint. Identical to `bp_pos` if `corrected` = `False`
 - `corrected_bp_dist_to_threep`: The genomic position of the closest 3' splice site that is downstream of the branchpoint, using the corrected position. Identical to `bp_dist_to_threep` if `corrected` = `False`
 - `corrected_genomic_bp_nt`: The nucleotide of the branchpoint in the genome, using the corrected position. Reverse-complemented if `strand` = `-`. Identical to `genomic_bp_nt` if `corrected` = `False`
 - `corrected_genomic_bp_context`: The genomic sequence from positions -8 to +8 of the branchpoint, using the corrected position. Reverse-complemented if `strand` = `-`. Identical to `genomic_bp_context` if `corrected` = `False`
 - `corrected_bp_mismatch `: `True` if `read_bp_nt` ≠ `corrected_genomic_bp_nt`, otherwise `False`. Identical to `bp_mismatch` if `corrected` = `False`
 - `corrected`: `True` if a more likely candidate for the branchpoint position was discovered, otherwise `False`

LariatMapper can use different 2 approaches to branchpoint correction:

**1. Position Weight Matrix (PWM)-based correction**

This uses a PWM of the expected branchpoint motif, including a marker for the location of the branchpoint within the PWM. The PWM is matched against the genomic sequence of the apparent branchpoint and each other position within the 3nt window. If a position in the window gets a higher match score than the apparent branchpoint position and the score is 0.8 or greater, it is deemed the correct branchpoint position. Multiple PWMs can be used, in which case the highest-scoring position across all PWMs is chosen. 

See https://doi.org/10.5281/zenodo.14735947 to download prebuilt PWMs for commonly-used reference genomes. See `scripts/pwm_build.R` to build a custom PWM from a FASTA file of branchpoint sequences.

Applied with `-P, --pwm_correction`

**2. Model-based correction**

This uses predictions from DeepEnsemble, a deep-learning-based branchpoint prediction model. A 3nt window downstream of the apparent branchpoint is checked, and if the apparent branchpoint is *not* predicted to be a branchpoint but a position within the 3nt window *is*, the later is deemed the correct branchpoint. 

See https://doi.org/10.5281/zenodo.14735947 to download precalculated predictions for commonly-used reference genomes. 

Applied with `-M, --model_correction`


### All optional arguments
<details>
<summary> Expand </summary>

```
-T TEMP_SWITCH_FILTER, --temp_switch_filter TEMP_SWITCH_FILTER
                     Set the parameters of the template-switching filter in the head-filtering step. Format = "N,M", where N is the number of downstream bases to check, and M is the minimum number of matches required to identify an alignment as template-switching. (Default = 2,2)
-m REF_REPEATMASKER, --ref_repeatmasker REF_REPEATMASKER
                     BED file of repetitive regions in the genome. Putative lariats that map to a repetitive region will be filtered out as false positives. May be gzip- compressed. (Default = REF_DIR/repeatmasker.bed if it's an existing file, otherwise skip repetitive region filtering)
-H REF_H2INDEX, --ref_h2index REF_H2INDEX
                     HISAT2 index of the reference genome. (Default = REF_DIR/hisat2_index)
-g REF_FASTA, --ref_fasta REF_FASTA
                     FASTA file of the reference genome. May be gzip- compressed. (Default = REF_DIR/genome.fa.gz)
-5 REF_5P_FASTA, --ref_5p_fasta REF_5P_FASTA
                     FASTA file of 5' splice site sequences, i.e. the first 20nt of all annotated introns. (Default = REF_DIR/fivep_sites.fa)
-e REF_EXONS, --ref_exons REF_EXONS
                     TSV file of all annotated introns. (Default = REF_DIR/exons.tsv.gz)
-i REF_INTRONS, --ref_introns REF_INTRONS
                     TSV file of all annotated introns. (Default = REF_DIR/introns.tsv.gz)
-P PWM_CORRECTION, --pwm_correction PWM_CORRECTION
                     RDS file with a position weight matrix to correct apparent branchpoint positions. Multiple files can be provided in comma-seperated format. Mutually exclusive with --model_correction. See https://doi.org/10.5281/zenodo.14735947 to download prebuilt PWMs. See scripts/pwm_build.R to build a custom matrix (Default = no correction)
-M MODEL_CORRECTION, --model_correction MODEL_CORRECTION
                     RDS file with predictions from DeepEnsemble, a deep- learning-based branchpoint prediction model. Mutually exclusive with --pwm_correction. See https://doi.org/10.5281/zenodo.14735947 to download predictions for specific reference genomes. (Default = no correction)
-p OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                     Add a prefix to output file names (-o OUT -p ABC -> OUT/ABC_lariat_reads.tsv). (Default = no prefix)
-u, --ucsc_track
                     Add an output file named "lariat_reads.bed". This can be used as a custom track in the UCSC Genome Browser to visualize lariat read alignments
-b, --keep_bam
                     Keep the BAM file produced in the initial linear mapping step (Default = delete)
-c, --keep_classes
                     Keep a file with per-read classification of non- linearly-aligned reads named "read_classes.tsv.gz" in the output (Default = delete)
-k, --keep_temp
                     Keep all temporary files created while running the pipeline. Forces -c and -b (Default = delete)
-t THREADS, --threads THREADS
                     Number of threads to use. (Default = 1)
-q, --quiet
                     Only print fatal error messages. Mutually exclusive with -w and -d
-w, --warning
                     Print warning messages and fatal error messages. Mutually exclusive with -q and -d
```

</details>


## Concepts
### Head and tail
![A diagram illustrating RNA-seq reads split into a head and a tail, and linear vs inverted gapped alignments](./resources/images/Linear%20vs%20Inverted%20alignments.svg)

The fundamental goal of LariatMapper is to identify lariat reads by their inverted gapped alignment to the genome. To find those alignments, we split RNA-seq reads into two fragments. For practicality, we name the 5'-end fragment the **head** and the 3'-end fragment the **tail**.

In the diagram above, we illustrate an mRNA read alignment as an example of a linear gapped alignment, with the head at the end of an exon and the tail at the start of a downstream exon. A lariat read alignment is used as an example of an inverted gapped alignment, with the head inside an intron and the tail at an upstream 5' splice site. 


### Steps


### Read classes
LariatMapper assigns a "read class" to all input reads which describes their RNA template. Some classes don't represent any particular type of RNA, but rather the extent of information that LariatMapper was able to obtain before the read got filtered out (e.g. **Fivep alignment**). Read classes can be divided into 2 groups based on whether or not the read aligned linearly to the genome in LariatMapper's first step:

1. The read *does* have a valid linear alignment to the genome, and...
	- **Exon-exon junction**: ... contains an exon-exon splice junction but no exon-intron junctions.
	- **Exon-intron junction**: ... contains an exon-intron junction.
	- **Exon only**: ... only overlaps exons.
	- **Intron only**: ... only overlaps introns.
	- **Intergenic or ambiguous**: ... falls within either no genes or multiple overlapping genes.
2. The read does *NOT* have a valid linear alignment to the genome, and...
	- **No alignment**: ... no 5' splice sites map to it.
	- **Fivep alignment**: ... at least one 5' splice site maps to it.
	- **Template-switching**: ... contains a template-switching event which occurred at the transition from tail to head during reverse transcription
	- **Circularized intron**: ... has an inverted gapped alignment that's characteristic of a lariat, but the apparent branchpoint is within 2nt of the 3' splice site
	- **Lariat**: ... has an inverted gapped alignment that's characteristic of a lariat



## Output
### Output files
All output will be written in the directory `OUT_DIR`. This includes:

- `lariat_reads.tsv`: A table of the reads classified as lariats
- `circularized_intron_reads.tsv`: A table of the reads classified as circularized introns
- `template_switching_reads.tsv`: A table of the reads classified as template-switching
- `output.bam_count.tsv`: A table of linearly-mapped read counts for each gene, stratified by read class
- `read_counts.tsv`: A table of counts for various read classes in machine-friendly format
- `summary.txt`: A collection of metadata and summary statistics for the run
- `plots/Branchpoint_base_composition.png`: A plot of the distribution of branchpoint nucleotides across all lariat reads
- `plots/Branchpoint_threep_distance.png`: A plot of the distribution of distances between the branchpoint and the 3' splice site across all lariat reads


### Output file columns 
<details>
<summary><code> lariat_reads.tsv </code></summary>

- `read_id`: The read's ID (unique)
- `gene_id`<sup>*</sup>: The gene ID of the intron that produced the lariat
- `chrom`: The chromosome of the intron that produced the lariat
- `strand`: The strand of the intron that produced the lariat. `+` for the forward strand and `-` for the reverse strand
- `fivep_pos`: The genomic position of the lariat's 5' splice site 
- `bp_pos`: The genomic position of the lariat's branchpoint 
- `threep_pos`: The genomic position of the closest 3' splice site that is downstream of the branchpoint
- `bp_dist_to_threep`: The distance from the branchpoint to the 3' splice site. Equal to `-|bp_pos - threep_pos|`
- `read_orient_to_gene`: `Forward` if the read's sequence matches the intron's sequence, `Reverse` if it matches the reverse-complement
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
- `gene_id`<sup>*</sup>: The gene ID of the intron
- `chrom`: The chromosome of the intron
- `strand`: The strand of the intron. `+` for the forward strand and `-` for the reverse strand
- `fivep_pos`: The genomic position of the intron's 5' splice site 
- `head_end_pos`: The genomic position of the head's end
- `threep_pos`: The genomic position of the closest 3' splice site that is downstream of the head's end
- `head_end_dist_to_threep`: The distance of the head's end to the 3' splice site. Equal to `-|head_end_pos - threep_pos|`
- `read_orient_to_gene`: `Forward` if the read's sequence matches the intron's sequence, `Reverse` if it matches the reverse-complement
- `read_seq_forward`: The read's DNA sequence. Reverse-complemented if `read_orient_to_gene` = `Reverse`
- `read_head_end_pos`: The end position of the head in the read
- `read_head_end_nt`: The nucleotide at the end position of the head in the read
- `genomic_head_end_nt`: The nucleotide at the end of the head alignment in the genome. Reverse-complemented if `strand` = `-`
- `genomic_head_end_context`: The genomic sequence from positions -8 to +8 of the end of the head alignment. Reverse-complemented if `strand` = `-`

<sup>*</sup> can be multiple comma-delimited values
</details>

<details>
<summary><code> template_switching_reads.tsv </code></summary>

- `read_id`: The read's ID (unique)
- `read_orient_to_gene`<sup>†</sup>: `Forward` if the read's sequence matches the intron's sequence, `Reverse` if it matches the reverse-complement
- `read_seq_forward`<sup>†</sup>: The read's DNA sequence. Reverse-complemented if `read_orient_to_gene` = `Reverse`
- `read_head_end_pos`<sup>†</sup>: The end position of the head in the read
- `fivep_seq`<sup>†</sup>: The 5' splice sites' DNA sequence. Reverse-complemented if `strand` = `-`
- `fivep_sites`<sup>*</sup><sup>†</sup>: The 5' splice site(s) that mapped to the read. Format is `CHROMOSOME;STRAND;POSITION;GENE_ID`, where `GENE_ID` may be multiple values delimited by ampersands (`&`)
- `head_chrom`<sup>†</sup>: The chromosome of the head alignment
- `head_strand`<sup>†</sup>: The strand of the gene overlapping the head alignment
- `head_gene_id`<sup>*</sup> <sup>†</sup>: The gene ID of the gene overlapping the head alignment
- `head_first_pos`<sup>†</sup>: The first position of the head alignment in the genome
- `head_last_pos`<sup>†</sup>: The last position of the head alignment in the genome. 0-based exclusive
- `head_end_pos`<sup>†</sup>: The genomic position of the head's end, where the reverse transcriptase attached during template-switching. Equal to `head_first_pos` if `head_strand` = `-`, or `head_last_pos - 1` if `head_strand` = `+`
- `threep_pos`<sup>†</sup>: The genomic position of the closest 3' splice site that is downstream of the head's end
- `head_end_dist_to_threep`<sup>†</sup>: The distance of the head's end to the 3' splice site. Equal to `|head_end_pos - threep_pos|`
- `genomic_head_end_context`<sup>†</sup>: The genomic sequence from positions -8 to +8 of the end of the head alignment. Reverse-complemented if `strand` = `-`

<sup>*</sup> can be multiple comma-delimited values

<sup>†</sup> can be multiple values delimited by vertical bars (`|`), for reads with multiple head alignments (of equal quality) that were caught by the template-switching filter. If the input read data is from paired-end sequencing, some reads may also have head alignments from both mates that were caught by the template-switching filter. In such cases, there may be 2 values for `read_orient_to_gene`, `read_seq_forward`, and/or `read_head_end_pos`.
</details>

<details>
<summary><code> output.bam_count.tsv </code></summary>

- `gene_id`: The gene's ID (unique)
- `gene`: The total number of reads mapped to the gene. 
- `exon_intron_junc`: The number of exon-intron junction reads mapped to the gene
- `exon_exon_junc`: The number of exon-exon junction reads mapped to the gene
- `exon_only`: The number of exon-only reads mapped to the gene
- `intron_only`: The number of intron-only reads mapped to the gene

`gene` = `exon_intron_junc`+`exon_only`+`exon_exon_junc`+`intron_only`

</details>

All position values are 0-based and inclusive, unless otherwise noted. 





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
