# LariatMapper (beta) | The Fairbrother Lab

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

### Options
#### Putative branchpoint correction
Content

#### Specifying refs
Content

#### All options
<details>
	<summary>Click me</summary>
</details>


## Output
All output will be written in the directory `OUT_DIR`. This includes:

- `lariat_reads.tsv`: A table of lariat reads and their alignments
- `circularized_intron_reads.tsv`: A table of circularized intron reads and their alignments
- `template_switching_reads.tsv`: A table of template-switching reads and their alignments
- `summary.txt`: A collection of metadata and summary statistics for the run

<details>
	<summary>Click me</summary>
`lariat_reads.tsv` and  columns:
</details>

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

<details>
	<summary>Click me</summary>
`circularized_intron_reads.tsv` columns:
</details>
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

<details>
	<summary>Click me</summary>
`template_switching_reads.tsv` columns:
</details>
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

## Contact


## Software attributions
<details>
	<summary>Click me</summary>
- **bedtoolsr**: Patwardhan, Mayura; Wenger, Craig D.; Davis, Eric S.; Phanstiel, Douglas H.: "Bedtoolsr: An R package for genomic data analysis and manipulation" (in preparation).
- **Biostrings**: Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). Biostrings: Efficient manipulation of biological strings. R package version 2.74.1, https://bioconductor.org/packages/Biostrings.
- **Bowtie2**: Langmead, B., Salzberg, S. Fast gapped-read alignment with Bowtie 2. Nat Methods 9, 357–359 (2012). https://doi.org/10.1038/nmeth.1923
- **fsspec**: https://github.com/fsspec/filesystem_spec
- **GenomicAlignments**, **GenomicFeatures**, **GenomicRanges**: Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan M, Carey V (2013). “Software for Computing and Annotating Genomic Ranges.” PLoS Computational Biology, 9. doi:10.1371/journal.pcbi.1003118, http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003118. 
- **ggplot2**: Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org.
- **gt**: Iannone R, Cheng J, Schloerke B, Hughes E, Lauer A, Seo J, Brevoort K, Roy O (2025). gt: Easily Create Presentation-Ready Display Tables. R package version 0.11.1.9000, https://github.com/rstudio/gt, https://gt.rstudio.com.
- **HISAT2**: Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019). https://doi.org/10.1038/s41587-019-0201-4
- **intervaltree**: https://github.com/chaimleib/intervaltree
- **magrittr**: Bache S, Wickham H (2022). magrittr: A Forward-Pipe Operator for R. https://magrittr.tidyverse.org, https://github.com/tidyverse/magrittr.  
- **NumPy**: Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2. (Publisher link).
- **optparse**: Davis T (2024). _optparse: Command Line Option Parser_. R package version 1.7.5, <https://CRAN.R-project.org/package=optparse>.
- **pandas**: The pandas development team. (2024). pandas-dev/pandas: Pandas (v2.2.3). Zenodo. https://doi.org/10.5281/zenodo.13819579
- **pyfaidx**: Shirley MD, Ma Z, Pedersen B, Wheelan S. Efficient "pythonic" access to FASTA files using pyfaidx. PeerJ PrePrints 3:e1196. 2015. 
- **pysam**: https://github.com/pysam-developers/pysam
- **rtracklayer**: Lawrence M, Gentleman R, Carey V (2009). “rtracklayer: an R package for interfacing with genome browsers.” Bioinformatics, 25, 1841-1842. doi:10.1093/bioinformatics/btp328, http://bioinformatics.oxfordjournals.org/content/25/14/1841.abstract.
- **tidyverse**:   Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” _Journal of Open Source Software_, *4*(43), 1686. doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>.
- **samtools**: Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li, Twelve years of SAMtools and BCFtools, GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008

LariatMapper was developed from an in-house analysis pipeline which was first publicized in "Large-scale mapping of branchpoints in human
pre-mRNA transcripts in vivo" by Taggart et al. (2012). 
</details>
## Additional information
See DEMO.md for a demonstration of a basic LariatMapper run.

See DESIGN.md for an overiview of LariatMapper's design and the theory behind it.