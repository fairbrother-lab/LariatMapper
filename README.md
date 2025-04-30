# LariatMapper (beta) | The Fairbrother Lab
## Overview

A pipeline for extracting lariats and their branchpoints from RNA-sequencing data.

NOTE: LariatMapper is currently in development. If you encounter any problems while using it, please let us know by [creating an issue on GitHub](https://github.com/fairbrother-lab/LariatMapper/issues/new?template=bug-report.md).



## Installation
LariatMapper can be installed in a Linux or macOS system. It will not work in Windows.

You can download it in a command line terminal with [git](https://git-scm.com/) by running

```
git clone https://github.com/fairbrother-lab/LariatMapper
```
and then 
```
cd LariatMapper
```
to enter the cloned directory.


## Setup
### Dependencies
LariatMapper's software dependencies are detailed in `requirements.txt`. Please note that the dependencies are named according to their package name in [conda](https://docs.conda.io/en/latest/), which may differ from their official name (for example, ggplot2 is named `r-ggplot2`). 
All the dependencies can be obtained from the conda channels `conda-forge` or `bioconda`. We **STRONGLY** recommend creating a dedicated environment for LariatMapper with [mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) to avoid dependency-related problems. 

With mamba installed, you can create a new environment named "larmap" by running
```
mamba create --name "larmap" --file requirements.txt --channel conda-forge --channel bioconda
```

The environment can then be activated by running
```
mamba activate larmap
```
before using LariatMapper.


### Reference files
LariatMapper needs a set of reference files to process RNA-seq data. You will have to build a "reference directory" that contains these files from an appropriate reference genome.

You will need:
- A FASTA file of the reference genome (`GENOME_FASTA`)
- A GTF or GFF file of annotations for the reference genome (`GENOME_ANNO`)
- A [HISAT2 index](https://daehwankimlab.github.io/hisat2/) of the reference genome (`HISAT2_INDEX`)
- An output path for the reference directory (`REF_DIR`) 

With these items determined, run
```
python build_references.py -f GENOME_FASTA -a GENOME_ANNO -i HISAT2_INDEX -o REF_DIR
```

You can then use `REF_DIR` as the reference directory when running LariatMapper.

`build_references.py` creates symbolic links to the input files by default. To copy the input files to `REF_DIR` instead, you can use the argument `-c, --copy`.

If you have a BED file of repetitive regions from RepeatMasker, you can use the argument `-r, --repeatmasker_bed` to include it in the reference directory. You can find RepeatMasker files for several reference genomes on the [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTables) in group "Repeats", track "RepeatMasker". If a RepeatMasker file is included in a run, LariatMapper will check putative lariat alignments for false positives which arise from repetitive regions. If the 5' splice site and branchpoint are both located in a reptitive region, the alignment will be filtered out.

To see all arguments for `build_references.py`, run
```
python build_references.py -h
```


## Running the Pipeline
### Required arguments 
You will need:
- A FASTQ file from single-end RNA-seq (`READ_FILE`) **OR** A pair of FASTQ files from paired-end RNA-seq (`READ_ONE`, `READ_TWO`)
- A reference directory (`REF_DIR`)
- An output path (`OUTPUT_DIR`)

For single-end RNA-seq data, run
```
python larmap.py -f READ_FILE -r REF_DIR -o OUTPUT_DIR
```
For paired-end RNA-seq data, run
```
python larmap.py -1 READ_ONE -2 READ_TWO -r REF_DIR -o OUTPUT_DIR
```
The FASTQ files can be uncompressed or gzip-compressed. The RNA-seq data should be preprocessed to remove low-quality reads, adapter sequences, and unique molecular identifiers for reliable results.

To see all arguments for `larmap.py`, run
```
python larmap.py -h
```


### Optional arguments
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
-d, --debug
                     Print extensive status messages. Mutually exclusive with -q and -w
```

</details>



## Additional information
See DEMO.md for a short demonstration of setting up and then running LariatMapper.

See DETAILS.md for more detailed information about LariatMapper's design and use.



## Contact us
You can reach us via email at jeremiah_buerer@brown.edu.

You can also [create an issue on GitHub](https://github.com/fairbrother-lab/LariatMapper/issues/new) to provide feedback, suggest changes, or report bugs.
