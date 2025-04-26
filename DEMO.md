# Demo
## Overview
Here, we look at a typical example of setting up and running LariatMapper. 

To follow this demonstration, open a Linux or MacOS command line terminal. As you read, copy the blocks of code into the terminal and run them. You will need to have [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) and [git](https://git-scm.com/) installed.

If you haven't downloaded LariatMapper yet, download it with git:

```
git clone https://github.com/fairbrother-lab/LariatMapper
```

## Running the Demo
First, download the demo data at https://github.com/fairbrother-lab/LariatMapper_aux with git:

```
git clone https://github.com/fairbrother-lab/LariatMapper_aux
```

Next, assign the path to the LariatMapper directory to the variable `larmap_dir`, and the path to the demo directory to the variable `demo_dir`. For example:

```
larmap_dir="/home/me/bioinformatics/LariatMapper"
demo_dir="/home/me/bioinformatics/LariatMapper_aux/demo"
```

This is the only code you'll have to edit before running – the rest of the code can be copied and pasted as-is.  

Next, create a mamba environment with LariatMapper's software dependencies:

```
mamba create -y -n LariatMapper --file "$larmap_dir/requirements.txt" -c conda-forge -c bioconda
```

and then activate it:

```
mamba activate LariatMapper
```

Next, create a directory with all of the neccessary reference files by calling `build_references.py`

```
genome_fasta="$demo_dir/reference_genome/genome.fa.gz"
genome_anno="$demo_dir/reference_genome/annotation.gtf.gz"
hisat2_index="$demo_dir/reference_genome/hisat2_index"
repeatmasker_bed="$demo_dir/reference_genome/repeatmasker.bed.gz"
ref_dir="./LariatMapper_ref_data"

python "$larmap_dir/build_references.py" -f "$genome_fasta" -a "$genome_anno" -i "$hisat2_index" -r "$repeatmasker_bed" -o "$ref_dir"
```

That should have created a directory named `LariatMapper_ref_data` in your working directory.

Finally, input the RNA-sequencing data and newly-created reference directory into LariatMapper by calling `larmap.py`

```
r1_reads="$demo_dir/sequencing_reads_R1.fastq.gz"
r2_reads="$demo_dir/sequencing_reads_R2.fastq.gz"
output_dir="./LariatMapper_output"

python "$larmap_dir/larmap.py" -r "$ref_dir" -1 "$r1_reads" -2 "$r2_reads" -o "$output_dir"
```

You should now have a directory named `LariatMapper_output` in your working directory containing the standard output files that LariatMapper produces. You can open the file named `summary.txt` to see metadata for the run and various summary statistics.

<br></br>

## The Demo Data
The files in `$demo_dir/reference_genome` were generated from the "Genome sequence, primary assembly (GRCh38)" FASTA file and "Comprehensive gene annotation" GTF file in [GENCODE Human release 44](https://www.gencodegenes.org/human/release_44.html). They only include chromosomes 20, 21, and 22 in order to reduce file sizes and processing time. 

The `sequencing_reads` files were generated from an in-house total RNA-seq experiment with HEK293T cells. They contain 100,000 reads from the original FASTQ file of ~53 million reads, which were semi-randomly sampled so that the distribution of output reads in the demo would be similar to the distribution observed when mapping the full FASTQ file against the full reference genome. At least 1 read was included for each output read class.