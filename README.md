# Lariat mapping pipeline

(to be modified)

## Environment

This pipeline depends on the command line tools and python packages contained in the conda environment at `/datasets2/lariat_mapping/larmap_env` and outline in the `environment.yaml` file in the `lariat_mapping` directory

## Workflow

0.  `pip install -r requirements.txt`

1.  Make a copy of a template info file from the `info_files` directory

2.  Fill out the info file. This file is tab-delimited and contains settings for the pipeline as well as a table describing the fastq files to be processed. The top of the file is populated with values for the variables below:

        fastq_dir - Path to the directory containing the fastq files
        scripts_dir - Path to directory were mapping scripts will be written (created if it does not exist)
        log_dir - Path to directory were mapping logs will be written (created if it does not exist)
        output_dir - Path to directory were mapping output will be written (created if it does not exist)
        results_path - Path to file were the final lariat mapping results will be written
        SB_cpus - Number of CPUs available for the mapping run
        SB_GBmem - Memory in GB available for the mapping run

    After these variables there is information on the reference genome data required to perform the mapping. This information will be constant across runs for a given species' genome build/gene annotation. The required files are:

        ref_b2index - Bowtie2 index of the whole genome
        ref_fasta - Fasta file of the whole genome
        ref_gtf - GTF file providing the gene annotation of the whole genome
        ref_5p_fasta - Fasta file containing the 5' sequences of all the annotated introns (first 20bp of the introns)
        ref_3p_b2index - Bowtie2 index of a custom 3' intronic genome (last 250bp of all annotated introns)
        ref_3p_lengths - The lengths of all the sequences in the custom 3' genome (not all introns are >=250bp, so some lengths will be shorter)
        ref_introns - BED file containing the intron annotation of the whole genome
        ref_repeatmasker - BED file containing the RepeatMasker annotation of the whole genome

    Following these settings is a tab-delimited table with a header that outlines all the fastq files which the run should process. Currently, the pipeline can handle processing of paired-end read data. Each line of the fastq file table must have values in the `read_one_file`, `read_two_file` and `output_base_name` columns. In addition, custom columns can be added after these in order to provide information on the samples being processed (eg. cell line, treatment, replicate number, etc.). Any additional columns will be propagated to the final lariat read results file.

3.  The `map_lariats_top.py` script prepares the directories and scripts for the lariat mapping run. It is called as follows:

        python /datasets2/lariat_mapping/scripts/map_lariats_top.py /path/to/info_file.txt

    This script will read the settings and read file information from the info file and generate scripts in the scripts*dir for mapping each of the sample's read one and read two files.
    Each read's script is titled larmap*[output_base_name]\_[R1/R2].sh. In addition it creates a \_sbatch_all.sh script which will submit all the created scripts for processing by SLURM. map_lariat_top.py outputs it's logs to top_log.out in the log_dir

4.  Run the read mapping scripts either through the individual `larmap*.sh` scripts or with the `\_sbatch_all.sh` script. Each `larmap*.sh` script will output logs to its own log file within the `log_dir/child_logs` directory.

5.  The `merge_filter_lariats.py` script combines the mapping results from each sample's read one and read two files and performs post-mapping filtering before outputting the final lariat mapping results. It is called as follows:

        python /datasets2/lariat_mapping/scripts/merge_filter_lariats.py /path/to/info_file.txt

    The candidate lariat reads from each sample are filtered based on the following criteria: - BP is within 2bp of a splice site (likely from an intron circle, not a lariat) - 5'SS and 3'SS are not in the correct order - Read maps to a Ubiquitin gene (likely false positive due to repetitive nature of gene) - There is a valid aligment for the 3' segment upstream of the 5' segment - Both the 5'SS and the BP overlap with repetitive regions from RepeatMasker (likely false positive) - NEEDS TO BE IMPLEMENTED: Filter out template switching reads (eg. 5'SS 6bp sequence matches 6bp sequence downstream of BP) - NEEDS TO BE IMPLEMENTED: Correct BP position to account for RT skipping

    After filtering, the final set of lariat reads for all the samples are merged and output to the `results_path` file. The initial columns of this tab-delimited file are any additional custom columns which were added to the info file. After these columns, the file contains the following data for each lariat read: - Gene name - Gene Ensembl ID - Gene type (e.g. protein-coding, lncRNA) - Read ID - Read sequence - Chromosome - Strand - 5' splice site coordinate - 3' splice site coordinate - Branchpoint coordinate - Read branchpoint nucleotide - Genomic branchpoint nucleotide - Genomic branchpoint sequence context (10bp window with the BP at position 5) - Branchpoint position relative to nearest 3' splice site - Total count of linearly-aligned reads for the sample the lariat read is from (useful for normalization when comparing lariat levels between samples)

## Pipeline process

### `map_lariats_top.py`

log file: `log_dir/top_log.out`

1. Read the info file and create a RunData object that holds the information
2. Make sure the required directories exist
3. For each specified fastq file, write a bash script in the scripts_dir directory that maps its reads. These scripts are designed to be run with SLURM through the sbatch command
4. Write a `\_sbatch_all.sh` script in the scripts_dir directory that runs an sbatch command on each mapping script

### `larmap*[output_base_name]*[R1/R2].sh` mapping script

log file: `log*dir/child_logs/larmap*[output_base_name]\_[R1/R2]-[SLURM JOB ID].out`

5. Activate the `larmap_env` conda virtual environment
6. Make a directory in the `output*dir` directory named `[output_base_name]*[R1/R2]`
7. Run the `map*lariats.sh` script on the read file. This will produce three files in the read output directory: `[output_base_name]*[R1/R2]_total_reads.txt` (one line file containing count of linearly-aligned reads from read file), `[output_base_name]_[R1/R2]_fivep_info_table.txt` (intermediate file containing info on the mapping of the 5'SS sequences to the unmapped reads), and `[output_base_name]_[R1/R2]\_final_info_table.txt` (results file containing candidate lariat reads obtained after mapping the 5'SS trimmed reads to the 3'SS region sequences)

### `map_lariats.sh`

8. Align reads to the reference genome with bowtie2 and create a BAM file of the unmapped reads with samtools
9. Convert the unmapped reads bam file to FASTA format with samtools
10. Build a bowtie2 index of the unmapped reads
11. Align a fasta file of 5'SS to the unmapped reads index
12. Trim reads with 5'SS alignments and write trimmed reads to fasta file
13. Align the trimmed reads to a Bowtie2 index of 3'SS regions
14. Take the mapped reads from step 6 and create an output file containing candidate lariat reads

### `merge_filter_lariats.py`

log file: `log_dir/merge_filter_log.out`

15. Read the info file and create a RunData object that holds the information
16. Load intron and gene information from provided annotation files
17. Merge the lariat mapping results from each pair of mate read files so that lariats discovered in both mates are only reported once
18. Filter the candidate lariat reads based on the criteria outlined above (Workflow step 5)
19. Output information on the lariat reads that pass filtering to the top level results file provided in the info file
