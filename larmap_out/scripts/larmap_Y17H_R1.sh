#!/bin/bash -l

SECONDS=0 
LOGFILE_PATH=larmap_out/logs/child_logs/larmap_Y17H_R1.out 
printf "Log file path = $LOGFILE_PATH\n"
source ~/anaconda3/etc/profile.d/conda.sh 
if [[ "$(conda config --show envs_dirs)" != *"- larmap_env"* ]]; then 
	conda config --append envs_dirs larmap_env 
fi 
conda activate larmap_env 
printf "conda environment = $CONDA_PREFIX\n"

mkdir -p larmap_out/output/Y17H_R1 

printf "Beginning Y17H read one lariat mapping...\n"
scripts/map_lariats.sh \
demo_files/demo_fastq_files_250k_bp/250k_cDBR1-Y17H_1.fq.gz \
larmap_out/output/Y17H_R1/Y17H_R1 \
4 \
demo_files/genomes/indices/bowtie2/mm39.fa \
demo_files/genomes/fasta_files/mm39.fa \
demo_files/genomes/annotations/mm39.gencode.basic.M32.annotation.gtf.gz \
demo_files/reference_files/mouse/mm39.gencode.basic.M32.fivep_sites.fa \
demo_files/reference_files/mouse/mm39.gencode.basic.M32.threep_sites.fa \
demo_files/reference_files/mouse/mm39.gencode.basic.M32.threep_seq_lens.txt \
scripts/filter_fivep_alignments.py \
scripts/filter_threep_alignments.py \
scripts/timestamp.sh 
printf "Mapping of Y17H read one file complete\n"

