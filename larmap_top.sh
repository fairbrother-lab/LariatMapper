#!/bin/bash

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
RESET='\033[0m' 

#=============================================================================#
#                                 INPUT SPECS                                 #
#=============================================================================#
fastq_dir=$1
read_1_a=$2
read_2_a=$3
read_1_b=$4
read_2_b=$5
output_base_name_a=$6
output_base_name_b=$7
num_cpus=$8

#=============================================================================#
#                               REFERENCE SPECS                               #
#=============================================================================#
ref_b2index=$9
ref_fasta=${10}
ref_gtf=${11}
ref_5p_fasta=${12}	
ref_3p_b2index=${13}
ref_3p_lengths=${14}
ref_introns=${15}
ref_repeatmasker=${16}	

#=============================================================================#
#                                OUTPUS SPECS                                 #
#=============================================================================#
results_filename=${17}

#=============================================================================#
#                              GENERATE INFO FILE                             #
#=============================================================================#
NOW="$(date +'%m.%d.%Y-%H.%M')"
LARMAP_OUTPUT_DIR="larmap_out"
INFO_FILE="$LARMAP_OUTPUT_DIR/info.$NOW.txt"
INFO="fastq_dir"$'\t'"$fastq_dir""
scripts_dir"$'\t'"$LARMAP_OUTPUT_DIR/scripts""
log_dir"$'\t'"$LARMAP_OUTPUT_DIR/logs""
output_dir"$'\t'"$LARMAP_OUTPUT_DIR/output""
results_path"$'\t'"$LARMAP_OUTPUT_DIR/$results_filename""
num_cpus"$'\t'"$num_cpus""

ref_b2index"$'\t'"$ref_b2index""
ref_fasta"$'\t'"$ref_fasta""
ref_gtf"$'\t'"$ref_gtf""
ref_5p_fasta"$'\t'"$ref_5p_fasta""
ref_3p_b2index"$'\t'"$ref_3p_b2index""
ref_3p_lengths"$'\t'"$ref_3p_lengths""
ref_introns"$'\t'"$ref_introns""
ref_repeatmasker"$'\t'"$ref_repeatmasker""

read_one_file"$'\t'"read_two_file"$'\t'"output_base_name""
$read_1_a"$'\t'"$read_2_a"$'\t'"$output_base_name_a""
$read_1_b"$'\t'"$read_2_b"$'\t'"$output_base_name_b"""
mkdir -p "$LARMAP_OUTPUT_DIR"
touch "$INFO_FILE"
echo "$INFO" > "$INFO_FILE"

#=============================================================================#
#                                    MAPPING                                  #
#=============================================================================#
START_TIME=$(date +"%Y-%m-%d %H:%M:%S")
echo -e "* ${RED}Started: ${GREEN}$START_TIME ${RESET}"

echo -e "* ${YELLOW}Generating scripts...${RESET}"
# prepares the directories and scripts for the lariat mapping run
python scripts/map_lariats_top_no_slurm.py $INFO_FILE
exit_code=$?
# Check the exit code and handle errors
if [ $exit_code -ne 0 ]; then
    echo "Error: Failed to execute map_lariats_top_no_slurm.py. Exit code: $exit_code"
    exit $exit_code
fi

echo -e "* ${YELLOW}Mapping...${RESET}"
# Run all read mapping scripts. Each larmap*.sh script will output logs to its own log file within the log_dir/child_logs directory.
./"$LARMAP_OUTPUT_DIR"/scripts/bash_all.sh &
exit_code=$?
# Check the exit code and handle errors
if [ $exit_code -ne 0 ]; then
    echo "Error: Failed to execute mapping shell scripts. Exit code: $exit_code"
    exit $exit_code
fi

wait

echo -e "* ${YELLOW}Integrating results...${RESET}"
# combines the mapping results from each sample's read one and read two files and performs post-mapping filtering before outputting the final lariat mapping results
python scripts/merge_filter_lariats_no_slurm.py $INFO_FILE
exit_code=$?
# Check the exit code and handle errors
if [ $exit_code -ne 0 ]; then
    echo "Error: Failed to execute merge_filter_lariats_no_slurm.py. Exit code: $exit_code"
    exit $exit_code
fi

END_TIME=$(date +"%Y-%m-%d %H:%M:%S")
echo -e "* ${RED}Finished: ${GREEN}$END_TIME ${RESET}"



