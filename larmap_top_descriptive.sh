#!/bin/bash

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
RESET='\033[0m' 

usage() {
  echo "Usage: $(basename "$0") -d <fastq_dir> -a <read_1_file> -b <read_2_file> -e <output_base_name> -c <num_cpus> -i <ref_b2index> -f <ref_fasta> -g <ref_gtf> -5 <ref_5p_fasta> -3 <ref_3p_b2index> -l <ref_3p_lengths> -n <ref_introns> -m <ref_repeatmasker> -s <results_filename>"
  exit 1
}

exit_abnormal() {                         
  usage
  exit 1
}


while getopts :d:a:b:e:c:i:f:g:5:3:l:n:m:s:-: opt; do        
    case $opt in                    
        d) fastq_dir=$OPTARG    ;;
        a) read_one_file=$OPTARG  ;;
        b) read_two_file=$OPTARG  ;;
        e) output_base_name=$OPTARG ;;
        s) results_path=$OPTARG ;;
        c) num_cpus=$OPTARG ;;
        i) ref_b2index=$OPTARG  ;;
        f) ref_fasta=$OPTARG    ;;
        g) ref_gtf=$OPTARG  ;;
        5) ref_5p_fasta=$OPTARG ;;
        3) ref_3p_b2index=$OPTARG   ;;
        l) ref_3p_lengths=$OPTARG   ;;
        n) ref_introns=$OPTARG  ;;
        m) ref_repeatmasker=$OPTARG ;;
        -)
            ;;
        :)                                   
            echo "Error: -${OPTARG} requires an argument."
            exit_abnormal   ;;
        *)                                   
            exit_abnormal   ;;
    esac
done

echo "fastq_dir: $fastq_dir"
echo "read_one_file: $read_one_file"
echo "read_two_file: $read_two_file"
echo "output_base_name: $output_base_name"
echo "results_path: $results_path"
echo "num_cpus: $num_cpus"
echo "ref_b2index: $ref_b2index"
echo "ref_fasta: $ref_fasta"
echo "ref_gtf: $ref_gtf"
echo "ref_5p_fasta: $ref_5p_fasta"
echo "ref_3p_b2index: $ref_3p_b2index"
echo "ref_3p_lengths: $ref_3p_lengths"
echo "ref_introns: $ref_introns"
echo "ref_repeatmasker: $ref_repeatmasker"


# Check if all required arguments are provided
if [[ -z $fastq_dir || -z $read_one_file || -z $read_two_file || -z $output_base_name || -z $results_path || -z $num_cpus || -z $ref_b2index || -z $ref_fasta || -z $ref_gtf || -z $ref_5p_fasta || -z $ref_3p_b2index || -z $ref_3p_lengths || -z $ref_introns || -z $ref_repeatmasker ]]; then
  echo "All arguments are required."
  exit_abnormal
fi

NOW="$(date +'%m.%d.%Y-%H.%M')"
LARMAP_OUTPUT_DIR="larmap_out"
# INFO_FILE="$LARMAP_OUTPUT_DIR/info.$NOW.txt"
# INFO="fastq_dir"$'\t'"$fastq_dir""
# scripts_dir"$'\t'"$LARMAP_OUTPUT_DIR/scripts""
# log_dir"$'\t'"$LARMAP_OUTPUT_DIR/logs""
# output_dir"$'\t'"$LARMAP_OUTPUT_DIR/output""
# results_path"$'\t'"$LARMAP_OUTPUT_DIR/$results_filename""
# num_cpus"$'\t'"$num_cpus""

# ref_b2index"$'\t'"$ref_b2index""
# ref_fasta"$'\t'"$ref_fasta""
# ref_gtf"$'\t'"$ref_gtf""
# ref_5p_fasta"$'\t'"$ref_5p_fasta""
# ref_3p_b2index"$'\t'"$ref_3p_b2index""
# ref_3p_lengths"$'\t'"$ref_3p_lengths""
# ref_introns"$'\t'"$ref_introns""
# ref_repeatmasker"$'\t'"$ref_repeatmasker""

# read_one_file"$'\t'"read_two_file"$'\t'"output_base_name""
# $read_1_file"$'\t'"$read_2_file"$'\t'"$output_base_name"""

# mkdir -p "$LARMAP_OUTPUT_DIR"
# touch "$INFO_FILE"
# echo "$INFO" > "$INFO_FILE"

#=============================================================================#
#                                    MAPPING                                  #
#=============================================================================#
START_TIME=$(date +"%Y-%m-%d %H:%M:%S")
echo -e "* ${RED}Started: ${GREEN}$START_TIME ${RESET}"

echo -e "* ${YELLOW}Generating scripts...${RESET}"
# prepares the directories and scripts for the lariat mapping run
python scripts/map_lariats_top_no_info.py $fastq_dir $read_one_file $read_two_file $output_base_name $results_path $num_cpus $ref_b2index $ref_fasta $ref_gtf $ref_5p_fasta $ref_3p_b2index $ref_3p_lengths $ref_introns $ref_repeatmasker

exit_code=$?
# Check the exit code and handle errors
if [ $exit_code -ne 0 ]; then
    echo "Error: Failed to execute map_lariats_top_no_info.py. Exit code: $exit_code"
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
python scripts/merge_filter_lariats_no_info.py $fastq_dir $read_one_file $read_two_file $output_base_name $results_path $num_cpus $ref_b2index $ref_fasta $ref_gtf $ref_5p_fasta $ref_3p_b2index $ref_3p_lengths $ref_introns $ref_repeatmasker
exit_code=$?
# Check the exit code and handle errors
if [ $exit_code -ne 0 ]; then
    echo "Error: Failed to execute merge_filter_lariats_no_info.py. Exit code: $exit_code"
    exit $exit_code
fi

END_TIME=$(date +"%Y-%m-%d %H:%M:%S")
echo -e "* ${RED}Finished: ${GREEN}$END_TIME ${RESET}"
