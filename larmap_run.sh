#!/bin/bash

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
RESET='\033[0m' 

usage() {
    echo ""
    echo "Options:"
    echo "  -d, --fastq_dir <fastq_dir>               directory of FASTQ files"
    echo "  -1, --read_1_file <read_1_file>           read 1 file"
    echo "  -2, --read_2_file <read_2_file>           read 2 file"
    echo "  -e, --output_base_name <output_base_name> output base name"
    echo "  -c, --num_cpus <num_cpus>                 number of CPUs to use"
    echo "  -i, --ref_b2index <ref_b2index>           reference B2 index file"
    echo "  -f, --ref_fasta <ref_fasta>               reference FASTA file"
    echo "  -g, --ref_gtf <ref_gtf>                   reference GTF file"
    echo "  -5, --ref_5p_fasta <ref_5p_fasta>         reference 5' FASTA file"
    echo "  -3, --ref_3p_b2index <ref_3p_b2index>     reference 3' B2 index file"
    echo "  -l, --ref_3p_lengths <ref_3p_lengths>     reference 3' lengths file"
    echo "  -n, --ref_introns <ref_introns>           reference introns file"
    echo "  -m, --ref_repeatmasker <ref_repeatmasker> reference repeatmasker file"
    echo "  -s, --results_filename <results_filename> file to store mapping results"
    echo ""
    exit 1
}

exit_abnormal() {                         
    usage
    exit 1
}


# https://stackoverflow.com/questions/402377/using-getopts-to-process-long-and-short-command-line-options

while getopts :d:1:2:e:s:c:i:f:g:5:3:l:n:m:-: opt; do        
    case $opt in                    
        d) fastq_dir=$OPTARG    ;;
        1) read_one_file=$OPTARG  ;;
        2) read_two_file=$OPTARG  ;;
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
            case "${OPTARG}" in
                fastq_dir)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    fastq_dir=$val  ;;
                read_1_file)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    read_one_file=$val  ;;
                read_2_file)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    read_two_file=$val  ;;
                output_base_name)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    output_base_name=$val  ;;
                results_path)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    results_path=$val  ;;
                num_cpus)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    num_cpus=$val  ;;
                ref_b2index)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_b2index=$val  ;;
                ref_fasta)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_fasta=$val ;;
                ref_gtf)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_gtf=$val ;;
                ref_5p_fasta)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_5p_fasta=$val ;;
                ref_3p_b2index)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_3p_b2index=$val ;;
                ref_3p_lengths)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_3p_lengths=$val ;;
                ref_introns)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_introns=$val ;;
                ref_repeatmasker)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_repeatmasker=$val ;;
                *)
                    exit_abnormal   ;;
            esac;;
        :)                                   
            echo "Error: -${OPTARG} requires an argument."
            exit_abnormal   ;;
        *)                                   
            exit_abnormal   ;;
    esac
done

echo ""
echo "fastq_dir: $fastq_dir"
echo "read_1_file: $read_one_file"
echo "read_2_file: $read_two_file"
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
echo ""

# Check if all required arguments are provided
if [[ -z $fastq_dir || -z $read_one_file || -z $read_two_file || -z $output_base_name || -z $results_path || -z $num_cpus || -z $ref_b2index || -z $ref_fasta || -z $ref_gtf || -z $ref_5p_fasta || -z $ref_3p_b2index || -z $ref_3p_lengths || -z $ref_introns || -z $ref_repeatmasker ]]; then
  echo "All arguments are required."
  exit_abnormal
fi

NOW="$(date +'%m.%d.%Y-%H.%M')"
LARMAP_OUTPUT_DIR="larmap_out"
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
