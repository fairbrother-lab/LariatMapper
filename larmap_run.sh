#!/bin/bash

usage() {
    echo ""
    echo "Options:"
    echo "  -d, --fastq_dir <fastq_dir>               directory of FASTQ files"
    echo "  -1, --read_1_file <read_1_file>           read 1 file"
    echo "  -2, --read_2_file <read_2_file>           read 2 file"
    echo "  -o, --output_dir <output_dir>             output directory"
    echo "  -e, --output_base_name <output_base_name> output base name"
    echo "  -c, --num_cpus <num_cpus>                 number of CPUs to use"
    echo "  -i, --ref_b2index <ref_b2index>           reference B2 index file"
    echo "  -f, --ref_fasta <ref_fasta>               reference FASTA file"
    echo "  -g, --ref_gtf <ref_gtf>                   reference GTF file"
    echo "  -5, --ref_5p_fasta <ref_5p_fasta>         FASTA file of first 20nt of 5'ss"
    echo "  -u, --ref_5p_upstream <ref_5p_upstream>   Custom file of sequences in 5nt window upstream of 5'ss"
    echo "  -3, --ref_3p_b2index <ref_3p_b2index>     reference 3' B2 index file"
    echo "  -l, --ref_3p_lengths <ref_3p_lengths>     reference 3' lengths file"
    echo "  -n, --ref_introns <ref_introns>           reference introns file"
    echo "  -m, --ref_repeatmasker <ref_repeatmasker> reference repeatmasker file"
    echo ""
    exit 1
}

exit_abnormal() {                         
    usage
    exit 1
}


# https://stackoverflow.com/questions/402377/using-getopts-to-process-long-and-short-command-line-options

while getopts :d:1:2:o:e:c:i:f:g:5:u:3:l:n:m:-: opt; do        
    case $opt in                    
        d) fastq_dir=$OPTARG ;;
        1) read_one_file=$OPTARG ;;
        2) read_two_file=$OPTARG ;;
        o) output_dir=$OPTARG ;;
        e) output_base_name=$OPTARG ;;
        c) num_cpus=$OPTARG ;;
        i) ref_b2index=$OPTARG ;;
        f) ref_fasta=$OPTARG ;;
        g) ref_gtf=$OPTARG ;;
        5) ref_5p_fasta=$OPTARG ;;
        u) ref_5p_upstream=$OPTARG ;;
        3) ref_3p_b2index=$OPTARG ;;
        l) ref_3p_lengths=$OPTARG ;;
        n) ref_introns=$OPTARG ;;
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
                output_dir)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    output_dir=$val  ;;
                output_base_name)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    output_base_name=$val  ;;
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
                ref_5p_upstream)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_5p_upstream=$val ;;
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
echo "output_dir: $output_dir"
echo "output_base_name: $output_base_name"
echo "num_cpus: $num_cpus"
echo "ref_b2index: $ref_b2index"
echo "ref_fasta: $ref_fasta"
echo "ref_gtf: $ref_gtf"
echo "ref_5p_fasta: $ref_5p_fasta"
echo "ref_5p_upstream: $ref_5p_upstream"
echo "ref_3p_b2index: $ref_3p_b2index"
echo "ref_3p_lengths: $ref_3p_lengths"
echo "ref_introns: $ref_introns"
echo "ref_repeatmasker: $ref_repeatmasker"
echo ""

# Check if all required arguments are provided
if [[ -z $fastq_dir || -z $read_one_file || -z $read_two_file || -z $output_dir || -z $output_base_name || -z $num_cpus || -z $ref_b2index || -z $ref_fasta || -z $ref_gtf || -z $ref_5p_fasta || -z ref_5p_upstream || -z $ref_3p_b2index || -z $ref_3p_lengths || -z $ref_introns || -z $ref_repeatmasker ]]; then
  printf "All arguments are required.\n"
  exit_abnormal
fi

printf "$(date +'%m/%d/%y - %H:%M:%S') | Starting lariat mapping run...\n"
script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $script_dir
#=============================================================================#
#                                    MAPPING                                  #
#=============================================================================#

# prepares the directories and scripts for the lariat mapping run
printf "$(date +'%m/%d/%y - %H:%M:%S') | Preparing directories...\n"
output_dir=$output_dir/$output_base_name"_lariat_mapping"
mkdir -p $output_dir/$output_base_name"_R1"
mkdir -p $output_dir/$output_base_name"_R2"

printf "$(date +'%m/%d/%y - %H:%M:%S') | Processing read one file...\n"
SECONDS=0
scripts/map_lariats.sh $fastq_dir/$read_one_file \
    $output_dir/$output_base_name"_R1/"$output_base_name"_R1" \
    $num_cpus $ref_b2index $ref_fasta $ref_gtf \
    $ref_5p_fasta $ref_5p_upstream \
    $ref_3p_b2index $ref_3p_lengths

printf "$(date +'%m/%d/%y - %H:%M:%S') | Processing read two file...\n"
SECONDS=0
scripts/map_lariats.sh $fastq_dir/$read_two_file \
    $output_dir/$output_base_name"_R2/"$output_base_name"_R2" \
    $num_cpus $ref_b2index $ref_fasta $ref_gtf \
    $ref_5p_fasta $ref_5p_upstream \
    $ref_3p_b2index $ref_3p_lengths

printf "$(date +'%m/%d/%y - %H:%M:%S') | Filtering results...\n"
# combines the mapping results from each sample's read one and read two files and performs post-mapping filtering before outputting the final lariat mapping results
python -u scripts/filter_lariats.py $fastq_dir $read_one_file $read_two_file $output_dir $output_base_name $num_cpus $ref_b2index $ref_fasta $ref_gtf $ref_5p_fasta $ref_3p_b2index $ref_3p_lengths $ref_introns $ref_repeatmasker
exit_code=$?
# Check the exit code and handle errors
if [ $exit_code -ne 0 ]; then
    printf "Error: Failed to execute merge_filter_lariats_no_info.py. Exit code: $exit_code"
    exit $exit_code
fi

printf "$(date +'%m/%d/%y - %H:%M:%S') | Finished.\n"
