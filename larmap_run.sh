#!/bin/bash
#=============================================================================#
#                                  Functions                                  #
#=============================================================================#
usage() {
    echo ""
    echo "Required arguments:"
    echo "  -r, --read_file <read_file>               FASTQ file"
    echo "  -o, --output_dir <output_dir>             Directory for output files"
    echo "  -e, --output_base_name <output_base_name> Prefix to add to output files"
    echo "  -c, --num_cpus <num_cpus>                 Number of CPUs available"
    echo "  -i, --ref_b2index <ref_b2index>           Bowtie2 index of the full reference genome"
    echo "  -f, --ref_fasta <ref_fasta>               FASTA file of the full reference genome"
    echo "  -g, --ref_gtf <ref_gtf>                   GTF file with gene annotation of the reference genome"
    echo "  -5, --ref_5p_fasta <ref_5p_fasta>         FASTA file with sequences of first 20nt from reference 5' splice sites (first 20nt of introns)"
    echo "  -u, --ref_5p_upstream <ref_5p_upstream>   Custom file of sequences in 5nt window upstream of 5' splice sites"
    echo "  -n, --ref_introns <ref_introns>           BED file of all introns in the reference genome"
    echo "  -x, --ref_exons <ref_exons>               BED file of all exons in the reference genome"
	echo "  -t, --ref_transcripts <ref_transcripts>   BED file of all transcripts in the reference genome, including the blockCount, blockSizes, and blockStarts columns with exon counts, exon lengths, and exon start positions, respectively"
    echo "  -m, --ref_repeatmasker <ref_repeatmasker> BED file of repetitive elements from repeatmasker"
	echo "Options:"
	echo "  -k, --keep_intermediates				  Don't delete the intermediate files created throughout the pipeline"
    echo ""
    exit 1
}


exit_abnormal() {                         
    usage
    exit 1
}



#=============================================================================#
#                                  Arguments                                  #
#=============================================================================#
# https://stackoverflow.com/questions/402377/using-getopts-to-process-long-and-short-command-line-options
keep_intermediates="False"
echo ""
while getopts :r:o:e:c:i:f:g:5:u:n:t:m:-: opt; do        
    case $opt in                    
        r) 
            read_file=$OPTARG 
            echo "read_file: $read_file" ;;
        o) 
            output_dir=$OPTARG 
            echo "output_dir: $output_dir" ;;
        e) 
            output_base_name=$OPTARG 
            echo "output_base_name: $output_base_name" ;;
        c) 
            num_cpus=$OPTARG 
            echo "num_cpus: $num_cpus" ;;
        i) 
            ref_b2index=$OPTARG 
            echo "ref_b2index: $ref_b2index" ;;
        f) 
            ref_fasta=$OPTARG
            echo "ref_fasta: $ref_fasta" ;; 
        g) 
            ref_gtf=$OPTARG 
            echo "ref_gtf: $ref_gtf" ;;
        5) 
            ref_5p_fasta=$OPTARG 
            echo "ref_5p_fasta: $ref_5p_fasta" ;;
        u) 
            ref_5p_upstream=$OPTARG 
            echo "ref_5p_upstream: $ref_5p_upstream" ;;
        n) 
            ref_introns=$OPTARG 
            echo "ref_introns: $ref_introns" ;;
        x) 
            ref_exons=$OPTARG 
            echo "ref_exons: $ref_exons" ;;
        t) 
            ref_transcripts=$OPTARG 
            echo "ref_transcripts: $ref_transcripts" ;;
        m) 
            ref_repeatmasker=$OPTARG 
            echo "ref_repeatmasker: $ref_repeatmasker" ;;
		k)
			keep_intermediates="True"
			echo "keep_intermediates" ;;
        -) 
            case "${OPTARG}" in
                read_file)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    read_file=$val  
                    echo "read_file: $read_file" ;;
                output_dir)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    output_dir=$val  
                    echo "output_dir: $output_dir" ;;
                output_base_name)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    output_base_name=$val  
                    echo "output_base_name: $output_base_name" ;;
                num_cpus)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    num_cpus=$val  
                    echo "num_cpus: $num_cpus" ;;
                ref_b2index)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_b2index=$val  
                    echo "ref_b2index: $ref_b2index" ;;
                ref_fasta)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_fasta=$val 
                    echo "ref_fasta: $ref_fasta" ;;
                ref_gtf)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_gtf=$val 
                    echo "ref_gtf: $ref_gtf" ;;
                ref_5p_fasta)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_5p_fasta=$val 
                    echo "ref_5p_fasta: $ref_5p_fasta" ;;
                ref_5p_upstream)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_5p_upstream=$val 
                    echo "ref_5p_upstream: $ref_5p_upstream" ;;
                ref_introns)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_introns=$val 
                    echo "ref_introns: $ref_introns" ;;
                ref_exons)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_exons=$val 
                    echo "ref_exons: $ref_exons" ;;
                ref_transcripts)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_transcripts=$val 
                    echo "ref_transcripts: $ref_transcripts" ;;
                ref_repeatmasker)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_repeatmasker=$val 
                    echo "ref_repeatmasker: $ref_repeatmasker" ;;
				keep_intermediates)
					keep_intermediates="True"
					echo "keep_intermediates" ;;
                *)
                    echo "Error: unrecognized option --${OPTARG}."
                    exit_abnormal   ;;
            esac;;
        :)  
            echo ""                                 
            echo "Error: -${OPTARG} requires an argument."
            exit_abnormal   ;;
        *)  
            echo ""
            echo "Error: unrecognized option -${OPTARG}."                              
            exit_abnormal   ;;
    esac
done


# Check if all required arguments are provided
if [[ -z $read_file || -z $output_dir || -z $output_base_name || -z $num_cpus || -z $ref_b2index || -z $ref_fasta || -z $ref_gtf || -z $ref_5p_fasta || -z $ref_5p_upstream || -z $ref_introns || -z $ref_exons || -z $ref_transcripts || -z $ref_repeatmasker ]]; then
  echo "Not all required arguments submitted"
  exit_abnormal
fi



#=============================================================================#
#                                    MAPPING                                  #
#=============================================================================#
printf "$(date +'%m/%d/%y - %H:%M:%S') | Starting lariat mapping run...\n"
script_dir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd $script_dir
printf "\n$(date +'%m/%d/%y - %H:%M:%S') | Preparing directories...\n"
output_dir=$output_dir/$output_base_name"_lariat_mapping"
mkdir -p $output_dir
output_base="$output_dir/"

printf "\n$(date +'%m/%d/%y - %H:%M:%S') | Processing read file...\n"
SECONDS=0
scripts/map_lariats.sh $read_file \
					$output_base \
					$num_cpus \
					$ref_b2index \
					$ref_fasta \
					$ref_gtf \
					$ref_5p_fasta \
					$ref_5p_upstream \
					$ref_transcripts \
					$ref_introns \
					$ref_exons \
					$ref_repeatmasker \
					$keep_intermediates
exit_code=$?
# Check the exit code and handle errors
if [ $exit_code -ne 0 ]; then
    printf "\nError: Failed to execute map_lariats.sh. Exit code: $exit_code"
    exit $exit_code
fi

printf "\n$(date +'%m/%d/%y - %H:%M:%S') | Finished.\n"