import sys
import os
import shutil
import subprocess



#=============================================================================#
#                                  Constants                                  #
#=============================================================================#



#=============================================================================#
#                                  Functions                                  #
#=============================================================================#
# def bowtie2_index(ref_fasta:str, output_dir:str) -> None:
# 	command = f'bowtie2 {ref_fasta} {output_dir}/ref_genome'.split(' ')
# 	subprocess.run(command, check=True)


# def ref_5p(ref_fasta:str, ref_gtf:str, output_dir:str) -> None:



# def ref_introns(ref_gtf:str, output_dir:str) -> None:



# def ref_transcripts(ref_gtf:str, output_dir:str) -> None:
	





#=============================================================================#
#                                    Main                                     #
#=============================================================================#
def main():
	ref_fasta, ref_gtf, ref_repeatmasker, output_dir = sys.argv[1:]
	ref_fasta = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.gencode.v44.primary_assembly.fa'
	ref_gtf = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.gencode.v44.basic.annotation.gtf'
	ref_repeatmasker = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references/hg38.repeat_masker.bed.gz'
	output_dir = '/Users/trumanmooney/Documents/GitHub/lariat_mapping/testing/references_test'

	# Make dir
	os.mkdir(output_dir)

	# Copy input reference files to dir
	shutil.copyfile(ref_fasta, f'{output_dir}/')
	shutil.copyfile(ref_gtf, f'{output_dir}/')
	shutil.copyfile(ref_repeatmasker, f'{output_dir}/')

	# Create remaining reference files
	bowtie2_index(ref_fasta, output_dir)
	ref_5p(ref_fasta, ref_gtf, output_dir)





if __name__ == '__main__':
	main()



# input        
	#     ref_fasta=$OPTARG
	#     ref_gtf=$OPTARG 
	#     ref_repeatmasker=$OPTARG 

# output
#    ref_b2index=$OPTARG 
        #     ref_5p_fasta=$OPTARG 
        #     ref_5p_upstream=$OPTARG 
        #     ref_introns=$OPTARG 
        #     ref_transcripts=$OPTARG 


