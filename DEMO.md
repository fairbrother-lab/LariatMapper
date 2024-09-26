# Demo
Here, we take a look at running LariatMapper via the command line. To follow this demonstration, open a Linux or MacOS command line terminal. As you read, copy the blocks of code into the terminal and run them.  

First, assign `larmap_dir` the path to the LariatMapper package directory:

	larmap_dir="<ENTER_PATH_HERE>"

For example, `larmap_dir="/home/me/bioinfo_packages/LariatMapper"`. This is the only code you'll have to edit before running – the rest of the code will work as-is.  

Next, set up a dedicated conda environment that contains all of LariatMapper's dependencies:

	demo_env="larmap_demo"

	conda create --name "$demo_env" --file "$larmap_dir/requirements.txt" --channel conda-forge --channel bioconda

and then activate it:

	conda activate "$demo_env"

Now create a directory with all the neccesary reference data by calling `build_references.py` 

	genome_fasta=$larmap_dir/resources/demo/hg38.demo.fa.gz
	genome_anno=$larmap_dir/resources/demo/hg38.demo.gtf.gz
	hisat2_index=$larmap_dir/resources/demo/hg38.demo.index
	repeatmasker_bed=$larmap_dir/resources/demo/hg38.demo.repeat_masker.bed.gz
	ref_dir=demo_references

	python "$larmap_dir/build_references.py" -f $genome_fasta -a $genome_anno -i $hisat2_index -r $repeatmasker_bed -o $ref_dir

Finally, input the example RNA-sequencing data into LariatMapper by calling `larmap.py`

	r1_reads=$larmap_dir/resources/demo/demo_reads_R1.fastq.gz
	r2_reads=$larmap_dir/resources/demo/demo_reads_R2.fastq.gz
	output_dir=./demo_output

	python "$larmap_dir/larmap.py" -r $ref_dir -1 $r1_reads -2 $r2_reads -o $output_dir

You should now have a directory named `demo_output` in your current working directory. It contains the standard output files that LariatMapper will produce. 