#!/usr/bin/env nextflow

params.pipeline_dir = '.'
params.read_file = false
params.read_one = false
params.read_two = false
params.strand = 'Unstranded'
params.ref_dir = false
params.output_dir = false
params.output_prefix = false
params.ucsc_track = false
params.keep_classes = false
params.keep_temp = false
params.skip_version_check = true
params.quiet = false
params.warning = false
params.debug = true
params.ref_h2index = false
params.ref_fasta = false
params.ref_5p_fasta = false
params.ref_introns = false
params.ref_repeatmasker = false


process validate_args { 

	input:
	path args_script
	val args_string
	path log_file

	output:
	val true

	script:
	"""
	python3 $args_script $args_string >> $log_file 2>&1
	"""
}


process map_single_end {

	input:
	val ready
	path read_file
	val genome_index
	path output_bam
	val hisat2_strand_arg
	path log_file

	output:
	path output_bam

	script:
	"""
	printf "\$(date +'%d/%b/%Y %H:%M:%S') | Mapping reads to genome...\n" >> $log_file
	hisat2 --no-softclip -k 1 --max-seeds 20 --pen-noncansplice 0 --n-ceil L,0,0.05 --score-min L,0,-0.24 --bowtie2-dp 1 \
	    	$hisat2_strand_arg --threads $task.cpus -x $genome_index -U $read_file \
		| samtools view --bam --with-header --add-flags PAIRED,READ1 \
		| samtools sort --threads $task.cpus --verbosity 0 \
		> $output_bam 2>>$log_file
	"""
}

process map_paired_end {

	input:
	val ready
	path read_one
	path read_two
	val genome_index
	path output_bam
	val hisat2_strand_arg
	path log_file

	output:
	path output_bam

	script:
	"""
	printf "\$(date +'%d/%b/%Y %H:%M:%S') | Mapping reads to genome...\n" >> $log_file
	hisat2 --no-softclip -k 1 --max-seeds 20 --pen-noncansplice 0 --n-ceil L,0,0.05 --score-min L,0,-0.24 --bowtie2-dp 1 \
		   $hisat2_strand_arg --threads $task.cpus -x $genome_index -1 $read_one -2 $read_two \
		| samtools view --bam --with-header \
		| samtools sort --threads $task.cpus --verbosity 0 \
		> $output_bam 2>>$log_file
	"""
}

process index_mapped_reads {

	input:
	path output_bam
	path log_file

	output:
	env unmapped_read_count

	script:
	index_cpus = task.cpus-1
	"""
	samtools index -@ $index_cpus $output_bam >> $log_file 2>&1
	unmapped_read_count=\$(samtools view --count --require-flags 4 $output_bam)
	if [ \$unmapped_read_count == 0 ];then
		printf "\$(date +'%d/%b/%Y %H:%M:%S') | No reads remaining.\n" >> $log_file
	fi
	"""
}

process build_unmapped_reads {

	input:
	val unmapped_read_count
	path output_bam
	path unmapped_fasta
	path log_file

	output:
	val continue_workflow

	script:
	continue_workflow = !unmapped_read_count.toInteger().equals(0)
	if (continue_workflow) {
		"""
		printf "\$(date +'%d/%b/%Y %H:%M:%S') | Creating fasta file of unmapped reads...\n" >> $log_file
		samtools fasta -N --require-flags 4 -o $unmapped_fasta $output_bam >> $log_file 2>&1
		printf "\$(date +'%d/%b/%Y %H:%M:%S') | Indexing unmapped reads fasta file...\n" >> $log_file
		samtools faidx $unmapped_fasta >> $log_file 2>&1
		"""
	} else {
		"""
		"""
	}
}

process index_unmapped_reads {

	input:
	val continue_workflow
	path unmapped_fasta
	val unmapped_index_base
	path log_file

	output:
	val continue_workflow

	script:
	if (continue_workflow) {
		"""
		printf "\$(date +'%d/%b/%Y %H:%M:%S') | Building bowtie2 index of unmapped fasta...\n" >> $log_file
		bowtie2-build --quiet --threads $task.cpus $unmapped_fasta $unmapped_index_base >> $log_file 2>&1
		"""
	} else {
		"""
		"""
	}
}

process align_fivep_sites {

	input:
	val continue_workflow
	val unmapped_index_base
	path fivep_fasta
	path fivep_to_reads
	path log_file

	output:
	val continue_workflow

	script:
	if (continue_workflow) {
		"""
		printf "\$(date +'%d/%b/%Y %H:%M:%S') | Mapping 5' splice sites to reads...\n" >> $log_file
		bowtie2 --end-to-end --sensitive --no-unal -f -k 10000 --score-min C,0,0 --threads $task.cpus -x $unmapped_index_base -U $fivep_fasta \
			| samtools sort --threads $task.cpus --verbosity 0 --output-fmt SAM -M \
			| samtools view \
			> $fivep_to_reads 2>>$log_file
		"""
	} else {
		"""
		"""
	}
}

process filter_fivep_aligns {

	input:
	val continue_workflow
	path fivep_to_reads
	path filter_fivep_script
	val output_base
	path log_file
	val log_level
	path genome_fasta
	path genome_faidx
	path fivep_fasta
	val strand

	output:
	env fivep_to_reads_count

	script:
	if (continue_workflow) {
		"""
		fivep_to_reads_count=\$(cat $fivep_to_reads | wc -l)
		if [ \$fivep_to_reads_count == 0 ];then
			printf "\$(date +'%d/%b/%Y %H:%M:%S') | No reads remaining.\n" >> $log_file
		else
			printf "\$(date +'%d/%b/%Y %H:%M:%S') | Finding 5' read alignments and trimming reads...\n" >> $log_file
			python3 -u $filter_fivep_script $output_base $log_file $log_level \$fivep_to_reads_count $genome_fasta $genome_faidx $fivep_fasta $strand $task.cpus >> $log_file 2>&1
		fi
		"""
	} else {
		"""
		fivep_to_reads_count=0
		"""
	}
}

process align_heads {

	input:
	val fivep_to_reads_count
	val genome_index
	path heads_fasta
	path heads_to_genome
	path log_file

	output:
	val continue_workflow

	script:
	continue_workflow = !fivep_to_reads_count.toInteger().equals(0)
	if (continue_workflow) {
		"""
		printf "\$(date +'%d/%b/%Y %H:%M:%S') | Mapping heads to genome...\n" >> $log_file
		hisat2 --no-softclip --no-spliced-alignment --very-sensitive -k 100 \
			--no-unal --threads $task.cpus -f -x $genome_index -U $heads_fasta \
			| samtools sort --threads $task.cpus --verbosity 0 --output-fmt SAM -n \
			| samtools view \
			> $heads_to_genome 2>>$log_file
		"""
	} else {
		"""
		"""
	}
}

process filter_head_aligns {

	input:
	val continue_workflow
	path heads_to_genome
	path filter_head_script
	path introns_tsv
	path genome_fasta
	val output_base
	path log_file
	val log_level

	output:
	env head_read_count

	script:
	if (continue_workflow) {
		"""
		head_read_count=\$(cat $heads_to_genome | wc -l)
		if [ \$head_read_count == 0 ];then
			printf "\$(date +'%d/%b/%Y %H:%M:%S') | No reads remaining.\n" >> $log_file
		else
			printf "\$(date +'%d/%b/%Y %H:%M:%S') | Analyzing head alignments and outputting lariat table...\n" >> $log_file
			python -u $filter_head_script $task.cpus \$head_read_count $introns_tsv $genome_fasta $output_base $log_file $log_level >> $log_file 2>&1
		fi
		"""
	} else {
		"""
		head_read_count=0
		"""
	}
}

process filter_lariats {

	input:
	val head_read_count
	path filter_lariats_script
	val output_base
	path log_file
	val log_level
	val seq_type
	val repeats_bed

	output:
	val continue_workflow

	script:
	continue_workflow = !head_read_count.toInteger().equals(0)
	if (continue_workflow) {
		"""
		printf "\$(date +'%d/%b/%Y %H:%M:%S') | Filtering putative lariat alignments...\n" >> $log_file
		python -u $filter_lariats_script $output_base $log_file $log_level $seq_type $repeats_bed >> $log_file 2>&1
		"""
	} else {
		"""
		"""
	}
}

process make_track {

	input:
	val continue_workflow
	val output_track
	path track_script
	val output_base
	path log_file
	val log_level

	output:
	val continue_workflow

	script:
	if (continue_workflow & output_track) {
		"""
		printf "\$(date +'%d/%b/%Y %H:%M:%S') | Making UCSC Genome Browser track...\n" >> $log_file
		python -u $track_script $output_base $log_file $log_level >> $log_file 2>&1
		"""
	} else {
		"""
		"""
	}
}

process end_run {

	input:
	val continue_workflow
	path end_script
	path ref_dir
	val output_base
	val keep_classes
	val keep_temp
	val seq_type
	path log_file
	val log_level
	path pipeline_dir

	script:
	"""
	bash $end_script $ref_dir $output_base $keep_classes $keep_temp $seq_type $log_file $log_level $pipeline_dir >> $log_file 2>&1
	"""
}

workflow { 

	// Setup pipeline scripts
	pipeline_dir = file(params.pipeline_dir)
	scripts_dir = pipeline_dir.resolve('scripts')
	args_script = scripts_dir.resolve('validate_args.py')
	filter_fivep_script = scripts_dir.resolve('filter_fivep_aligns.py')
	filter_head_script = scripts_dir.resolve('filter_head_aligns.py')
	filter_lariats_script = scripts_dir.resolve('filter_lariats.py')
	track_script = scripts_dir.resolve('make_track.py')
	end_script = scripts_dir.resolve('end_run.sh')
	
	// Collect user-supplied and default options
	paired_end = params.read_file == false
	validate_args_string = ''
	if (paired_end) {
		validate_args_string = '-1 ' + params.read_one + ' -2 ' + params.read_two
		read_one = file(params.read_one)
		read_two = file(params.read_two)
		seq_type = 'paired'
	} else {
		validate_args_string = '-f ' + params.read_file
		read_file = file(params.read_file)
		seq_type = 'single'
	}
	
	validate_args_string = validate_args_string + ' -r ' + params.ref_dir + ' -o ' + params.output_dir
	
	binary_params = [params.ucsc_track, params.keep_classes, params.keep_temp, params.skip_version_check]
	binary_options = ['-u', '-c', '-k', '-x']
	for (int i = 0; i < binary_options.size(); i++) {
		if (binary_params[i]) {
			validate_args_string = validate_args_string + ' ' + binary_options[i]
		}
	}

	other_params = [params.output_prefix, params.ref_h2index, params.ref_fasta, params.ref_5p_fasta, params.ref_introns, params.ref_repeatmasker]
	other_options = ['-p', '-i', '-g', '-5', '-n', '-m']
	for (int i = 0; i < other_options.size(); i++) {
		if (other_params[i]) {
			validate_args_string = validate_args_string + ' ' + other_options[i] + ' ' + other_params[i]
		}
	}

	hisat2_strand_arg = ''
	if (params.strand) {
		validate_args_string = validate_args_string + ' -s ' + params.strand
		if (params.strand == 'Forward') {
			hisat2_strand_arg = '--rna-strandness F'
			if (paired_end) {
				hisat2_strand_arg = hisat2_strand_arg + 'R'
			}
		} else {
			hisat2_strand_arg = '--rna-strandness R'
			if (paired_end) {
				hisat2_strand_arg = hisat2_strand_arg + 'F'
			}
		}
	}
	
	if (params.output_prefix) {
		output_prefix = params.output_prefix + '_'
	} else {
		output_prefix = ''
	}
	ref_dir = file(params.ref_dir)
	output_dir = file(params.output_dir)
	output_base = params.output_dir + '/' + output_prefix

	log_file = output_dir.resolve(output_prefix + 'mapping_log.txt')
	log_level = 'INFO'
	if (params.debug) {
		validate_args_string = validate_args_string + ' -d'
		log_level = 'DEBUG'
	} else if (params.warning){
		validate_args_string = validate_args_string + ' -w'
		log_level = 'WARNING'
	} else if (params.quiet) {
		validate_args_string = validate_args_string + ' -q'
		log_level = 'ERROR'
	}
	
	// Check that all input files, reference files and other options are valid
	args_done = validate_args(args_script, validate_args_string, log_file)
	
	output_bam = output_dir.resolve(output_prefix + 'output.bam')
	if (params.ref_h2index) {
		genome_index = file(params.ref_h2index)
	} else {
		genome_index = ref_dir.resolve('hisat2_index')
	}
	// Map reads to genome and keep unmapped reads. Lariat reads crossing the brachpoint will not be able to map to the gene they're from
	if (paired_end) {
		output_bam = map_paired_end(args_done, read_one, read_two, genome_index, output_bam, hisat2_strand_arg, log_file)
	} else {
		output_bam = map_single_end(args_done, read_file, genome_index, output_bam, hisat2_strand_arg, log_file)
	}

	unmapped_fasta = output_dir.resolve(output_prefix + 'unmapped_reads.fa')
	// Index mapped reads BAM
	unmapped_read_count = index_mapped_reads(output_bam, log_file) 
	
	// Create and index fasta file of unmapped reads 
	continue_workflow = build_unmapped_reads(unmapped_read_count, output_bam, unmapped_fasta, log_file)

	unmapped_index_base = output_base + 'unmapped_reads.fa'
	// Build a bowtie2 index of unmapped reads
	continue_workflow = index_unmapped_reads(continue_workflow, unmapped_fasta, unmapped_index_base, log_file)

	if (params.ref_5p_fasta) {
		fivep_fasta = file(params.ref_5p_fasta)
	} else {
		fivep_fasta = ref_dir.resolve('fivep_sites.fa')
	}
	fivep_to_reads = output_dir.resolve(output_prefix + 'fivep_to_reads.sam')
	// Align fasta file of all 5' splice sites (first 20nts of introns) to unmapped reads index
	continue_workflow = align_fivep_sites(continue_workflow, unmapped_index_base, fivep_fasta, fivep_to_reads, log_file)

	if (params.ref_fasta) {
		genome_fasta = file(params.ref_fasta)
		genome_faidx = file(params.ref_fasta + '.fai')
	} else {
		genome_fasta = ref_dir.resolve('genome.fa')
		genome_faidx = ref_dir.resolve('genome.fa.fai')
	}
	// Extract reads with a mapped 5' splice site and trim it off
	fivep_to_reads_count = filter_fivep_aligns(continue_workflow, fivep_to_reads, filter_fivep_script, output_base, log_file, log_level, genome_fasta, genome_faidx, fivep_fasta, params.strand)

	heads_fasta = output_dir.resolve(output_prefix + 'heads.fa')
	heads_to_genome = output_dir.resolve(output_prefix + 'heads_to_genome.sam')
	// Map trimmed read heads to genome
	continue_workflow = align_heads(fivep_to_reads_count, genome_index, heads_fasta, heads_to_genome, log_file)

	if (params.ref_introns) {
		introns_tsv = file(params.ref_introns)
	} else {
		introns_tsv = ref_dir.resolve('introns.tsv.gz')
	}
	// Filter head alignments
	head_read_count = filter_head_aligns(continue_workflow, heads_to_genome, filter_head_script, introns_tsv, genome_fasta, output_base, log_file, log_level)
	
	if (params.ref_repeatmasker) {
		repeats_bed = params.ref_repeatmasker
	} else {
		repeats_bed = 'null'
	}
	// Filter lariat mappings and choose 1 for each read
	continue_workflow = filter_lariats(head_read_count, filter_lariats_script, output_base, log_file, log_level, seq_type, repeats_bed)

	// Make a custom track BED file of identified lariats 
	continue_workflow = make_track(continue_workflow, params.ucsc_track, track_script, output_base, log_file, log_level)

	// Classify reads, summarise results, and create empty output files if they don't exist
	end_run(continue_workflow, end_script, ref_dir, output_base, params.keep_classes, params.keep_temp, seq_type, log_file, log_level, pipeline_dir)

	
} 