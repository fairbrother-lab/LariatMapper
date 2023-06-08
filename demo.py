from larmap import (LarMap, map_lariats)

# Alternatively, call the map_lariats() function on a LarMap object

larmap_runObj = LarMap(
    fastq_dir='demo_files/demo_fastq_files_250k_bp', 
    read_1_a='250k_cWT_1.fq.gz',
    read_2_a='250k_cWT_2.fq.gz',
    read_1_b='250k_cDBR1-Y17H_1.fq.gz',
    read_2_b='250k_cDBR1-Y17H_2.fq.gz',
    output_base_name_a='WT',
    output_base_name_b='Y17H',
    num_cpus=4,
    ref_b2index='demo_files/genomes/indices/bowtie2/mm39.fa',
    ref_fasta='demo_files/genomes/fasta_files/mm39.fa ',
    ref_gtf='demo_files/genomes/annotations/mm39.gencode.basic.M32.annotation.gtf.gz',
    ref_5p_fasta='demo_files/reference_files/mouse/mm39.gencode.basic.M32.fivep_sites.fa',
    ref_3p_b2index='demo_files/reference_files/mouse/mm39.gencode.basic.M32.threep_sites.fa',
    ref_3p_lengths='demo_files/reference_files/mouse/mm39.gencode.basic.M32.threep_seq_lens.txt',
    ref_introns='demo_files/genomes/annotations/mm39.gencode.basic.M32.introns.bed.gz',
    ref_repeatmasker='demo_files/genomes/annotations/mm39.repeat_masker.bed.gz',
    results_filename='demo_mapped_reads_250k.txt'
    )

map_lariats(larmap_runObj)