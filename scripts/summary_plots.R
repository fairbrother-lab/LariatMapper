suppressPackageStartupMessages(require(ggplot2))


args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL, 
                help = "input bam file", metavar = "BAM_PATH"),
    make_option(c("-f", "--file"), type = "character", default = NULL, 
                help = "Path to utils R file", metavar = "FILE_PATH"),
    make_option(c("-r", "--ref_dir"), type = "character", default = NULL, 
                help = "Path to reference dir", metavar = "DIR_PATH"),
    make_option(c("-o", "--output_base"), type = "character", default = NULL, 
                help = "Path to output", metavar = "OUT_PATH"),
    make_option(c("-l", "--read_layout"), type = "character", default = "paired", 
                help = "single = single-end sequencing data, paired = paired-end sequencing data", metavar = "READ_LAYOUT")
)


output_base <- opts$output_base