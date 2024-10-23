suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(GenomicAlignments))
suppressPackageStartupMessages(require(magrittr))
suppressPackageStartupMessages(require(stringr))

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

parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)

input_bam <- opts$input
utils_R <- opts$file
ref_dir <- opts$ref_dir
output_base <- opts$output_base
read_layout <- opts$read_layout

if (is.null(read_layout) || !read_layout %in% c("single", "paired")) {
    stop("--read_layout must be 'single' or 'paired'.")
}

### Source utils
source(utils_R)
###

### Prepare annotation file
gene_gr <- readRDS(file.path(ref_dir, "gene_gr.rds"))
exon_gr <- readRDS(file.path(ref_dir, "exon_gr.rds"))
intron_gr <- readRDS(file.path(ref_dir, "intron_gr.rds"))
###

### Run counting step
singleEnd_mode <- ifelse(read_layout == "paired", F, T)
bam_res <- linear_map_integrate(input_bam,
                                gene_gr = gene_gr, exon_gr = exon_gr, intron_gr = intron_gr,
                                singleEnd = singleEnd_mode)
###

### Output counts table

write.table(bam_res, paste0(output_base, "output.bam_count.tsv"),
            row.names = F, col.names = T, quote = F, sep = "\t")
###