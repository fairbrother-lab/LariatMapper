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
    make_option(c("-g", "--gtf"), type = "character", default = NULL, 
                help = "Path to GTF dir", metavar = "DIR_PATH"),
    make_option(c("-o", "--output"), type = "character", default = NULL, 
                help = "Path to output dir", metavar = "OUT_PATH"),
    make_option(c("-r", "--read_layout"), type = "character", default = "pairEnd", 
                help = "singleEnd or pairEnd datasets", metavar = "READ_LAYOUT")
)

parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)

input_bam <- opts$input
utils_R <- opts$file
gtf_dir <- opts$gtf
output_dir <- opts$output
read_layout <- opts$read_layout

if (is.null(read_layout) || !read_layout %in% c("singleEnd", "pairEnd")) {
    stop("--read_layout must be singleEnd or pairEnd.")
}

### Source utils
source(utils_R)
###

### Prepare annotation file
gene_gr <- readRDS(file.path(gtf_dir, "gene_gr.rds"))
exon_gr <- readRDS(file.path(gtf_dir, "exon_gr.rds"))
intron_gr <- readRDS(file.path(gtf_dir, "intron_gr.rds"))
###

### Run counting step
singleEnd_mode <- ifelse(read_layout == "pairEnd", F, T)
bam_res <- linear_map_integrate(input_bam,
                                gene_gr = gene_gr, exon_gr = exon_gr, intron_gr = intron_gr,
                                singleEnd = singleEnd_mode)
###

### Output counts table
# Check if the directory exists
if (!file.exists(output_dir)) {
  # Create the directory
  dir.create(output_dir, recursive = TRUE)
  message("Directory created: ", output_dir)
} else {
  message("Directory already exists: ", output_dir)
}

write.table(bam_res, file.path(output_dir, 
                               paste0(basename(input_bam), "_", "count.tsv")),
            row.names = F, col.names = T, quote = F, sep = "\t")
###
