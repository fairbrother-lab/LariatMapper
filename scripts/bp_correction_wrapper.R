suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(magrittr))
suppressPackageStartupMessages(require(GenomicRanges))
suppressPackageStartupMessages(require(Biostrings))

args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL, 
                help = "Input file: lariat_read.tsv", metavar = "INPUT_PATH"),
    make_option(c("-f", "--file"), type = "character", default = NULL, 
                help = "Path to utils R file", metavar = "FILE_PATH"),
    make_option(c("-s", "--window_size"), type = "numeric", default = 2, 
                help = "Window size for correction. i.e. If the value is set to 2, it will only search for the optimal BPs with 2-nt distance relevant to putative BPs", 
                metavar = "WINDOW_SIZE"),
    make_option(c("-d", "--both_upstream_downstream"), action = "store_true", default = FALSE,
                help = "Correcting both upstream and downstream of putative BPs (default: FALSE)"),
    make_option(c("-m", "--method"), type = "character", default = "PWM", 
                help = "The method used for branchpoint correction. Only support PWM method for now", metavar = "CORRECTION_METHOD"),
    make_option(c("-w", "--PWM_path"), type = "character", default = NULL, 
                help = "Path to PWM matrix file. Multiple paths can be provided separated by commas.", metavar = "PWM_PATH"),
    make_option(c("-o", "--output_base"), type = "character", default = NULL, 
                help = "Path to output", metavar = "OUT_PATH")
)

parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)

input_lariat <- opts$input
utils_R <- opts$file
offset <- opts$window_size
correct_upstream <- opts$both_upstream_downstream
correction_method <- opts$method
pwm_path <- opts$PWM_path
output_dir <- opts$output_base

# Process paths
pwm_l <- list()
pwm_l_c <- 1
if (!is.null(pwm_path)) {
  # Split the input paths by commas
  paths <- unlist(strsplit(pwm_path, ","))
  
  # Validate each path
  for (path in paths) {
    if (!file.exists(path)) {
      stop(paste("Path to PWM does not exist:", path))
    } else{
        pwm_l[[pwm_l_c]] <- readRDS(path)
        pwm_l_c <- pwm_l_c + 1
    }
  }
  
  # Print valid paths
  cat("Valid PWM provided:\n")
  print(paths)
} else {
  cat("No PWM provided.\n")
}

file <- read.csv(input_lariat, sep = "\t")
gr <- GRanges(seqnames = file$chrom,
              IRanges(file$bp_pos + 1),
              strand = file$strand)
genome <- file$genomic_bp_context

### Source utils
source(utils_R)
###

### Running the step
if(length(pwm_l) == 1){
    corrected_gr <- pwm_search(gr, pwm_l[[1]], offset, genome, correct_upstream, debug = T)
} else{
    corrected_gr <- pwm_l_search(gr, pwm_l, offset, genome, correct_upstream, debug = T)
}
###

### modify lariat_reads.tsv
corrected_bp_pos <- start(corrected_gr$bp_pos) - 1
file <- as.data.frame(append(file, list("corrected_bp_pos" = corrected_bp_pos), after = 6))

corrected_bp_dist_to_threep <- ifelse(file$strand == "+", 
                                      file$corrected_bp_pos - file$threep_pos, 
                                      file$threep_pos - file$corrected_bp_pos)
file <- as.data.frame(append(file, list("corrected_bp_dist_to_threep" = corrected_bp_dist_to_threep), after = 9))

context_seq <- file$genomic_bp_context %>% DNAStringSet() %>% as.matrix()
shift_loc <- corrected_gr$adjust_loc
c_loc <- ((ncol(context_seq) - 1)/2)+1
corrected_genomic_bp_nt <- sapply(seq_along(shift_loc), function(x){
  context_seq[x, c_loc + shift_loc[x]]
})
file <- as.data.frame(append(file, list("corrected_genomic_bp_nt" = corrected_genomic_bp_nt), after = 15))

file$correction_status <- corrected_gr$status

write.table(file, file.path(output_dir, "lariat_reads_corrected.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
###