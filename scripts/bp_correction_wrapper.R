suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(magrittr))
suppressPackageStartupMessages(require(GenomicRanges))
suppressPackageStartupMessages(require(Biostrings))
suppressPackageStartupMessages(require(GenomicFeatures))


args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL, 
                help = "Input file: lariat_read.tsv", metavar = "INPUT_PATH"),
    make_option(c("-r", "--ref_fasta"), type = "character", default = NULL, 
                help = "m"),
    make_option(c("-f", "--file"), type = "character", default = NULL, 
                help = "Path to utils R file", metavar = "FILE_PATH"),
    make_option(c("-s", "--window_size"), type = "numeric", default = 3, 
                help = "Window size for correction. i.e. If the value is set to 3, it will only search for the optimal BPs with 3-nt distance relevant to putative BPs", 
                metavar = "WINDOW_SIZE"),
    make_option(c("-d", "--both_upstream_downstream"), action = "store_true", default = FALSE,
                help = "Correcting both upstream and downstream of putative BPs (default: FALSE)"),
    make_option(c("-m", "--method"), type = "character", default = "PWM", 
                help = "The method used for branchpoint correction. PWM (Position Weight Matrix) or Model-based. Model-based method only applies for the human genome", metavar = "CORRECTION_METHOD"),
    make_option(c("-w", "--PWM_path"), type = "character", default = NULL, 
                help = "Path to PWM matrix file, when PWM is enabled for --method. Multiple paths can be provided separated by commas", metavar = "PWM_PATH"),
    make_option(c("-e", "--model_path"), type = "character", default = NULL, 
                help = "Path to pre-computed model file, when Model-based is enabled for --method", metavar = "MODEL_PATH"),
    make_option(c("-c", "--out_context_size"), type = "numeric", default = 8, 
                help = "k"),
    make_option(c("-o", "--output_base"), type = "character", default = NULL, 
                help = "Path to output", metavar = "OUT_PATH"),
	make_option(c("-l", "--log_level"), type = "character", default = "INFO", 
                help = "Log level (DEBUG, INFO, WARNING, or ERROR)")
)

parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)

input_lariat <- opts$input
ref_fasta <- opts$ref_fasta
utils_R <- opts$file
offset <- opts$window_size
correct_upstream <- opts$both_upstream_downstream
correction_method <- opts$method
pwm_path <- opts$PWM_path
model_path <- opts$model_path
out_context_size <- opts$out_context_size
output_dir <- opts$output_base
log_level <- opts$log_level

### Source utils
source(utils_R)
###

# Process paths
if(correction_method == "PWM"){
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
    if (log_level == "DEBUG") {
      cat("Valid PWM provided:\n")
      print(paths)
    }
  } else {
	if (log_level == "DEBUG") {
      cat("No PWM provided.\n")
	}
  }
} else if (correction_method == "Model-based") {

    if (!file.exists(model_path)) {
          stop(paste("Path to pre-computed model file does not exist:", model_path))
        } else{
            cbp_prob <- readRDS(model_path)
        } 
	if (log_level == "DEBUG") {
	  cat("Valid pre-computed model file provided:\n")
	  print(model_path)
	}

} else{
  stop("The correction method can either be PWM or Model-based")
}

file <- read.csv(input_lariat, sep = "\t")
gr <- GRanges(seqnames = file$chrom,
              IRanges(file$bp_pos + 1, width = 1),
              strand = file$strand)
genome <- get_context_seq(file, ref_fasta, out_context_size, offset, correction_method, pwm_l)

### Running the step
debug = ifelse(log_level %in% c("DEBUG", "INFO"), T, F)
if(correction_method == "PWM"){
  corrected_gr <- pwm_l_search(gr, pwm_l, offset, genome, correct_upstream, debug, make_plot = F)
} else if (correction_method == "Model-based") {
  corrected_gr <- model_based_search(gr, cbp_prob, offset, correct_upstream, debug, make_plot = F)
}
###

### Add <corrected_XXX> columns to lariat_reads.tsv
context_seq <- genome %>% DNAStringSet() %>% as.matrix()
shift_loc <- corrected_gr$adjust_loc
c_loc <- ((ncol(context_seq) - 1)/2)+1

corrected_bp_pos <- start(corrected_gr$bp_pos) - 1
file <- as.data.frame(append(file, list("corrected_bp_pos" = corrected_bp_pos), after = 6))
corrected_bp_dist_to_threep <- ifelse(file$strand == "+", 
                                      file$corrected_bp_pos - file$threep_pos, 
                                      file$threep_pos - file$corrected_bp_pos)
file <- as.data.frame(append(file, list("corrected_bp_dist_to_threep" = corrected_bp_dist_to_threep), after = 9))
corrected_genomic_bp_nt <- sapply(seq_along(shift_loc), function(x){
  context_seq[x, c_loc+shift_loc[x]]
})
file <- as.data.frame(append(file, list("corrected_genomic_bp_nt" = corrected_genomic_bp_nt), after = 15))
corrected_genomic_bp_context <- sapply(seq_along(shift_loc), function(x){
  start <- c_loc + shift_loc[x] - out_context_size
  end <- c_loc + shift_loc[x] + out_context_size
  trimmed_seq <- context_seq[x, start:end]
  paste(trimmed_seq, collapse="")
})
file <- as.data.frame(append(file, list("corrected_genomic_bp_context" = corrected_genomic_bp_context), after = 18))
corrected_bp_mismatch = ifelse(file$read_bp_nt != file$corrected_genomic_bp_nt, "True", "False")
file <- as.data.frame(append(file, list("corrected_bp_mismatch" = corrected_bp_mismatch), after = 17))
file$corrected <- ifelse(corrected_gr$status == "corrected", "True", "False")

write.table(file, file.path(output_dir, "lariat_reads.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
###