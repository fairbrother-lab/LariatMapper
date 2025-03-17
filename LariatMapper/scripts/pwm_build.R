suppressPackageStartupMessages(require(magrittr))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(Biostrings))

args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL, 
                help = "Input fasta file for building PWM", metavar = "INPUT_PATH"),
    make_option(c("-d", "--prior_distribution"), type = "character", default = "A=0.25;C=0.25;G=0.25;T=0.25", 
                help = 'Dirichlet prior distribution of nucleoides as for building PWM. Default: "A=0.25;C=0.25;G=0.25;T=0.25"', metavar = "PRIOR"),            
    make_option(c("-p", "--bp_pos"), type = "numeric", default = NULL, 
                help = "Branchpoint position in fasta file. 1-based indexing", metavar = "BP_POS"),
    make_option(c("-o", "--output_base"), type = "character", default = NULL, 
                help = "Path to output", metavar = "OUT_PATH")
)

parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)

seq_path <- opts$input
prior <- opts$prior_distribution
bp_pos <- opts$bp_pos
save_path <- opts$output_base

## check parameters
prior <- strsplit(prior, ";")[[1]] %>% lapply(strsplit, "=") %>% data.frame() %>% t %>% as.data.frame()
colnames(prior) <- c("nt", "prob")
rownames(prior) <- prior$nt
prior$prob <- as.numeric(prior$prob)
stopifnot(sum(prior$prob) == 1)

seq <- readDNAStringSet(seq_path)
stopifnot(length(unique(width(seq))) == 1)
stopifnot(bp_pos > 0 & bp_pos <= unique(width(seq)))
##

bp_pfm <- consensusMatrix(seq)
bp_pwm <- PWM(bp_pfm, prior.params = c(A = prior["A", "prob"], C = prior["C", "prob"], G = prior["G", "prob"], T = prior["T", "prob"]))
pwm_l <- list("pwm" = bp_pwm, "bp_pos" = bp_pos)

saveRDS(pwm_l, file.path(save_path, "pwm.rds"))
