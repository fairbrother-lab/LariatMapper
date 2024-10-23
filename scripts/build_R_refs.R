suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))

# Get args
args <- commandArgs(trailingOnly = TRUE)

# Define arg options
option_list <- list(
    make_option(c("-a", "--anno"), type = "character", default = NULL, 
                help = "input annotation file (GTF or GFF)", metavar = "ANNO_PATH"),
	make_option(c("-g", "--g_attr"), type = "character", default = "gene_id", 
                help = "The name of the attribute that encodes gene ID (default=gene_id)", metavar = "OUT_PATH"),
	make_option(c("-t", "--t_attr"), type = "character", default = "transcript_id", 
                help = "The name of the attribute that encodes transcript ID (default=transcript_id)", metavar = "OUT_PATH"),
	make_option(c("-o", "--output"), type = "character", default = NULL, 
                help = "Path to output dir", metavar = "OUT_PATH"),
	make_option(c("-q", "--quiet"), action = "store_true", default = FALSE, 
                help = "Don't print any status messages", metavar = "QUIET"),
	make_option(c("-d", "--debug"), action = "store_true", default = FALSE, 
                help = "Print extensive status messages", metavar = "DEBUG")
)

# Parse args
parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)

anno_file <- opts$anno
g_attr <- opts$g_attr
t_attr <- opts$t_attr
output_dir <- opts$output
quiet <- opts$quiet
debug <- opts$debug

# Import annotations
anno <- import(anno_file)

# Assign gene ID and transcript ID
mcols(anno)['gene_id'] = mcols(anno)[[g_attr]] %>% sapply(., function(x) if (length(x) == 0) NA else x)
mcols(anno)['tx_name'] = mcols(anno)[[t_attr]] %>% sapply(., function(x) if (length(x) == 0) NA else x)

# Make TxDb object
txdb <- makeTxDbFromGRanges(anno)

# Split into exon/intron regions
exon_tx <- exonsBy(txdb, by = "tx")
intron_tx <- intronsByTranscript(txdb)

# Get gene-transcript mapping
meta_info <- transcripts(txdb, column=c("tx_name", "gene_id"))
meta_info$gene_id <- as.character(meta_info$gene_id)

# Add gene ID to exon/intron regions
exon_tx <- exon_tx %>% unlist %>% granges()
exon_tx$gene_id <- meta_info[names(exon_tx) %>% as.numeric()]$gene_id

intron_tx <- intron_tx %>% unlist %>% granges()
intron_tx$gene_id <- meta_info[names(intron_tx) %>% as.numeric()]$gene_id

# Combine exon/intron regions and group by gene ID
exon_intron <- c(exon_tx, intron_tx)
exon_intron <- exon_intron %>% split(exon_intron$gene_id) %>% GRangesList()

# Disjoin fragments per each gene
exon_intron <- exon_intron %>% disjoin() %>% unlist()
exon_intron$class <- "intron"
exon_intron$class[exon_intron %over% exon_tx] <- "exon"

# Extract exon/intron regions
exon_gr <- exon_intron[exon_intron$class == "exon"]
intron_gr <- exon_intron[exon_intron$class == "intron"]

# Group and collapse exon/intron regions by gene ID
exon_gr <- exon_gr %>% split(names(exon_gr)) %>% GRangesList() %>% reduce()
names(exon_gr) <- paste(names(exon_gr), "exon", sep = "_")
intron_gr <- intron_gr %>% split(names(intron_gr)) %>% GRangesList() %>% reduce()
names(intron_gr) <- paste(names(intron_gr), "intron", sep = "_")

if (debug) {
	print(paste(" ", mean(intron_gr %over% exon_tx)))
	print(paste(" ", mean(intron_tx %over% exon_tx)))
	print(paste(" ", mean(exon_gr %over% intron_gr)))
	print(paste(" ", mean(intron_gr %over% exon_gr)))
}

# Make gene regions object
gene_gr <- exon_intron %>% split(names(exon_intron)) %>% GRangesList() %>% reduce()

# Write GRanges objects to file
saveRDS(exon_gr, file.path(output_dir, "exon_gr.rds"))
saveRDS(intron_gr, file.path(output_dir, "intron_gr.rds"))
saveRDS(gene_gr, file.path(output_dir, "gene_gr.rds"))