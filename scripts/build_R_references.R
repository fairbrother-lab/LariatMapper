
library(magrittr)
library(rtracklayer)
library(GenomicFeatures)

# dir_path <- "/users/yzhong36/data_rsingh47/yzhong36/BP_prediction/lariat_simulation/linear_mapping/ref"

# # Check if the directory exists
# if (!file.exists(dir_path)) {
#   # Create the directory
#   dir.create(dir_path, recursive = TRUE)
#   message("Directory created: ", dir_path)
# } else {
#   message("Directory already exists: ", dir_path)
# }

# ### prepare annotation file
# gtf <- import("/users/yzhong36/data_rsingh47/yzhong36/BP_prediction/datasets/annotations/gencode/gencode.v44.basic.annotation.gtf.gz") 
# gtf
# table(gtf$gene_type)

# # BAD_GENE_TYPES = c('transcribed_unprocessed_pseudogene',
# #                    'unprocessed_pseudogene',
# #                    'processed_pseudogene' ,
# #                    'transcribed_processed_pseudogene' ,
# #                    'TEC',
# #                    'transcribed_unitary_pseudogene' ,
# #                    'rRNA_pseudogene',
# #                    'unitary_pseudogene' ,
# #                    'pseudogene' ,
# #                    'IG_V_pseudogene' ,
# #                    'TR_V_pseudogene' ,
# #                    'IG_C_pseudogene',
# #                    'TR_J_pseudogene' ,
# #                    'IG_J_pseudogene' ,
# #                    'IG_pseudogene',
# #                    'artifact' 
# #                    )

# BAD_GENE_TYPES <- c()

# filtered_gtf <- gtf[!gtf$gene_type %in% BAD_GENE_TYPES]
# table(filtered_gtf$gene_type)

# txdb <- makeTxDbFromGRanges(filtered_gtf)
# ###

# Get args
args <- commandArgs(trailingOnly = TRUE)

# Define arg options
option_list <- list(
    make_option(c("-a", "--anno"), type = "character", default = NULL, 
                help = "input annotation file (GTF or GFF)", metavar = "ANNO_PATH"),
	make_option(c("-o", "--output"), type = "character", default = NULL, 
                help = "Path to output dir", metavar = "OUT_PATH"),
	make_option(c("-d", "--quiet"), action = "store_true", default = FALSE, 
                help = "Don't print any status messages", metavar = "QUIET"),
	make_option(c("-d", "--debug"), action = "store_true", default = FALSE, 
                help = "Print extensive status messages", metavar = "DEBUG"),
)

# Parse args
parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)

anno_file <- opts$anno
output_dir <- opts$output

# Import annotations
anno = import(anno_file)
txdb = makeTxDbFromGRanges(anno)

### split into gene/exon/intron regions
exon_tx <- exonsBy(txdb, by = "tx")
intron_tx <- intronsByTranscript(txdb)

meta_info <- transcripts(txdb)
meta_gtf <- filtered_gtf[filtered_gtf$type == "transcript"]
names(meta_gtf) <- meta_gtf$transcript_id
meta_info$gene_id <- meta_gtf[meta_info$tx_name]$gene_id

exon_tx <- exon_tx %>% unlist %>% granges()
exon_tx$gene_id <- meta_info[names(exon_tx) %>% as.numeric()]$gene_id

intron_tx <- intron_tx %>% unlist %>% granges()
intron_tx$gene_id <- meta_info[names(intron_tx) %>% as.numeric()]$gene_id

exon_intron <- c(exon_tx, intron_tx)
exon_intron <- exon_intron %>% split(exon_intron$gene_id) %>% GRangesList()

exon_intron <- exon_intron %>% disjoin() %>% unlist()
exon_intron$class <- "intron"
exon_intron$class[exon_intron %over% exon_tx] <- "exon"

exon_gr <- exon_intron[exon_intron$class == "exon"]
intron_gr <- exon_intron[exon_intron$class == "intron"]

exon_gr <- exon_gr %>% split(names(exon_gr)) %>% GRangesList() %>% reduce()
intron_gr <- intron_gr %>% split(names(intron_gr)) %>% GRangesList() %>% reduce()

mean(intron_gr %over% exon_tx)
mean(intron_tx %over% exon_tx)
mean(exon_gr %over% intron_gr)
mean(intron_gr %over% exon_gr)

gene_gr <- exon_intron %>% split(names(exon_intron)) %>% GRangesList() %>% reduce()

names(exon_gr) <- paste(names(exon_gr), "exon", sep = "_")
names(intron_gr) <- paste(names(intron_gr), "intron", sep = "_")
###

saveRDS(gene_gr, file.path(output_dir, "gene_gr.rds"))
saveRDS(exon_gr, file.path(output_dir, "exon_gr.rds"))
saveRDS(intron_gr, file.path(output_dir, "intron_gr.rds"))