require(GenomicAlignments)
require(magrittr)
require(stringr)

total_reads_bam <- function(file, ...){

    param <- ScanBamParam(flag = scanBamFlag(...), what = c("qname"))
    bam <- scanBam(file, param = param)
    return(length(bam[[1]]$qname))

}

total_reads_tsv <- function(file){

    res <- read.csv(file, sep = "\t", header = T)
    stopifnot(mean(duplicated(res$read_id)) == 0)
    return(nrow(res))

}

linear_map <- function(file, gene_gr, exon_gr, intron_gr, ...){

    gene_counts <- summarizeOverlaps(features = gene_gr, reads = file, mode = "IntersectionStrict", ...) %>% assay()
    colnames(gene_counts) <- "count"

    exon_intron_gr <- c(exon_gr, intron_gr)
    exon_intron_counts <- summarizeOverlaps(features = exon_intron_gr, reads = file, mode = "IntersectionStrict", ...) %>% assay() 
    colnames(exon_intron_counts) <- "count"

    exon_only_counts <- summarizeOverlaps(features = exon_gr %>% unlist, reads = file, mode = "IntersectionStrict", ...) %>% assay() 
    colnames(exon_only_counts) <- "count"

    exon_counts = exon_intron_counts[names(exon_gr),, drop = F]
    intron_counts = exon_intron_counts[names(intron_gr),, drop = F]

    exon_only_counts <- cbind(exon_only_counts, "gene_id" = rownames(exon_only_counts))
    exon_only_counts <- exon_only_counts %>% as.data.frame()
    exon_only_counts$count <- exon_only_counts$count %>% as.numeric()

    exon_only_counts_by_gene <- aggregate(count ~ gene_id, data = exon_only_counts, FUN = sum)
    colnames(exon_only_counts_by_gene) <- c("gene", "count")
    rownames(exon_only_counts_by_gene) <- exon_only_counts_by_gene$gene
    exon_only_counts_by_gene$gene <- NULL
    exon_only_counts_by_gene <- exon_only_counts_by_gene %>% as.matrix()
    exon_only_counts_by_gene <- exon_only_counts_by_gene[rownames(exon_counts),, drop = F]

    exon_intron_counts <- exon_intron_counts %>% as.data.frame()
    exon_intron_counts$gene_id = (rownames(exon_intron_counts) %>% str_split("_", simplify = T))[,1]

    exon_intron_counts_by_gene <- aggregate(count ~ gene_id, data = exon_intron_counts, FUN = sum)
    colnames(exon_intron_counts_by_gene) <- c("gene", "count")
    rownames(exon_intron_counts_by_gene) <- exon_intron_counts_by_gene$gene
    exon_intron_counts_by_gene$gene <- NULL
    exon_intron_counts_by_gene <- exon_intron_counts_by_gene %>% as.matrix()
    exon_intron_counts_by_gene <- exon_intron_counts_by_gene[rownames(gene_counts),, drop = F]

    premrna_counts <- gene_counts - exon_intron_counts_by_gene
    exon_junc_counts <- exon_counts - exon_only_counts_by_gene ## not accurate; see leafcutter or regtools extract
    premrna_counts[premrna_counts < 0] <- 0
    exon_junc_counts[exon_junc_counts < 0] <- 0

    return(list(gene_counts = gene_counts,
                premrna_counts = premrna_counts,
                exon_counts = exon_only_counts_by_gene, 
                exon_junc_counts = exon_junc_counts,
                intron_counts = intron_counts))

}

linear_map_integrate <- function(file, gene_gr, exon_gr, intron_gr, ...){

    map_res <- linear_map(file, gene_gr = gene_gr, exon_gr = exon_gr, intron_gr = intron_gr,
                          ...)

    gene_counts <- map_res$gene_counts %>% data.frame()
    premrna_counts <- map_res$premrna_counts %>% data.frame()
    exon_counts <- map_res$exon_counts %>% data.frame()
    exon_junc_counts <- map_res$exon_junc_counts %>% data.frame()
    intron_counts <- map_res$intron_counts %>% data.frame()

    colnames(gene_counts) <- "gene"
    colnames(premrna_counts) <- "exon_intron_junc"
    colnames(exon_counts) <- "exon_only"
    colnames(exon_junc_counts) <- "exon_exon_junc"
    colnames(intron_counts) <- "intron_only"

    gene_counts$gene_id <- (str_split(rownames(gene_counts), "_", simplify = T))[,1]
    premrna_counts$gene_id <- (str_split(rownames(premrna_counts), "_", simplify = T))[,1]
    exon_counts$gene_id <- (str_split(rownames(exon_counts), "_", simplify = T))[,1]
    exon_junc_counts$gene_id <- (str_split(rownames(exon_junc_counts), "_", simplify = T))[,1]
    intron_counts$gene_id <- (str_split(rownames(intron_counts), "_", simplify = T))[,1]

    count_df <- merge(gene_counts, premrna_counts, by = "gene_id", all = T) %>% 
    merge(exon_counts, by = "gene_id", all = T) %>% 
    merge(exon_junc_counts, by = "gene_id", all = T) %>% 
    merge(intron_counts, by = "gene_id", all = T) 

    count_df[is.na(count_df)] <- 0

    return(count_df)

}

linear_map_integrate_stats <- function(file, gene_gr, exon_gr, intron_gr, ...){

    map_res <- linear_map_integrate(file, gene_gr = gene_gr, exon_gr = exon_gr, intron_gr = intron_gr,
                                    ...)
    map_res_count <- colSums(map_res[, -1])   

    return(map_res_count)                 

}


linear_map_stats <- function(file, gene_gr, exon_gr, intron_gr, ...){

    map_res <- linear_map(file, gene_gr = gene_gr, exon_gr = exon_gr, intron_gr = intron_gr,
                          ...)
    map_res_count <- sapply(map_res, colSums)
    names(map_res_count) <- c("Gene", "pre-mRNA", "Exonic", "Exonic-junc", "Intronic")

    return(map_res_count)

}