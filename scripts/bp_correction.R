require(Biostrings)
require(magrittr)
require(GenomicRanges)
require(GenomicFeatures)
require(bedtoolsr)

shift_matching <- function(bp_gr, pattern_l, offset, genome, correct_upstream = T, debug = T, make_plot = T){

    stopifnot(offset > 0 & offset <= 3)

    pattern <- pattern_l[["pattern"]]
    size <- pattern_l[["pattern"]] %>% nchar()
    bp_pos <- pattern_l[["bp_pos"]]

    bp_gr_extend <- bp_gr %>% resize(bp_pos, fix = "end") %>% resize(size, fix = "start")
    bp_gr$bp_pos <- bp_gr
    bp_gr$status <- "uncorrected"

    do_index <- seq_along(bp_gr_extend)

    if(correct_upstream){
        offset_l <- c(0:-offset, 1:offset)
        if(debug){
            adjusted_pos <- vector(length=offset*2+1, mode = "numeric")
            names(adjusted_pos) <- offset_l %>% as.character()
        }
    } else{
        offset_l <- 0:offset
        if(debug){
            adjusted_pos <- vector(length=offset+1, mode = "numeric")
            names(adjusted_pos) <- offset_l %>% as.character()
        }
    }
    for(i in offset_l){
        matched_index <- pattern_index(bp_gr_extend, do_index, i, pattern, genome) 
        if(length(matched_index) > 0){
            bp_gr$bp_pos[matched_index] <- shiftStranded(bp_gr$bp_pos[matched_index], i) 
            bp_gr$status[matched_index] <- ifelse(i == 0, "uncorrected", "corrected")

            if(debug){
                adjusted_pos[as.character(i)] <- length(matched_index)
            }
            
            do_index <- setdiff(do_index, matched_index)
            if(length(do_index) == 0){
                break
            }
        } 
    }

    if(debug){
        adjusted_pos <- adjusted_pos[order(names(adjusted_pos) %>% as.numeric())]
        adjusted_pos <- adjusted_pos[-which(names(adjusted_pos) == "0")]
		
        cat("Count of correction for each position relvant to the potential BP:", "\n")
        print(adjusted_pos)

        if(sum(adjusted_pos) > 0){
            cat("Ratio of correction for each position relvant to the potential BP:", "\n")
            adjusted_pos_ratio <- adjusted_pos %>% prop.table()
            print(adjusted_pos_ratio)
            if(make_plot){
                barplot(adjusted_pos_ratio, main = paste("Number of reads being corrected:",sum(adjusted_pos)))  
            }
        }
    }

    return(bp_gr)
}

pattern_index <- function(gr, index, offset, pattern, genome){

    gr_shift <- shiftStranded(gr[index], offset)
    gr_shift_seq <- getSeq(genome, gr_shift)
    gr_shift_seq_match <- pattern_matching(gr_shift_seq, pattern)
    match_index <- sapply(gr_shift_seq_match, length) == 1
    match_index <- index[match_index]
    return(match_index)

}

pattern_matching <- function(dna_string, pattern){
    
    res <- vmatchPattern(pattern, dna_string, fixed = F)
    return(res)

}

pwm_l_search <- function(bp_gr, pwm_ll, offset, genome, correct_upstream = T, debug = T, make_plot = T){

    stopifnot(offset > 0 & offset <= 3)

    bp_gr$bp_pos <- bp_gr
    bp_gr$status <- "uncorrected"

    adjust_loc <- rep(0, length(bp_gr))
    names(adjust_loc) <- rep(-Inf, length(bp_gr))
    for(pwm_l in pwm_ll){

        pwm <- pwm_l[["pwm"]]
        size <- pwm_l[["pwm"]] %>% ncol()
        bp_pos <- pwm_l[["bp_pos"]]

        if(class(genome) == "character"){
        bp_gr_extend_seq <- subset_seq(string = genome, 
                                       win_size = ((genome[1] %>% nchar) - 1)/2, 
                                       pwm_len = size, 
                                       bp_pos = bp_pos, 
                                       offset = offset, 
                                       correct_upstream = correct_upstream)
        } else{
            bp_gr_extend <- bp_gr %>% resize(bp_pos, fix = "end") %>% resize(size, fix = "start") 
            if(correct_upstream){
                bp_gr_extend <- bp_gr_extend + offset
            } else{
                bp_gr_extend <- (bp_gr_extend + offset) %>% resize(size + offset, fix = "end")
            }
            bp_gr_extend_seq <- getSeq(genome, bp_gr_extend)
        }
        
        bp_gr_extend_pwm <- sapply(bp_gr_extend_seq, pwm_matching, pwm, offset, bp_pos, correct_upstream)
        
        bp_gr_extend_pwm_max_index <- apply(bp_gr_extend_pwm, 2, which.max)
        bp_gr_extend_pwm_max_value <- apply(bp_gr_extend_pwm, 2, max)
        
        max_index_value_df <- cbind(index = bp_gr_extend_pwm_max_index, value = bp_gr_extend_pwm_max_value) %>% data.frame()
        max_index_value_df$index[max_index_value_df$value < 0.8] <- 1
        pwm_loc <- rownames(bp_gr_extend_pwm)[max_index_value_df$index] %>% as.numeric()
        adjust_loc_tmp <- pwm_loc - bp_pos
        names(adjust_loc_tmp) <- max_index_value_df$value

        adjust_loc[as.numeric(names(adjust_loc)) < as.numeric(names(adjust_loc_tmp))] <- adjust_loc_tmp[as.numeric(names(adjust_loc)) < as.numeric(names(adjust_loc_tmp))]
        names(adjust_loc)[as.numeric(names(adjust_loc)) < as.numeric(names(adjust_loc_tmp))] <- names(adjust_loc_tmp)[as.numeric(names(adjust_loc)) < as.numeric(names(adjust_loc_tmp))]

    }
    
    bp_gr$bp_pos <- shiftStranded(bp_gr$bp_pos, adjust_loc)
    bp_gr$status[adjust_loc != 0] <- "corrected"
    bp_gr$adjust_loc <- adjust_loc

    if(debug){
        if(mean(adjust_loc) == 0){
            if(correct_upstream){
                offset_l <- c(0:-offset, 1:offset)
                adjusted_pos <- vector(length=offset*2+1, mode = "numeric")
                names(adjusted_pos) <- offset_l %>% as.character()
            } else{
                offset_l <- 0:offset
                adjusted_pos <- vector(length=offset+1, mode = "numeric")
                names(adjusted_pos) <- offset_l %>% as.character()
            }
        } else{
            adjusted_pos <- adjust_loc %>% table 
            adjusted_pos <- as.vector(adjusted_pos)
            names(adjusted_pos) <- names(adjust_loc %>% table)
            
            if(length(adjusted_pos) != ifelse(correct_upstream, 2*offset+1, offset+1)){
                if(correct_upstream){
                    add_name <- c(0:-offset, 1:offset) %>% as.character()
                } else{
                    add_name <- c(0:offset) %>% as.character()
                }
                add_name <- add_name[!add_name %in% c(names(adjusted_pos))]
                for (i in seq_along(add_name)) {
                    adjusted_pos[add_name[i]] <- 0
                }
            }
        }
        adjusted_pos <- adjusted_pos[order(names(adjusted_pos) %>% as.numeric())]
        adjusted_pos <- adjusted_pos[-which(names(adjusted_pos) == "0")]

        cat("Count of correction for each position relvant to the potential BP:", "\n")
        print(adjusted_pos)

        if(sum(adjusted_pos) > 0){
            cat("Ratio of correction for each position relvant to the potential BP:", "\n")
            adjusted_pos_ratio <- adjusted_pos %>% prop.table()
            print(adjusted_pos_ratio)
            if(make_plot){
                barplot(adjusted_pos_ratio, main = paste("Number of reads being corrected:",sum(adjusted_pos)))
            }
        }
    }

    return(bp_gr)

}

# pwm_search <- function(bp_gr, pwm_l, offset, genome, correct_upstream = T, debug = T){

#     stopifnot(offset > 0 & offset <= 3)

#     bp_gr$bp_pos <- bp_gr
#     bp_gr$status <- "uncorrected"

#     pwm <- pwm_l[["pwm"]]
#     size <- pwm_l[["pwm"]] %>% ncol()
#     bp_pos <- pwm_l[["bp_pos"]]

#     if(class(genome) == "character"){
#         bp_gr_extend_seq <- subset_seq(string = genome, 
#                                        win_size = ((genome[1] %>% nchar) - 1)/2, 
#                                        pwm_len = size, 
#                                        bp_pos = bp_pos, 
#                                        offset = offset, 
#                                        correct_upstream = correct_upstream)
#     } else{
#         bp_gr_extend <- bp_gr %>% resize(bp_pos, fix = "end") %>% resize(size, fix = "start") 
#         if(correct_upstream){
#             bp_gr_extend <- bp_gr_extend + offset
#         } else{
#             bp_gr_extend <- (bp_gr_extend + offset) %>% resize(size + offset, fix = "end")
#         }
#         bp_gr_extend_seq <- getSeq(genome, bp_gr_extend)
#     }

#     bp_gr_extend_pwm <- sapply(bp_gr_extend_seq, pwm_matching, pwm, offset, bp_pos, correct_upstream)
#     bp_gr_extend_pwm_max_index <- apply(bp_gr_extend_pwm, 2, function(x) {
#         index <- intersect(which(x >= 0.8), which.max(x))
#         return(ifelse(length(index) == 0, 1, index))
#     })
#     pwm_loc <- rownames(bp_gr_extend_pwm)[bp_gr_extend_pwm_max_index] %>% as.numeric()
#     adjust_loc <- pwm_loc - bp_pos

#     bp_gr$bp_pos <- shiftStranded(bp_gr$bp_pos, adjust_loc)
#     bp_gr$status[adjust_loc != 0] <- "corrected"
#     bp_gr$adjust_loc <- adjust_loc

#     if(debug){
#         if(mean(adjust_loc) == 0){
#             if(correct_upstream){
#                 offset_l <- c(0:-offset, 1:offset)
#                 adjusted_pos <- vector(length=offset*2+1, mode = "numeric")
#                 names(adjusted_pos) <- offset_l %>% as.character()
#             } else{
#                 offset_l <- 0:offset
#                 adjusted_pos <- vector(length=offset+1, mode = "numeric")
#                 names(adjusted_pos) <- offset_l %>% as.character()
#             }
#         } else{
#             adjusted_pos <- adjust_loc %>% table 
#             adjusted_pos <- as.vector(adjusted_pos)
#             names(adjusted_pos) <- names(adjust_loc %>% table)
            
#             if(length(adjusted_pos) != ifelse(correct_upstream, 2*offset+1, offset+1)){
#                 if(correct_upstream){
#                     add_name <- c(0:-offset, 1:offset) %>% as.character()
#                 } else{
#                     add_name <- c(0:offset) %>% as.character()
#                 }
#                 add_name <- add_name[!add_name %in% c(names(adjusted_pos))]
#                 for (i in seq_along(add_name)) {
#                     adjusted_pos[add_name[i]] <- 0
#                 }
#             }
#         }
#         adjusted_pos <- adjusted_pos[order(names(adjusted_pos) %>% as.numeric())]
#         adjusted_pos <- adjusted_pos[-which(names(adjusted_pos) == "0")]

#         cat("Count of correction for each position relvant to the potential BP:", "\n")
#         print(adjusted_pos)

#         if(sum(adjusted_pos) > 0){
#             cat("Ratio of correction for each position relvant to the potential BP:", "\n")
#             adjusted_pos_ratio <- adjusted_pos %>% prop.table()
#             print(adjusted_pos_ratio)
#             barplot(adjusted_pos_ratio, main = paste("Number of reads being corrected:",sum(adjusted_pos)))
#         }
#     }

#     return(bp_gr)

# }

pwm_matching <- function(dna_string, pwm, offset, ref_pos, correct_upstream){

    if(correct_upstream){
        pos_index <- c(ref_pos:(ref_pos - offset), (ref_pos + 1):(ref_pos + offset))
        pos_start_index <- c(offset + 1, (offset):1, (offset + 2):(length(dna_string) - ncol(pwm) + 1))
    } else{
        pos_index <- ref_pos:(ref_pos + length(dna_string) - ncol(pwm))
        pos_start_index <- pos_index - (ref_pos - 1)
    }

    pwm_score <- PWMscoreStartingAt(pwm, dna_string, 
                                    starting.at = pos_start_index)
    names(pwm_score) <- pos_index
    return(pwm_score)

}

model_based_search <- function(bp_gr, cbp_prob, offset, correct_upstream = T, debug = T, make_plot = T){

    stopifnot(offset > 0 & offset <= 3)

    og_bp_gr <- bp_gr
    toDO_index <- !bp_gr %over% cbp_prob
    bp_gr <- bp_gr[toDO_index]

    if(correct_upstream){
        bp_gr_extend <- bp_gr + offset
    } else{
        bp_gr_extend <- (bp_gr + offset) %>% resize(1 + offset, fix = "end")
    }

    og_bp_gr$bp_pos <- og_bp_gr
    og_bp_gr$status <- "uncorrected"
    og_bp_gr$adjust_loc <- 0
    bp_gr$bp_pos <- bp_gr
    bp_gr$status <- "uncorrected"
    bp_gr$adjust_loc <- 0

    mapped_res <- bp_map2_cbp(cbp_prob, bp_gr_extend)

    if(correct_upstream){
        shift_pos <- mapped_res$start - (offset + 1)
    } else{
        shift_pos <- mapped_res$start - 1
    }

    bp_gr$bp_pos[mapped_res$transcriptsHits] <- shiftStranded(bp_gr$bp_pos[mapped_res$transcriptsHits], shift_pos)
    bp_gr$status[start(bp_gr) != start(bp_gr$bp_pos)] <- "corrected"
    bp_gr$adjust_loc[mapped_res$transcriptsHits] <- shift_pos

    if(debug){
        if(is.null(mapped_res)){
            if(correct_upstream){
                offset_l <- c(0:-offset, 1:offset)
                adjusted_pos <- vector(length=offset*2+1, mode = "numeric")
                names(adjusted_pos) <- offset_l %>% as.character()
            } else{
                offset_l <- 0:offset
                adjusted_pos <- vector(length=offset+1, mode = "numeric")
                names(adjusted_pos) <- offset_l %>% as.character()
            }
        } else{
            adjust_loc <- shift_pos
            adjusted_pos <- adjust_loc %>% table 
            adjusted_pos <- as.vector(adjusted_pos)
            names(adjusted_pos) <- names(adjust_loc %>% table)

            if(length(adjusted_pos) != ifelse(correct_upstream, 2*offset+1, offset+1)){
                if(correct_upstream){
                    add_name <- c(0:-offset, 1:offset) %>% as.character()
                } else{
                    add_name <- c(0:offset) %>% as.character()
                }
                add_name <- add_name[!add_name %in% c(names(adjusted_pos))]
                for (i in seq_along(add_name)) {
                    adjusted_pos[add_name[i]] <- 0
                }
            }
        }
        adjusted_pos <- adjusted_pos[order(names(adjusted_pos) %>% as.numeric())]
        adjusted_pos <- adjusted_pos[-which(names(adjusted_pos) == "0")]

        cat("Count of correction for each position relvant to the potential BP:", "\n")
        print(adjusted_pos)

        if(sum(adjusted_pos) > 0){
            cat("Ratio of correction for each position relvant to the potential BP:", "\n")
            adjusted_pos_ratio <- adjusted_pos %>% prop.table()
            print(adjusted_pos_ratio)
            if(make_plot){
                barplot(adjusted_pos_ratio, main = paste("Number of reads being corrected:",sum(adjusted_pos)))
            }
        }
    }

    og_bp_gr$bp_pos[toDO_index] <- bp_gr$bp_pos
    og_bp_gr$status[toDO_index] <- bp_gr$status
    og_bp_gr$adjust_loc[toDO_index] <- bp_gr$adjust_loc
    return(og_bp_gr)

}

bp_map2_cbp <- function(bp_gr, bp_extend){

    names(bp_extend) <- seq_along(bp_extend)
    map_bp2_bp_extend <- mapToTranscripts(bp_gr, bp_extend) 
    
    if(length(map_bp2_bp_extend) == 0){
        return(NULL)
    }

    mcols(map_bp2_bp_extend)$bp_prob <- bp_gr$BP_prob[map_bp2_bp_extend$xHits]
    agg_result <- aggregate(bp_prob ~ transcriptsHits, data = map_bp2_bp_extend %>% as.data.frame, FUN = max) %>% 
    merge(map_bp2_bp_extend %>% as.data.frame, by = c("transcriptsHits", "bp_prob"), sort = F)

    return(agg_result)

}

shiftStranded <- function(x, shift=0L,...) GenomicRanges::shift(x ,shift=shift*ifelse('-'==strand(x),-1,1),...)

subset_seq <- function(string, win_size, pwm_len, bp_pos, offset, correct_upstream = T){

    genome <- DNAStringSet(string)
    n_len <- 2 * win_size + 1
    bp_loc <- ((n_len - 1) / 2) + 1

    stopifnot(bp_loc - bp_pos + 1 - offset > 0)

    if(correct_upstream){
        genome_seq <- subseq(genome, start = bp_loc - bp_pos + 1 - offset, 
                             width = pwm_len + 2*offset)
    } else{
        genome_seq <- subseq(genome, start = bp_loc - bp_pos + 1, 
                             width = pwm_len + offset)
    }
    return(genome_seq)

}

get_context_seq <- function(file, ref_fasta, correction_context_size){
	context_bed = data.frame(
		'chrom' = file$chrom,
		'start' = c(file$bp_pos - correction_context_size),
		'end' = c(file$bp_pos + correction_context_size + 1),
		'id' = file$read_id,
		'score' = rep(0, nrow(file)),
		'strand' = file$strand
	)

	context_seq <- bt.getfasta(ref_fasta, 
							context_bed,
							bedOut=T,
							s=T)
	
	context_seq <- context_seq[[7]] %>% DNAStringSet() %>% as.matrix()
	
	return(context_seq)
}

trim_context <- function(context_seq, shift_loc, out_context_size){


}