library(InteractionSet)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library("R.utils")
library(GenomicInteractions)
library(dplyr)

library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(plyranges)

library(tibble)

## functions for inclusive hierarchical annotation
# annotates a Granges object based on the provided objects in hierarchical order
annotate_with_priority <- function(query_gr, annotation_layers) {
  query_gr$annotation <- NA_character_
  
  for (i in seq_along(annotation_layers)) {
    layer_name <- names(annotation_layers)[i]
    
    # find overlaps between current feature and Granges
    overlaps <- findOverlaps(query_gr, annotation_layers[[i]])
    # annotate only if not already annotated
    query_gr[queryHits(overlaps)]$annotation <- 
      ifelse(is.na(query_gr[queryHits(overlaps)]$annotation),
             layer_name,
             query_gr[queryHits(overlaps)]$annotation)
  }
  query_gr$annotation[is.na(query_gr$annotation)] <- "Null"
  return(query_gr)
}

# annotates Ginteractions (loops)
annotate_ginteractions <- function(gi_obj, annotation_layers, order=c('E','P','CTCF','Null')) {
  # annotate each of the anchors 
  anchors1 <- annotate_with_priority(anchors(gi_obj, type="first"), annotation_layers)
  anchors2 <- annotate_with_priority(anchors(gi_obj, type="second"), annotation_layers)
  
  # set the null annotation
  anchors1$annotation[is.na(anchors1$annotation)] <- "Null"
  anchors2$annotation[is.na(anchors2$annotation)] <- "Null"
  
  # TODO: if the metadata of the original GInteractions is not none, set the columns to be tibble columns 
  
  # re-arrange the anchors such that the pairs are symmetrical (i.e. E-P = P-E)
  bind_cols(
    as_tibble(anchors1) %>% 
      dplyr::select(seqnames, start, end, annotation) %>%
      rename_with(~paste0(.,"1")),
    
    as_tibble(anchors2) %>% 
      dplyr::select(seqnames, start, end, annotation) %>% 
      rename_with(~paste0(.,"2")) %>%
      # somehow need to make this happen in custom order so that symmetry is preserved
      mutate(
        rank1 = match(anchors1$annotation, order),
        rank2 = match(anchors2$annotation, order),
        anno_low  = if_else(rank1 <= rank2, anchors1$annotation, anchors2$annotation),
        anno_high = if_else(rank1 > rank2, anchors1$annotation, anchors2$annotation),
        pair = paste(anno_low, anno_high, sep = "-")
      ) %>%
      select(-rank1, -rank2, -anno_low, -anno_high)
   ) 
}


## functions for exclusive annotation 
create_annotation_matrix <- function(query_gr, annotation_layers, parallel=FALSE, mcols_use=list()) {
  # Convert to named GRangesList if not already
  if (!is(annotation_layers, "GRangesList")) {
    annotation_layers <- GRangesList(annotation_layers)
  }
  
  if (is.null(names(annotation_layers))) {
    names(annotation_layers) <- paste0("layer", seq_along(annotation_layers))
  }
  
  # Parallel or serial processing
  if (parallel) {
    library(BiocParallel)
    register(MulticoreParam())
    map_fun <- bplapply
  } else {
    map_fun <- lapply
  }
  
  # Create logical overlap matrix and collect metadata columns
  results <- map_fun(seq_along(annotation_layers), function(i) {
    annot_gr <- annotation_layers[[i]]
    name <- names(annotation_layers)[i]
    ov_logical <- overlapsAny(query_gr, annot_gr, maxgap = -1L, type = "any")
    
    # Start with logical column
    out <- tibble(!!name := ov_logical)
    
    # If metadata exists, extract and align to query_gr
    for (colname in mcols_use[[name]]) {
      if (colname %in% colnames(mcols(annot_gr))) {
        hits <- findOverlaps(query_gr, annot_gr, maxgap = -1L, type = "any")
        values <- rep(NA, length(query_gr))
        values[queryHits(hits)] <- mcols(annot_gr)[[colname]][subjectHits(hits)]
        out[[paste0(name, '_', colname)]] <- values
      }
    }
    out
  })
  
  # Combine results into a single tibble
  bind_cols(
    as_tibble(query_gr) %>%
      select(seqnames, start, end),
    do.call(bind_cols, results)
  )
}

collapse_gr<-function(query_gr, hierarchy_group) {
  # get names of query gr and match them to the names to collapse
  for (group in hierarchy_group) {
    group_cols <- group
    print(group_cols)
    
    mat <- as.matrix(mcols(query_gr)[, group, drop=FALSE])
    true_idx <- max.col(mat, ties.method = "first")
    
    # Create a logical matrix with only the highest-priority TRUE per row
    collapsed <- matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat))
    collapsed[cbind(seq_len(nrow(mat)), true_idx)] <- rowSums(mat) > 0
    
    # Replace the relevant columns in query_gr
    mcols(query_gr)[, group] <- collapsed
  }
  return(query_gr)
}

collapse_matrix<-function(query_mat, hierarchy_group) {
  mat<-query_mat[, 4:ncol(query_mat)]
  # get names of query gr and match them to the names to collapse
  for (group in hierarchy_group) {
    group_cols <- group
    
    print(group_cols)
        true_idx <- max.col(mat, ties.method = "first")
    
    # Create a logical matrix with only the highest-priority TRUE per row
    collapsed <- matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat))
    collapsed[cbind(seq_len(nrow(mat)), true_idx)] <- rowSums(mat) > 0
    
    # Replace the relevant columns in query_gr
    query_mat[, 4:ncol(query_mat)] <- collapsed
  }
  return(query_mat)
}

get_exclusive_annotations <- function(query_gr) {
  logical_cols <- sapply(mcols(query_gr), is.logical)
  logic_df <- as.data.frame(mcols(query_gr)[, logical_cols, drop = FALSE])
  
  # Keep only rows that have 0 or 1 TRUE (drop rows with >1 TRUE)
  keep_mask <- rowSums(logic_df) <= 1
  gr_subset <- query_gr[keep_mask]
  
  # Build annotation vector for the kept rows:
  # - rows with 1 TRUE -> the column name with TRUE
  # - rows with 0 TRUE -> "Null"
  kept_logic <- logic_df[keep_mask, , drop = FALSE]
  annotation <- rep("Null", nrow(kept_logic))
  one_true_mask <- rowSums(kept_logic) == 1
  if (any(one_true_mask)) {
    annotation[one_true_mask] <- names(kept_logic)[max.col(kept_logic[one_true_mask, , drop = FALSE], ties.method = "first")]
  }
  
  # Replace metadata with a single column "annotation"
  mcols(gr_subset) <- DataFrame(annotation = annotation)
  
  return(gr_subset)
}

process_excl<-function(anchor_gr, annotation_layers, collapse){
  annotated_anchors<-create_annotation_matrix(anchor_gr, annotation_layers)
  anchor_gr<-makeGRangesFromDataFrame(annotated_anchors,seqnames = 'seqnames', start.field = 'start', end.field = 'end',     
                                      keep.extra.columns = TRUE)
  anchors_collapsed<-collapse_gr(anchor_gr, collapse)
  exclusive_anchors<-get_exclusive_annotations(anchors_collapsed)
  
  return(exclusive_anchors)
}

get_exclusive_loops<-function(gi, annotation_layers,
                              order=c('E','P','CTCF','Null'), collapse=list(c('P','E'))){
  
  exclusive_anchors1<-process_excl(anchorOne(gi), annotation_layers, collapse)
  exclusive_anchors2<-process_excl(anchorTwo(gi), annotation_layers, collapse)
  
  # Match exclusive anchors back to original GInteractions
  hits1 <- findMatches(anchorOne(gi), exclusive_anchors1)
  hits2 <- findMatches(anchorTwo(gi), exclusive_anchors2)
  
  # Keep only loops where both anchors are exclusive
  common_idxs <- intersect(queryHits(hits1), queryHits(hits2))
  gi_exclusive <- gi[common_idxs]
  
  # attach exclusive annotations from both anchors
  anchor1_annot <- mcols(exclusive_anchors1)$annotation[subjectHits(hits1)[match(common_idxs, queryHits(hits1))]]
  anchor2_annot <- mcols(exclusive_anchors2)$annotation[subjectHits(hits2)[match(common_idxs, queryHits(hits2))]]
  
  # Resolve pair annotation with order if provided
  rank1 <- match(anchor1_annot, order)
  rank2 <- match(anchor2_annot, order)
  
  anno_low <- ifelse(rank1 <= rank2, anchor1_annot, anchor2_annot)
  anno_high <- ifelse(rank1 > rank2, anchor1_annot, anchor2_annot)
  pair <- paste(anno_low, anno_high, sep = "-")
  
  # Add all annotation info to gi
  mcols(gi_exclusive)$anchor1_annotation <- anchor1_annot
  mcols(gi_exclusive)$anchor2_annotation <- anchor2_annot
  mcols(gi_exclusive)$class <- pair
  
  return(gi_exclusive)
}
