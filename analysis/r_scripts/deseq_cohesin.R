library("biomaRt")
library(tidyverse)
library("R.utils")
library("DESeq2")
library(dplyr)
library("pheatmap")
library(tximport)

library(AnnotationDbi)
library(org.Hs.eg.db) 
library(EnhancedVolcano)

DATAPATH <- '/Users/varshiniramanathan/Documents/HansenLab/blood_data'
filename="rmdup_counts_v1.txt"

# featureCounts outputs a comment on the first line, so we ignore it
fc_data <- read.table(file.path(DATAPATH, filename), header=TRUE, row.names=1, comment.char="#", stringsAsFactors=FALSE)

## process counts matrix. all the conditions were sequenced together 
inds<-c(10, 11, 13, 14, 15, 17, 18, 20, 21, 22, 23, 24)
SAMPLES_TO_INCLUDE=colnames(fc_data)[inds]
count_matrix <- fc_data[, 6:ncol(fc_data)]

DONOR_REPLICATES <- c("D1", "C2", "C3", "C2", "D2", "C3", "C3", "D2","C2","D1", "C2", "C3")

CONDITIONS <- c("AAVS1_P3","AAVS1","STAG2","STAG2",
                "AAVS1_P3","NIPBL",
                "STAG1","NFE2_P3",
                "NIPBL","NFE2_P3", "STAG1",
                "AAVS1")


CONTRASTS <- list(
  c("condition", "NIPBL", "AAVS1") 
)

count_matrix <- fc_data[, SAMPLES_TO_INCLUDE]
contrast_nipbl<-c("condition", "NIPBL", "AAVS1") 
is_exp1 <- DONOR_REPLICATES %in% c("D1", "D2")

## split by cohesin vs NFE2 so that the replicate variable isn't confounded
cohesin_reps <- DONOR_REPLICATES[!is_exp1]
cohesin_conditions <- CONDITIONS[!is_exp1]
cohesin_samples <- SAMPLES_TO_INCLUDE[!is_exp1]
cohesin_counts <- count_matrix[, cohesin_samples]

col_data_cohesin <- data.frame(
  row.names = cohesin_samples,
  donor = factor(cohesin_reps), # donor is just for consistency here, they weren't actually donors 
  condition = factor(cohesin_conditions)
)

col_data_cohesin$donor <- droplevels(col_data_cohesin$donor)
col_data_cohesin$condition <- droplevels(col_data_cohesin$condition)

dds_cohesin <- DESeqDataSetFromMatrix(countData = cohesin_counts,
                                      colData = col_data_cohesin,
                                      design = ~ donor + condition)
smallestGroupSize <- 3
keep <- rowSums(counts(dds_cohesin) >= 10) >= smallestGroupSize
dds_cohesin <- dds_cohesin[keep,]
dds_cohesin <- DESeq(dds_cohesin)

# NIPBL
contrast_curr<-contrast_nipbl

res_volcano <- results(dds_cohesin, contrast = contrast_curr, alpha = 0.05)
res_volcano <- lfcShrink(dds_cohesin, contrast=contrast_curr, type="ashr")
rownames(res_volcano) <- gsub("\\..*$", "", rownames(res_volcano))

down<-res_volcano[which(res_volcano$log2FoldChange < 0 & res_volcano$padj<0.05), ]
up<-res_volcano[which(res_volcano$log2FoldChange > 0 & res_volcano$padj<0.05), ]

down<-down[order(down$padj), ]
up<-up[order(up$padj), ]

write.table(down, file=file.path(DATAPATH, sprintf("%s_%s_down_lfcShrink_sepExp.txt", contrast_curr[2], contrast_curr[3])), sep='\t', quote=F)
write.table(up, file=file.path(DATAPATH, sprintf("%s_%s_up_lfcShrink_sepExp.txt", contrast_curr[2], contrast_curr[3])), sep='\t', quote=F)
