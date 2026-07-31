library("biomaRt")
library(tidyverse)
library("R.utils")
library("DESeq2")
library(dplyr)
library("pheatmap")
library(tximport)

library(AnnotationDbi)
library(org.Hs.eg.db) # Human annotation database
library(EnhancedVolcano)

DATAPATH <- '/Users/varshiniramanathan/Documents/HansenLab/blood_data'
filename="rmdup_counts_v2.txt"

# featureCounts outputs a comment on the first line, so we ignore it
fc_data <- read.table(file.path(DATAPATH, filename), header=TRUE, row.names=1, comment.char="#", stringsAsFactors=FALSE)
inds<-c(9, 12, 15, 16)
SAMPLES_TO_INCLUDE=colnames(fc_data)[inds]
count_matrix <- fc_data[, 6:ncol(fc_data)]

DONOR_REPLICATES <- c("D1", "D2", "D2", "D1")

CONDITIONS <- c("AAVS1_P3",
                "AAVS1_P3",
                "NFE2_P3",
                "NFE2_P3")


CONTRASTS <- list(
  c("condition", "NFE2_P3", "AAVS1_P3")
)

count_matrix <- fc_data[, SAMPLES_TO_INCLUDE]


contrast_nfe2<-c("condition", "NFE2_P3", "AAVS1_P3")
contrast_nipbl<-c("condition", "NIPBL", "AAVS1") 
is_exp1 <- DONOR_REPLICATES %in% c("D1", "D2")

nfe2_reps <- DONOR_REPLICATES[is_exp1]
nfe2_conditions <- CONDITIONS[is_exp1]
nfe2_samples <- SAMPLES_TO_INCLUDE[is_exp1]
nfe2_counts <- count_matrix[, nfe2_samples]

col_data_nfe2 <- data.frame(
  row.names = nfe2_samples,
  donor = factor(nfe2_reps),
  condition = factor(nfe2_conditions)
)

col_data_nfe2$donor <- droplevels(col_data_nfe2$donor)
col_data_nfe2$condition <- droplevels(col_data_nfe2$condition)

dds_nfe2 <- DESeqDataSetFromMatrix(countData = nfe2_counts,
                                   colData = col_data_nfe2,
                                   design = ~ donor + condition)
smallestGroupSize <- 3
keep <- rowSums(counts(dds_nfe2) >= 10) >= smallestGroupSize
dds_nfe2 <- dds_nfe2[keep,]
dds_nfe2 <- DESeq(dds_nfe2)

contrast_curr<-contrast_nfe2

res_volcano <- results(dds_nfe2, contrast = contrast_curr, alpha = 0.05)
res_volcano <- lfcShrink(dds_nfe2, contrast=contrast_curr, type="ashr")
rownames(res_volcano) <- gsub("\\..*$", "", rownames(res_volcano))

down<-res_volcano[which(res_volcano$log2FoldChange<0 & res_volcano$padj<0.05), ]
up<-res_volcano[which(res_volcano$log2FoldChange>0 & res_volcano$padj<0.05), ]

down<-down[order(down$padj), ]
up<-up[order(up$padj), ]

write.table(down, file=file.path(DATAPATH, sprintf("%s_%s_down_lfcShrink_sepExp2.txt", contrast_curr[2], contrast_curr[3])), sep='\t', quote=F)
write.table(up, file=file.path(DATAPATH, sprintf("%s_%s_up_lfcShrink_sepExp2.txt", contrast_curr[2], contrast_curr[3])), sep='\t', quote=F)

