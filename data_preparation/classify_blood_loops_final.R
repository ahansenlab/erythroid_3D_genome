source("loop_anchor_classification.R")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
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

hg38<-TxDb.Hsapiens.UCSC.hg38.knownGene
exons <- exons(TxDb.Hsapiens.UCSC.hg38.knownGene)
cds <- cds(TxDb.Hsapiens.UCSC.hg38.knownGene)
transcripts <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)

DATAPATH='/Users/varshiniramanathan/Documents/HansenLab/blood_data/remapped'
ENHANCER_FNAME='merged_enhancers_500bp.bed'
PROMOTER_FNAME='promoters_hg38.bed'
CTCF_FNAME='merged_CTCF_P123.bed'

enhancer = read_tsv(file.path(DATAPATH, ENHANCER_FNAME), col_names = c('chrom', 'start', 'end'), show_col_types = FALSE) |> 
  makeGRangesFromDataFrame()

promoter = read_tsv(file.path(DATAPATH, PROMOTER_FNAME), col_names = c('chrom', 'start', 'end', 'name', 'alt_tss', 'strand'), show_col_types = FALSE) |> 
  makeGRangesFromDataFrame()

ctcf = read_tsv(file.path(DATAPATH, CTCF_FNAME), col_names = c('chrom', 'start', 'end'), show_col_types = FALSE) |> 
  makeGRangesFromDataFrame()

annotation_layers<-GRangesList(P=promoter, E=enhancer, CTCF=ctcf)

# THIS GENERATES THE LOOP CLASSES FILE
###
LOOPS_FNAME='merged_calls_125kb_q0.01.bedpe'
ANCHORS_FNAME='merged_calls_2kb_st0.88_q0.02_merged.bed'
###

loops = import(file.path(DATAPATH, LOOPS_FNAME)) |> makeGInteractionsFromGRangesPairs()

excl_loops <- get_exclusive_loops(loops, annotation_layers)


loopnames<-c(
  "032524_combined_DAP_JAN_028_ATAC_v1_both.bedpe", 
  "erythroid_associated_ACRs_v1_both.bedpe",
  "acrs_unique_v1_both.bedpe",
  "erythroid_tf_acrs_unique_v1_both.bedpe")

loopnames<-c('merged_calls_125kb_q0.01.bedpe')

for (lc in loopnames){
  print(lc)
  loops_curr=import(file.path(DATAPATH, lc)) |> makeGInteractionsFromGRangesPairs()
  bn<-strsplit(lc, '.', fixed=TRUE)[[1]][1]
  annotated_gi<-get_exclusive_loops(loops_curr, annotation_layers)
  write.table(annotated_gi,
                           file = file.path(DATAPATH, sprintf("%s_excl_classes.txt", bn)),
                           sep = "\t",
                           quote = FALSE,
                           row.names = FALSE
              )
  }
#anchs = import(file.path(DATAPATH, ANCHORS_FNAME)) |> GRanges()
#anch_matrix = create_annotation_matrix(anchs, annotation_layers, mcols_use = list('transcript'='tx_name'))
#annotated_granges<-annotate_with_priority(anchs, GRangesList(P=promoter, E=enhancer, CTCF=ctcf))
#anch_matrix$annotation<-annotated_granges$annotation
# annotation_layers<-GRangesList(P=promoter, E=enhancer, CTCF=ctcf,
#                                exon=exons, 
#                                cds=cds, transcript=transcripts)

#annotated_gi<-annotate_ginteractions(loops, annotation_layers)
# write.table(annotated_gi,
#             file = file.path(DATAPATH, "merged_calls_125kb_q0.01_classes_new.txt"),
#             sep = "\t",
#             quote = FALSE,
#             row.names = FALSE
# )


# write.table(excl_loops,
#             file=file.path(DATAPATH,"merged_calls_125kb_q0.01_classes_excl.txt"),
#             sep="\t", quote=FALSE, row.names=FALSE)




#write.table(anch_matrix,
#  file = ANCHOR_OUTNAME,
#  sep = "\t",
#  quote = FALSE,
#  row.names = FALSE
#)