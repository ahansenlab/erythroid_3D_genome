library("GenomicFeatures")
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(ChIPseeker)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(plyranges)
library(tidyverse)
library("R.utils")

hg38_txdb <-TxDb.Hsapiens.UCSC.hg38.refGene  # shorthand (for convenience)

DATAPATH='~/Documents/HansenLab/blood_data'
FNAME='quant3_100kb_merge10.0kb_clus.bed'

curr_file = read_tsv(file.path(DATAPATH, FNAME), col_names = c('chrom', 'start', 'end'), show_col_types = FALSE) |> 
  makeGRangesFromDataFrame()
peakAnno <- annotatePeak(curr_file, tssRegion=c(-2000, 1000),
                         TxDb=hg38_txdb, annoDb="org.Hs.eg.db")

pie<-plotAnnoPie(peakAnno)
pie

peakprof<-plotPeakProf2(curr_file, upstream = rel(1), downstream = rel(1),
                       conf = 0.95, by = "gene", type = "body", nbin = 800,
                       TxDb = hg38_txdb, ignore_strand = T)
peakprof

promoter <- getPromoters(TxDb=hg38_txdb, upstream=2000, downstream=1000)
tagMatrix<-getTagMatrix(curr_file, windows=promoter)
plotAvgProf(tagMatrix, xlim=c(-2000, 1000), conf = 0.95, resample = 500)