#library(scater)
library("Seurat")
library(tidyverse)
library(cowplot)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)

library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)

DATAPATH='~/Documents/HansenLab/blood_data/'

seurat <- readRDS(file.path(DATAPATH, 'ChromVar_NT_only_JAN_028_annotated.rds'))
metadata <- seurat@meta.data
counts <- seurat@assays$RNA@counts 

Idents(seurat) <- seurat@meta.data$new_CellType

new_annotations <- c("MEP" = "P1", "BFU-E" = "P1", "CFU-E" = "P1", 
                     "Pro-Erythroblast" = "P2", "Basophilic Erythroblast" = "P2",
                     "Orthochromatic Erythroblast" = "P3", "Polychromatic Erythroblast" = "P3")

seurat <- RenameIdents(seurat, new_annotations)
seurat$phase_groups <- Idents(seurat)
DefaultAssay(seurat) <- "RNA"

### bulk expression 

bulk_expr <- AggregateExpression(
  seurat,
  return.seurat = TRUE,
  assays = "RNA",
  group.by = c("phase_groups")
)

#write.table(bulk_expr[["RNA"]]$data,  file.path(DATAPATH, "all_bulk_expr.tsv"), sep="\t", quote=F)

## then, do aggregated by replicate 
meta_columns <- c("replicate", "phase_groups")

bulk <- AggregateExpression(
  seurat,
  return.seurat = TRUE,
  assays = "RNA",
  group.by = c("phase_groups", "replicate")
)
Idents(bulk) <- "phase_groups"

n_cells <- seurat@meta.data %>% 
  dplyr::count(phase_groups, replicate) 

meta_bulk <- left_join(bulk@meta.data, n_cells)
rownames(meta_bulk) <- meta_bulk$orig.ident
bulk@meta.data <- meta_bulk

# Turn condition into a factor
bulk$condition <- factor(bulk$phase_groups, levels=c("P1", "P2", "P3"))
bulk@meta.data %>% head()

bulk_phase <- subset(bulk, subset=  (condition %in% c("P1", "P2", "P3")))


ggplot(bulk_phase@meta.data, aes(x=orig.ident, y=n, fill=condition)) +
  geom_bar(stat="identity", color="black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x="Sample name", y="Number of cells") +
  geom_text(aes(label=n), vjust=-0.5)


bulk_phase_filtered <- subset(
  bulk_phase,
  subset = n >= 100
)
cluster_counts <- FetchData(bulk_phase_filtered, layer="counts", vars=rownames(bulk_phase))


# Create DESeq2 object

# transpose it to get genes as rows
dds <- DESeqDataSetFromMatrix(t(cluster_counts),
                              colData = bulk_phase_filtered@meta.data,
                              design = ~ condition)

## make pca of replicate data for sanity check
rld <- rlog(dds, blind=TRUE)
pca_data_condition <- plotPCA(rld, intgroup=c("condition"), returnData = TRUE) 

ggplot(pca_data_condition, aes(x = PC1, y = PC2, color = condition, label = name_parsed)) +
  geom_point() + 
  geom_text(vjust = 1.5, hjust = 0.5, show.legend = FALSE) +
  theme_classic() +
  xlab(paste0("PC1: ", round(attr(pca_data_condition, "percentVar")[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(attr(pca_data_condition, "percentVar")[2] * 100), "% variance")) 


## now run differential expression
dds <- DESeq(dds)

CONTRASTS <- list(c("condition", "P3", "P1")
)
for (contrast_curr in CONTRASTS) {
  res <- results(dds, 
                 contrast=contrast_curr,
                 alpha = 0.05)
  
  padj.cutoff <- 0.05
  
  # Turn the results object into a tibble for use with tidyverse functions
  dge_deseq2 <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble()
  
  # Subset the significant results
  sig_res <- dplyr::filter(dge_deseq2, 
                           padj < padj.cutoff)


}



### findMarkers single cell 
seurat <- subset(seurat,idents=c("P1", "P2", "P3", "EoBasoMast Precursor"))

p3_markers <- FindMarkers(seurat,
                          ident.1="P3",
                          only.pos = TRUE
)
p3_markers <- p3_markers %>% subset(p_val_adj < 0.05)

ebm_markers  <- FindMarkers(seurat,
                            ident.1="EoBasoMast Precursor",
                            only.pos = TRUE
)
ebm_markers <- ebm_markers %>% subset(p_val_adj < 0.05)



dge_p32 <- FindMarkers(seurat_norm,
                       ident.1="P3",
                       ident.2="P2",
                       only.pos = TRUE
)
dge_p32_sig <- dge_p32 %>% subset(p_val_adj < 0.05)

dge_p21 <- FindMarkers(seurat_norm,
                       ident.1="P2",
                       ident.2="P1",
                       only.pos = TRUE
)
dge_p21_sig <- dge_p21 %>% subset(p_val_adj < 0.05)


## save significant results 
#write.csv(dge_p21_sig, file.path(DATAPATH, "findmarkers_p21.csv"))
#write.csv(dge_p31_sig, file.path(DATAPATH, "findmarkers_p31.csv"))
#write.csv(dge_p32_sig, file.path(DATAPATH, "findmarkers_p32.csv"))
#write.csv(p3_markers, file.path(DATAPATH, "findmarkers_p3.csv"))
#write.table(ebm_markers, file.path(DATAPATH, "findmarkers_ebm.tsv"), sep="\t", quote=F)
#write.table(p3_markers, file.path(DATAPATH, "findmarkers_p3_ebm.tsv"), sep="\t", quote=F)

