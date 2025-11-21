library(Seurat)
library(tidyverse)
library(limma)
library(clustree)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(scCustomize)
library(SCpubr)
library(patchwork)
library(scCustomize)
library(reticulate)
library(harmony)
library(msigdbr)
library(monocle3)

set.seed(123)
emt_tfs <- c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2")
# Loading Colors
load("Proteoglycan_Project_Colors.RData")
load("Tyser_processed_files.RData")
load("Interesting_GAG_Genes.RData")

## Separating the endoderm and mesoderm lineages

# Endoderm
cds_endo <- choose_graph_segments(cds, clear_cds = F)
root_plot <- plot_cells(cds_endo, color_cells_by="cluster",
                        group_label_size=4, graph_label_size=3.5,
                        label_cell_groups=T, label_principal_points=TRUE,
                        label_groups_by_cluster=FALSE) 
root_plot
fData(cds_endo)$gene_short_name <- rownames(fData(cds_endo))

cds_endo <- order_cells(cds_endo, root_pr_nodes = "Y_2")
cds_endo$Pseudotime_Endo <- pseudotime(cds_endo)
endo_cells <- colnames(cds_endo)

pt_endo <- matrix(monocle3::pseudotime(cds_endo), ncol=1,
                  dimnames=list(endo_cells, "Endo"))
cw_endo <- matrix(1, nrow=length(endo_cells), ncol=1,
                  dimnames=list(endo_cells, "Endo"))

## Mesoderm 
cds_meso <- choose_graph_segments(cds, clear_cds = F)
root_plot <- plot_cells(cds_meso, color_cells_by="cluster",
                        group_label_size=4, graph_label_size=3.5,
                        label_cell_groups=T, label_principal_points=TRUE,
                        label_groups_by_cluster=FALSE) 
root_plot
cds_meso <- order_cells(cds_meso, root_pr_nodes = "Y_2")
cds_meso$Pseudotime_Meso<- pseudotime(cds_meso)
fData(cds_meso)$gene_short_name <- rownames(fData(cds_meso))

meso_cells <- colnames(cds_meso)
pt_meso <- matrix(monocle3::pseudotime(cds_meso), ncol=1,
                  dimnames=list(meso_cells, "Meso"))
cw_meso <- matrix(1, nrow=length(meso_cells), ncol=1,
                  dimnames=list(meso_cells, "Meso"))
## --- Build lineage-specific count matrices (UMI counts!) ------------------
counts_endo <- seurat[["RNA"]]$counts[, endo_cells, drop=FALSE]
counts_meso <- seurat[["RNA"]]$counts[, meso_cells, drop=FALSE]

## restrict to genes present
goi <- intersect(unique(c(emt_tfs, proteoglycans$Name, gag$genes)), rownames(counts_meso))
genes <- goi

library(tradeSeq)

## Fit tradeSeq per lineage

sce_endo <- fitGAM(counts = counts_endo,
                   pseudotime = pt_endo,
                   cellWeights = cw_endo,
                   nknots = 6,
                   verbose = TRUE,
                   gene = genes,
                   BPPARAM = BiocParallel::MulticoreParam(workers = 1))

sce_meso <- fitGAM(counts = counts_meso,
                   pseudotime = pt_meso,
                   cellWeights = cw_meso,
                   nknots = 6,
                   verbose = TRUE,
                   gene = genes,
                   BPPARAM = BiocParallel::MulticoreParam(workers = 1))

## Association tests from Tradeseq

asso_endo <- associationTest(sce_endo) %>% as_tibble(rownames="gene") %>%
  mutate(padj = p.adjust(pvalue, "fdr"))

asso_endo$Type <- ifelse(asso_endo$gene %in% proteoglycans$Name, "Proteoglycan",
                         ifelse(asso_endo$gene %in% gag$genes, "GAG", "EMT TFs"))

writexl::write_xlsx(asso_endo, "Association_towards_Endoderm.xlsx")

asso_meso <- associationTest(sce_meso) %>% as_tibble(rownames="gene") %>%
  mutate(padj = p.adjust(pvalue, "fdr"))

asso_meso$Type <- ifelse(asso_meso$gene %in% proteoglycans$Name, "Proteoglycan",
                         ifelse(asso_meso$gene %in% gag$genes, "GAG", "EMT TFs"))

writexl::write_xlsx(asso_meso, "Association_towards_Mesoderm.xlsx")

save(sce_endo, cds_endo, sce_meso, cds_meso, asso_meso, asso_endo, counts_endo, counts_meso,
     proteoglycans, gag, emt_tfs, file = "Separate_CDS.RData")