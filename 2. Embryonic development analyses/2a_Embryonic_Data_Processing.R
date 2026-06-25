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

# Loading raw reads downloaded from the authors github repository
load("Tyser_et_al_Gastrula.RData")
seurat <- gastrula
rm(gastrula)
gc()
subset_embryo_colors <- c("Ectoderm" = "#66C2A5",
                          "Endoderm" = "#1F78B4",
                          "Mesoderm" = "#F22233",
                          "Primitive Streak" = "purple",
                          "Epiblast" = "skyblue")



seurat$Cell_Type <- seurat$cluster_id

seurat$Cell_Type <- gsub("ExE Mesoderm", "Extraembryonic Mesoderm", seurat$Cell_Type)

# Mapping to Broader Cell Types
Subset_Cell_Types_mapping <- list(
  "Endoderm" = "Endoderm",
  "Ectoderm" = "Non-Neural Ectoderm",
  "Epiblast" = "Epiblast",
  "Primitive Streak" = "Primitive Streak",
  "Hematopoietic Precursors" = c("Hemogenic Endothelial Progenitors", "Erythroblasts"),
  "Mesoderm" = c("Emergent Mesoderm", "Nascent Mesoderm", "Advanced Mesoderm", "Axial Mesoderm", "Extraembryonic Mesoderm")
)

assign_subset_cell_type <- function(cluster) {
  for (type in names(Subset_Cell_Types_mapping)) {
    if (cluster %in% Subset_Cell_Types_mapping[[type]]) {
      return(type)
    }
  }
  return(NA) # Return NA if the cluster doesn't match any category
}

# Applying the function to the 'Subset Cell_Type' column to create the 'Major Cell Type' column
seurat$Subset_Cell_Types <- sapply(seurat$Cell_Type, assign_subset_cell_type)

seurat <- SCTransform(seurat)
seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, dims = 1:20, reduction.name = "UMAP")

DimPlot_scCustom(seurat, group.by = "Subset_Cell_Types", colors_use = subset_embryo_colors)


# Only focusing on the Primitive Streak, Epiblasts and Three Germ layers for the trajectory analysis
important <- c("Endoderm", "Mesoderm", "Primitive Streak", 
               "Epiblast", "Ectoderm")

seurat <- subset(seurat, subset = Subset_Cell_Types %in% important)
seurat <- SCTransform(seurat)
seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, dims = 1:20, reduction.name = "UMAP")

DimPlot_scCustom(seurat, group.by = "Subset_Cell_Types", colors_use = subset_embryo_colors)

save(seurat, file="Subsetted_Tyser_Processed.RData")
load("Subsetted_Tyser_Processed.RData")

# Using monocle3 to create the trajectory
cds <- SeuratWrappers::as.cell_data_set(seurat,assay = "SCT")

set.seed(42) # Setting seeds

# Only making clusters
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds_clusters <- clusters(cds)
seurat$Monocle3_Clusters <- cds_clusters[Cells(seurat)]


cds <- learn_graph(cds, use_partition=FALSE, close_loop=FALSE)
root_plot <- plot_cells(cds, color_cells_by="cluster",
                        group_label_size=4, graph_label_size=3.5,
                        label_cell_groups=T, label_principal_points=TRUE,
                        label_groups_by_cluster=FALSE) 

cell_ids <- colnames(cds)[seurat$Subset_Cell_Types ==  "Epiblast"]

closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
closest_vertex <- closest_vertex[cell_ids, ]
closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
mst <- principal_graph(cds)$UMAP
root_pr_nodes <- igraph::V(mst)$name[closest_vertex]

cds <- order_cells(cds, root_pr_nodes = "Y_2")

##Visualization
library(viridis)
trajec_plot <- plot_cells(cds, color_cells_by = "cluster",label_groups_by_cluster = T,
                          label_roots = F,
                          label_leaves = F,
                          label_principal_points = F,
                          label_branch_points = F,
                          group_label_size = 5) 
trajec_plot

pseudo_plot <- plot_cells(cds,
                          color_cells_by = "pseudotime",
                          group_cells_by = "cluster",
                          label_cell_groups = T,
                          label_groups_by_cluster=T,
                          label_leaves=FALSE,
                          label_branch_points=FALSE,
                          label_roots = FALSE,
                          trajectory_graph_color = "red",
                          group_label_size = 5) 
pseudo_plot
ggsave("PseudoPlot.svg", height = 4, width = 6, pseudo_plot, bg = "transparent", dpi = 300)

cds$monocle3_pseudotime <- pseudotime(cds)

seurat <- AddMetaData(
  object = seurat,
  metadata = cds$monocle3_pseudotime,
  col.name = "Monocle3_Pseudotime"
)
embryo_colors <- c(lineage_colors, cell_type_colors, broad_cell_type_colors, Major_Cell_Type_colors,  Sample_ID_colors, subset_embryo_colors)

library(SCpubr)

p1 <- do_BoxPlot(seurat, feature = "Monocle3_Pseudotime", group.by = "Sample_ID", order = T, legend.position = "none")
p2 <- do_BoxPlot(seurat, feature = "Monocle3_Pseudotime", group.by = "Monocle3_Clusters", order = T, legend.position = "none")
p3 <- DimPlot_scCustom(seurat, group.by = "Monocle3_Clusters", ggplot_default_colors = T, figure_plot = T, reduction = "UMAP")
p4 <- FeaturePlot_scCustom(seurat, features = "Monocle3_Pseudotime",colors_use = rev(viridis_plasma_dark_high), figure_plot = T, reduction = "UMAP")
p5 <- DimPlot_scCustom(seurat, group.by = "Subset_Cell_Types", reduction = "UMAP", colors_use = subset_embryo_colors, figure_plot = F) + theme(legend.position = "none")
p6 <- DimPlot_scCustom(seurat, group.by = "Subset_Cell_Types", reduction = "UMAP", colors_use = subset_embryo_colors, figure_plot = T)

pdf("Monocle3_Tyser_Dim_Plots.pdf", height = 6, width = 8)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(pseudo_plot)
print(trajec_plot)
dev.off()


### Only focusing on three germ layers to check the genes expression
sub_seurat <- subset(seurat, subset = Subset_Cell_Types %in% c("Ectoderm", "Endoderm",
                                                               "Mesoderm"))

sub_seurat <- SetupForWGCNA(
  sub_seurat,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "layers" # the name of the hdWGCNA experiment
)

sub_seurat <- MetacellsByGroups(
  seurat_obj = sub_seurat,
  group.by = c("Subset_Cell_Types"), # specify the columns in seurat_obj@meta.data to group by
  min_cells = 20,
  k = 5, # nearest-neighbors parameter
  max_shared = 5, # maximum number of shared cells between two metacells
  ident.group = 'Subset_Cell_Types' # set the Idents of the metacell seurat object
)

seurat_metacell <- GetMetacellObject(sub_seurat)

save(seurat, subset_embryo_colors, sub_seurat, cds, seurat_metacell, file="Subsetted_Tyser_Processed_Seurat.RData")


### Boxplots of expression
seurat_metacell@meta.data$Subset_Cell_Types <- factor(
  seurat_metacell@meta.data$Subset_Cell_Types, 
  levels = c( "Endoderm", "Mesoderm", "Ectoderm")
)

proteo_genes <- intersect(proteoglycans$Gene, rownames(sub_seurat))
gags <- intersect(gag$Gene, rownames(sub_seurat))

features_list <- list(
  EMT_TFs = c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2"),
  Proteoglycans = sort(proteo_genes),
  Gags = sort(gags)
)
library(gridExtra)

pdf("S2G_Metacells_Germ_Layers_BoxPlots_Paged.pdf", width = 10, height = 10)
for (category in names(features_list)) {
  
  genes <- features_list[[category]]
  genes <- genes[genes %in% rownames(seurat_metacell)]
  
  if (length(genes) == 0) next
  
  expr_data <- FetchData(seurat_metacell, vars = genes)
  expr_data$Subset_Cell_Types <- seurat_metacell$Subset_Cell_Types
  
  long_df <- expr_data %>%
    pivot_longer(
      cols = all_of(genes),
      names_to = "Gene",
      values_to = "Expression"
    )
  
  long_df$Gene <- factor(long_df$Gene, levels = genes)
  
  # one full category per page
  p <- ggplot(long_df, aes(
    x = Subset_Cell_Types,
    y = Expression,
    fill = Subset_Cell_Types
  )) +
    geom_boxplot(
      outlier.shape = 16,
      outlier.size = 0.35,
      linewidth = 0.35,
      width = 0.65,
      color = "black"
    ) + stat_compare_means(
      comparisons = my_comparisons,
      method = "t.test",
      label = "p.signif",
      hide.ns = TRUE,
      size = 3
    ) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 7) +
    scale_fill_manual(values = subset_embryo_colors) +
    labs(
      title = category,
      x = "Germ Layer",
      y = "Expression"
    ) +
    theme_pubr(base_size = 9) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9),
      strip.text = element_text(face = "bold", size = 7),
      strip.background = element_rect(fill = "grey90", color = "grey75"),
      panel.background = element_rect(fill = "white", color = "white"),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.3),
      axis.line = element_line(linewidth = 0.4, color = "black")
    )
  
  print(p)
}

dev.off()


