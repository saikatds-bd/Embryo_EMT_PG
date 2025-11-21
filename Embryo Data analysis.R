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
# Loading Colors
load("Proteoglycan_Project_Colors.RData")

# Loading raw reads downloaded from the authors github repository
load("Tyser_et_al_Gastrula.RData")
seurat <- gastrula
subset_embryo_colors <- c("Ectoderm" = "#F22233",
                          "Endoderm" = "#66C2A5",
                          "Mesoderm" = "sienna1",
                          "Epiblast" = "#1F78B4",
                          "Primitive Streak" = "#32CD32",
                          "Hematopoietic Precursors" ="#FB9A99")


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

dput(unique(seurat$Subset_Cell_Types))

seurat <- SCTransform(seurat)
seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, dims = 1:20, reduction.name = "UMAP")

DimPlot_scCustom(seurat, group.by = "Subset_Cell_Types", colors_use = subset_embryo_colors)
save(seurat, file="Processed_Tyler_Data.RData")
load(file="Processed_Tyler_Data.RData")
dput(unique(seurat$Subset_Cell_Types))
p6 <- DimPlot_scCustom(seurat, group.by = "Subset_Cell_Types", reduction = "UMAP", colors_use = subset_embryo_colors, figure_plot = F) + theme(legend.position = "none")
p7 <- DimPlot_scCustom(seurat, group.by = "Subset_Cell_Types", reduction = "UMAP", colors_use = subset_embryo_colors, figure_plot = T)
p7

pdf("Tyser_Dim_Plots.pdf", height = 6, width = 8)
print(p6)
print(p7)
dev.off()

dput(unique(seurat$Subset_Cell_Types))

# Only focusing on the Primitive Streak, Epiblasts and Three Germ layers
important <- c("Endoderm", "Mesoderm", "Primitive Streak", 
               "Epiblast", "Ectoderm")

seurat <- subset(seurat, subset = Subset_Cell_Types %in% important)
seurat <- SCTransform(seurat)
seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, dims = 1:20, reduction.name = "UMAP")

DimPlot_scCustom(seurat, group.by = "Subset_Cell_Types", colors_use = subset_embryo_colors)
save(seurat, file="Subsetted_Tyler_Processed.RData")


cds <- SeuratWrappers::as.cell_data_set(seurat,assay = "SCT")
set.seed(42)

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
p3 <- DimPlot_scCustom(seurat, group.by = "Sample_ID", colors_use = Sample_ID_colors, figure_plot = T, reduction = "UMAP")
p4 <- DimPlot_scCustom(seurat, group.by = "Monocle3_Clusters", ggplot_default_colors = T, figure_plot = T, reduction = "UMAP")
p5 <- FeaturePlot_scCustom(seurat, features = "Monocle3_Pseudotime",colors_use = rev(viridis_plasma_dark_high), figure_plot = T, reduction = "UMAP")
p6 <- DimPlot_scCustom(seurat, group.by = "Subset_Cell_Types", reduction = "UMAP", colors_use = embryo_colors, figure_plot = F) + theme(legend.position = "none")
p7 <- DimPlot_scCustom(seurat, group.by = "Subset_Cell_Types", reduction = "UMAP", colors_use = embryo_colors, figure_plot = T)

pdf("Subset_Monocle3_Tyser_Dim_Plots.pdf", height = 6, width = 8)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(pseudo_plot)
print(trajec_plot)
dev.off()

### Creating Metacell
library(hdWGCNA)

DefaultAssay(seurat) <- "RNA"
Idents(seurat) <- "Subset_Cell_Types"


seurat <- SetupForWGCNA(
  seurat,
  gene_select = "fraction", 
  fraction = 0.05, 
  wgcna_name = "gastrulation" 
)

seurat <- MetacellsByGroups(
  seurat_obj = seurat,
  group.by = "Subset_Cell_Types", 
  k = 10, 
  max_shared = 5, 
  ident.group = 'Subset_Cell_Types',
  min_cells = 20
)

# Working with Metacell object for visualization and further processing
seurat_metacell <- GetMetacellObject(seurat)

seurat_metacell <- SCTransform(seurat_metacell)

seurat_metacell <- RunPCA(seurat_metacell)
seurat_metacell <- RunUMAP(seurat_metacell, reduction.name = "UMAP", dims = 1:20)
DimPlot_scCustom(seurat_metacell, group.by = "Subset_Cell_Types", colors_use = subset_embryo_colors)
dim(seurat)
dim(seurat_metacell)
table(seurat_metacell$Subset_Cell_Types)
seurat_metacell@meta.data$Subset_Cell_Types <- factor(
  seurat_metacell@meta.data$Subset_Cell_Types, 
  levels = c("Epiblast", "Primitive Streak", "Endoderm", "Mesoderm", "Ectoderm")
)

# Loading the proteoglycan and gags
load("Interesting_GAG_Genes.RData")

# To make sure that all the genes are already present in seurat object
proteo_genes <- intersect(proteo$Erik.et.al, rownames(seurat))
gags <- intersect(gag$genes, rownames(seurat))

features_list <- list(
  EMT_TFs = c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2"),
  Proteoglycans = sort(proteo_genes),
  Gags = sort(gags)
)
library(gridExtra)

pdf("Metacells_Germ_Layers_BoxPlots_Paged.pdf", width = 12, height = 8)

for (category in names(features_list)) {
  genes <- features_list[[category]]
  
  # Filter only those present in the data
  genes <- genes[genes %in% rownames(seurat_metacell)]
  
  if (length(genes) == 0) next  # skip if no genes are found
  
  # Prepare expression data
  expr_data <- FetchData(seurat_metacell, vars = genes)
  expr_data$Subset_Cell_Types <- seurat_metacell$Subset_Cell_Types
  long_df <- expr_data %>%
    pivot_longer(
      cols = all_of(genes),
      names_to = "Gene",
      values_to = "Expression"
    )
  
  long_df$Gene <- factor(long_df$Gene, levels = genes)
  
  # Split genes into chunks of 6
  gene_chunks <- split(genes, ceiling(seq_along(genes) / 6))
  
  # Loop over chunks to create 6-gene-per-page plots
  for (chunk in gene_chunks) {
    chunk_df <- long_df %>% filter(Gene %in% chunk)
    chunk_df$Gene <- factor(chunk_df$Gene, levels = chunk)
    
    p <- ggplot(chunk_df, aes(x = Subset_Cell_Types, y = Expression, fill = Subset_Cell_Types)) +
      geom_boxplot(outlier.shape = NA, width = 0.15, color = "Black") +
      facet_wrap(~ Gene, scales = "free_y") +
      scale_fill_manual(values = lineage_colors) +
      theme_bw(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.position = "none"
      ) +
      labs(title = category, x = "Germ Layer", y = "Expression")
    
    print(p)
  }
}

dev.off()

save(seurat, seurat_metacell, cds, file= "Tyser_processed_files.RData")
