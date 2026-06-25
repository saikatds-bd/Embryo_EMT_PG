library(tidyverse)
library(data.table)
library(ggthemes)
library(ComplexHeatmap)
library(dendsort)
library(circlize)
library(RColorBrewer)
library(qs2)

proteoglycan_df <- readxl::read_xlsx("Proteoglycans_GAGs_List.xlsx", sheet = "Proteoglycans")
gag_df <- readxl::read_xlsx("Proteoglycans_GAGs_List.xlsx", sheet = "GAG")

emt_tf <- c("SNAI1", "SNAI2", "ZEB1", "ZEB2","TWIST1", "TWIST2")

data <- fread("rna_single_cell_type_germ_layer.csv",
              drop = 1) # Adult tissue data downloaded. File in long format. 

length(unique(data$`Cell type`))


#Subset the data based on the gene of interest 
genes <- proteoglycan_df$Gene

genes_1 <- c(proteoglycan_df$Gene, emt_tf) # Combining with the emt tfs

data.subset <- data %>% dplyr::filter(`Gene name` %in% genes_1) # Making a subset of the expression data with our genes of interest.
length(unique(data.subset$`Gene name`)) # 48
length(unique(data.subset$`Cell type`)) # 81


####################
## log2 transform ##
####################

data.subset <- data.subset %>%
  group_by(`Cell type`,`Gene name`) %>%
  summarize(Expression = nTPM) 

data.subset <- data.subset %>% select(all_of(c("Gene name", "Cell type", "Expression")))

data.subset <- pivot_wider(data.subset, 
                           names_from = `Cell type`, 
                           values_from = Expression) # Converting to wide format for heatmap


data.subset <- column_to_rownames(data.subset, "Gene name")

data.subset <- log2(data.subset+1)


# Expression of the EMT TFs in cell types
ann1 <- subset(data.subset, rownames(data.subset) %in% emt_tf)
ann1 <- as.data.frame(t(ann1))
ann1 <- rownames_to_column(ann1, "Cell type")

# Making annotation for the Cell types with Origin 

ann2 <- data %>% select(all_of(c("Cell type", "origin"))) %>% distinct(`Cell type`, .keep_all = T)

data.subset <- data.subset %>% dplyr::select(ann2$`Cell type`)
all(names(data.subset) == ann2$`Cell type`) ## Confirming the orders match

ann2 <- column_to_rownames(ann2, "Cell type")
ann1 <- column_to_rownames(ann1, "Cell type")
ann1 <- ann1[colnames(expr),]

#Heatmap will only include the proteoglycans
expr <- subset(data.subset, rownames(data.subset) %in% genes)

ann3 <- proteoglycan_df
ann3$Category <- "Proteoglycan"
ann3 <- column_to_rownames(ann3, "Gene")

ann3 <- ann3[rownames(expr),]


all(rownames(ann1) == colnames(expr))
all(rownames(ann2) == colnames(expr))
all(rownames(ann3) == rownames(expr))


#####################################################

library(clValid)

internal_cols <- clValid(t(expr),
                         clMethods = "hierarchical",
                         nClust = 2:10,
                         validation = "internal",
                         method = "ward",
                         metric = "correlation")

## Cluster GENES (original rows)
internal_rows <- clValid(expr,
                         clMethods = "hierarchical",
                         nClust = 2:10,
                         validation = "internal",
                         method = "ward",
                         metric = "correlation")

pdf("Ward_Proteoglycan_HPA_Stability_CellType.pdf", height = 3, width = 4.5)
plot(internal_cols, legend = FALSE)
dev.off()

pdf("Ward_Protoglycan_HPA_Stability_Genes.pdf", height = 3, width = 4.5)
plot(internal_rows, legend = FALSE)
dev.off()

### Making a dataframe for later use
expr_df <- merge(ann3, expr, by = 0)
writexl::write_xlsx(expr_df, "Proteoglycans_in_adult_tissue.xlsx")

h <- as.matrix(expr)

## To check if we have any zero variance genes
# Compute row variances
row_variances <- apply(h, 1, var)

# Identify rows with zero variance
zero_variance_rows <- which(row_variances == 0)

# Display zero variance rows
print(zero_variance_rows)

column_variances <- apply(h, 2, var)

# Identify columns with zero variance
zero_variance_columns <- which(column_variances == 0)

# Display zero variance columns
print(zero_variance_columns)


h_breaks <- seq(min(h), max(h), length.out = 11)
h_breaks

ann1$SNAI1 <- as.numeric(ann1$SNAI1)
a <- ann1$SNAI1 
a_breaks <- seq(min(a), max(a), length.out = 10)
a_breaks

ann1$SNAI2 <- as.numeric(ann1$SNAI2)
b <- ann1$SNAI2 
b_breaks <- seq(min(b), max(b), length.out = 10)
b_breaks

ann1$ZEB1 <- as.numeric(ann1$ZEB1)
d <- ann1$ZEB1 
d_breaks <- seq(min(d), max(d), length.out = 10)
d_breaks

ann1$ZEB2 <- as.numeric(ann1$ZEB2)
e <- ann1$ZEB2
e_breaks <- seq(min(e), max(e), length.out = 10)
e_breaks

ann1$TWIST1 <- as.numeric(ann1$TWIST1)
f <- ann1$TWIST1 
f_breaks <- seq(min(f), max(f), length.out = 10)
f_breaks

ann1$TWIST2 <- as.numeric(ann1$TWIST2)
g <- ann1$TWIST2
g_breaks <- seq(min(g), max(g), length.out = 10)
g_breaks

#Reposition the breaks in quantile positions
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

h_breaks <- quantile_breaks(h, n = 11)
h_breaks

a_breaks <- quantile_breaks(a, n = 11)
a_breaks

b_breaks <- quantile_breaks(b, n = 11)
b_breaks

d_breaks <- quantile_breaks(d, n = 11)
d_breaks

e_breaks <- quantile_breaks(e, n = 11)
e_breaks

f_breaks <- quantile_breaks(f, n = 11)
f_breaks

g_breaks <- quantile_breaks(g, n = 11)
g_breaks

#Make clusters based on their pearson correlation
row_dist <- as.dist(1-cor(t(h), method = "pearson"))
col_dist <- as.dist(1-cor(h, method = "pearson"))
col_hc <- hclust(col_dist, method = "ward.D2")
row_hc <- hclust(row_dist, method = "ward.D2")


#We can use the dendsort package and reorder the clustering
Rowv=dendsort(as.dendrogram(row_hc), isRevers=TRUE, type = "average") 
Colv=dendsort(as.dendrogram(col_hc), type = "average")

#Using custom colors
cl <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)

RdBu_ramp <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))

RdBu_colors <- RdBu_ramp(9)

col_fun <- colorRamp2(quantile_breaks(h, n = 11), RdBu_colors)

ha = HeatmapAnnotation(foo = 1:10, col = list(foo = col_fun))

ha <- HeatmapAnnotation("SNAI1"= ann1$SNAI1,
                        "SNAI2" = ann1$SNAI2,
                        "ZEB1" = ann1$ZEB1,
                        "ZEB2" = ann1$ZEB2,
                        "TWIST1" = ann1$TWIST1,
                        "TWIST2" = ann1$TWIST2,
                        "Germ Layer" = ann2$origin,
                        col = list("SNAI1"= col_fun,
                                   "SNAI2" = col_fun,
                                   "ZEB1" = col_fun,
                                   "ZEB2" = col_fun,
                                   "TWIST1" = col_fun,
                                   "TWIST2" = col_fun,
                                   "Germ Layer" = c(
                                     "Ectoderm" = "#66C2A5",
                                     "Endoderm" = "#1F78B4",
                                     "Mesoderm" = "#F22233",
                                     "Trophoblast" = "purple3", 
                                     "Varies" = "orange3"
                                   )),
                        simple_anno_size = unit(0.50, "cm"),
                        annotation_name_side = "right",
                        border = T,
                        annotation_name_gp= gpar(fontsize = 8,fontface="bold"),
                        show_legend = T,
                        annotation_legend_param = list(title_gp = gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=8)))

ha_2 = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#1F78B4","#F22233","#66C2A5")),
                                          labels = c("Endoderm", "Mesoderm", "Ectoderm"), 
                                          labels_gp = gpar(col = "white", fontsize = 18)))

ann3[is.na(ann3)] = "NA"

row_ha <- rowAnnotation(
  "GAG Chain" = ann3$Chain,
  "Location" = ann3$Location,
  col = list(  "GAG Chain" = c(
    "CS" = "orange3",
    "KS" = "#F22233",
    "CS/DS" = "#1F78B4", 
    "CS/KS" = "purple3", 
    "CS/HS" = "#66C2A5",
    "NA" = "grey"),
    "Location" = c("Pericellular" = "orange3", 
                   "Extracellular" = "#F22233", 
                   "Cell surface" = "#66C2A5", 
                   "Intracellular" = "#1F78B4"
    )),
  simple_anno_size = unit(0.50, "cm"),
  border = T,
  annotation_name_gp= gpar(fontsize = 8,fontface="bold"),
  show_legend = T,
  annotation_legend_param = list(title_gp = gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=8)))



ht <- Heatmap(h, col = col_fun,
              cluster_columns = Colv,
              name = "Expression Values",
              show_heatmap_legend = T,
              top_annotation = ha,
              right_annotation = row_ha,
              show_column_names = F,
              show_row_names = T,cluster_rows = Rowv,
              row_title_gp = gpar(fontsize=8),
              column_title_gp = gpar(fontsize=8),
              height = unit(13, "cm"),
              width = unit(13, "cm"),
              column_dend_reorder = T,
              row_dend_reorder = T,
              show_row_dend = T,
              border = T,
              column_dend_height = unit(2, "cm"),
              column_names_rot = 90,
              row_names_gp = grid::gpar(fontsize= 9),
              column_names_gp = grid::gpar(fontsize = 8),
              heatmap_legend_param = list(title="Log2 + 1 nTPM", legend_height=unit(3, "cm"),title_gp=gpar(fontsize=9, fontface="bold"),labels_gp = gpar(fontsize=8)),
              row_split = 3, column_split = 3,
              column_gap = unit(c(3,3,3), "mm"),
              row_gap = unit(c(3,3,3), "mm"))

draw(ht,column_title = "Proteoglycans in Normal Tissue Cell Types", column_title_gp = gpar(fontsize = 8, "bold"), 
     merge_legend=TRUE, padding = unit(c(2, 2, 2, 2), "mm"))



ht <- draw(ht, column_title_gp = gpar(fontsize = 8), 
           merge_legend=FALSE, padding = unit(c(2, 2, 2, 2), "mm"))

pdf("Ward_HPA_Proteoglycans.pdf",  width=20,height=18)
print(ht)
dev.off()

#For getting the gene orders
for (i in 1:length(row_order(ht))){   if (i == 1) {
  clu <- t(t(row.names(expr[row_order(ht)[[i]],])))
  out <- cbind(clu, paste("Cluster ", i, sep=""))
  colnames(out) <- c("GeneID", "Heatmap Cluster")   } else {
    clu <- t(t(row.names(expr[row_order(ht)[[i]],])))
    clu <- cbind(clu, paste("Cluster ", i, sep=""))
    out <- rbind(out, clu)   } 
}
out

out <- as.data.frame(out)
rownames(out) <- out$GeneID
ann3_merged <- merge(ann3, out, by = 0)
#ann3_merged <- ann3_merged[,c(10,3,4:7, 11)]
write.csv2(ann3_merged, file = "WARD_HPA_Proteo_Gene_Order.csv")

#For getting the cell type orders

for (i in 1:length(column_order(ht)))   if (i == 1) {
  clu <- t(t(colnames(expr[,column_order(ht)[[i]]])))
  out <- cbind(clu, paste("Cluster ", i, sep=""))
  colnames(out) <- c("Cell Type", "Heatmap Cluster")   } else {
    clu <- t(t(colnames(expr[,column_order(ht)[[i]]])))
    clu <- cbind(clu, paste("Cluster ", i, sep=""))
    out <- rbind(out, clu)   } 

out <- as.data.frame(out)
rownames(out) <- out$`Cell Type`

ann2_merged <- merge(ann2, out, by = 0)


write.csv2(ann2_merged, "20260408_WARD_HPA_Proteo_Celltype_Order.csv")


merged <- ann2_merged[,-1]
cell_table <- merged  %>%
  distinct(`Cell Type`, `Heatmap Cluster`, origin) %>%  # Keep unique combinations
  group_by(`Heatmap Cluster`, origin) %>%
  summarise(Unique_Cell_Count = n())

#Subset only for the three major germ layers
cell_table <- subset(cell_table, 
                     cell_table$origin %in% c("Endoderm", "Ectoderm", "Mesoderm"))

names(cell_table) <- gsub("Heatmap Cluster", "Cluster", names(cell_table))
cell_table

#Create count for endoderm in cluster 3 so to fill out the tileplot
add <- data.frame(
  Cluster = "Cluster 3",
  origin = "Endoderm", 
  Unique_Cell_Count = 0
)

#row bind the two df
cell_table <- rbind(cell_table, add)

add <- data.frame(
  Cluster = "Cluster 2",
  origin = "Endoderm", 
  Unique_Cell_Count = 0
)
cell_table <- rbind(cell_table, add)

cell_table$origin <- factor(cell_table$origin,
                            levels = c("Ectoderm", "Mesoderm", "Endoderm"))

#create the first tileplot 
library(ggpubr)
tileplot <- ggplot(cell_table,aes(x = Cluster, y = origin))+ 
  geom_tile(aes(fill = Unique_Cell_Count), color = "white") +
  scale_fill_distiller(palette = "RdYlBu") +
  #scale_fill_gradient(low = "blue", high = "red")+ 
  geom_text(aes(label = Unique_Cell_Count), color= "white", size = 12) +
  scale_x_discrete(position = "top")+
  theme_transparent()+
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))

tileplot

###############
# Fisher Test #
###############
fish <- pivot_wider(cell_table, 
                    names_from = origin, values_from = Unique_Cell_Count, values_fill = 0)

fish <- column_to_rownames(fish, "Cluster")

# Perform Chi-Square Test
fish_result <- fisher.test(fish)
print(fish_result)

p_value <- fish_result$p.value
p_value <- format(p_value, digits = 3, scientific = TRUE) 


# Add the Fisher results to the tileplot 
tileplot <- tileplot + 
  ggtitle(paste("Fisher's Exact Test p-value:", p_value)) +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

#tileplot + theme(text = element_text(size = 9))

pdf("Fisher test_Ward_HPA_Proteo.pdf",
    width = 6, height = 5)
print(tileplot)
dev.off()

### EMT TF Boxplots
library(ggpubr)
tf_data <- data.subset[emt_tf,]
tf_data <- as.data.frame(t(tf_data))
tf_data <- rownames_to_column(tf_data, var = "CellType")

merged_ann <- merge(ann2_merged, tf_data, by.x = "Cell Type", by.y = "CellType")
merged_ann <- pivot_longer(merged_ann, cols = 5:10,
                           names_to = "gene",
                           values_to = "expression")

names(merged_ann) <- gsub("Heatmap Cluster", "Cluster", names(merged_ann))

merged_ann$gene <- factor(
  merged_ann$gene,
  levels = c("SNAI1", "SNAI2",  "TWIST1", "TWIST2", "ZEB1", "ZEB2")
)

library(ggpubr)
library(stats)
library(rstatix)

pw <- merged_ann %>%
  group_by(gene) %>%
  pairwise_t_test(expression ~ Cluster, p.adjust.method = "fdr") %>%
  add_significance("p.adj") %>%     
  add_xy_position(x = "Cluster")

p <- ggplot(merged_ann, aes(x = Cluster, y = expression, fill = Cluster)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA
  ) +
  labs(
    y = "log2(nTPM+1)",
    title = "EMT TF Expression: Heatmap Clusters",
    x = ""
  ) +
  theme_pubr() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 9),
    axis.text.y = element_text(size = 8),
    panel.background = element_rect(fill = "white", colour = "white"),
    axis.line = element_line(linewidth = 0.7, colour = "black"),
    strip.text = element_text(size = 8)
  ) +
  ggsci::scale_fill_lancet() +
  facet_wrap(
    ~gene,
    nrow = 2,
    ncol = 3,
    scales = "free"
  ) +
  stat_pvalue_manual(
    pw,
    label = "p.signif",
    tip.length = 0.05,
    hide.ns = TRUE
  )

p

pdf("Ward_Proteo_GAG_EMT_Genes_Saikat", height = 8, width = 6)
print(p)
dev.off()


