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
library(gridExtra)
library(tidyverse)
library(tradeSeq)
library(ggplot2)
library(ggpubr)
library(ggsci)

set.seed(123)
emt_tfs <- c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2")

proteoglycans <- readxl::read_xlsx("Proteoglycans_GAGs_List.xlsx", sheet = "Proteoglycans")
gag <- readxl::read_xlsx("Proteoglycans_GAGs_List.xlsx", sheet = "GAG")

source("Pseudotime_Helpers_All_Genes.R")

load(file= "Embryo_Separate_CDS.RData")

subset_embryo_colors <- c("Ectoderm" = "#66C2A5",
                          "Endoderm" = "#1F78B4",
                          "Mesoderm" = "#F22233")

my_comparisons <- list(
  c("Endoderm", "Mesoderm"),
  c("Endoderm", "Ectoderm"),
  c("Mesoderm", "Ectoderm")
)

set.seed(123)

proteo_genes <- intersect(proteoglycans$Gene, rownames(sce_endo))
gags <- intersect(gag$Gene, rownames(sce_endo))

features_list <- list(
  EMT_TFs = c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2"),
  Proteoglycans = sort(proteo_genes),
  Gags = sort(gags)
)

##Endoderm
pg_pool  <- intersect(proteoglycans$Gene, rownames(counts_endo))
emt_pool <- intersect(emt_tfs,              rownames(counts_endo))
gag_pool <- intersect(gag$Gene,            rownames(counts_endo))

endo_pools <- list(
  `EMT TFs` = emt_pool,
  Proteoglycans = pg_pool,
  GAGs = gag_pool
)
res_endo <- run_layer_plots(
  sce = sce_endo, asso_df = asso_endo, counts_mat = counts_endo,
  pools = endo_pools, layer_name = "Endoderm",
  alpha = 0.05, out_dir = "pseudotime_plots/Endoderm",
  nPoints = 200, zscore = TRUE, genes_per_page = 16
)


##Mesoderm
pg_pool  <- intersect(proteoglycans$Name, rownames(counts_meso))
emt_pool <- intersect(emt_tfs,              rownames(counts_meso))
gag_pool <- intersect(gag$genes,            rownames(counts_meso))
meso_pools <- list(
  `EMT TFs` = emt_pool,
  Proteoglycans = pg_pool,
  GAGs = gag_pool
)
res_meso <- run_layer_plots(
  sce = sce_meso, asso_df = asso_meso, counts_mat = counts_meso,
  pools = meso_pools, layer_name = "Mesoderm",
  alpha = 0.05, out_dir = "pseudotime_plots/Mesoderm",
  nPoints = 200, zscore = TRUE, genes_per_page = 16
)

gene_sets <- list(
  emt = emt_tfs,
  pg  = proteoglycans$Gene,
  gags = gag$Gene  
)

gags = gag$Gene

grad_endo <- c("#1F78B4", "#32CD32", "#66C2A5")
library(ggpubr)

# Following will create the pairwise module correlation plots as well. But this plot may not be easy to understand
res_endo_pg  <- run_endo_proteoglycans()
res_endo_gag  <- run_endo_gags()
res_meso_pg  <- run_meso_proteoglycans()
res_meso_gag  <- run_meso_gags()

cor_dataframe <- list("Endo_Proteo" = res_endo_pg[["cor_pairs_wide"]],
                      "Endo_GAG" = res_endo_gag[["cor_pairs_wide"]],
                      "Meso_Proteo" = res_meso_pg[["cor_pairs_wide"]],
                      "Meso_GAG" = res_meso_gag[["cor_pairs_wide"]])

writexl::write_xlsx(cor_dataframe, "Only Significant_Correlation_with_Pseudotime_Dynamic.xlsx")


endo_corr_excel <- make_layer_emt_correlation_excel(
  sce = sce_endo,
  asso_df = asso_endo,
  counts_mat = counts_endo,
  proteo_genes = proteoglycans$Gene,
  gag_genes = gag$Gene,
  emt_genes = emt_tfs,
  layer_name = "Endoderm",
  layer_short = "Endo",
  out_dir = "pseudotime_plots/Endoderm",
  nPoints = 200,
  cor_method = "pearson",
  zscore = TRUE,
  cor_threshold = 0.4
)

meso_corr_excel <- make_layer_emt_correlation_excel(
  sce = sce_meso,
  asso_df = asso_meso,
  counts_mat = counts_meso,
  proteo_genes = proteoglycans$Gene,
  gag_genes = gag$Gene,
  emt_genes = emt_tfs,
  layer_name = "Mesoderm",
  layer_short = "Meso",
  out_dir = "pseudotime_plots/Mesoderm",
  nPoints = 200,
  cor_method = "pearson",
  zscore = TRUE,
  cor_threshold = 0.4
)


## Getting a summary of the correlations
summarise_tf_cor_counts <- function(
    df,
    tf_cols = c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2"),
    cor_threshold = 0.4
) {
  tf_cols <- intersect(tf_cols, colnames(df))
  
  df %>%
    dplyr::select(dplyr::all_of(tf_cols)) %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "TF",
      values_to = "correlation"
    ) %>%
    dplyr::group_by(.data$TF) %>%
    dplyr::summarise(
      n_positive = sum(.data$correlation >= cor_threshold, na.rm = TRUE),
      n_negative = sum(.data$correlation <= -cor_threshold, na.rm = TRUE),
      n_total_tested = sum(!is.na(.data$correlation)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      n_above_threshold = .data$n_positive + .data$n_negative,
      positive_fraction = .data$n_positive / .data$n_total_tested,
      negative_fraction = .data$n_negative / .data$n_total_tested
    )
}

endo_proteo_tf_counts <- summarise_tf_cor_counts(
  df = endo_corr_excel$Proteo,
  cor_threshold = 0.4
)

endo_gag_tf_counts <- summarise_tf_cor_counts(
  df = endo_corr_excel$GAG,
  cor_threshold = 0.4
)

meso_proteo_tf_counts <- summarise_tf_cor_counts(
  df = meso_corr_excel$Proteo,
  cor_threshold = 0.4
)

meso_gag_tf_counts <- summarise_tf_cor_counts(
  df = meso_corr_excel$GAG,
  cor_threshold = 0.4
)


## Visualization
# Endoderm ------------------------------------------------------------
endo_pools <- list(
  `EMT TFs` = intersect(emt_tfs, rownames(counts_endo)),
  Proteoglycans = intersect(proteoglycans$Gene, rownames(counts_endo)),
  GAGs = intersect(gag$Gene, rownames(counts_endo))
)

three_endo <- make_three_module_trajectory_final(
  sce = sce_endo,
  asso_df = asso_endo,
  counts_mat = counts_endo,
  pools = endo_pools,
  layer_name = "Endoderm",
  alpha = 0.05,
  nPoints = 200,
  zscore_gene_for_plot = FALSE,
  use_only_sig = TRUE,
  cor_method = "pearson",
  out_dir = "pseudotime_plots/Endoderm",
  show_ribbon = FALSE
)

# Mesoderm ------------------------------------------------------------
meso_pools <- list(
  `EMT TFs` = intersect(emt_tfs, rownames(counts_meso)),
  Proteoglycans = intersect(proteoglycans$Gene, rownames(counts_meso)),
  GAGs = intersect(gag$Gene, rownames(counts_meso))
)

three_meso <- make_three_module_trajectory_final(
  sce = sce_meso,
  asso_df = asso_meso,
  counts_mat = counts_meso,
  pools = meso_pools,
  layer_name = "Mesoderm",
  alpha = 0.05,
  nPoints = 200,
  zscore_gene_for_plot = FALSE,
  use_only_sig = TRUE,
  cor_method = "pearson",
  out_dir = "pseudotime_plots/Mesoderm",
  show_ribbon = FALSE
)

### Facet plot
facet_df <- bind_rows(
  three_endo$module_df_raw %>%
    mutate(Germ_layer = "Endoderm"),
  three_meso$module_df_raw %>%
    mutate(Germ_layer = "Mesoderm")
) %>%
  mutate(
    Germ_layer = factor(Germ_layer, levels = c("Endoderm", "Mesoderm")),
    Module = factor(Module, levels = c("EMT TFs", "Proteoglycans", "GAGs"))
  )

# ------------------------------------------------------------
# Combine only EMT TFs vs PG and EMT TFs vs GAG R labels
# ------------------------------------------------------------

r_facet_df <- bind_rows(
  three_endo$r_inside_df %>%
    mutate(Germ_layer = "Endoderm"),
  three_meso$r_inside_df %>%
    mutate(Germ_layer = "Mesoderm")
) %>%
  mutate(
    Germ_layer = factor(Germ_layer, levels = c("Endoderm", "Mesoderm"))
  ) %>%
  group_by(Germ_layer) %>%
  summarise(
    r_label = paste(label, collapse = "\n"),
    .groups = "drop"
  )

label_pos_df <- facet_df %>%
  group_by(Germ_layer) %>%
  summarise(
    x = min(pt, na.rm = TRUE) + 0.03 * diff(range(pt, na.rm = TRUE)),
    y = max(mean_log, na.rm = TRUE) - 0.05 * diff(range(mean_log, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(r_facet_df, by = "Germ_layer")

# ------------------------------------------------------------
# Colors
# ------------------------------------------------------------

module_cols <- c(
  "EMT TFs" = "#FFD92F",
  "Proteoglycans" = "#00468BFF",
  "GAGs" = "#ED0000FF"
)

# ------------------------------------------------------------
# Faceted plot
# ------------------------------------------------------------

p_facet_log <- ggplot(
  facet_df,
  aes(x = pt, y = mean_log, color = Module)
) +
  geom_line(linewidth = 1.25, na.rm = TRUE) +
  geom_text(
    data = label_pos_df,
    aes(x = x, y = y, label = r_label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.3,
    color = "black"
  ) +
  facet_wrap(Germ_layer ~ ., scales = "free_y", ncol = 1) +
  scale_color_manual(values = module_cols, drop = FALSE) +
  labs(
    title = "Module dynamics along pseudotime",
    subtitle = "Exact log-transformed module trajectories used for Pearson correlation",
    x = "Scaled pseudotime",
    y = "log2(1 + mean smoothed module expression)",
    color = NULL
  ) +
  theme_pubr(base_size = 11) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 12
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 9
    ),
    axis.title = element_text(
      face = "bold",
      size = 10
    ),
    axis.text = element_text(size = 9),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    strip.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.3
    ),
    legend.position = "top",
    legend.text = element_text(size = 9),
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.4
    )
  )

p_facet_log
ggsave(
  filename = "pseudotime_plots/2C.pdf",
  plot = p_facet_log,
  width = 5.5,
  height = 6.7
)
ggsave(
  filename = "pseudotime_plots/2C.svg",
  plot = p_facet_log,
  width = 6.8,
  height = 7.2
)

