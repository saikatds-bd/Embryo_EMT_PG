library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(stringr)
library(grid)
library(writexl)

# ============================================================
# 1. Paths and global settings
# ============================================================

sup_dir <- "Supplementary/"

emt_tfs <- c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2")

r_cutoff <- 0.4

# Output folders
emb_bubble_out <- file.path(sup_dir, "Embryo_TF_Correlation_Bubble_Plots")
adult_bubble_out <- file.path(sup_dir, "Adult_TF_Correlation_Bubble_Plots")

dir.create(emb_bubble_out, recursive = TRUE, showWarnings = FALSE)
dir.create(adult_bubble_out, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# General helper: convert correlation table to long format
# ============================================================

make_cor_long <- function(df,
                          layer_name,
                          gene_category,
                          gene_col = NULL,
                          tf_order = emt_tfs) {
  
  if (is.null(gene_col)) {
    gene_col <- colnames(df)[1]
  }
  
  df %>%
    dplyr::rename(Module_gene = all_of(gene_col)) %>%
    tidyr::pivot_longer(
      cols = all_of(tf_order),
      names_to = "EMT_gene",
      values_to = "r"
    ) %>%
    dplyr::mutate(
      Layer = layer_name,
      Gene_category = gene_category,
      abs_r = abs(r),
      Direction = dplyr::case_when(
        r > 0 ~ "Positive",
        r < 0 ~ "Negative",
        TRUE ~ "Zero"
      )
    )
}


# ============================================================
# Functions
# ============================================================

plot_tf_cor_bubbles <- function(df_plot,
                                title,
                                gene_order = NULL,
                                tf_order = emt_tfs,
                                layer_order = NULL,
                                r_cutoff = 0.4,
                                base_size = 8,
                                dot_range = c(1.5, 5.5),
                                show_facet = TRUE) {
  
  if (nrow(df_plot) == 0) {
    stop("No correlations passed the cutoff for: ", title)
  }
  
  if (is.null(gene_order)) {
    gene_order <- unique(df_plot$Module_gene)
  }
  
  if (is.null(layer_order)) {
    layer_order <- unique(df_plot$Layer)
  }
  
  df_plot <- df_plot %>%
    mutate(
      EMT_gene = factor(EMT_gene, levels = tf_order),
      Module_gene = factor(Module_gene, levels = rev(gene_order)),
      Layer = factor(Layer, levels = layer_order),
      Direction = factor(Direction, levels = c("Positive", "Negative"))
    )
  
  p <- ggplot(df_plot, aes(x = EMT_gene, y = Module_gene)) +
    geom_point(
      aes(size = abs_r, fill = r),
      shape = 21,
      color = "white",
      stroke = 0.15,
      alpha = 0.95
    ) +
    scale_fill_gradient2(
      low = "#4575b4",
      mid = "white",
      high = "#d73027",
      midpoint = 0,
      limits = c(-1, 1),
      name = "r"
    ) +
    scale_size_continuous(
      range = dot_range,
      limits = c(r_cutoff, 1),
      breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
      name = "|r|",
      guide = guide_legend(
        override.aes = list(
          shape = 21,
          fill = "black",
          color = "black",
          alpha = 1,
          stroke = 0.2
        )
      )
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = title
    ) +
    theme_pubr(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = base_size),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.key.size = grid::unit(0.45, "cm"),
      legend.spacing.y = grid::unit(0.12, "cm")
    )
  
  if (show_facet) {
    p <- p +
      facet_grid(
        . ~ Layer,
        scales = "free_x",
        space = "free_x"
      ) +
      theme(
        strip.background = element_rect(
          fill = "white",
          color = "black",
          linewidth = 0.5
        ),
        strip.text = element_text(face = "bold", size = base_size)
      )
  }
  
  return(p)
}


save_cor_bubble_plot <- function(plot,
                                 df_plot,
                                 filename,
                                 out_dir,
                                 width = 5.2,
                                 min_height = 2.5,
                                 height_per_gene = 0.22,
                                 bg = "transparent") {
  
  n_genes <- length(unique(df_plot$Module_gene))
  plot_height <- max(min_height, n_genes * height_per_gene)
  
  ggsave(
    filename = file.path(out_dir, filename),
    plot = plot,
    width = width,
    height = plot_height,
    dpi = 300,
    bg = bg
  )
  
  message("Saved: ", file.path(out_dir, filename))
  message("Number of plotted genes: ", n_genes)
  message("Height used: ", round(plot_height, 2), " inches")
}


# ============================================================
# Embryonic data: read tables
# ============================================================

emb_endo_pg <- readxl::read_xlsx(
  paste0(sup_dir, "Table S10.xlsx"),
  sheet = "Endo_Proteo"
)
emb_endo_pg <- emb_endo_pg[1:37, 1:7]

emb_meso_pg <- readxl::read_xlsx(
  paste0(sup_dir, "Table S10.xlsx"),
  sheet = "Meso_Proteo"
)
emb_meso_pg <- emb_meso_pg[1:37, 1:7]

emb_endo_gag <- readxl::read_xlsx(
  paste0(sup_dir, "Table S10.xlsx"),
  sheet = "Endo_GAG"
)
emb_endo_gag <- emb_endo_gag[1:46, 1:7]

emb_meso_gag <- readxl::read_xlsx(
  paste0(sup_dir, "Table S10.xlsx"),
  sheet = "Meso_GAG"
)
emb_meso_gag <- emb_meso_gag[1:46, 1:7]


# ============================================================
# Embryonic data: long format
# ============================================================

emb_pg_long <- bind_rows(
  make_cor_long(
    df = emb_endo_pg,
    layer_name = "Endoderm",
    gene_category = "Proteoglycan"
  ),
  make_cor_long(
    df = emb_meso_pg,
    layer_name = "Mesoderm",
    gene_category = "Proteoglycan"
  )
)

emb_gag_long <- bind_rows(
  make_cor_long(
    df = emb_endo_gag,
    layer_name = "Endoderm",
    gene_category = "GAG"
  ),
  make_cor_long(
    df = emb_meso_gag,
    layer_name = "Mesoderm",
    gene_category = "GAG"
  )
)


# ============================================================
#  Embryonic data: filter abs(r) >= 0.4
# ============================================================

emb_pg_plot_df <- emb_pg_long %>%
  filter(!is.na(r), abs_r >= r_cutoff)

emb_gag_plot_df <- emb_gag_long %>%
  filter(!is.na(r), abs_r >= r_cutoff)


# ============================================================
#  Embryonic data: gene order
# ============================================================

emb_pg_gene_order <- emb_pg_long %>%
  pull(Module_gene) %>%
  unique()

emb_gag_gene_order <- emb_gag_long %>%
  pull(Module_gene) %>%
  unique()


p_emb_pg <- plot_tf_cor_bubbles(
  df_plot = emb_pg_plot_df,
  title = "Proteoglycan genes correlated with EMT transcription factors",
  gene_order = emb_pg_gene_order,
  layer_order = c("Endoderm", "Mesoderm"),
  r_cutoff = r_cutoff,
  dot_range = c(1.5, 5.5),
  show_facet = TRUE
)

p_emb_gag <- plot_tf_cor_bubbles(
  df_plot = emb_gag_plot_df,
  title = "GAG-related genes correlated with EMT transcription factors",
  gene_order = emb_gag_gene_order,
  layer_order = c("Endoderm", "Mesoderm"),
  r_cutoff = r_cutoff,
  dot_range = c(1.5, 5.5),
  show_facet = TRUE
)

p_emb_pg
p_emb_gag


save_cor_bubble_plot(
  plot = p_emb_pg,
  df_plot = emb_pg_plot_df,
  filename = "20260528_Embryo_Proteoglycan_TF_Correlation_BubblePlot.pdf",
  out_dir = emb_bubble_out,
  width = 5.2,
  min_height = 2.5,
  height_per_gene = 0.22,
  bg = "transparent"
)

save_cor_bubble_plot(
  plot = p_emb_gag,
  df_plot = emb_gag_plot_df,
  filename = "20260528_Embryo_GAG_TF_Correlation_BubblePlot.pdf",
  out_dir = emb_bubble_out,
  width = 5.2,
  min_height = 2.5,
  height_per_gene = 0.22,
  bg = "transparent"
)


# ============================================================
# Adult data: read tables
# ============================================================

adult_pg <- readxl::read_xlsx(
  paste0(sup_dir, "Paper3_Table S11.xlsx"),
  sheet = "Proteoglycans"
)
adult_pg <- adult_pg[1:42, 1:7]

adult_gag <- readxl::read_xlsx(
  paste0(sup_dir, "Paper3_Table S11.xlsx"),
  sheet = "GAG_enzymes"
)
adult_gag <- adult_gag[1:48, 1:7]


adult_pg_long <- make_cor_long(
  df = adult_pg,
  layer_name = "Adult",
  gene_category = "Proteoglycan",
  gene_col = "Gene"
)

adult_gag_long <- make_cor_long(
  df = adult_gag,
  layer_name = "Adult",
  gene_category = "GAG",
  gene_col = "Gene"
)


adult_pg_plot_df <- adult_pg_long %>%
  filter(!is.na(r), abs_r >= r_cutoff)

adult_gag_plot_df <- adult_gag_long %>%
  filter(!is.na(r), abs_r >= r_cutoff)

adult_pg_gene_order <- adult_pg_long %>%
  pull(Module_gene) %>%
  unique()

adult_gag_gene_order <- adult_gag_long %>%
  pull(Module_gene) %>%
  unique()

p_adult_pg <- plot_tf_cor_bubbles(
  df_plot = adult_pg_plot_df,
  title = "Adult proteoglycan genes correlated with EMT transcription factors",
  gene_order = adult_pg_gene_order,
  layer_order = "Adult",
  r_cutoff = r_cutoff,
  dot_range = c(1.5, 5.5),
  show_facet = FALSE
)

p_adult_gag <- plot_tf_cor_bubbles(
  df_plot = adult_gag_plot_df,
  title = "Adult GAG-related genes correlated with EMT transcription factors",
  gene_order = adult_gag_gene_order,
  layer_order = "Adult",
  r_cutoff = r_cutoff,
  dot_range = c(1.5, 5.5),
  show_facet = FALSE
)

p_adult_pg
p_adult_gag

save_cor_bubble_plot(
  plot = p_adult_pg,
  df_plot = adult_pg_plot_df,
  filename = "20260528_Adult_Proteoglycan_TF_Correlation_BubblePlot.pdf",
  out_dir = adult_bubble_out,
  width = 5.2,
  min_height = 2.5,
  height_per_gene = 0.22,
  bg = "transparent"
)

save_cor_bubble_plot(
  plot = p_adult_gag,
  df_plot = adult_gag_plot_df,
  filename = "20260528_Adult_GAG_TF_Correlation_BubblePlot.pdf",
  out_dir = adult_bubble_out,
  width = 5.2,
  min_height = 2.5,
  height_per_gene = 0.22,
  bg = "transparent"
)

# ============================================================
# Summary plots: adult + embryo consistency
#   Keep genes with same-direction correlation between Adult
#   and at least one germ layer for at least one EMT TF
# ============================================================

make_consistent_adult_embryo_genes <- function(adult_long,
                                               embryo_long,
                                               r_cutoff = 0.4) {
  
  adult_hits <- adult_long %>%
    filter(!is.na(r), abs_r >= r_cutoff) %>%
    mutate(
      Adult_r = r,
      Adult_direction = case_when(
        r > 0 ~ "Positive",
        r < 0 ~ "Negative",
        TRUE ~ "Zero"
      )
    ) %>%
    select(
      Module_gene,
      EMT_gene,
      Gene_category,
      Adult_r,
      Adult_direction
    )
  
  embryo_hits <- embryo_long %>%
    filter(!is.na(r), abs_r >= r_cutoff) %>%
    mutate(
      Embryo_r = r,
      Embryo_direction = case_when(
        r > 0 ~ "Positive",
        r < 0 ~ "Negative",
        TRUE ~ "Zero"
      )
    ) %>%
    select(
      Module_gene,
      EMT_gene,
      Gene_category,
      Layer,
      Embryo_r,
      Embryo_direction
    )
  
  consistent_pairs <- adult_hits %>%
    inner_join(
      embryo_hits,
      by = c("Module_gene", "EMT_gene", "Gene_category")
    ) %>%
    filter(
      Adult_direction == Embryo_direction
    ) %>%
    mutate(
      Consistent_with = Layer,
      Consistency_summary = paste0(
        EMT_gene, " (", Adult_direction, " in Adult and ", Layer, ")"
      )
    )
  
  consistent_genes <- consistent_pairs %>%
    distinct(Module_gene) %>%
    pull(Module_gene)
  
  return(
    list(
      genes = consistent_genes,
      pairs = consistent_pairs
    )
  )
}

# ============================================================
# Identify consistent proteoglycan and GAG genes
# ============================================================

pg_consistency <- make_consistent_adult_embryo_genes(
  adult_long = adult_pg_long,
  embryo_long = emb_pg_long,
  r_cutoff = r_cutoff
)

gag_consistency <- make_consistent_adult_embryo_genes(
  adult_long = adult_gag_long,
  embryo_long = emb_gag_long,
  r_cutoff = r_cutoff
)

pg_consistent_genes <- pg_consistency$genes
gag_consistent_genes <- gag_consistency$genes

pg_consistent_pairs <- pg_consistency$pairs
gag_consistent_pairs <- gag_consistency$pairs

# ============================================================
# Build summary plot data
# ============================================================

summary_pg_long <- bind_rows(
  adult_pg_long,
  emb_pg_long
) %>%
  filter(
    Module_gene %in% pg_consistent_genes,
    !is.na(r),
    abs_r >= r_cutoff
  ) %>%
  mutate(
    Layer = factor(
      Layer,
      levels = c("Adult", "Endoderm", "Mesoderm")
    )
  )

summary_gag_long <- bind_rows(
  adult_gag_long,
  emb_gag_long
) %>%
  filter(
    Module_gene %in% gag_consistent_genes,
    !is.na(r),
    abs_r >= r_cutoff
  ) %>%
  mutate(
    Layer = factor(
      Layer,
      levels = c("Adult", "Endoderm", "Mesoderm")
    )
  )


summary_pg_gene_order <- bind_rows(adult_pg_long, emb_pg_long) %>%
  filter(Module_gene %in% pg_consistent_genes) %>%
  pull(Module_gene) %>%
  unique()

summary_gag_gene_order <- bind_rows(adult_gag_long, emb_gag_long) %>%
  filter(Module_gene %in% gag_consistent_genes) %>%
  pull(Module_gene) %>%
  unique()


p_summary_pg <- plot_tf_cor_bubbles(
  df_plot = summary_pg_long,
  title = "Proteoglycan genes with conserved adult–embryonic EMT-TF correlations",
  gene_order = summary_pg_gene_order,
  layer_order = c("Adult", "Endoderm", "Mesoderm"),
  r_cutoff = r_cutoff,
  dot_range = c(1.5, 5.5),
  show_facet = TRUE
)

p_summary_gag <- plot_tf_cor_bubbles(
  df_plot = summary_gag_long,
  title = "GAG-related genes with conserved adult–embryonic EMT-TF correlations",
  gene_order = summary_gag_gene_order,
  layer_order = c("Adult", "Endoderm", "Mesoderm"),
  r_cutoff = r_cutoff,
  dot_range = c(1.5, 5.5),
  show_facet = TRUE
)

p_summary_pg
p_summary_gag

# ============================================================
# Save Figure 3E and 3F summary plots
# ============================================================

summary_bubble_out <- file.path(sup_dir, "Summary_TF_Correlation_Bubble_Plots")
dir.create(summary_bubble_out, recursive = TRUE, showWarnings = FALSE)

save_cor_bubble_plot(
  plot = p_summary_pg,
  df_plot = summary_pg_long,
  filename = "Summary_Proteoglycan_Adult_Embryo_Consistent_TF_Correlation_BubblePlot.pdf",
  out_dir = summary_bubble_out,
  width = 5.5,
  min_height = 2.5,
  height_per_gene = 0.22,
  bg = "transparent"
)

save_cor_bubble_plot(
  plot = p_summary_gag,
  df_plot = summary_gag_long,
  filename = "Summary_GAG_Adult_Embryo_Consistent_TF_Correlation_BubblePlot.pdf",
  out_dir = summary_bubble_out,
  width = 5.5,
  min_height = 2.5,
  height_per_gene = 0.22,
  bg = "transparent"
)