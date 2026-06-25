library(Seurat)
library(tidyverse)
library(limma)
library(clustree)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(patchwork)
conflicted::conflicts_prefer(
  dplyr::filter,
  dplyr::select,
  dplyr::rename,
  dplyr::mutate,
  dplyr::arrange,
  dplyr::summarise,
  dplyr::count,
  dplyr::transmute,
  dplyr::bind_rows
)
subset_embryo_colors <- c("Ectoderm" = "#66C2A5",
                          "Endoderm" = "#1F78B4",
                          "Mesoderm" = "#F22233")

#Subsetted ones

load("Subsetted_Tyser_Processed_Seurat.RData")

proteoglycans <- readxl::read_xlsx("Proteoglycans_GAGs_List.xlsx", sheet = "Proteoglycans")
gag <- readxl::read_xlsx("Proteoglycans_GAGs_List.xlsx", sheet = "GAG")



#load(file="C:/Users/ssa214/UiT Office 365/O365-Serglycin - General/Core protein asociation towards EMT/Saikat/New/Subsetted_Processed_Germ_Layers.RData")
#load("C:/Users/ssa214/UiT Office 365/O365-Serglycin - General/Core protein asociation towards EMT/Saikat/New/Proteoglycan_Project_Colors.RData")
#load("C:/Users/ssa214/UiT Office 365/O365-Serglycin - General/Core protein asociation towards EMT/Saikat/Anne_figures/Interesting_GAG_Genes.RData")
table(sub_seurat$Subset_Cell_Types)

# Set your interesting gene set
emt_tfs <- c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2")
interesting_genes <- unique(c(proteoglycans$Gene, gag$Gene, emt_tfs))

features_list <- list(
  Proteoglycans = sort(proteoglycans$Gene),
  GAG_genes = sort(gag$Gene),
 emt_tfs = emt_tfs
)

ordered_genes <- c(sort(proteoglycans$Gene), sort(gag$Gene), emt_tfs)
proteo_genes <- proteoglycans$Gene
gags <- gag$Gene


my_comparisons <- list(
  c("Endoderm", "Mesoderm"),
  c("Endoderm", "Ectoderm"),
  c("Mesoderm", "Ectoderm")
)
library(ggpubr)

##Based on third quantile
library(dplyr)
library(tidyr)
library(rlang)

layer_levels <- c("Endoderm","Mesoderm","Ectoderm")

make_layer_quantile_seurat <- function(
    obj,
    gene_set,
    prob = 0.75,                 # 3rd quantile
    group_var = "Subset_Cell_Types",
    layers = layer_levels,
    nd_threshold = -Inf          # set to 0 to tag all-zero Q3 as "Not Detected"
) {
  gene_set <- sort(unique(gene_set))
  genes_present <- intersect(gene_set, rownames(obj))
  if (length(genes_present) == 0) stop("None of the requested genes are present in the Seurat object.")
  
  # Pull expression + group
  df <- FetchData(obj, vars = c(genes_present, group_var))
  names(df)[ncol(df)] <- group_var
  
  long_df <- df %>%
    pivot_longer(cols = all_of(genes_present), names_to = "Gene", values_to = "Expression") %>%
    filter(.data[[group_var]] %in% layers)
  
  qtab <- long_df %>%
    group_by(Gene, !!sym(group_var)) %>%
    summarise(q_expr = quantile(Expression, probs = prob, na.rm = TRUE), .groups = "drop")
  
  qtab %>%
    complete(
      Gene = gene_set,
      !!sym(group_var) := factor(layers, levels = layers),
      fill = list(q_expr = NA_real_)
    ) %>%
    pivot_wider(names_from = !!sym(group_var), values_from = q_expr) %>%
    rowwise() %>%
    mutate(
      Top_Layer = {
        vals <- as.numeric(c_across(all_of(layers)))
        if (all(is.na(vals))) {
          NA_character_
        } else if (max(vals, na.rm = TRUE) <= nd_threshold) {
          "Not Detected"
        } else {
          winners <- layers[vals == max(vals, na.rm = TRUE)]
          paste(winners, collapse = "+")
        }
      }
    ) %>%
    ungroup() %>%
    select(Gene, all_of(layers), Top_Layer) %>%
    arrange(Gene)
}

# Q3 summaries for embryonic data
proteoglycan_q3_emb <- make_layer_quantile_seurat(
  seurat_metacell, proteoglycans$Gene, prob = 0.75, layers = layer_levels
)
gag_q3_emb <- make_layer_quantile_seurat(
  seurat_metacell, gag$Gene, prob = 0.75, layers = layer_levels
)

dataframe <- list("Proteoglycans Q3" = proteoglycan_q3_emb,
                  "GAG Q3" = gag_q3_emb)
library(openxlsx)
write.xlsx(dataframe,file = "20260504_Embryonic_Tissue_ECMs_Q3.xlsx")

## Adult tissue
library(dplyr)
library(tidyr)
library(data.table)

# === Load and Prepare Data ===

all_data <- fread("C:/Users/ssa214/UiT Office 365/O365-Serglycin - General/Anne/Normal tissue expression/rna_single_cell_type_germ_layer.csv")

subset_data <- all_data %>% filter(`Gene name` %in% interesting_genes)

layer_levels <- c("Endoderm","Mesoderm","Ectoderm")


emt_tfs <- c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2")

features_list <- list(
  Proteoglycans = sort(proteoglycans$Gene),
  GAG_genes = sort(gag$Gene)
)

proteo_genes <- intersect(proteoglycans$Gene, subset_data$`Gene name`)
gags <- intersect(gag$Gene, subset_data$`Gene name`)

features_list <- list(
  Proteoglycans = sort(proteo_genes),
  Gags = sort(gags)
)
library(gridExtra)

subset_data$Gene <- subset_data$`Gene name`
subset_data$Expression <- log2(1+subset_data$nTPM)
subset_embryo_colors <- c("Ectoderm" = "#66C2A5",
                          "Endoderm" = "#1F78B4",
                          "Mesoderm" = "#F22233")

subset_data <- subset_data %>% dplyr::filter(origin %in% c("Endoderm", "Mesoderm", "Ectoderm"))
subset_data$origin <- factor(subset_data$origin, levels = c("Endoderm", "Mesoderm", "Ectoderm"))

pdf("S2H_Adult_Germ_Layers_BoxPlots_Paged.pdf", width = 10, height = 10)
for (category in names(features_list)) {
  
  genes <- features_list[[category]]
  genes <- genes[genes %in% subset_data$Gene]
  
  if (length(genes) == 0) next
  
  plot_df <- subset_data %>%
    filter(Gene %in% genes)
  
  plot_df$Gene <- factor(plot_df$Gene, levels = genes)
  
  p <- ggplot(plot_df, aes(
    x = origin,
    y = Expression,
    fill = origin
  )) +
    geom_boxplot(
      outlier.shape = 16,
      outlier.size = 0.35,
      linewidth = 0.35,
      width = 0.65,
      color = "black",
      na.rm = TRUE
    ) +
    stat_compare_means(
      comparisons = my_comparisons,
      method = "t.test",
      label = "p.signif",
      hide.ns = TRUE,
      size = 3
    ) +
    facet_wrap(
      ~ Gene,
      scales = "free_y",
      ncol = 7
    ) +
    scale_fill_manual(values = subset_embryo_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
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
    ) +
    coord_cartesian(clip = "off")
  
  print(p)
}

dev.off()


make_layer_q3 <- function(df, gene_set, prob = 0.75) {
  gene_set <- sort(unique(gene_set))
  
  df %>%
    filter(Gene %in% gene_set, origin %in% layer_levels) %>%
    group_by(Gene, origin) %>%
    summarise(q_expr = quantile(Expression, probs = prob, na.rm = TRUE), .groups = "drop") %>%
    # ensure every gene has all three layers present (even if NA)
    complete(
      Gene = gene_set,
      origin = factor(layer_levels, levels = layer_levels),
      fill = list(q_expr = NA_real_)
    ) %>%
    pivot_wider(names_from = origin, values_from = q_expr) %>%
    # pick the layer(s) with the highest Q3 (ties shown as "Endoderm+Mesoderm", etc.)
    rowwise() %>%
    mutate(
      Top_Layer = {
        vals <- c(Endoderm, Mesoderm, Ectoderm)
        if (all(is.na(vals))) NA_character_ else {
          winners <- layer_levels[which(vals == max(vals, na.rm = TRUE))]
          paste(winners, collapse = "+")
        }
      }
    ) %>%
    ungroup() %>%
    select(Gene, all_of(layer_levels), Top_Layer) %>%
    arrange(Gene)
}


# Two dataframes (Q3-based):
proteoglycan_q3 <- make_layer_q3(subset_data, proteoglycans$Gene)
gag_q3          <- make_layer_q3(subset_data, gag$Gene)

dataframe <- list("Proteoglycans Q3" = proteoglycan_q3,
                  "GAG Q3" = gag_q3)
library(openxlsx)
write.xlsx(dataframe,file = "Adult_Tissue_ECMs_Q3.xlsx")

## Merging and preparing for the Sankey Plot
# ============================================================
# Merge embryonic and adult Q3 summaries directly from objects
# ============================================================

library(dplyr)
library(openxlsx)

# Helper function to standardize column names
rename_q3_table <- function(df, suffix) {
  
  new_names <- c(
    "Gene",
    paste0("Endoderm_", suffix, "_Q3"),
    paste0("Mesoderm_", suffix, "_Q3"),
    paste0("Ectoderm_", suffix, "_Q3"),
    paste0("Top_Layer_", suffix)
  )
  
  df <- df[, c("Gene", "Endoderm", "Mesoderm", "Ectoderm", "Top_Layer")]
  names(df) <- new_names
  
  return(df)
}

# Helper function to clean Top_Layer annotation
clean_top_layer <- function(x) {
  case_when(
    is.na(x) ~ "Not Detected",
    x == "Endoderm+Mesoderm+Ectoderm" ~ "Sparse",
    TRUE ~ x
  )
}

# Helper function to merge embryonic and adult Q3 tables
merge_embryonic_adult_q3 <- function(emb_df, adult_df) {
  
  emb_clean <- emb_df %>%
    rename_q3_table("Embryonic")
  
  adult_clean <- adult_df %>%
    rename_q3_table("Adult")
  
  merged_df <- emb_clean %>%
    full_join(adult_clean, by = "Gene") %>%
    mutate(
      Top_Layer_Embryonic = clean_top_layer(Top_Layer_Embryonic),
      Top_Layer_Adult = clean_top_layer(Top_Layer_Adult),
      Overlap = ifelse(
        Top_Layer_Adult == Top_Layer_Embryonic,
        "Yes",
        "No"
      )
    ) %>%
    arrange(Gene)
  
  return(merged_df)
}

# Merge proteoglycans and GAG genes
merged_pg <- merge_embryonic_adult_q3(
  emb_df = proteoglycan_q3_emb,
  adult_df = proteoglycan_q3
)

merged_gag <- merge_embryonic_adult_q3(
  emb_df = gag_q3_emb,
  adult_df = gag_q3
)

# Check overlap summaries
table(merged_pg$Overlap)
table(merged_gag$Overlap)

table(merged_pg$Top_Layer_Embryonic)
table(merged_pg$Top_Layer_Adult)

table(merged_gag$Top_Layer_Embryonic)
table(merged_gag$Top_Layer_Adult)

# Final supplementary table
table_s8 <- list(
  "Proteoglycans Q3" = merged_pg,
  "GAG Q3" = merged_gag
)

write.xlsx(
  table_s8,
  file = "Table S8.xlsx",
  overwrite = TRUE
)

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggalluvial)

make_germlayer_sankey <- function(
    merged_df,
    output_file,
    plot_width = 8,
    plot_height = 5,
    label_genes = TRUE,
    remove_categories = c("Sparse", "Not Detected")
) {
  
  # -----------------------------
  # Prepare Sankey input
  # -----------------------------
  
  df <- merged_df %>%
    dplyr::select(Gene, Top_Layer_Embryonic, Top_Layer_Adult) %>%
    dplyr::rename(
      Embryonic = Top_Layer_Embryonic,
      Adult = Top_Layer_Adult
    )
  
  df_use <- df %>%
    dplyr::filter(
      !(Embryonic %in% remove_categories),
      !(Adult %in% remove_categories),
      !is.na(Embryonic),
      !is.na(Adult)
    )
  
  # Count Adult-Embryonic pairs
  pair_counts <- df_use %>%
    dplyr::count(Adult, Embryonic, name = "n") %>%
    dplyr::arrange(dplyr::desc(n))
  
  # Long alluvial format
  long_alluvial <- pair_counts %>%
    dplyr::mutate(
      alluvium = paste(Adult, Embryonic, sep = "/")
    ) %>%
    tidyr::pivot_longer(
      cols = c(Embryonic, Adult),
      names_to = "Stage",
      values_to = "Layer"
    ) %>%
    dplyr::select(alluvium, Stage, Layer, n)
  
  # Define ordered layer levels
  layer_levels <- c(
    "Endoderm",
    "Mesoderm",
    "Ectoderm",
    "Endoderm+Mesoderm",
    "Endoderm+Ectoderm",
    "Mesoderm+Ectoderm",
    "Endoderm+Mesoderm+Ectoderm"
  )
  
  layer_levels <- layer_levels[layer_levels %in% unique(long_alluvial$Layer)]
  
  long_alluvial <- long_alluvial %>%
    dplyr::mutate(
      Stage = factor(Stage, levels = c("Embryonic", "Adult")),
      Layer = factor(Layer, levels = layer_levels)
    )
  
  # Labels for stratum counts
  label_df <- long_alluvial %>%
    dplyr::group_by(Stage, Layer) %>%
    dplyr::summarise(
      n_total = sum(n),
      .groups = "drop"
    )
  
  stage_totals <- long_alluvial %>%
    dplyr::group_by(Stage) %>%
    dplyr::summarise(
      total = sum(n),
      .groups = "drop"
    )
  
  # Colour vector
  layer_cols <- c(
    "Ectoderm" = "#66C2A5",
    "Endoderm" = "#1F78B4",
    "Mesoderm" = "#F22233",
    "Endoderm+Ectoderm" = "purple",
    "Endoderm+Mesoderm" = "sienna1",
    "Mesoderm+Ectoderm" = "goldenrod",
    "Endoderm+Mesoderm+Ectoderm" = "skyblue",
    "Sparse" = "skyblue",
    "Not Detected" = "grey70"
  )
  
  layer_cols_use <- layer_cols[layer_levels]
  
  # -----------------------------
  # Base Sankey plot
  # -----------------------------
  
  p <- ggplot(
    long_alluvial,
    aes(
      x = Stage,
      stratum = Layer,
      alluvium = alluvium,
      y = n,
      fill = Layer
    )
  ) +
    geom_flow(
      alpha = 0.75,
      aes.flow = "forward",
      knot.pos = 0.4,
      color = "white",
      linewidth = 0.25
    ) +
    geom_stratum(
      width = 0.22,
      color = "white",
      linewidth = 0.35
    ) +
    geom_text(
      data = label_df,
      stat = "stratum",
      aes(
        x = Stage,
        stratum = Layer,
        y = n_total,
        label = paste0(Layer, "\n", n_total)
      ),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold",
      lineheight = 0.95,
      vjust = 0.5
    ) +
    geom_text(
      data = stage_totals,
      aes(
        x = Stage,
        y = total * 1.06,
        label = Stage
      ),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 3.8
    ) +
    scale_fill_manual(
      values = layer_cols_use,
      drop = FALSE
    ) +
    scale_x_discrete(expand = c(0.02, 0.02)) +
    scale_y_continuous(expand = c(0.02, 0.10)) +
    labs(fill = "Germ layer") +
    theme_void(base_size = 12) +
    theme(
      legend.position = "none",
      plot.margin = margin(6, 8, 6, 8)
    )
  
  # -----------------------------
  # Optional gene labels
  # -----------------------------
  
  if (label_genes) {
    
    collapse_lines <- function(x) {
      paste(sort(unique(x)), collapse = "\n")
    }
    
    labels_both <- dplyr::bind_rows(
      df_use %>%
        dplyr::transmute(
          Stage = "Embryonic",
          Layer = Embryonic,
          Gene = Gene
        ) %>%
        dplyr::filter(!is.na(Layer)) %>%
        dplyr::group_by(Stage, Layer) %>%
        dplyr::summarise(
          n_total = dplyr::n(),
          label = collapse_lines(Gene),
          .groups = "drop"
        ),
      
      df_use %>%
        dplyr::transmute(
          Stage = "Adult",
          Layer = Adult,
          Gene = Gene
        ) %>%
        dplyr::filter(!is.na(Layer)) %>%
        dplyr::group_by(Stage, Layer) %>%
        dplyr::summarise(
          n_total = dplyr::n(),
          label = collapse_lines(Gene),
          .groups = "drop"
        )
    ) %>%
      dplyr::mutate(
        Stage = factor(Stage, levels = c("Embryonic", "Adult")),
        Layer = factor(Layer, levels = layer_levels)
      )
    
    left_labels <- labels_both %>%
      dplyr::filter(Stage == "Embryonic")
    
    right_labels <- labels_both %>%
      dplyr::filter(Stage == "Adult")
    
    p <- p +
      geom_text(
        data = left_labels,
        stat = "stratum",
        aes(
          x = Stage,
          stratum = Layer,
          y = n_total,
          label = label
        ),
        inherit.aes = FALSE,
        hjust = 3,
        size = 3.0,
        lineheight = 0.92
      ) +
      geom_text(
        data = right_labels,
        stat = "stratum",
        aes(
          x = Stage,
          stratum = Layer,
          y = n_total,
          label = label
        ),
        inherit.aes = FALSE,
        hjust = -2,
        size = 3.0,
        lineheight = 0.92
      ) +
      coord_cartesian(clip = "off") +
      theme(
        plot.margin = margin(6, 90, 6, 90)
      )
  }
  
  # Save plot
  ggsave(
    filename = output_file,
    plot = p,
    width = plot_width,
    height = plot_height,
    units = "in",
    bg = "transparent"
  )
  
  return(p)
}

p_sankey_gag_labeled <- make_germlayer_sankey(
  merged_df = merged_gag,
  output_file = "GAG_Germlayer_specificity.svg",
  plot_width = 8,
  plot_height = 5,
  label_genes = TRUE
)

p_sankey_gag_labeled

p_sankey_pg_labeled <- make_germlayer_sankey(
  merged_df = merged_pg,
  output_file = "Proteoglycan_Germlayer_specificity.svg",
  plot_width = 10,
  plot_height = 7,
  label_genes = TRUE
)

p_sankey_pg_labeled