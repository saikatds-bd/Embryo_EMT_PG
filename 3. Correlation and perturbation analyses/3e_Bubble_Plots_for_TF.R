library(readxl)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(writexl)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("reduce", "purrr")



excel_dir <- "TF_RNA_Excel"

gene_file <- "PG_GAG_Genes.RData"

out_dir <- file.path(excel_dir, "Bubble_Plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. Load gene categories
# ============================================================

load(gene_file)

emt_tfs <- c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2")

proteo_genes <- unique(proteoglycans$Name)
gags <- unique(gag$genes)

interesting_genes <- unique(c(
  proteo_genes,
  gags,
  emt_tfs
))

features_list <- list(
  Proteoglycans = sort(proteo_genes),
  GAG_genes = sort(gags)
)


excel_files_all <- list.files(
  path = excel_dir,
  pattern = "\\.xlsx$",
  full.names = TRUE
)

excel_files_all

get_tf_name_from_file <- function(file) {
  fname <- basename(file)
  
  tf <- str_extract(fname, "SNAI1|SNAI2|TWIST1|TWIST2|ZEB1|ZEB2")
  
  if (is.na(tf)) {
    warning("Could not detect TF name from file: ", fname)
  }
  
  return(tf)
}

tf_file_table <- tibble(
  file = excel_files_all,
  TF = map_chr(excel_files_all, get_tf_name_from_file)
) %>%
  filter(!is.na(TF)) %>%
  distinct(TF, .keep_all = TRUE) %>%
  arrange(match(TF, emt_tfs))

tf_file_table


read_one_tf_excel <- function(file, tf_name, interesting_genes) {
  
  sheet_name <- paste0(tf_name, "_DE_Converted")
  
  available_sheets <- readxl::excel_sheets(file)
  
  if (!sheet_name %in% available_sheets) {
    stop(
      "Sheet not found: ", sheet_name,
      "\nFile: ", file,
      "\nAvailable sheets: ", paste(available_sheets, collapse = ", ")
    )
  }
  
  df <- readxl::read_excel(file, sheet = sheet_name) %>%
    filter(Symbol %in% interesting_genes) %>%
    pivot_longer(
      cols = -Symbol,
      names_to = "CellLine",
      values_to = "logFC"
    ) %>%
    mutate(
      TF = tf_name,
      Perturbation = case_when(
        str_detect(CellLine, "\\bKD\\b|\\bKO\\b") ~ "KD/KO",
        str_detect(CellLine, "\\bOE\\b") ~ "OE",
        TRUE ~ "Unknown"
      )
    )
  
  return(df)
}

long_logfc <- purrr::pmap_dfr(
  list(
    file = tf_file_table$file,
    tf_name = tf_file_table$TF
  ),
  ~ read_one_tf_excel(
    file = ..1,
    tf_name = ..2,
    interesting_genes = interesting_genes
  )
)

long_logfc <- long_logfc %>%
  mutate(
    TF = factor(TF, levels = emt_tfs),
    Gene_Category = case_when(
      Symbol %in% proteo_genes ~ "Proteoglycan",
      Symbol %in% gags ~ "GAG",
      Symbol %in% emt_tfs ~ "EMT_TF",
      TRUE ~ "Other"
    )
  )

head(long_logfc)



gene_tf_qtiles <- long_logfc %>%
  filter(!is.na(logFC)) %>%
  group_by(Symbol, TF) %>%
  summarise(
    n = sum(!is.na(logFC)),
    Q1 = quantile(logFC, 0.25, na.rm = TRUE),
    Q3 = quantile(logFC, 0.75, na.rm = TRUE),
    Median = median(logFC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n >= 3)

gene_tf_qtiles


make_bubble_pairs <- function(gene_list,
                              qtiles_df = gene_tf_qtiles,
                              up_thr = 1,
                              down_thr = -1) {
  
  df <- qtiles_df %>%
    filter(Symbol %in% gene_list) %>%
    mutate(
      TF = as.character(TF),
      Regulation = case_when(
        Q3 > up_thr ~ "Upregulated by EMT TF",
        Q1 < down_thr ~ "Downregulated by EMT TF",
        TRUE ~ NA_character_
      ),
      Effect = case_when(
        Q3 > up_thr ~ Q3,
        Q1 < down_thr ~ abs(Q1),
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(Regulation)) %>%
    transmute(
      Module_gene = Symbol,
      EMT_gene = TF,
      Regulation,
      Effect,
      n,
      Q1,
      Q3,
      Median
    )
  
  return(df)
}

proteo_pairs <- make_bubble_pairs(proteo_genes)
gag_pairs <- make_bubble_pairs(gags)

proteo_pairs
gag_pairs

plot_tf_bubbles <- function(df_pairs,
                            title,
                            tf_order = emt_tfs,
                            gene_order = NULL,
                            effect_range = c(1.8, 7.5)) {
  
  if (nrow(df_pairs) == 0) {
    stop("No regulated gene-TF pairs found for: ", title)
  }
  
  if (is.null(gene_order)) {
    gene_order <- sort(unique(df_pairs$Module_gene))
  }
  
  df_plot <- df_pairs %>%
    mutate(
      EMT_gene = factor(EMT_gene, levels = tf_order),
      Module_gene = factor(Module_gene, levels = rev(gene_order)),
      Regulation = factor(
        Regulation,
        levels = c(
          "Upregulated by EMT TF",
          "Downregulated by EMT TF"
        )
      )
    )
  
  p <- ggplot(df_plot, aes(x = EMT_gene, y = Module_gene)) +
    geom_point(
      aes(size = Effect, fill = Regulation),
      shape = 21,
      color = "black",
      stroke = 0.25,
      alpha = 0.9
    ) +
    scale_fill_manual(
      values = c(
        "Upregulated by EMT TF" = "#d73027",
        "Downregulated by EMT TF" = "#4575b4"
      ),
      name = NULL
    ) +
    scale_size(
      range = effect_range,
      name = expression(paste("|Log"[2], "FC| effect"))
    ) +
    labs(
      x = "EMT transcription factor",
      y = NULL,
      title = title
    ) +
    theme_pubr(base_size = 11) +
    theme(
      #panel.grid.major = element_line(linewidth = 0.25, color = "grey85"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 8),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
  
  return(p)
}

# ============================================================
# 8. Make plots
# ============================================================

proteo_gene_order <- sort(unique(proteo_pairs$Module_gene))
gag_gene_order <- sort(unique(gag_pairs$Module_gene))

p_proteo <- plot_tf_bubbles(
  df_pairs = proteo_pairs,
  title = "Proteoglycans regulated by EMT transcription factors",
  gene_order = proteo_gene_order
)

p_gag <- plot_tf_bubbles(
  df_pairs = gag_pairs,
  title = "GAG-related genes regulated by EMT transcription factors",
  gene_order = gag_gene_order
)

p_proteo
p_gag

#Extract
# ============================================================
# 9. Save plots
# ============================================================

save_bubble_plot <- function(plot, df_pairs, filename, width = 7) {
  
  n_genes <- length(unique(df_pairs$Module_gene))
  plot_height <- max(2.5, n_genes * 0.22)
  
  ggsave(
    filename = file.path(out_dir, filename),
    plot = plot,
    width = width,
    height = plot_height,
    dpi = 300,
    bg = "transparent"
  )
}

save_bubble_plot(
  plot = p_proteo,
  df_pairs = proteo_pairs,
  filename = "20260528_Proteoglycans_TF_BubblePlot.pdf",
  width = 6
)

save_bubble_plot(
  plot = p_gag,
  df_pairs = gag_pairs,
  filename = "20260528_GAG_genes_TF_BubblePlot.pdf",
  width = 6
)


## Excel

table_s12_out <- file.path(
  out_dir,
  "20260528_Paper3_Table_S12_TF_RNA_Bubble_Source.xlsx"
)

logfc_detailed_sheet <- long_logfc %>%
  mutate(
    TF = as.character(TF),
    Category = Gene_Category,
    Column_Name = paste(TF, CellLine, sep = "_")
  ) %>%
  select(
    Symbol,
    Category,
    Column_Name,
    logFC
  ) %>% dplyr::filter(!(Symbol %in% emt_tfs)) %>%
  distinct() %>%
  pivot_wider(
    names_from = Column_Name,
    values_from = logFC
  ) %>%
  mutate(
    across(
      where(is.numeric),
      ~ tidyr::replace_na(.x, 0)
    )
  ) %>%
  arrange(
    factor(
      Category,
      levels = c("Proteoglycan", "GAG", "EMT_TF", "Other")
    ),
    Symbol
  )

proteoglycan_summary_sheet <- proteo_pairs %>%
  arrange(EMT_gene, Regulation, Module_gene) %>%
  select(
    Module_gene,
    EMT_gene,
    Regulation,
    Effect,
    n,
    Q1,
    Q3,
    Median
  )

gag_summary_sheet <- gag_pairs %>%
  arrange(EMT_gene, Regulation, Module_gene) %>%
  select(
    Module_gene,
    EMT_gene,
    Regulation,
    Effect,
    n,
    Q1,
    Q3,
    Median
  )


make_tf_summary_wide <- function(gene_list,
                                 qtiles_df = gene_tf_qtiles,
                                 tf_order = emt_tfs,
                                 up_thr = 1,
                                 down_thr = -1) {
  
  qtiles_sub <- qtiles_df %>%
    filter(Symbol %in% gene_list) %>%
    mutate(
      TF = as.character(TF),
      Up_call = Q3 > up_thr,
      Down_call = Q1 < down_thr
    )
  
  q3_wide <- qtiles_sub %>%
    select(Symbol, TF, Q3) %>%
    mutate(TF_col = paste0(TF, "_Q3_LogFC")) %>%
    select(Symbol, TF_col, Q3) %>%
    pivot_wider(
      names_from = TF_col,
      values_from = Q3
    )
  
  call_summary <- qtiles_sub %>%
    group_by(Symbol) %>%
    summarise(
      Upregulated_by = paste(TF[Up_call], collapse = ", "),
      Downregulated_by = paste(TF[Down_call], collapse = ", "),
      .groups = "drop"
    ) %>%
    mutate(
      Upregulated_by = na_if(Upregulated_by, ""),
      Downregulated_by = na_if(Downregulated_by, "")
    )
  
  regulation_summary <- qtiles_sub %>%
    filter(Up_call | Down_call) %>%
    mutate(
      call_text = case_when(
        Up_call ~ paste0(TF, " (Upregulated)"),
        Down_call ~ paste0(TF, " (Downregulated)")
      )
    ) %>%
    group_by(Symbol) %>%
    summarise(
      Regulation_summary = paste(call_text, collapse = ", "),
      .groups = "drop"
    )
  
  out <- tibble(Symbol = gene_list) %>%
    left_join(q3_wide, by = "Symbol") %>%
    left_join(call_summary, by = "Symbol") %>%
    left_join(regulation_summary, by = "Symbol") %>%
    mutate(
      across(
        ends_with("_Q3_LogFC"),
        ~ tidyr::replace_na(.x, 0)
      ),
      Upregulated_by = tidyr::replace_na(Upregulated_by, ""),
      Downregulated_by = tidyr::replace_na(Downregulated_by, ""),
      Regulation_summary = tidyr::replace_na(Regulation_summary, "")
    )
  
  wanted_q3_cols <- paste0(tf_order, "_Q3_LogFC")
  existing_q3_cols <- intersect(wanted_q3_cols, colnames(out))
  
  out <- out %>%
    select(
      Symbol,
      all_of(existing_q3_cols),
      Upregulated_by,
      Downregulated_by,
      Regulation_summary
    ) %>%
    arrange(Symbol)
  
  return(out)
}


proteoglycan_summary_sheet <- make_tf_summary_wide(
  gene_list = sort(proteo_genes),
  qtiles_df = gene_tf_qtiles,
  tf_order = emt_tfs,
  up_thr = 1,
  down_thr = -1
)

gag_summary_sheet <- make_tf_summary_wide(
  gene_list = sort(gags),
  qtiles_df = gene_tf_qtiles,
  tf_order = emt_tfs,
  up_thr = 1,
  down_thr = -1
)


writexl::write_xlsx(
  list(
    LogFC_Detailed = logfc_detailed_sheet,
    Proteoglycan_TF_Summary = proteoglycan_summary_sheet,
    GAG_TF_Summary = gag_summary_sheet
  ),
  path = table_s12_out
)

cat("Saved Excel file:\n", table_s12_out, "\n")
