suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readxl)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})


# ============================================================
# 1. Paths and settings
# ============================================================

raw_dir <- "raw_Agilent_files"
metadata_file <- "TF_perturbation_Agilent_metadata.xlsx"
metadata_sheet <- 1

out_dir <- "Agilent_single_channel_limma_outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dataset_id  <- "GSE155107"
platform_id <- "Agilent_single_channel_array"

sample_col    <- "sample_id"
file_col      <- "raw_file"
condition_col <- "condition"

contrast_table <- data.frame(
  contrast_name = c("KD_vs_WT"),
  reference_condition = c("WT"),
  case_condition = c("KD"),
  stringsAsFactors = FALSE
)

annotation_priority <- c(
  "GeneName",
  "GeneSymbol",
  "SYMBOL",
  "GENE_SYMBOL",
  "REFSEQ",
  "ENSEMBLTRANS",
  "SystematicName",
  "ProbeName"
)


# ============================================================
# 2. Helper functions
# ============================================================

make_safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}


normalize_condition <- function(x) {
  x <- as.character(x)
  
  x <- str_replace_all(x, regex("^control$|^ctrl$|^cntrl$|^untreated$|^mock$|^wt$", ignore_case = TRUE), "WT")
  x <- str_replace_all(x, regex("overexpression|over_exp|overexp|^oe$", ignore_case = TRUE), "OE")
  x <- str_replace_all(x, regex("knockdown|^kd$|^sh$|shRNA|siRNA", ignore_case = TRUE), "KD")
  x <- str_replace_all(x, regex("knockout|^ko$|CRISPR", ignore_case = TRUE), "KO")
  
  x
}


collapse_probe_to_symbol <- function(df, symbol_col = "Symbol") {
  
  df %>%
    filter(!is.na(.data[[symbol_col]]), .data[[symbol_col]] != "") %>%
    group_by(.data[[symbol_col]]) %>%
    summarise(
      across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    rename(Symbol = all_of(symbol_col))
}


extract_gene_symbol <- function(genes_df) {
  
  available_cols <- intersect(annotation_priority, colnames(genes_df))
  
  if (length(available_cols) == 0) {
    stop(
      "No recognized annotation column found. Available columns are: ",
      paste(colnames(genes_df), collapse = ", ")
    )
  }
  
  symbol_source <- available_cols[1]
  
  symbol <- genes_df[[symbol_source]]
  
  if (symbol_source == "REFSEQ") {
    symbol <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = as.character(symbol),
      column = "SYMBOL",
      keytype = "REFSEQ",
      multiVals = "first"
    )
  }
  
  if (symbol_source == "ENSEMBLTRANS") {
    symbol <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = as.character(symbol),
      column = "SYMBOL",
      keytype = "ENSEMBLTRANS",
      multiVals = "first"
    )
  }
  
  data.frame(
    ProbeID = rownames(genes_df),
    Symbol = as.character(symbol),
    annotation_source = symbol_source,
    stringsAsFactors = FALSE
  )
}


standardize_de_output <- function(de_table,
                                  dataset_id,
                                  platform_id,
                                  contrast_name,
                                  ref_condition,
                                  case_condition,
                                  metadata_subset,
                                  analysis_method) {
  
  tf_value <- if ("tf" %in% colnames(metadata_subset)) unique(metadata_subset$tf) else NA
  perturbation_value <- if ("perturbation" %in% colnames(metadata_subset)) unique(metadata_subset$perturbation) else NA
  cell_model_value <- if ("cell_model" %in% colnames(metadata_subset)) unique(metadata_subset$cell_model) else NA
  
  de_table %>%
    mutate(
      dataset_id = dataset_id,
      platform = platform_id,
      TF = paste(na.omit(tf_value), collapse = ";"),
      perturbation = paste(na.omit(perturbation_value), collapse = ";"),
      cell_model = paste(na.omit(cell_model_value), collapse = ";"),
      contrast = contrast_name,
      reference_condition = ref_condition,
      case_condition = case_condition,
      analysis_method = analysis_method
    ) %>%
    relocate(
      dataset_id,
      platform,
      TF,
      perturbation,
      cell_model,
      contrast,
      reference_condition,
      case_condition,
      analysis_method,
      Symbol
    )
}


make_mds_plot <- function(expr_symbol, sample_data, title_text) {
  
  mat <- expr_symbol %>%
    column_to_rownames("Symbol") %>%
    as.matrix()
  
  mds <- limma::plotMDS(mat, plot = FALSE)
  
  plot_df <- data.frame(
    Dim1 = mds$x,
    Dim2 = mds$y,
    sample_id = colnames(mat),
    stringsAsFactors = FALSE
  ) %>%
    left_join(sample_data, by = "sample_id")
  
  ggplot(plot_df, aes(x = Dim1, y = Dim2)) +
    geom_point(aes(colour = condition), size = 3) +
    geom_text(aes(label = sample_id), vjust = -0.8, size = 3, show.legend = FALSE) +
    theme_bw(base_size = 11) +
    labs(
      title = title_text,
      x = "Leading logFC dimension 1",
      y = "Leading logFC dimension 2",
      colour = "Condition"
    )
}


# ============================================================
# 3. Load metadata
# ============================================================

metadata <- readxl::read_xlsx(metadata_file, sheet = metadata_sheet) %>%
  as.data.frame() %>%
  rename(
    sample_id = all_of(sample_col),
    raw_file = all_of(file_col),
    condition = all_of(condition_col)
  ) %>%
  mutate(
    sample_id = as.character(sample_id),
    raw_file = as.character(raw_file),
    condition_raw = condition,
    condition = normalize_condition(condition)
  )

required_cols <- c("sample_id", "raw_file", "condition")
missing_cols <- setdiff(required_cols, colnames(metadata))

if (length(missing_cols) > 0) {
  stop("Missing metadata column(s): ", paste(missing_cols, collapse = ", "))
}

raw_files <- file.path(raw_dir, metadata$raw_file)

if (!all(file.exists(raw_files))) {
  missing_files <- raw_files[!file.exists(raw_files)]
  stop("Missing Agilent raw file(s):\n", paste(missing_files, collapse = "\n"))
}


# ============================================================
# 4. Read and normalize Agilent files
# ============================================================

targets <- data.frame(
  FileName = raw_files,
  SampleID = metadata$sample_id,
  stringsAsFactors = FALSE
)

raw_agilent <- limma::read.maimages(
  files = targets$FileName,
  source = "agilent",
  green.only = TRUE
)

colnames(raw_agilent$E) <- metadata$sample_id

bg_corrected <- limma::backgroundCorrect(
  raw_agilent,
  method = "normexp",
  offset = 50
)

norm_agilent <- limma::normalizeBetweenArrays(
  bg_corrected,
  method = "quantile"
)

expr_probe <- norm_agilent$E
colnames(expr_probe) <- metadata$sample_id

gene_mapping <- extract_gene_symbol(norm_agilent$genes)

expr_probe_df <- as.data.frame(expr_probe) %>%
  rownames_to_column("ProbeID") %>%
  left_join(gene_mapping, by = "ProbeID")

expr_symbol <- expr_probe_df %>%
  select(Symbol, all_of(metadata$sample_id)) %>%
  collapse_probe_to_symbol(symbol_col = "Symbol")

save(
  metadata,
  raw_agilent,
  bg_corrected,
  norm_agilent,
  gene_mapping,
  expr_probe_df,
  expr_symbol,
  file = file.path(out_dir, paste0(dataset_id, "_Agilent_normalized_expression.RData"))
)


# ============================================================
# 5. Limma DE function
# ============================================================

run_agilent_limma <- function(expr_symbol,
                              metadata,
                              contrast_name,
                              ref_condition,
                              case_condition) {
  
  message("Running contrast: ", contrast_name)
  
  md <- metadata %>%
    filter(condition %in% c(ref_condition, case_condition))
  
  if (nrow(md) < 2) {
    message("Skipping ", contrast_name, ": fewer than two samples.")
    return(NULL)
  }
  
  replicate_table <- table(md$condition)
  
  if (!all(replicate_table[c(ref_condition, case_condition)] >= 2)) {
    message("Skipping ", contrast_name, ": insufficient biological replication.")
    return(NULL)
  }
  
  expr <- expr_symbol %>%
    select(Symbol, all_of(md$sample_id))
  
  mat <- expr %>%
    column_to_rownames("Symbol") %>%
    as.matrix()
  
  md <- md[match(colnames(mat), md$sample_id), ]
  
  group <- factor(md$condition)
  group <- relevel(group, ref = ref_condition)
  
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  fit <- limma::lmFit(mat, design)
  
  contrast_string <- paste0(case_condition, " - ", ref_condition)
  
  contrast_matrix <- limma::makeContrasts(
    contrasts = contrast_string,
    levels = design
  )
  
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2, robust = TRUE)
  
  de <- limma::topTable(
    fit2,
    number = Inf,
    sort.by = "P",
    adjust.method = "BH"
  ) %>%
    rownames_to_column("Symbol")
  
  expr_out <- as.data.frame(mat) %>%
    rownames_to_column("Symbol")
  
  de_out <- de %>%
    left_join(expr_out, by = "Symbol")
  
  de_out <- standardize_de_output(
    de_table = de_out,
    dataset_id = dataset_id,
    platform_id = platform_id,
    contrast_name = contrast_name,
    ref_condition = ref_condition,
    case_condition = case_condition,
    metadata_subset = md,
    analysis_method = "Agilent_single_channel_normexp_quantile_limma"
  )
  
  p_mds <- make_mds_plot(
    expr_symbol = expr,
    sample_data = md,
    title_text = paste0(dataset_id, ": ", contrast_name)
  )
  
  out_prefix <- file.path(
    out_dir,
    paste0(make_safe_filename(dataset_id), "_", make_safe_filename(contrast_name))
  )
  
  save(
    de_out,
    md,
    fit2,
    p_mds,
    file = paste0(out_prefix, "_DE_object.RData")
  )
  
  write.csv(
    de_out,
    file = paste0(out_prefix, "_DE_results.csv"),
    row.names = FALSE
  )
  
  ggsave(
    filename = paste0(out_prefix, "_MDS_plot.pdf"),
    plot = p_mds,
    width = 6,
    height = 5
  )
  
  de_out
}


# ============================================================
# 6. Run contrasts
# ============================================================

results_list <- lapply(seq_len(nrow(contrast_table)), function(i) {
  
  run_agilent_limma(
    expr_symbol = expr_symbol,
    metadata = metadata,
    contrast_name = contrast_table$contrast_name[i],
    ref_condition = contrast_table$reference_condition[i],
    case_condition = contrast_table$case_condition[i]
  )
})

names(results_list) <- contrast_table$contrast_name
results_list <- results_list[!vapply(results_list, is.null, logical(1))]

merged_de_results <- bind_rows(results_list)

write.csv(
  merged_de_results,
  file = file.path(out_dir, paste0(dataset_id, "_Agilent_DE_results_merged.csv")),
  row.names = FALSE
)

save(
  results_list,
  merged_de_results,
  file = file.path(out_dir, paste0(dataset_id, "_Agilent_DE_results_merged.RData"))
)


# ============================================================
# 7. Session info
# ============================================================

writeLines(
  capture.output(sessionInfo()),
  con = file.path(out_dir, "sessionInfo_Agilent_single_channel_limma.txt")
)

message("Agilent single-channel representative workflow completed.")