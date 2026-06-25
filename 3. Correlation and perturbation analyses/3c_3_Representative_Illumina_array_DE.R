suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readxl)
  library(openxlsx)
  library(ggplot2)
})


# ============================================================
# 1. Paths and settings
# ============================================================

expression_file <- "GSE56164_Illumina_expression_matrix.xlsx"
metadata_file   <- "TF_perturbation_Illumina_metadata.xlsx"
metadata_sheet  <- 1

out_dir <- "Illumina_array_limma_outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dataset_id  <- "GSE56164"
platform_id <- "Illumina_expression_array"

probe_col  <- "ProbeID"
symbol_col <- "Symbol"

sample_col    <- "sample_id"
condition_col <- "condition"

contrast_table <- data.frame(
  contrast_name = c("KD_vs_Control"),
  reference_condition = c("WT"),
  case_condition = c("KD"),
  stringsAsFactors = FALSE
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
    condition = all_of(condition_col)
  ) %>%
  mutate(
    sample_id = as.character(sample_id),
    condition_raw = condition,
    condition = normalize_condition(condition)
  )

required_cols <- c("sample_id", "condition")
missing_cols <- setdiff(required_cols, colnames(metadata))

if (length(missing_cols) > 0) {
  stop("Missing metadata column(s): ", paste(missing_cols, collapse = ", "))
}


# ============================================================
# 4. Load Illumina expression matrix
# ============================================================

if (grepl("\\.xlsx$", expression_file, ignore.case = TRUE)) {
  expr_raw <- readxl::read_xlsx(expression_file) %>% as.data.frame()
} else if (grepl("\\.csv$", expression_file, ignore.case = TRUE)) {
  expr_raw <- read.csv(expression_file, check.names = FALSE)
} else if (grepl("\\.tsv$|\\.txt$", expression_file, ignore.case = TRUE)) {
  expr_raw <- read.delim(expression_file, check.names = FALSE)
} else {
  stop("Unsupported expression file format.")
}

required_expr_cols <- c(probe_col, symbol_col)
missing_expr_cols <- setdiff(required_expr_cols, colnames(expr_raw))

if (length(missing_expr_cols) > 0) {
  stop("Missing expression column(s): ", paste(missing_expr_cols, collapse = ", "))
}

available_samples <- intersect(metadata$sample_id, colnames(expr_raw))

if (length(available_samples) < 2) {
  stop("Fewer than two metadata samples matched the expression matrix.")
}

metadata <- metadata %>%
  filter(sample_id %in% available_samples) %>%
  arrange(match(sample_id, available_samples))

expr_raw <- expr_raw %>%
  select(
    ProbeID = all_of(probe_col),
    Symbol = all_of(symbol_col),
    all_of(metadata$sample_id)
  )

expr_values <- expr_raw %>%
  select(all_of(metadata$sample_id)) %>%
  as.matrix()

mode(expr_values) <- "numeric"

if (max(expr_values, na.rm = TRUE) > 100) {
  expr_values <- log2(expr_values + 1)
}

expr_values <- limma::normalizeBetweenArrays(
  expr_values,
  method = "quantile"
)

expr_norm_probe <- data.frame(
  ProbeID = expr_raw$ProbeID,
  Symbol = expr_raw$Symbol,
  expr_values,
  check.names = FALSE
)

expr_symbol <- expr_norm_probe %>%
  select(Symbol, all_of(metadata$sample_id)) %>%
  collapse_probe_to_symbol(symbol_col = "Symbol")

save(
  metadata,
  expr_raw,
  expr_norm_probe,
  expr_symbol,
  file = file.path(out_dir, paste0(dataset_id, "_Illumina_normalized_expression.RData"))
)


# ============================================================
# 5. Limma DE function
# ============================================================

run_illumina_limma <- function(expr_symbol,
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
    analysis_method = "Illumina_expression_array_quantile_limma"
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
  
  run_illumina_limma(
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
  file = file.path(out_dir, paste0(dataset_id, "_Illumina_DE_results_merged.csv")),
  row.names = FALSE
)

save(
  results_list,
  merged_de_results,
  file = file.path(out_dir, paste0(dataset_id, "_Illumina_DE_results_merged.RData"))
)


# ============================================================
# 7. Session info
# ============================================================

writeLines(
  capture.output(sessionInfo()),
  con = file.path(out_dir, "sessionInfo_Illumina_array_limma.txt")
)

message("Illumina expression-array representative workflow completed.")