# ============================================================
# 01_RNAseq_counts_edgeR_limma_voom.R
# ============================================================
#
# Representative workflow for TF perturbation datasets analyzed
# from RNA-seq/raw-count matrices.
#
# This script illustrates the count-based differential-expression
# workflow used for TF perturbation datasets where raw counts or
# count-like matrices were available.
#


suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(edgeR)
  library(limma)
  library(ggplot2)
})




# Folder containing raw count files.
# Each file should contain one gene column and one or more count columns.
count_dir <- "raw_counts"

# Metadata file describing datasets, samples, conditions and contrasts.
metadata_file <- "TF_perturbation_RNAseq_metadata.xlsx"

# Output directory
out_dir <- "RNAseq_edgeR_limma_voom_outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Optional: name of metadata sheet
metadata_sheet <- 1

# Column names expected in the metadata file
dataset_col   <- "dataset_id"
sample_col    <- "sample_id"
condition_col <- "condition"

# Name of the gene column in count files.
# Change this based on the files
gene_col <- "Name"


normalize_condition <- function(x) {
  x <- as.character(x)
  
  x <- str_replace_all(x, regex("^control$|^ctrl$|^cntrl$|^untreated$|^mock$", ignore_case = TRUE), "WT")
  x <- str_replace_all(x, regex("overexpression|over_exp|overexp|^oe$", ignore_case = TRUE), "OE")
  x <- str_replace_all(x, regex("knockdown|^kd$|^sh$|shRNA|siRNA", ignore_case = TRUE), "KD")
  x <- str_replace_all(x, regex("knockout|^ko$|CRISPR", ignore_case = TRUE), "KO")
  
  x
}


clean_count_colnames <- function(x) {
  x <- as.character(x)
  
  # Remove extra text after spaces
  x <- gsub(" .*", "", x)
  
  # Remove common suffixes introduced during import
  x <- gsub("_\\d+$", "", x)
  
  x
}


make_safe_filename <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}


read_count_file <- function(file, gene_col = "Name") {
  
  message("Reading count file: ", basename(file))
  
  if (grepl("\\.xlsx$", file, ignore.case = TRUE)) {
    df <- readxl::read_xlsx(file)
  } else if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    df <- read.csv(file, check.names = FALSE)
  } else if (grepl("\\.tsv$|\\.txt$", file, ignore.case = TRUE)) {
    df <- read.delim(file, check.names = FALSE)
  } else {
    stop("Unsupported count file format: ", file)
  }
  
  if (!gene_col %in% colnames(df)) {
    stop("Gene column '", gene_col, "' was not found in ", basename(file))
  }
  
  df <- df %>%
    dplyr::rename(Symbol = all_of(gene_col)) %>%
    dplyr::filter(!is.na(Symbol), Symbol != "")
  
  # Keep only gene column and numeric count columns
  count_cols <- setdiff(colnames(df), "Symbol")
  
  df[count_cols] <- lapply(df[count_cols], function(z) {
    suppressWarnings(as.numeric(z))
  })
  
  df <- df %>%
    dplyr::select(Symbol, where(is.numeric))
  
  colnames(df) <- clean_count_colnames(colnames(df))
  
  df
}


collapse_duplicate_genes <- function(count_df) {
  
  count_df %>%
    dplyr::group_by(Symbol) %>%
    dplyr::summarise(
      dplyr::across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )
}


make_mds_plot <- function(lcpm, sample_data, colour_by = "condition", label_by = "sample_id",
                          title = "MDS plot") {
  
  mds <- limma::plotMDS(lcpm, plot = FALSE)
  
  plot_df <- data.frame(
    Dim1 = mds$x,
    Dim2 = mds$y,
    sample_id = colnames(lcpm),
    stringsAsFactors = FALSE
  )
  
  sample_data <- as.data.frame(sample_data)
  sample_data$sample_id <- as.character(sample_data$sample_id)
  
  plot_df <- plot_df %>%
    dplyr::left_join(sample_data, by = "sample_id")
  
  p <- ggplot(plot_df, aes(x = Dim1, y = Dim2)) +
    geom_point(aes(colour = .data[[colour_by]]), size = 3) +
    geom_text(aes(label = .data[[label_by]]), vjust = -0.8, size = 3, show.legend = FALSE) +
    theme_bw(base_size = 11) +
    labs(
      title = title,
      x = "Leading logFC dimension 1",
      y = "Leading logFC dimension 2",
      colour = colour_by
    )
  
  p
}


standardize_de_output <- function(de_table,
                                  dataset_id,
                                  contrast_name,
                                  ref_condition,
                                  case_condition,
                                  metadata_subset,
                                  analysis_method) {
  
  tf_value <- if ("tf" %in% colnames(metadata_subset)) unique(metadata_subset$tf) else NA
  perturbation_value <- if ("perturbation" %in% colnames(metadata_subset)) unique(metadata_subset$perturbation) else NA
  cell_model_value <- if ("cell_model" %in% colnames(metadata_subset)) unique(metadata_subset$cell_model) else NA
  
  tf_value <- paste(na.omit(tf_value), collapse = ";")
  perturbation_value <- paste(na.omit(perturbation_value), collapse = ";")
  cell_model_value <- paste(na.omit(cell_model_value), collapse = ";")
  
  de_table %>%
    dplyr::mutate(
      dataset_id = dataset_id,
      TF = tf_value,
      perturbation = perturbation_value,
      cell_model = cell_model_value,
      contrast = contrast_name,
      reference_condition = ref_condition,
      case_condition = case_condition,
      analysis_method = analysis_method
    ) %>%
    dplyr::relocate(
      dataset_id,
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


metadata <- readxl::read_xlsx(metadata_file, sheet = metadata_sheet) %>%
  as.data.frame()

required_cols <- c(dataset_col, sample_col, condition_col)

missing_cols <- setdiff(required_cols, colnames(metadata))

if (length(missing_cols) > 0) {
  stop(
    "The metadata file is missing required column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

metadata <- metadata %>%
  dplyr::rename(
    dataset_id = all_of(dataset_col),
    sample_id = all_of(sample_col),
    condition = all_of(condition_col)
  ) %>%
  dplyr::mutate(
    dataset_id = as.character(dataset_id),
    sample_id = clean_count_colnames(sample_id),
    condition_raw = condition,
    condition = normalize_condition(condition)
  )

if (!"replicate" %in% colnames(metadata)) {
  metadata$replicate <- metadata$sample_id
}



count_files <- list.files(
  path = count_dir,
  pattern = "\\.xlsx$|\\.csv$|\\.tsv$|\\.txt$",
  full.names = TRUE
)

if (length(count_files) == 0) {
  stop("No count files found in: ", count_dir)
}

count_list <- lapply(count_files, read_count_file, gene_col = gene_col)

names(count_list) <- tools::file_path_sans_ext(basename(count_files))

merged_counts <- purrr::reduce(count_list, dplyr::full_join, by = "Symbol")

merged_counts <- collapse_duplicate_genes(merged_counts)

count_mat <- merged_counts %>%
  tibble::column_to_rownames("Symbol") %>%
  as.matrix()

mode(count_mat) <- "numeric"

colnames(count_mat) <- clean_count_colnames(colnames(count_mat))

# Replace missing values introduced by merging with zero
count_mat[is.na(count_mat)] <- 0



run_rnaseq_de_for_dataset <- function(dataset_id, metadata, count_mat) {
  
  message("\n============================================================")
  message("Dataset: ", dataset_id)
  message("============================================================")
  
  md <- metadata %>%
    dplyr::filter(dataset_id == !!dataset_id)
  
  matched_samples <- intersect(md$sample_id, colnames(count_mat))
  
  if (length(matched_samples) < 2) {
    message(dataset_id, ": fewer than two matched samples. Skipping.")
    return(NULL)
  }
  
  cm <- count_mat[, matched_samples, drop = FALSE]
  
  md <- md %>%
    dplyr::filter(sample_id %in% matched_samples) %>%
    dplyr::arrange(match(sample_id, colnames(cm)))
  
  cm <- cm[, md$sample_id, drop = FALSE]
  
  group <- factor(md$condition)
  
  if (nlevels(group) < 2) {
    message(dataset_id, ": only one condition present. Skipping.")
    return(NULL)
  }
  
  # ------------------------------------------------------------
  # Determine contrast
  # ------------------------------------------------------------
  # If contrast_ref and contrast_case are supplied in metadata,
  # use them. Otherwise use the first two condition levels.
  # For publication use, it is best to provide contrast_ref and
  # contrast_case explicitly in the metadata table.
  
  if (all(c("contrast_ref", "contrast_case") %in% colnames(md))) {
    
    ref_condition <- unique(na.omit(md$contrast_ref))
    case_condition <- unique(na.omit(md$contrast_case))
    
    ref_condition <- normalize_condition(ref_condition)
    case_condition <- normalize_condition(case_condition)
    
    if (length(ref_condition) != 1 || length(case_condition) != 1) {
      stop(dataset_id, ": contrast_ref and contrast_case must each contain one unique value.")
    }
    
  } else {
    
    levs <- levels(group)
    
    if ("WT" %in% levs) {
      ref_condition <- "WT"
      case_condition <- setdiff(levs, ref_condition)[1]
    } else {
      ref_condition <- levs[1]
      case_condition <- levs[2]
    }
  }
  
  if (!all(c(ref_condition, case_condition) %in% levels(group))) {
    stop(
      dataset_id,
      ": contrast conditions not found in group levels. ",
      "Reference = ", ref_condition,
      "; Case = ", case_condition,
      "; Available = ", paste(levels(group), collapse = ", ")
    )
  }
  
  group <- relevel(group, ref = ref_condition)
  
  contrast_name <- paste0(case_condition, "_vs_", ref_condition)
  contrast_name_safe <- make_safe_filename(contrast_name)
  dataset_id_safe <- make_safe_filename(dataset_id)
  
  message("Contrast: ", contrast_name)
  

  replicate_table <- table(group)
  
  message("Sample numbers:")
  print(replicate_table)
  
  has_replicates <- all(replicate_table[c(ref_condition, case_condition)] >= 2)
  
  # Keep only samples from the target contrast
  keep_samples <- group %in% c(ref_condition, case_condition)
  
  cm <- cm[, keep_samples, drop = FALSE]
  md <- md[keep_samples, , drop = FALSE]
  group <- droplevels(group[keep_samples])
  
  # ------------------------------------------------------------
  # Replicated analysis: edgeR + limma-voom
  # ------------------------------------------------------------
  
  if (has_replicates) {
    
    y <- edgeR::DGEList(counts = cm, group = group)
    
    y$samples$sample_id <- colnames(cm)
    y$samples$condition <- group
    y$samples$replicate <- md$replicate
    
    keep_genes <- edgeR::filterByExpr(y, group = group)
    
    y <- y[keep_genes, , keep.lib.sizes = FALSE]
    
    y <- edgeR::calcNormFactors(y, method = "TMM")
    
    lcpm <- edgeR::cpm(y, log = TRUE, prior.count = 0.25)
    
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    
    v <- limma::voomWithQualityWeights(y, design, plot = FALSE)
    
    fit <- limma::lmFit(v, design)
    
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
      tibble::rownames_to_column("Symbol")
    
    lcpm_df <- as.data.frame(lcpm) %>%
      tibble::rownames_to_column("Symbol")
    
    de_out <- de %>%
      dplyr::left_join(lcpm_df, by = "Symbol")
    
    de_out <- standardize_de_output(
      de_table = de_out,
      dataset_id = dataset_id,
      contrast_name = contrast_name,
      ref_condition = ref_condition,
      case_condition = case_condition,
      metadata_subset = md,
      analysis_method = "edgeR_TMM_limma_voomWithQualityWeights"
    )
    
    p_mds <- make_mds_plot(
      lcpm = lcpm,
      sample_data = md,
      colour_by = "condition",
      label_by = "sample_id",
      title = paste0(dataset_id, ": ", contrast_name)
    )
    
    # Save outputs
    out_prefix <- file.path(out_dir, paste0(dataset_id_safe, "_", contrast_name_safe))
    
    save(
      de_out,
      md,
      y,
      lcpm,
      v,
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
    
    message("Saved replicated DE output: ", out_prefix)
    
    return(de_out)
  }
  
  # ------------------------------------------------------------
  # Non-replicated analysis: descriptive log2 fold-change only
  # ------------------------------------------------------------
  #
  # This branch is included because some public perturbation
  # datasets lack sufficient biological replication. No moderated
  # p-values are estimated here. The output is retained only as a
  # descriptive log2 fold-change table.
  
  if (!has_replicates) {
    
    message(dataset_id, ": insufficient replication. Calculating descriptive logFC only.")
    
    cm_pc <- cm + 0.1
    
    ref_samples <- which(group == ref_condition)
    case_samples <- which(group == case_condition)
    
    log_expr <- log2(cm_pc)
    
    logFC <- rowMeans(log_expr[, case_samples, drop = FALSE]) -
      rowMeans(log_expr[, ref_samples, drop = FALSE])
    
    AveExpr <- rowMeans(log_expr)
    
    de_out <- data.frame(
      Symbol = rownames(cm_pc),
      logFC = logFC,
      AveExpr = AveExpr,
      P.Value = NA_real_,
      adj.P.Val = NA_real_,
      stringsAsFactors = FALSE
    )
    
    raw_count_df <- as.data.frame(cm) %>%
      tibble::rownames_to_column("Symbol")
    
    de_out <- de_out %>%
      dplyr::left_join(raw_count_df, by = "Symbol")
    
    de_out <- standardize_de_output(
      de_table = de_out,
      dataset_id = dataset_id,
      contrast_name = contrast_name,
      ref_condition = ref_condition,
      case_condition = case_condition,
      metadata_subset = md,
      analysis_method = "descriptive_log2FC_no_replicated_DE_test"
    )
    
    out_prefix <- file.path(out_dir, paste0(dataset_id_safe, "_", contrast_name_safe))
    
    save(
      de_out,
      md,
      file = paste0(out_prefix, "_descriptive_logFC_object.RData")
    )
    
    write.csv(
      de_out,
      file = paste0(out_prefix, "_descriptive_logFC_results.csv"),
      row.names = FALSE
    )
    
    message("Saved descriptive logFC output: ", out_prefix)
    
    return(de_out)
  }
}


# ============================================================
# 6. Run analysis over all RNA-seq TF perturbation datasets
# ============================================================

all_datasets <- unique(metadata$dataset_id)

results_list <- lapply(
  all_datasets,
  run_rnaseq_de_for_dataset,
  metadata = metadata,
  count_mat = count_mat
)

names(results_list) <- all_datasets

results_list <- results_list[!vapply(results_list, is.null, logical(1))]

save(
  results_list,
  file = file.path(out_dir, "RNAseq_TF_perturbation_DE_results_list.RData")
)



merged_de_results <- dplyr::bind_rows(results_list)

write.csv(
  merged_de_results,
  file = file.path(out_dir, "RNAseq_TF_perturbation_DE_results_merged.csv"),
  row.names = FALSE
)

save(
  merged_de_results,
  file = file.path(out_dir, "RNAseq_TF_perturbation_DE_results_merged.RData")
)
