library(readxl)
library(dplyr)
library(purrr)
library(limma)
library(edgeR)
library(ggpubr)

load("Interesting_GAG_Genes.RData")

## We will merge all the excel files or raw counts together first

excel_files <- list.files(pattern = "\\.xlsx$", full.names = TRUE)

expr_list <- lapply(excel_files, function(file) {
  df <- read_excel(file)
  df <- dplyr::select(df, contains(c("Name", "Total")))
  return(df)
})


names(expr_list) <- tools::file_path_sans_ext(basename(excel_files))

expr_list[["Metadata"]] <- NULL

all_metadata <- readxl::read_xlsx("Metadata.xlsx")
merged_expr <- reduce(expr_list, full_join, by = "Name")
colnames(merged_expr) <- gsub(" .*", "", colnames(merged_expr))
colnames(merged_expr) <- gsub("_\\d+$", "", gsub(" .*", "", colnames(merged_expr)))

project_col <- "Bioproject"

normalize_condition <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "(?i)control|ctrl|cntrl", "WT")
  x <- str_replace_all(x, "(?i)oe|overexp|overexpression", "OE")
  x <- str_replace_all(x, "(?i)knockout|ko", "KO")
  x <- str_replace_all(x, "(?i)knockdown|kd|sh", "KD")
  x
}
all_metadata <- all_metadata %>% mutate(Condition = normalize_condition(Condition))


analyze_project <- function(project_id, metadata, expr_mat) {
  md <- metadata %>% filter(.data[[project_col]] == project_id)
  
  runs <- intersect(md$Run, colnames(expr_mat))
  if (length(runs) < 2) {
    message(project_id, ": <2 samples after matching runs — skipping.")
    return(NULL)
  }
  
  # subset & order counts by metadata order
  cm <- expr_mat[, runs, drop = FALSE]
  md <- md[match(colnames(cm), md$Run), , drop = FALSE]
  
  # groups
  group <- factor(md$Condition)
  if (nlevels(group) < 2) {
    message(project_id, ": only one condition present — skipping.")
    return(NULL)
  }
  
  # Fixing the contrasts
  levs <- levels(group)
  cond1 <- levs[1]
  cond2 <- levs[2]
  contrast_name <- paste0(cond2, "vs", cond1)
  
  # Check replicates per condition
  tab <- table(group)
  has_reps <- all(tab >= 2)
  
  if (has_reps) {
    ## limma pipeline
    y <- DGEList(counts = cm, group = group)
    y$samples$Sample_Names <- colnames(cm)
    
    # if metadta has replicates info:
    if ("Replicate" %in% colnames(md)) {
      y$samples$replicate <- md$Replicate
    } else {
      y$samples$replicate <- factor(seq_len(nrow(y$samples)))
    }
    
    keep <- filterByExpr(y, group = y$samples$group)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y, method = "TMM")
    
    lcpm <- edgeR::cpm(y, log = TRUE, prior.count = 0.25)
    
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    
    v <- voomWithQualityWeights(y, design, plot = FALSE)
    vfit <- lmFit(v, design)
    
    # build contrast like "OE - WT"
    ct_string <- paste0(cond2, " - ", cond1)
    contr <- makeContrasts(contrasts = ct_string, levels = colnames(vfit$coefficients))
    cfit <- contrasts.fit(vfit, contr)
    efit <- eBayes(cfit, robust = TRUE)
    
    top <- topTable(efit, sort.by = "P", n = Inf, adjust.method = "BH")
    top$Symbol <- rownames(top)
    
    lcpm_df <- as.data.frame(lcpm)
    lcpm_df$Symbol <- rownames(lcpm_df)
    res <- merge(top, lcpm_df, by = "Symbol")
    
    # MDS plots will be checked again to find batch effects
    p_mds <- plotMDS_ggpubr(
      lcpm = lcpm, 
      sample_data = y$samples, 
      color_var = "replicate", # Use 'Condition' column for colors
      label_var = "Sample_Names", # Use 'Replicate' column for labels
      palette = "Set2",         # Change color palette
      title = "Original MDS Plot",
      dims=c(1,2)
    )
    
    # Save
    out_rdata <- paste0(project_id, "_DE_", contrast_name, ".RData")
    meta_copy <- md
    save(res, meta_copy, p_mds, file = out_rdata)
    message(project_id, ": replicated analysis saved -> ", out_rdata)
    return(invisible(res))
    
  } else {
    ## Single-rep: manual logFC
    # Pseudocount
    df <- cm + 0.1
    
    # If there are multiple samples in a condition, using mean
    g1 <- which(group == cond1)
    g2 <- which(group == cond2)
    if (length(g1) == 0 || length(g2) == 0) {
      message(project_id, ": missing one condition — skipping.")
      return(NULL)
    }
    
    logFC <- log2(rowMeans(df[, g2, drop = FALSE])) - log2(rowMeans(df[, g1, drop = FALSE]))
    res <- data.frame(
      Symbol = rownames(df),
      logFC = logFC,
      AveExpr = rowMeans(log2(df)),   
      stringsAsFactors = FALSE
    )
    
    # attach the raw columns for traceability
    res <- cbind(res, as.data.frame(df))
    
    out_rdata <- paste0(project_id, "_SingleRep_", cond2, "vs", cond1, ".RData")
    meta_copy <- md
    save(res, meta_copy, file = out_rdata)
    message(project_id, ": single-rep analysis saved -> ", out_rdata)
    return(invisible(res))
  }
}

# running over all projects that appear in the metadata
all_projects <- unique(all_metadata[[project_col]])
results_list <- lapply(all_projects, analyze_project, metadata = all_metadata, expr_mat = expr_mat)
names(results_list) <- all_projects
save(results_list, file = "TF_DE.RData")
