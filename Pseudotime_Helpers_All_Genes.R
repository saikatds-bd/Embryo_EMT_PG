analyze_emt_vs_module <- function(
    sce, counts, asso, 
    emt_tfs, 
    module_genes,           
    module_label = c("Proteoglycans", "GAGs"),
    lineage_label = NULL,   
    alpha = 0.05,
    nPoints = 200
) {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(tradeSeq)
    library(ggplot2)
    library(colorspace)
    library(RColorBrewer)
  })
  filter    <- dplyr::filter
  select    <- dplyr::select
  mutate    <- dplyr::mutate
  summarise <- dplyr::summarise
  arrange   <- dplyr::arrange
  pull      <- dplyr::pull
  
  module_label <- match.arg(module_label)
  
  # Coerce module genes
  if (is.data.frame(module_genes) && "Name" %in% colnames(module_genes)) {
    module_genes <- module_genes$Name
  } else if (!is.vector(module_genes)) {
    stop("`module_genes` must be a vector of gene IDs or a data.frame with a 'Name' column.")
  }
  
  get_smooths_df <- function(sce, genes, nPoints = 200, group_label) {
    pt_grid <- seq(0, 1, length.out = nPoints)
    mats <- lapply(genes, function(g) {
      sm <- as.numeric(predictSmooth(sce, gene = g, nPoints = nPoints, tidy = FALSE))
      tibble(gene = g, pt = pt_grid, value = sm)
    })
    bind_rows(mats) %>% mutate(group = group_label)
  }
  
  # Pools limited to available genes
  module_pool <- intersect(module_genes, rownames(counts))
  emt_pool    <- intersect(emt_tfs,       rownames(counts))
  
  # Significant sets
  module_sig <- asso %>%
    filter( gene %in% module_pool) %>%
    arrange(padj) %>% pull(gene) %>% unique()
  
  emt_sig <- asso %>%
    filter( gene %in% emt_pool) %>%
    arrange(padj) %>% pull(gene) %>% unique()
  
  # Display picks
  emt_display <- intersect(c("SNAI1","SNAI2","TWIST1","ZEB1","ZEB2"), emt_sig)
  module_display <- head(module_sig, 10)
  
  # Smoothed per-gene z-scored curves for display sets
  df_gene <- bind_rows(
    get_smooths_df(sce, emt_display,   nPoints = nPoints, group_label = "EMT TFs"),
    get_smooths_df(sce, module_display, nPoints = nPoints, group_label = module_label)
  ) %>%
    group_by(gene) %>%
    mutate(z = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)) %>%
    ungroup()
  
  # Colors (unchanged logic)
  cols <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Set2"))
  
  title_suffix <- if (is.null(lineage_label)) "" else paste0(lineage_label, ": ")
  pairwise_gene_plot <-
    ggplot(df_gene, aes(pt, z, color = gene)) +
    geom_line(size = 0.9, alpha = 0.95) +
    scale_color_manual(values = cols) +
    facet_wrap(~ group, ncol = 1, scales = "free_y") +
    labs(
      title = paste0(title_suffix, "All genes smoothed expression (z-score)"),
      x = "Scaled pseudotime (0-1)", y = "Smoothed expression (z)"
    ) +
    theme_pubr() +
    theme(legend.position = "right")
  
  print(pairwise_gene_plot)
  
  # Module means across ALL significant genes (for correlation/peak comparisons)
  df_all <- bind_rows(
    get_smooths_df(sce, emt_sig,    nPoints = nPoints, group_label = "EMT TFs"),
    get_smooths_df(sce, module_sig, nPoints = nPoints, group_label = module_label)
  )
  
  module_means <- df_all %>%
    group_by(group, pt) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = mean_value)
  
  mean_module_plot <-
    ggplot(module_means %>% pivot_longer(-pt, names_to = "group", values_to = "mean_value"),
           aes(pt, mean_value, color = group)) +
    geom_line(size = 1.2) +
    labs(
      title = paste0(title_suffix, "EMT vs ", module_label, " module means (FDR-selected)"),
      x = "Scaled pseudotime (0-1)", y = "Mean smoothed log-expression"
    ) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "bottom")
  
  print(mean_module_plot)
  
  # Quick summaries
  peak_stats <- module_means %>%
    summarise(
      EMT_peak_pt = pt[which.max(`EMT TFs`)],
      MODULE_peak_pt  = pt[which.max(!!sym(module_label))],
      EMT_max     = max(`EMT TFs`, na.rm = TRUE),
      MODULE_max  = max(!!sym(module_label), na.rm = TRUE)
    ) %>%
    mutate(peak_lag = MODULE_peak_pt - EMT_peak_pt)
  
  scatter_df <- module_means %>%
    transmute(
      pt,
      EMT = `EMT TFs`,
      MODULE = !!sym(module_label),
      EMT_log = log2(1 + EMT),
      MODULE_log = log2(1 + MODULE)
    )
  
  
  emt_module_cor <- suppressWarnings(
    cor(scatter_df$EMT_log, scatter_df$MODULE_log,
        method = "pearson", use = "complete.obs")
  )
  
  # Pivot all smooths wide: one column per gene
  wide_smooths <- df_all %>%
    select(-group) %>%                      # <-- remove 'group' so id is only pt
    arrange(pt) %>%
    distinct(pt, gene, .keep_all = TRUE) %>%# guard against dupes
    pivot_wider(
      id_cols = pt,                         # <-- single time base
      names_from = gene,
      values_from = value
    )
  
  emt_genes_present    <- intersect(emt_sig,    colnames(wide_smooths))
  module_genes_present <- intersect(module_sig, colnames(wide_smooths))
  
  # Drop genes that are entirely NA after smoothing
  all_na <- function(v) all(is.na(v))
  emt_genes_present    <- emt_genes_present[!vapply(wide_smooths[emt_genes_present], all_na, TRUE)]
  module_genes_present <- module_genes_present[!vapply(wide_smooths[module_genes_present], all_na, TRUE)]
  
  if (length(emt_genes_present) == 0 || length(module_genes_present) == 0) {
    cor_pairs <- tibble(EMT_gene = character(), Module_gene = character(), correlation = double())
  } else {
    X <- as.matrix(wide_smooths[, emt_genes_present, drop = FALSE])
    Y <- as.matrix(wide_smooths[, module_genes_present, drop = FALSE])
    
    C <- suppressWarnings(cor(X, Y, use = "pairwise.complete.obs", method = "pearson"))
    
    cor_pairs <- C %>%
      as.data.frame(check.names = FALSE) %>%
      tibble::rownames_to_column("EMT_gene") %>%
      pivot_longer(-EMT_gene, names_to = "Module_gene", values_to = "correlation")
  }
  
  if (nrow(cor_pairs) == 0) {
    cor_pairs_wide <- tibble(Module_gene = character())
  } else {
    cor_pairs_wide <- cor_pairs %>%
      arrange(Module_gene, EMT_gene) %>%
      pivot_wider(
        id_cols   = Module_gene,
        names_from = EMT_gene,
        values_from = correlation
      )
  }
  
  library(viridis)
  phase_scatter_plot <- ggscatter(
    scatter_df, x = "EMT_log", y = "MODULE_log", color = "pt",
    add = "reg.line", conf.int = TRUE
  ) +
    stat_cor(method = "pearson") +
    theme(legend.position = "right") +
    geom_path(arrow = arrow(type = "closed", length = unit(0.12, "cm"))) +
    geom_point(size = 0.6) +
    theme_pubr() +
    scale_color_viridis(option = "magma") +
    labs(
      title = paste0(lineage_label, " phase plot"),
      x = "EMT module", y = paste(module_label, "module")
    )
  
  print(phase_scatter_plot)
  
  invisible(list(
    pairwise_gene_plot = pairwise_gene_plot,
    mean_module_plot = mean_module_plot,
    df_gene = df_gene,
    df_all = df_all,
    module_means = module_means,
    peak_stats = peak_stats,
    correlation = emt_module_cor,
    emt_display = emt_display,
    module_display = module_display,
    module_sig = module_sig,
    emt_sig = emt_sig,
    cor_pairs = cor_pairs,
    module_label = module_label,
    lineage_label = lineage_label,
    phase_scatter_plot = phase_scatter_plot,
    cor_pairs_wide = cor_pairs_wide
  ))
}

# Convenience runners for your two lineages and two module types ----

run_endo_proteoglycans <- function() {
  analyze_emt_vs_module(
    sce = sce_endo, counts = counts_endo, asso = asso_endo,
    emt_tfs = emt_tfs,
    module_genes = proteoglycans,   # accepts proteoglycans$Name or a vector
    module_label = "Proteoglycans",
    lineage_label = "Endoderm"
  )
}

run_endo_gags <- function() {
  analyze_emt_vs_module(
    sce = sce_endo, counts = counts_endo, asso = asso_endo,
    emt_tfs = emt_tfs,
    module_genes = gags,            # accepts gags$Name or a vector
    module_label = "GAGs",
    lineage_label = "Endoderm"
  )
}

run_meso_proteoglycans <- function() {
  analyze_emt_vs_module(
    sce = sce_meso, counts = counts_meso, asso = asso_meso,
    emt_tfs = emt_tfs,
    module_genes = proteoglycans,
    module_label = "Proteoglycans",
    lineage_label = "Mesoderm"
  )
}

run_meso_gags <- function() {
  analyze_emt_vs_module(
    sce = sce_meso, counts = counts_meso, asso = asso_meso,
    emt_tfs = emt_tfs,
    module_genes = gags,
    module_label = "GAGs",
    lineage_label = "Mesoderm"
  )
}


phase_scatter_from_res <- function(res, tf_gene, module_gene,
                                   group_tf = "EMT TFs",
                                   group_module = "Proteoglycans",
                                   zscore = FALSE,
                                   title_prefix = NULL) {
  df <- res$df_all
  if (is.null(df) || !"pt" %in% names(df)) stop("res$df_all is missing.")
  
  # Pull smoothed series for the two genes
  xdat <- df %>% filter(gene == tf_gene, group == group_tf) %>% select(pt, x = value)
  ydat <- df %>% filter(gene == module_gene, group == group_module) %>% select(pt, y = value)
  
  if (nrow(xdat) == 0) stop(sprintf("TF '%s' not found in df_all (likely not FDR-significant in this lineage).", tf_gene))
  if (nrow(ydat) == 0) stop(sprintf("Module gene '%s' not found in df_all (likely not FDR-significant in this lineage).", module_gene))
  
  dat <- left_join(xdat, ydat, by = "pt")
  
  # optional z-score (shape-only comparison)
  if (zscore) {
    dat <- dat %>%
      mutate(x = (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE),
             y = (y - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE))
    x_lab <- paste0(tf_gene, " (z)")
    y_lab <- paste0(module_gene, " (z)")
  } else {
    x_lab <- tf_gene
    y_lab <- module_gene
  }
  
  rP <- suppressWarnings(cor(dat$x, dat$y, method = "pearson",  use = "complete.obs"))
  rS <- suppressWarnings(cor(dat$x, dat$y, method = "spearman", use = "complete.obs"))
  
  ttl <- paste0(ifelse(is.null(title_prefix), "", paste0(title_prefix, ": ")),
                module_gene, " vs ", tf_gene, " (smoothed)")
  
  p <- ggscatter(dat, x = "x", y = "y", color = "pt",
                 add = "reg.line", conf.int = TRUE) +
    scale_color_viridis(option = "magma", name = "Pseudotime") +
    labs(title = ttl, x = x_lab, y = y_lab,
         subtitle = sprintf("Pearson r = %.2f", rP)) +
    theme_pubr() +
    theme(legend.position = "right") +
    geom_path(data = dat %>% arrange(pt),
              aes(x = x, y = y), inherit.aes = FALSE,
              arrow = arrow(type = "closed", length = unit(0.12, "cm")),
              linewidth = 0.3, alpha = 0.7)
  
  print(p)
  invisible(list(plot = p, data = dat, r_pearson = rP, r_spearman = rS))
}
