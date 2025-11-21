suppressPackageStartupMessages({
  library(tidyverse)
  library(tradeSeq)
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
})

load("Separate_CDS.RData")
load("Interesting_GAG_Genes.RData")
set.seed(123)
emt_tfs <- c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2")
# Loading Colors
load("Proteoglycan_Project_Colors.RData")
interpolate_lancet <- function(n) {
  base <- ggsci::pal_lancet()(9)
  grDevices::colorRampPalette(base)(n)
}


# tradeSeq smoother function

get_smooths_df <- function(sce, genes, nPoints = 200, group_label = NULL) {
  if (length(genes) == 0) {
    return(tibble(gene = character(), pt = numeric(), value = numeric(),
                  group = character()))
  }
  pt_grid <- seq(0, 1, length.out = nPoints)
  mats <- lapply(genes, function(g) {
    sm <- tryCatch(
      as.numeric(predictSmooth(sce, gene = g, nPoints = nPoints, tidy = FALSE)),
      error = function(e) rep(NA_real_, nPoints)
    )
    tibble(gene = g, pt = pt_grid, value = sm)
  })
  bind_rows(mats) %>% mutate(group = group_label %||% "group")
}

# Returns a list with ggplot objects and the vectors of significant genes used

plot_sig_category <- function(
    sce, asso_df, counts_mat, gene_pool, layer_name,
    category_name, alpha = 0.05,
    nPoints = 200, zscore = TRUE,
    out_dir = ".", overlay_file_prefix = NULL,
    genes_per_page = 16L
) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Keep only genes present
  gene_pool <- intersect(gene_pool, rownames(counts_mat))
  
  # Significant genes by FDR
  sig_genes <- asso_df %>%
    filter(padj < alpha, gene %in% gene_pool) %>%
    arrange(padj) %>%
    pull(gene) %>%
    unique()
  
  if (length(sig_genes) == 0) {
    message(sprintf("[%s/%s] No significant genes at FDR < %.2f",
                    layer_name, category_name, alpha))
    return(list(
      overlay = NULL, faceted = NULL, sig_genes = character(0)
    ))
  }
  
  df <- get_smooths_df(sce, sig_genes, nPoints = nPoints, group_label = category_name)
  
  if (zscore) {
    df <- df %>%
      group_by(gene) %>%
      mutate(
        mu = mean(value, na.rm = TRUE),
        sdv = sd(value, na.rm = TRUE),
        z = ifelse(is.finite(sdv) & sdv > 0, (value - mu)/sdv, 0)
      ) %>%
      ungroup()
    yvar <- "z"
    ylab <- "Smoothed expression (z)"
  } else {
    yvar <- "value"
    ylab <- "Smoothed log-expression"
  }
  
  # Colors
  pg_cols <- setNames(interpolate_lancet(length(sig_genes)), sig_genes)
  
  # Overlay
  p_overlay <-
    ggplot(df, aes(.data$pt, .data[[yvar]], color = gene)) +
    geom_line(size = 1.0, alpha = 0.95, na.rm = TRUE) +
    labs(
      title = sprintf("Significant Genes (FDR<%.2f)", layer_name, category_name, alpha),
      x = "Scaled pseudotime (0->1)", y = ylab
    ) +
    theme_pubr() +
    theme(legend.position = "right") +
    scale_color_manual(values = pg_cols, name = category_name)
  
  # Faceted multi-page
  genes <- unique(df$gene)
  chunks <- split(genes, ceiling(seq_along(genes)/genes_per_page))
  
  faceted_pdf <- file.path(
    out_dir,
    sprintf("%s_%s_sig_faceted.pdf", layer_name, category_name)
  )
  grDevices::pdf(faceted_pdf, width = 10, height = 8)
  for (chunk in chunks) {
    dpage <- df %>% filter(gene %in% chunk)
    p_faceted <-
      ggplot(dpage, aes(.data$pt, .data[[yvar]], color = gene)) +
      geom_line(size = 0.9, na.rm = TRUE) +
      facet_wrap(~ gene, scales = "free_y") +
      labs(
        title = sprintf("per-gene smoothed expression (FDR<%.2f)", layer_name, category_name, alpha),
        x = "Scaled pseudotime (0->1)", y = ylab
      ) +
      theme_pubr() +
      theme(legend.position = "none") +
      scale_color_manual(values = pg_cols)
    print(p_faceted)
  }
  grDevices::dev.off()
  
  # Save overlay as SVG/PDF
  if (is.null(overlay_file_prefix)) {
    overlay_file_prefix <- sprintf("%s_%s_sig_overlay", layer_name, category_name)
  }
  overlay_svg <- file.path(out_dir, paste0(overlay_file_prefix, ".svg"))
  overlay_pdf <- file.path(out_dir, paste0(overlay_file_prefix, ".pdf"))
  ggsave(overlay_svg, p_overlay, width = 7, height = 5)
  ggsave(overlay_pdf, p_overlay, width = 7, height = 5)
  
  invisible(list(
    overlay = p_overlay,
    faceted_pdf = faceted_pdf,
    overlay_svg = overlay_svg,
    overlay_pdf = overlay_pdf,
    sig_genes = sig_genes,
    colors = pg_cols
  ))
}


run_layer_plots <- function(
    sce, asso_df, counts_mat, pools, layer_name,
    alpha = 0.05, out_dir = file.path("pseudotime_plots", layer_name),
    nPoints = 200, zscore = TRUE, genes_per_page = 16L
) {
  res <- list()
  for (nm in names(pools)) {
    res[[nm]] <- plot_sig_category(
      sce = sce,
      asso_df = asso_df,
      counts_mat = counts_mat,
      gene_pool = pools[[nm]],
      layer_name = layer_name,
      category_name = nm,
      alpha = alpha,
      nPoints = nPoints,
      zscore = zscore,
      out_dir = out_dir,
      genes_per_page = genes_per_page
    )
  }
  invisible(res)
}

##Endoderm
pg_pool  <- intersect(proteoglycans$Name, rownames(counts_endo))
emt_pool <- intersect(emt_tfs,              rownames(counts_endo))
gag_pool <- intersect(gag$genes,            rownames(counts_endo))

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

# Finally the correlation plots

source("Pseudotime_Helpers_All_Genes.R")

gene_sets <- list(
  emt = emt_tfs,
  pg  = proteoglycans$Name,
  gags = gag$genes
)

gags = gag$genes

grad_endo <- c("#1F78B4", "#32CD32", "#66C2A5")
library(ggpubr)

# Following will create the pairwise module correlation plots as well
res_endo_pg  <- run_endo_proteoglycans()
res_endo_gag  <- run_endo_gags()
res_meso_pg  <- run_meso_proteoglycans()
res_meso_gag  <- run_meso_gags()

cor_dataframe <- list("Endo_Proteo" = res_endo_pg[["cor_pairs_wide"]],
                      "Endo_GAG" = res_endo_gag[["cor_pairs_wide"]],
                      "Meso_Proteo" = res_meso_pg[["cor_pairs_wide"]],
                      "Meso_GAG" = res_meso_gag[["cor_pairs_wide"]])

writexl::write_xlsx(cor_dataframe, "Correlation_with_Pseudotime_Dynamic.xlsx")
