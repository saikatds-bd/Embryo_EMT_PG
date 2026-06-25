# === Load Libraries ===
library(tidyverse)
library(data.table)
library(ggpubr)
library(patchwork)

# === Load and Prepare Data ===
raw_data <- fread("rna_single_cell_type_germ_layer.csv",
                  drop = 1)
raw_data <- as.data.frame(raw_data)

# Use gene symbols
raw_data$Gene <- toupper(raw_data$`Gene name`) 
raw_data$nTPM <- log2(1+raw_data$nTPM)

raw_data <- raw_data %>%
  group_by(Gene, `Cell type`) %>%
  summarise(nTPM = mean(nTPM), .groups = "drop")  # Collapse duplicates

# === Load Gene Sets ===
proteoglycans <- readxl::read_xlsx("Proteoglycans_GAGs_List.xlsx", sheet = "Proteoglycans")

gag <- readxl::read_xlsx("Proteoglycans_GAGs_List.xlsx", sheet = "GAG")

# Ensure gene symbol format
emt_tfs <- toupper(c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2"))
proteoglycan_genes <- toupper(proteoglycans$Gene)
gag_genes <- toupper(gag$Gene)

interesting_genes <- unique(c(proteoglycan_genes, gag_genes, emt_tfs))

# === Build Feature List by Category ===
features_list <- list(
  Proteoglycans = sort(intersect(proteoglycan_genes, raw_data$Gene)),
  GAG_genes = sort(intersect(gag_genes, raw_data$Gene)),
  emt_tfs = emt_tfs
)

# === Create Expression Matrix (genes ?? cell types) ===
expr_mat <- raw_data %>%
  filter(Gene %in% interesting_genes) %>%
  pivot_wider(names_from = `Cell type`, values_from = nTPM) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

expr_mat[is.na(expr_mat)] <- 0

emt_tfs_in_data <- intersect(emt_tfs, rownames(expr_mat))
other_genes <- setdiff(rownames(expr_mat), emt_tfs_in_data)

# === Compute Correlation Matrix ===
make_normal_tissue_emt_correlation_excel <- function(
    expr_mat,
    proteoglycan_genes,
    gag_genes,
    emt_genes = c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2"),
    out_file = "Normal_tissue_EMT_pairwise_correlations.xlsx",
    cor_method = "pearson",
    cor_threshold = 0.4
) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required. Install with install.packages('openxlsx').")
  }
  
  # Make everything uppercase to match data
  rownames(expr_mat) <- toupper(rownames(expr_mat))
  proteoglycan_genes <- toupper(proteoglycan_genes)
  gag_genes <- toupper(gag_genes)
  emt_genes <- toupper(emt_genes)
  
  # Keep only genes present in expression matrix
  proteoglycan_genes <- sort(unique(intersect(proteoglycan_genes, rownames(expr_mat))))
  gag_genes <- sort(unique(intersect(gag_genes, rownames(expr_mat))))
  
  # Keep EMT TFs in preferred biological order
  emt_order <- c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2")
  emt_genes <- intersect(emt_order, emt_genes)
  emt_genes <- emt_genes[emt_genes %in% rownames(expr_mat)]
  
  if (length(emt_genes) == 0) {
    stop("No EMT TFs were found in expr_mat.")
  }
  
  if (length(proteoglycan_genes) == 0 && length(gag_genes) == 0) {
    stop("No proteoglycan or GAG genes were found in expr_mat.")
  }
  
  make_one_sheet <- function(module_genes) {
    
    module_genes <- sort(unique(intersect(module_genes, rownames(expr_mat))))
    
    if (length(module_genes) == 0) {
      return(tibble::tibble(
        module_gene = character(),
        `# positive` = integer(),
        `# negative` = integer()
      ))
    }
    
    out_df <- purrr::map_dfr(module_genes, function(module_gene) {
      
      cors <- purrr::map_dbl(emt_genes, function(tf) {
        
        x <- as.numeric(expr_mat[module_gene, ])
        y <- as.numeric(expr_mat[tf, ])
        
        ok <- is.finite(x) & is.finite(y)
        
        if (sum(ok) < 3) {
          return(NA_real_)
        }
        
        if (sd(x[ok], na.rm = TRUE) == 0 || sd(y[ok], na.rm = TRUE) == 0) {
          return(NA_real_)
        }
        
        suppressWarnings(
          stats::cor(x[ok], y[ok], method = cor_method)
        )
      })
      
      names(cors) <- emt_genes
      
      tibble::tibble(
        module_gene = module_gene
      ) %>%
        dplyr::bind_cols(tibble::as_tibble_row(cors))
    }) %>%
      dplyr::mutate(
        `# positive` = rowSums(
          dplyr::across(
            dplyr::all_of(emt_genes),
            ~ !is.na(.x) & .x >= cor_threshold
          )
        ),
        `# negative` = rowSums(
          dplyr::across(
            dplyr::all_of(emt_genes),
            ~ !is.na(.x) & .x <= -cor_threshold
          )
        )
      ) %>%
      dplyr::select(
        "module_gene",
        dplyr::all_of(emt_genes),
        "# positive",
        "# negative"
      ) %>%
      dplyr::arrange(.data$module_gene)
    
    out_df
  }
  
  proteo_df <- make_one_sheet(proteoglycan_genes)
  gag_df    <- make_one_sheet(gag_genes)
  
  wb <- openxlsx::createWorkbook()
  
  openxlsx::addWorksheet(wb, "Normal_Proteo")
  openxlsx::writeData(wb, "Normal_Proteo", proteo_df)
  openxlsx::freezePane(wb, "Normal_Proteo", firstRow = TRUE)
  openxlsx::setColWidths(wb, "Normal_Proteo", cols = 1:ncol(proteo_df), widths = "auto")
  
  openxlsx::addWorksheet(wb, "Normal_GAG")
  openxlsx::writeData(wb, "Normal_GAG", gag_df)
  openxlsx::freezePane(wb, "Normal_GAG", firstRow = TRUE)
  openxlsx::setColWidths(wb, "Normal_GAG", cols = 1:ncol(gag_df), widths = "auto")
  
  openxlsx::saveWorkbook(wb, out_file, overwrite = TRUE)
  
  message("Saved Excel file: ", out_file)
  
  invisible(list(
    Proteo = proteo_df,
    GAG = gag_df,
    excel_file = out_file
  ))
}


normal_corr_excel <- make_normal_tissue_emt_correlation_excel(
  expr_mat = expr_mat,
  proteoglycan_genes = proteoglycans$Gene,
  gag_genes = gag$Gene,
  emt_genes = emt_tfs,
  out_file = "Normal_tissue_EMT_pairwise_correlations.xlsx",
  cor_method = "pearson",
  cor_threshold = 0.4
)


print_total_pos_neg <- function(
    df,
    tf_cols = c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2"),
    cor_threshold = 0.4,
    label = NULL
) {
  tf_cols <- intersect(tf_cols, colnames(df))
  
  if (length(tf_cols) == 0) {
    stop("None of the EMT TF columns were found in the dataframe.")
  }
  
  vals <- df %>%
    dplyr::select(dplyr::all_of(tf_cols)) %>%
    unlist(use.names = FALSE)
  
  total_positive <- sum(vals >= cor_threshold, na.rm = TRUE)
  total_negative <- sum(vals <= -cor_threshold, na.rm = TRUE)
  total_above_threshold <- total_positive + total_negative
  total_tested <- sum(!is.na(vals))
  
  if (!is.null(label)) {
    cat("\n", label, "\n", sep = "")
    cat(strrep("-", nchar(label)), "\n", sep = "")
  }
  
  cat("Overall counts\n")
  cat("Positive correlations >= ", cor_threshold, ": ", total_positive, "\n", sep = "")
  cat("Negative correlations <= -", cor_threshold, ": ", total_negative, "\n", sep = "")
  cat("Total above threshold: ", total_above_threshold, "\n", sep = "")
  cat("Total tested: ", total_tested, "\n\n", sep = "")
  
  cat("Per EMT TF counts\n")
  
  for (tf in tf_cols) {
    x <- df[[tf]]
    
    tf_positive <- sum(x >= cor_threshold, na.rm = TRUE)
    tf_negative <- sum(x <= -cor_threshold, na.rm = TRUE)
    tf_above_threshold <- tf_positive + tf_negative
    tf_tested <- sum(!is.na(x))
    
    cat(
      tf, 
      " | positive: ", tf_positive,
      " | negative: ", tf_negative,
      " | total above threshold: ", tf_above_threshold,
      " | tested: ", tf_tested,
      "\n",
      sep = ""
    )
  }
  
  invisible(NULL)
}

print_total_pos_neg(
  normal_corr_excel$Proteo,
  cor_threshold = 0.4,
  label = "Normal Proteoglycans"
)

print_total_pos_neg(
  normal_corr_excel$GAG,
  cor_threshold = 0.4,
  label = "Normal GAGs"
)


# For visualization
cor_matrix <- cor(t(expr_mat[other_genes, ]), t(expr_mat[emt_tfs_in_data, ]), method = "pearson")
cor_df <- as.data.frame(as.table(cor_matrix))
colnames(cor_df) <- c("Gene", "EMT_TF", "Correlation")

# === Rank Genes by Max Absolute Correlation ===
cor_df_summary <- cor_df %>%
  group_by(Gene) %>%
  summarise(max_cor = max(abs(Correlation))) %>%
  arrange(desc(max_cor))
emt_tfs_ordered <- c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2")
emt_tfs_in_data <- emt_tfs_ordered[emt_tfs_ordered %in% rownames(expr_mat)]

# Ordered gene categories
proteoglycan_genes_in_data <- sort(intersect(proteoglycan_genes, rownames(expr_mat)))
gag_genes_in_data          <- sort(intersect(gag_genes, rownames(expr_mat)))


# ============================================================
# Normal tissue correlation plots
# Proteoglycans first, then GAGs
# Genes alphabetically ordered
# EMT TFs ordered as SNAI1, SNAI2, TWIST1, TWIST2, ZEB1, ZEB2
# ============================================================

# === Define ordered categories ===
category_gene_list <- list(
  Proteoglycans = proteoglycan_genes_in_data,
  GAGs = gag_genes_in_data
)

category_titles <- c(
  Proteoglycans = "Correlation between proteoglycan core enzymes and EMT TFs",
  GAGs = "Correlation between GAG biosynthetic enzymes and EMT TFs"
)

# === Build expression table for plotting ===
genes_to_keep <- unique(c(
  proteoglycan_genes_in_data,
  gag_genes_in_data,
  emt_tfs_in_data
))

long_expr <- raw_data %>%
  dplyr::filter(.data$Gene %in% genes_to_keep) %>%
  dplyr::select(`Cell type`, Gene, nTPM) %>%
  tidyr::pivot_wider(
    id_cols = `Cell type`,
    names_from = Gene,
    values_from = nTPM
  )

# ============================================================
# Plot function
# ============================================================

make_gene_tf_plots <- function(
    gene,
    tf_list,
    expr_long_df,
    base_family = "Times New Roman",
    text_size = 8
) {
  plots <- list()
  
  for (tf in tf_list) {
    
    # If either gene or TF is missing, add empty spacer
    if (!(gene %in% colnames(expr_long_df)) || !(tf %in% colnames(expr_long_df))) {
      plots[[length(plots) + 1]] <- patchwork::plot_spacer()
      next
    }
    
    plot_df <- expr_long_df %>%
      dplyr::select(
        x = dplyr::all_of(tf),
        y = dplyr::all_of(gene)
      ) %>%
      dplyr::filter(
        is.finite(.data$x),
        is.finite(.data$y)
      )
    
    # Only plot when enough cell types are available
    if (nrow(plot_df) > 3 &&
        stats::sd(plot_df$x, na.rm = TRUE) > 0 &&
        stats::sd(plot_df$y, na.rm = TRUE) > 0) {
      
      p <- ggplot(plot_df, aes(x = x, y = y)) +
        geom_point(
          size = 1.4,
          alpha = 0.8,
          color = "#1F78B4"
        ) +
        geom_smooth(
          method = "lm",
          se = TRUE,
          color = "#FF0000",
          linewidth = 0.7
        ) +
        ggpubr::stat_cor(
          method = "pearson",
          cor.coef.name = "r",
          size = 2.5,
          family = base_family,
          label.x.npc = "left",
          label.y.npc = "top"
        ) +
        labs(
          x = tf,
          y = gene,
          title = paste0(gene, " vs ", tf)
        ) +
        ggpubr::theme_pubr(
          base_family = base_family,
          base_size = text_size
        ) +
        theme(
          text = element_text(
            family = base_family,
            size = text_size
          ),
          plot.title = element_text(
            family = base_family,
            size = text_size,
            face = "bold",
            hjust = 0.5
          ),
          axis.title = element_text(
            family = base_family,
            size = text_size
          ),
          axis.text = element_text(
            family = base_family,
            size = text_size
          ),
          axis.line = element_line(
            linewidth = 0.4,
            colour = "black"
          ),
          plot.margin = margin(3, 3, 3, 3)
        )
      
    } else {
      p <- patchwork::plot_spacer()
    }
    
    plots[[length(plots) + 1]] <- p
  }
  
  plots
}


# ============================================================
# Export PDF
# ============================================================

out_pdf <- "Normal_Tissue_Log_Correlation_Plots_Proteo_GAG_A4.pdf"

cairo_pdf(
  filename = out_pdf,
  height = 11.69,   # A4 landscape width in inches
  width = 8.27,   # A4 landscape height in inches
  family = "Times New Roman"
)

for (category in names(category_gene_list)) {
  
  genes_in_cat <- category_gene_list[[category]]
  genes_in_cat <- sort(genes_in_cat)
  
  if (length(genes_in_cat) == 0) {
    next
  }
  
  # 3 genes per page, 6 EMT TF plots per gene
  gene_chunks <- split(
    genes_in_cat,
    ceiling(seq_along(genes_in_cat) / 6)
  )
  
  for (page_i in seq_along(gene_chunks)) {
    
    chunk <- gene_chunks[[page_i]]
    
    all_plots <- list()
    
    for (gene in chunk) {
      
      plots <- make_gene_tf_plots(
        gene = gene,
        tf_list = emt_tfs_in_data,
        expr_long_df = long_expr,
        base_family = "Times New Roman",
        text_size = 8
      )
      
      all_plots <- append(all_plots, plots)
    }
    
    # Pad to 36 panels: 6 rows x 6 columns
    while (length(all_plots) < 36) {
      all_plots[[length(all_plots) + 1]] <- patchwork::plot_spacer()
    }
    
    p_page <- patchwork::wrap_plots(
      all_plots,
      ncol = 6,
      nrow = 6
    )
    
    # Category title only on first page of each category
    if (page_i == 1) {
      p_page <- p_page +
        patchwork::plot_annotation(
          title = category_titles[[category]],
          theme = theme(
            plot.title = element_text(
              family = "Times New Roman",
              size = 9,
              face = "bold",
              hjust = 0
            )
          )
        )
    }
    
    print(p_page)
  }
}

dev.off()




