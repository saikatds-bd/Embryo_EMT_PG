
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(writexl)
  library(ggplot2)
  library(ggpubr)
  library(scales)
})



base_dir <- "Overlap"

gene_list_rdata <- "PG_GAG_Genes.RData"

ref_tbl_path <- "Chipseq_Annotation.xlsx"

tf_pattern <- "(SNAI1|SNAI2|ZEB1|ZEB2)"

tf_order <- c("SNAI1", "SNAI2", "ZEB1", "ZEB2")

PROP_ANY  <- 0.50
PROP_SAME <- 0.50

out_xlsx <- "20260529_Detailed_Chipseq_Result.xlsx"

out_pg_plot  <- "20260529_Proteoglycans_Strict_Regulation_BubblePlot.pdf"
out_gag_plot <- "20260529_GAGs_Strict_Regulation_BubblePlot.pdf"


pick_first_idx <- function(nm_low, cands) {
  idx <- which(nm_low %in% cands)
  if (length(idx)) idx[1] else NA_integer_
}



read_overlap_sheet <- function(fp, tf_pattern) {
  base <- basename(fp)
  
  m <- stringr::str_match(
    base,
    paste0("^(.*)_", tf_pattern, "_Chipseq_Overlap\\.xlsx$")
  )
  
  if (any(is.na(m))) {
    stop("Bad filename: ", base)
  }
  
  cell_line <- m[, 2]
  TF        <- m[, 3]
  
  raw <- read_xlsx(fp, sheet = "overlaps")
  
  if (!nrow(raw)) {
    return(
      tibble(
        TF = TF,
        cell_line = cell_line,
        gene = character(),
        genehancer_id = character()
      )
    )
  }
  
  ln <- tolower(names(raw))
  
  i_gene <- pick_first_idx(
    ln,
    c("gene", "symbol", "genes", "gene_symbol")
  )
  
  i_ghid <- pick_first_idx(
    ln,
    c("genehancer_id", "gh_id", "genehancer")
  )
  
  stopifnot(!is.na(i_gene), !is.na(i_ghid))
  
  tibble(
    TF = TF,
    cell_line = cell_line,
    gene = str_squish(as.character(raw[[i_gene]])),
    genehancer_id = str_squish(as.character(raw[[i_ghid]]))
  ) %>%
    filter(!is.na(gene), nzchar(gene)) %>%
    separate_rows(genehancer_id, sep = "[,;|]") %>%
    mutate(genehancer_id = na_if(str_squish(genehancer_id), "")) %>%
    distinct()
}


# Map comma-separated GeneHancer IDs to comma-separated annotation names.
map_ids_to_names <- function(id_string, ref_tbl) {
  if (is.na(id_string) || id_string == "") {
    return(NA_character_)
  }
  
  ids <- str_split(id_string, ",\\s*")[[1]]
  
  name_tbl <- ref_tbl %>%
    filter(genehancer_id %in% ids) %>%
    select(genehancer_id, `Reannotated name`)
  
  mapped <- vapply(
    ids,
    function(id) {
      nm <- name_tbl$`Reannotated name`[name_tbl$genehancer_id == id]
      if (length(nm) == 0) NA_character_ else nm[1]
    },
    character(1)
  )
  
  paste(mapped, collapse = ", ")
}


# Map comma-separated GeneHancer IDs to comma-separated GeneHancer types.
ids_to_types <- function(id_string, annot_tbl) {
  if (is.na(id_string) || id_string == "") {
    return(NA_character_)
  }
  
  ids <- str_split(id_string, ",\\s*")[[1]]
  
  mapped <- vapply(
    ids,
    function(id) {
      t <- annot_tbl$Type[annot_tbl$genehancer_id == id]
      if (length(t) == 0) NA_character_ else t[1]
    },
    character(1)
  )
  
  paste(mapped, collapse = ", ")
}


# Collapse a comma-separated type string into one plotting category.
types_to_category <- function(types_str) {
  if (is.na(types_str) || types_str == "") {
    return("Unknown")
  }
  
  u <- unique(str_trim(unlist(str_split(types_str, ",\\s*"))))
  u <- u[!is.na(u)]
  
  if (!length(u)) {
    return("Unknown")
  }
  
  u_low <- tolower(u)
  
  only_prom <- all(u_low == "promoter")
  only_enh  <- all(u_low == "enhancer")
  
  if (only_prom) {
    "Promoter"
  } else if (only_enh) {
    "Enhancer"
  } else {
    "Mixed"
  }
}


# Make strict summary table.
# Strict means:
#   Regulated_in_Above_50pct_Lines == TRUE
#   Same_GH_in_Above_50pct_Lines   == TRUE
mk_strict_summary <- function(df, genes, tf_order) {
  dat <- df %>%
    filter(Gene %in% genes) %>%
    mutate(.tf_rank = match(TF, tf_order)) %>%
    arrange(.tf_rank, TF)
  
  reg_by <- dat %>%
    group_by(Gene) %>%
    summarise(
      Regulated_By = {
        tfs <- unique(TF)
        tfs <- intersect(tf_order, tfs)
        
        if (length(tfs) == 0) {
          NA_character_
        } else {
          paste(tfs, collapse = " + ")
        }
      },
      .groups = "drop"
    )
  
  tf_to_ids <- dat %>%
    group_by(Gene) %>%
    summarise(
      Consensus_GH_IDs_By_TF = paste(
        sprintf("%s: %s", TF, Consensus_GeneHancer_IDs),
        collapse = " | "
      ),
      .groups = "drop"
    )
  
  tf_to_types <- dat %>%
    group_by(Gene) %>%
    summarise(
      Consensus_GH_Types_By_TF = paste(
        sprintf("%s: %s", TF, Consensus_GH_Types),
        collapse = " | "
      ),
      .groups = "drop"
    )
  
  reg_by %>%
    left_join(tf_to_ids, by = "Gene") %>%
    left_join(tf_to_types, by = "Gene") %>%
    arrange(Gene)
}


# Make bubble plot from strict bubble data.
make_bubble_plot <- function(plot_df, title_text, subtitle_text) {
  category_colors <- c(
    "Promoter" = "#d73027",
    "Enhancer" = "#4575b4",
    "Mixed"    = "#9467bd",
    "Unknown"  = "grey80"
  )
  
  ggplot(plot_df, aes(x = TF, y = Gene)) +
    geom_point(
      aes(
        size = Fraction_with_ConsensusGH,
        fill = Type_Category
      ),
      shape = 21,
      color = "black",
      alpha = 0.9
    ) +
    scale_fill_manual(
      values = category_colors,
      drop = FALSE
    ) +
    scale_size_continuous(
      name = "Consensus lines / assayed",
      limits = c(0.50, 1.00),
      breaks = c(0.50, 0.75, 1.00),
      labels = percent_format(accuracy = 1),
      range = c(3, 10)
    ) +
    labs(
      x = "EMT TF",
      y = "Gene",
      title = title_text,
      subtitle = subtitle_text
    ) +
    theme_pubr(base_size = 11) +
    theme(
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      )
    )
}


files_overlap <- list.files(
  path = base_dir,
  pattern = paste0("_", tf_pattern, "_Chipseq_Overlap\\.xlsx$"),
  full.names = TRUE,
  recursive = TRUE
)

if (!length(files_overlap)) {
  stop("No *_Chipseq_Overlap.xlsx files found in: ", base_dir)
}


# This records how many cell lines were assayed per TF.

assay_design <- tibble(file = basename(files_overlap)) %>%
  mutate(
    TF = stringr::str_match(
      file,
      paste0("_", tf_pattern, "_Chipseq_Overlap\\.xlsx$")
    )[, 2],
    cell_line = stringr::str_match(
      file,
      paste0("^(.*)_", tf_pattern, "_Chipseq_Overlap\\.xlsx$")
    )[, 2]
  ) %>%
  distinct(TF, cell_line) %>%
  filter(!is.na(TF), !is.na(cell_line))

lines_per_TF <- assay_design %>%
  count(TF, name = "n_lines_assayed")

stopifnot(nrow(lines_per_TF) > 0)



hits_all <- purrr::map_dfr(
  files_overlap,
  read_overlap_sheet,
  tf_pattern = tf_pattern
)


gh_sets <- hits_all %>%
  group_by(TF, gene, cell_line) %>%
  summarise(
    gh_ids = list(sort(unique(na.omit(genehancer_id)))),
    .groups = "drop"
  )

per_line <- gh_sets %>%
  mutate(has_any = lengths(gh_ids) > 0L) %>%
  select(TF, gene, cell_line, has_any)

per_gh <- gh_sets %>%
  unnest_longer(gh_ids, values_to = "gh_id") %>%
  filter(!is.na(gh_id)) %>%
  distinct(TF, gene, cell_line, gh_id) %>%
  count(TF, gene, gh_id, name = "n_lines_with_gh")



sum_any <- per_line %>%
  group_by(TF, gene) %>%
  summarise(
    n_lines_detected = sum(has_any),
    .groups = "drop"
  )

sum_gh <- per_gh %>%
  group_by(TF, gene) %>%
  summarise(
    max_gh_lines = if (n() == 0L) {
      0L
    } else {
      max(n_lines_with_gh)
    },
    consensus_gh_ids = if (n() == 0L) {
      NA_character_
    } else {
      ids <- gh_id[n_lines_with_gh == max(n_lines_with_gh)]
      paste(sort(unique(ids)), collapse = ", ")
    },
    .groups = "drop"
  )

summary_TF_gene <- sum_any %>%
  full_join(sum_gh, by = c("TF", "gene")) %>%
  mutate(
    n_lines_detected = coalesce(as.integer(n_lines_detected), 0L),
    max_gh_lines = coalesce(as.integer(max_gh_lines), 0L)
  ) %>%
  left_join(lines_per_TF, by = "TF") %>%
  mutate(
    n_lines_assayed = coalesce(as.integer(n_lines_assayed), 0L),
    
    prop_any = ifelse(
      n_lines_assayed > 0,
      n_lines_detected / n_lines_assayed,
      NA_real_
    ),
    
    consensus_gh_prop = ifelse(
      n_lines_assayed > 0,
      max_gh_lines / n_lines_assayed,
      NA_real_
    ),
    
    regulated_any_50_raw = is.finite(prop_any) &
      prop_any >= PROP_ANY,
    
    regulated_sameGH_50_raw = is.finite(consensus_gh_prop) &
      consensus_gh_prop >= PROP_SAME,
    
    regulated_any_50 = regulated_any_50_raw,
    
    regulated_sameGH_50 = regulated_any_50_raw &
      regulated_sameGH_50_raw
  ) %>%
  arrange(
    TF,
    desc(regulated_sameGH_50),
    desc(consensus_gh_prop),
    desc(prop_any),
    gene
  )


# ============================================================
# Sanity checks
# ============================================================

viol <- summary_TF_gene %>%
  filter(max_gh_lines > n_lines_detected)

if (nrow(viol)) {
  message(
    "WARNING: ",
    nrow(viol),
    " rows where max_gh_lines > n_lines_detected. Inspect object: viol"
  )
}

audit <- summary_TF_gene %>%
  filter(regulated_sameGH_50 & !regulated_any_50)

if (nrow(audit)) {
  message(
    "WARNING: ",
    nrow(audit),
    " rows where sameGH_50 TRUE but any_50 FALSE. Inspect object: audit"
  )
}


export_df <- summary_TF_gene %>%
  select(
    TF,
    gene,
    n_lines_detected,
    max_gh_lines,
    consensus_gh_ids,
    n_lines_assayed,
    prop_any,
    consensus_gh_prop,
    regulated_any_50,
    regulated_sameGH_50
  )



load("PG_GAG_Genes.RData")




export_df <- export_df %>%
  mutate(
    Category = case_when(
      gene %in% proteoglycans$Name ~ "Proteoglycan",
      gene %in% gag$genes ~ "GAG",
      TRUE ~ "GAG"
    )
  )


export_df <- export_df %>%
  dplyr::rename(
    Gene = gene,
    Cell_Lines_with_Any_Hit = n_lines_detected,
    Cell_Lines_with_ConsensusGeneHancer = max_gh_lines,
    Consensus_GeneHancer_IDs = consensus_gh_ids,
    Total_Cell_Lines = n_lines_assayed,
    Fraction_with_Any_Hit = prop_any,
    Fraction_with_ConsensusGH = consensus_gh_prop,
    Regulated_in_Above_50pct_Lines = regulated_any_50,
    Same_GH_in_Above_50pct_Lines = regulated_sameGH_50
  )



ref_df <- read_xlsx(ref_tbl_path)
ref_df$`Reannotated name` <- ref_df$Category
ref_df_names <- ref_df %>%
  select(genehancer_id, `Reannotated name`) %>%
  distinct(genehancer_id, .keep_all = TRUE)




export_df <- export_df %>%
  mutate(
    Type_of_Genehancer_ID = vapply(
      Consensus_GeneHancer_IDs,
      map_ids_to_names,
      character(1),
      ref_tbl = ref_df_names
    )
  )



pg_df <- export_df %>%
  filter(Gene %in% proteoglycans$Name)

gag_df <- export_df %>%
  filter(Gene %in% gag$genes)




type_col <- c("Reannotated name")[
  c("Reannotated name") %in% names(ref_df)
][1]

if (is.na(type_col)) {
  stop("ref_df must contain either 'Category' or 'Reannotated name'.")
}

gh_annot <- ref_df %>%
  select(genehancer_id, Type = all_of(type_col)) %>%
  distinct(genehancer_id, .keep_all = TRUE)



strict_df <- export_df %>%
  filter(
    Regulated_in_Above_50pct_Lines,
    Same_GH_in_Above_50pct_Lines
  ) %>%
  mutate(
    Consensus_GH_Types = map_chr(
      Consensus_GeneHancer_IDs,
      ~ ids_to_types(.x, gh_annot)
    )
  )


pg_genes <- unique(strict_df$Gene[strict_df$Gene %in% proteoglycans$Name])
gag_genes <- unique(strict_df$Gene[strict_df$Gene %in% gag$genes])

pg_summary_strict <- mk_strict_summary(
  df = strict_df,
  genes = pg_genes,
  tf_order = tf_order
)

gag_summary_strict <- mk_strict_summary(
  df = strict_df,
  genes = gag_genes,
  tf_order = tf_order
)


write_xlsx(
  list(
    
    Proteoglycans_Filtered = pg_summary_strict,
    GAG_Filtered = gag_summary_strict,
    Proteoglycans_Detailed = pg_df,
    GAG_Detailed = gag_df
    
  ),
  out_xlsx
)


bubble_df <- strict_df %>%
  mutate(
    TF = factor(TF, levels = tf_order),
    Type_Category = vapply(
      Consensus_GH_Types,
      types_to_category,
      character(1)
    )
  ) %>%
  mutate(
    Gene = factor(Gene, levels = rev(sort(unique(Gene))))
  ) %>%
  select(
    Gene,
    TF,
    Fraction_with_ConsensusGH,
    Consensus_GeneHancer_IDs,
    Consensus_GH_Types,
    Type_Category
  )



pg_bubble_df <- bubble_df %>%
  filter(Gene %in% pg_summary_strict$Gene)

pg_bubble <- make_bubble_plot(
  plot_df = pg_bubble_df,
  title_text = "Proteoglycans (Any > 50% AND Same-GH > 50%)",
  subtitle_text = "Bubble size = fraction of assayed lines using consensus GH; color = consensus GH type"
)

pg_bubble

ggsave(
  filename = out_pg_plot,
  plot = pg_bubble,
  width = 7.5,
  height = 6,
  dpi = 300
)



gag_bubble_df <- bubble_df %>%
  filter(Gene %in% gag_summary_strict$Gene)

gag_bubble <- make_bubble_plot(
  plot_df = gag_bubble_df,
  title_text = "GAGs (Any > 50% AND Same-GH > 50%)",
  subtitle_text = "Bubble size = fraction of assayed lines using consensus GH; color = consensus GH type"
)

gag_bubble

ggsave(
  filename = out_gag_plot,
  plot = gag_bubble,
  width = 7.5,
  height = 6,
  dpi = 300
)
