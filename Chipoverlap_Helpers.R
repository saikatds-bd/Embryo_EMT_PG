suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
  library(readr)
  library(readxl)
  library(writexl)
  library(ggplot2)
  library(ggpubr)
  library(stringr)
})


`%||%` <- function(a, b) if (!is.null(a)) a else b
norm_sym  <- function(x) str_squish(toupper(as.character(x)))
norm_syms <- function(x) unique(norm_sym(x))

pick_col <- function(df, choices, required = TRUE) {
  nm <- intersect(choices, names(df))
  if (!length(nm) && required)
    stop("Missing one of columns: ", paste(choices, collapse = ", "))
  if (!length(nm)) return(NULL)
  nm[1]
}

read_table_auto <- function(x) {
  if (is.data.frame(x)) return(x)
  stopifnot(is.character(x), length(x) == 1)
  ext <- tolower(tools::file_ext(x))
  if (ext %in% c("xlsx","xls")) {
    readxl::read_excel(x)
  } else if (ext %in% c("csv")) {
    readr::read_csv(x, show_col_types = FALSE)
  } else {
    readr::read_tsv(x, show_col_types = FALSE)
  }
}

# Path hygiene
.clean_path <- function(x) {
  p <- as.character(x)[1]
  p <- trimws(p)
  p <- sub("^file://", "", p)
  p <- gsub("\r|\n", "", p)
  p <- iconv(p, from = "", to = "UTF-8")
  path.expand(p)
}
.file_exists <- function(p) {
  p0 <- .clean_path(p)
  if (file.exists(p0)) return(TRUE)
  suppressWarnings({
    pn <- try(normalizePath(p0, winslash = "/", mustWork = FALSE), silent = TRUE)
  })
  isTRUE(file.exists(pn))
}
.dir_glimpse <- function(p) {
  d <- dirname(.clean_path(p))
  if (!dir.exists(d)) return(sprintf("Parent dir does NOT exist: %s", d))
  kids <- tryCatch(list.files(d, all.files = TRUE, no.. = TRUE), error = function(e) character(0))
  paste0("Parent dir: ", d, "\nContains: ", paste(head(kids, 50), collapse = ", "))
}

# =========================
# Peaks: narrowPeak/BED -> GRanges
# =========================
parse_narrowpeak <- function(path) {
  p <- .clean_path(path)
  if (!.file_exists(p)) {
    stop(sprintf("Peak file not found:\n  %s\ngetwd(): %s\n%s",
                 p, getwd(), .dir_glimpse(p)))
  }
  p_abs <- tryCatch(normalizePath(p, winslash = "/", mustWork = TRUE), error = function(e) p)
  
  df <- readr::read_tsv(
    p_abs, col_names = FALSE, comment = "#",
    col_types = cols(.default = col_character()), progress = FALSE
  )
  if (nrow(df) == 0) stop("Peak file has 0 rows: ", p_abs)
  if (ncol(df) < 3)  stop("Need b	%3 columns (chrom,start,end) in: ", p_abs)
  
  chrom  <- df[[1]]
  start0 <- suppressWarnings(as.numeric(df[[2]]))   # 0-based
  end1   <- suppressWarnings(as.numeric(df[[3]]))   # 1-based
  if (anyNA(start0) || anyNA(end1)) stop("Non-numeric start/end in: ", p_abs)
  
  has_strand <- ncol(df) >= 6 && all(df[[6]] %in% c("+","-",".","*"), na.rm = TRUE)
  strand <- if (has_strand) { x <- df[[6]]; x[is.na(x) | x == "."] <- "*"; x } else rep("*", nrow(df))
  
  gr <- GRanges(seqnames = chrom,
                ranges   = IRanges(start = start0 + 1L, end = end1),
                strand   = strand)
  
  mcols(gr)$name  <- if (ncol(df) >= 4) df[[4]] else paste0("peak_", seq_len(nrow(df)))
  mcols(gr)$score <- if (ncol(df) >= 5) suppressWarnings(as.integer(df[[5]])) else NA_integer_
  
  off <- if (has_strand) 1L else 0L
  mcols(gr)$signalValue <- if (ncol(df) >= 6 + off) suppressWarnings(as.numeric(df[[6 + off]])) else NA_real_
  mcols(gr)$pValue      <- if (ncol(df) >= 7 + off) suppressWarnings(as.numeric(df[[7 + off]])) else NA_real_
  mcols(gr)$qValue      <- if (ncol(df) >= 8 + off) suppressWarnings(as.numeric(df[[8 + off]])) else NA_real_
  mcols(gr)$peak        <- if (ncol(df) >= 9 + off) suppressWarnings(as.integer(df[[9 + off]])) else NA_integer_
  gr
}

# =========================
# Reference -> GRanges with standardized Label
# =========================
ref_to_granges_safe <- function(ref_tbl) {
  ref <- read_table_auto(ref_tbl)
  
  chr_col    <- pick_col(ref, c("Chromosome","chrom","chr","seqnames"))
  region_col <- pick_col(ref, c("Region","region"), required = FALSE)
  
  if (!is.null(region_col)) {
    region_parsed <- stringr::str_split_fixed(ref[[region_col]], "\\.\\.", 2)
    start_num <- as.numeric(gsub("[^0-9]", "", region_parsed[,1]))
    end_num   <- as.numeric(gsub("[^0-9]", "", region_parsed[,2]))
  } else {
    start_col <- pick_col(ref, c("start","Start","tss_start","TSS_start","prom_start","Promoter_start"))
    end_col   <- pick_col(ref, c("end","End","tss_end","TSS_end","prom_end","Promoter_end"))
    start_num <- as.numeric(ref[[start_col]])
    end_num   <- as.numeric(ref[[end_col]])
  }
  
  seqs <- as.character(ref[[chr_col]])
  if (!any(grepl("^chr", seqs))) seqs <- paste0("chr", seqs)
  
  gene_col <- pick_col(ref, c("Target","Target_upper","connected_gene","Connected gene",
                              "connected_genes","Gene","genes","Genes"))
  gh_col   <- pick_col(ref, c("genehancer_id","GeneHancer_ID","genehancer","GH_ID"), required = FALSE)
  type_col <- pick_col(ref, c("Type","Category","class"), required = FALSE)
  
  gr <- GRanges(seqnames = seqs, ranges = IRanges(start = start_num, end = end_num))
  
  # coalesce label columns in preference order
  lab_cols <- intersect(c("Category","Consensus_Label","Final_Label","Manual_Consensus","Type"),
                        names(ref))
  coalesce_chars <- function(vecs) {
    out <- as.character(vecs[[1]]); out[!nzchar(out)] <- NA_character_
    if (length(vecs) > 1) {
      for (i in 2:length(vecs)) {
        v <- as.character(vecs[[i]]); v[!nzchar(v)] <- NA_character_
        idx <- is.na(out); out[idx] <- v[idx]
      }
    }
    out
  }
  labels_raw <- if (length(lab_cols)) coalesce_chars(lapply(lab_cols, function(nm) ref[[nm]])) else rep(NA_character_, nrow(ref))
  labels0 <- tolower(labels_raw)
  
  Label <- dplyr::case_when(
    grepl("both", labels0)              ~ "Both",
    grepl("prom|tss|prox", labels0)     ~ "Promoter",
    grepl("enh|distal", labels0)        ~ "Enhancer",
    TRUE                                ~ NA_character_
  )
  
  mcols(gr)$Label         <- Label
  mcols(gr)$genes_raw     <- as.character(ref[[gene_col]])
  mcols(gr)$genehancer_id <- if (!is.null(gh_col)) as.character(ref[[gh_col]]) else NA_character_
  mcols(gr)$Type          <- if (!is.null(type_col)) as.character(ref[[type_col]]) else NA_character_
  mcols(gr)$region_id     <- seq_len(nrow(ref))
  mcols(gr)$Category_raw  <- labels_raw
  
  keep_idx <- !is.na(mcols(gr)$Label)
  if (!any(keep_idx)) {
    bad <- sort(table(labels_raw), decreasing = TRUE)
    stop("ref_to_granges_safe: 0 rows passed label filter. Top raw labels: ",
         paste(head(sprintf("%s=%d", names(bad), bad), 10), collapse="; "),
         ". Adjust mapping.")
  }
  gr[keep_idx]
}

# =========================
# Harmonize seqlevels (prune to common in-use)
# =========================
prune_to_common_seqlevels <- function(peaks, ref_gr) {
  # normalize style robustly
  suppressWarnings({
    seqlevelsStyle(peaks)  <- "UCSC"
    seqlevelsStyle(ref_gr) <- "UCSC"
  })
  common <- union(
    intersect(seqlevelsInUse(peaks),  seqlevels(ref_gr)),
    intersect(seqlevelsInUse(ref_gr), seqlevels(peaks))
  )
  if (!length(common)) stop("No common seqlevels between peaks and reference. Check genome build/naming.")
  peaks2  <- keepSeqlevels(peaks,  common, pruning.mode = "coarse")
  ref2    <- keepSeqlevels(ref_gr, common, pruning.mode = "coarse")
  list(peaks = peaks2, ref = ref2)
}

# =========================
# Core analysis (single sample)
# =========================
analyze_chip_pg_gag <- function(
    peaks_file,
    ref_tbl,
    proteoglycans,   # df with column 'Name' OR character vector
    gag,             # list/df with 'genes' OR character vector
    sample_id    = "Sample",
    tf_name      = "TF",
    out_dir      = NULL,
    make_plots   = TRUE,
    save_outputs = TRUE
) {
  peaks  <- parse_narrowpeak(peaks_file)
  ref_gr <- ref_to_granges_safe(ref_tbl)
  
  # prune seqlevels to common to avoid Seqinfo warnings
  hz <- prune_to_common_seqlevels(peaks, ref_gr)
  peaks  <- hz$peaks
  ref_gr <- hz$ref
  
  # Target gene sets
  pg_vec  <- if (is.data.frame(proteoglycans) && "Name" %in% names(proteoglycans)) proteoglycans$Name else as.character(proteoglycans)
  gag_vec <- if (is.list(gag) && "genes" %in% names(gag)) gag$genes else as.character(gag)
  target_genes <- norm_syms(c(pg_vec, gag_vec))
  
  # Overlaps
  hits <- findOverlaps(peaks, ref_gr, ignore.strand = TRUE)
  
  # Expand overlaps (use normalized Label carried in ref_gr)
  expand_hits <- function(subject_idx, gr_obj) {
    tibble(
      region_id       = mcols(gr_obj)$region_id[subject_idx],
      genes_raw       = mcols(gr_obj)$genes_raw[subject_idx],
      Label           = mcols(gr_obj)$Label[subject_idx],         # Promoter/Enhancer/Both
      chr             = as.character(seqnames(gr_obj))[subject_idx],
      start           = start(gr_obj)[subject_idx],
      end             = end(gr_obj)[subject_idx],
      genehancer_id   = mcols(gr_obj)$genehancer_id[subject_idx],
      Type            = mcols(gr_obj)$Type[subject_idx]
    ) %>%
      tidyr::separate_rows(genes_raw, sep = "[,;|]") %>%
      tidyr::separate_rows(genehancer_id, sep = "[,;|]") %>%     # <-- NEW: split GH IDs if multiple
      mutate(
        genehancer_id = stringr::str_squish(genehancer_id),
        gene          = stringr::str_squish(genes_raw)
      ) %>%
      filter(nzchar(gene)) %>%
      distinct(region_id, chr, start, end, Label, Type, gene, genehancer_id, .keep_all = TRUE)
  }
  
  ol_df <- if (length(hits)) expand_hits(subjectHits(hits), ref_gr) else
    tibble(region_id = integer(), genes_raw = character(), Label = character(),
           chr = character(), start = integer(), end = integer(),
           genehancer_id = character(), Type = character(), gene = character())
  
  # standardized feature
  ol_df <- ol_df %>%
    mutate(feature = factor(Label, levels = c("Enhancer","Promoter","Both")))
  
  # Targets only (element-wise normalization!)
  ol_targets <- ol_df %>%
    mutate(gene_norm = norm_sym(gene)) %>%
    filter(gene_norm %in% target_genes)
  
  # Per-gene summary
  per_gene <- ol_targets %>%
    count(gene = gene_norm, feature, name = "hits") %>%
    tidyr::pivot_wider(names_from = feature, values_from = hits, values_fill = 0)
  for (cn in c("Promoter","Enhancer","Both")) if (!cn %in% names(per_gene)) per_gene[[cn]] <- 0L
  per_gene <- per_gene %>%
    mutate(
      n_total_hits = Promoter + Enhancer + Both,
      !!paste0("regulated_by_", tf_name) := n_total_hits > 0
    ) %>%
    arrange(desc(Promoter > 0), desc(n_total_hits), gene)
  
  # Summary counts
  n_prom <- if ("Promoter" %in% names(per_gene)) sum(per_gene$Promoter > 0) else 0L
  n_enh  <- if ("Enhancer" %in% names(per_gene)) sum(per_gene$Enhancer > 0) else 0L
  n_both <- if ("Both"     %in% names(per_gene)) sum(per_gene$Both     > 0) else 0L
  n_any  <- if (nrow(per_gene)) sum(per_gene[[paste0("regulated_by_", tf_name)]]) else 0L
  summary_tbl <- tibble(
    Metric = c("Called peaks", "Target genes (PG+GAG)",
               "Promoter-bound targets", "Enhancer-bound targets", "Both-labeled targets",
               "Any bound (promoter/enhancer/both)"),
    N = c(length(peaks), length(target_genes), n_prom, n_enh, n_both, n_any)
  )
  
  # Per-region summary
  per_region <- ol_targets %>%
    group_by(region_id, chr, start, end, feature) %>%
    summarise(
      genes   = paste0(sort(unique(gene)), collapse = ", "),
      n_genes = n_distinct(gene),
      .groups = "drop"
    ) %>%
    arrange(feature, chr, start) %>%
    mutate(region_uid = paste0(chr, ":", start, "-", end, "|", feature))
  
  # Output dir
  if (is.null(out_dir)) out_dir <- paste0(tf_name, "_", sample_id, "_bound_outputs")
  if (save_outputs) dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save tables & BEDs
  if (save_outputs) {
    ol_targets_export <- ol_targets %>%
      transmute(gene = gene_norm, chr, start, end, feature, region_id,
                Label, genehancer_id, Type)
    
    readr::write_csv(ol_targets_export, file.path(out_dir, sprintf("%s_%s_overlap_PG_GAG_prom_enh.csv", sample_id, tf_name)))
    readr::write_csv(per_gene,          file.path(out_dir, sprintf("%s_%s_per_gene_hits.csv",   sample_id, tf_name)))
    readr::write_csv(summary_tbl,       file.path(out_dir, sprintf("%s_%s_summary_counts.csv",  sample_id, tf_name)))
    readr::write_csv(per_region,        file.path(out_dir, sprintf("%s_%s_overlapped_regions_with_ids.csv", sample_id, tf_name)))
    
    to_bed <- function(df) tibble(
      chrom = df$chr,
      start = as.integer(df$start) - 1L,  # BED is 0-based
      end   = as.integer(df$end),
      name  = paste0(df$feature, "_", df$region_id),
      score = 0L, strand = "."
    )
    unique_regions <- per_region %>% distinct(region_id, .keep_all = TRUE)
    readr::write_tsv(to_bed(unique_regions),
                     file.path(out_dir, sprintf("%s_%s_overlapped_regions.bed", sample_id, tf_name)),
                     col_names = FALSE)
    readr::write_tsv(to_bed(unique_regions %>% filter(feature == "Promoter")),
                     file.path(out_dir, sprintf("%s_%s_overlapped_promoters.bed", sample_id, tf_name)),
                     col_names = FALSE)
    readr::write_tsv(to_bed(unique_regions %>% filter(feature == "Enhancer")),
                     file.path(out_dir, sprintf("%s_%s_overlapped_enhancers.bed", sample_id, tf_name)),
                     col_names = FALSE)
    
    writexl::write_xlsx(
      list(
        overlaps   = ol_targets_export,
        per_gene   = per_gene,
        per_region = per_region,
        summary    = summary_tbl
      ),
      path = file.path(out_dir, sprintf("%s_%s_Chipseq_Overlap.xlsx", sample_id, tf_name))
    )
    
    save(ol_targets, per_gene, per_region, summary_tbl, ref_gr, peaks,
         file = file.path(out_dir, sprintf("%s_%s_Chipseq.RData", sample_id, tf_name)))
  }
  
  # Plots
  plots <- list()
  if (make_plots) {
    category_colors <- c("Promoter"="#d73027","Enhancer"="#4575b4","Both"="#fdae61")
    
    sum_tbl <- tibble(Label = c("Promoter","Enhancer","Both"),
                      N = c(n_prom, n_enh, n_both)) %>%
      filter(N > 0) %>%
      mutate(Label = factor(Label, levels = c("Enhancer","Promoter","Both")))
    
    if (nrow(sum_tbl)) {
      plots$bar <- ggpubr::ggbarplot(sum_tbl, x = "Label", y = "N",
                                     fill = "Label") +
        scale_fill_manual(values = category_colors[levels(sum_tbl$Label)]) +
        labs(title = sprintf("%s binding at PG/GAG regulatory regions (%s)", tf_name, sample_id),
             x = NULL, y = "# genes")
    }
    
    make_dot <- function(target_set, title) {
      if (!nrow(per_gene)) return(NULL)
      keep <- norm_syms(target_set)
      df <- per_gene %>%
        filter(gene %in% keep) %>%
        select(gene, Promoter, Enhancer, Both) %>%
        tidyr::pivot_longer(c(Promoter, Enhancer, Both),
                            names_to = "Feature", values_to = "hits") %>%
        mutate(present = hits > 0)
      if (!nrow(df)) return(NULL)
      ggplot(df, aes(Feature, gene)) +
        geom_point(aes(size = pmin(hits, 5), alpha = present, shape = Feature, colour = Feature)) +
        scale_alpha_manual(values = c(`FALSE`=0.15, `TRUE`=1)) +
        scale_color_manual(values = category_colors) +
        theme_pubr(base_size = 9) +
        theme(axis.text.y = element_text(size = 6), legend.position = "right") +
        labs(title = title, x = NULL, y = NULL)
    }
    
    plots$pg  <- make_dot(if (is.data.frame(proteoglycans)) proteoglycans$Name else proteoglycans,
                          sprintf("Proteoglycans %s hits - %s", tf_name, sample_id))
    plots$gag <- make_dot(if (is.list(gag)) gag$genes else gag,
                          sprintf("GAG enzymes %s hits - %s", tf_name, sample_id))
    
    if (save_outputs && length(plots)) {
      pdf(file.path(out_dir, sprintf("%s_%s_Chipseq_Plots.pdf", sample_id, tf_name)),
          width = 6.5, height = 6)
      lapply(plots, function(p) if (!is.null(p)) print(p))
      dev.off()
    }
  }
  
  list(
    peaks = peaks,
    ref_gr = ref_gr,
    overlaps_all = ol_df,
    overlaps_targets = ol_targets,
    per_gene = per_gene,
    per_region = per_region,
    summary = summary_tbl,
    plots = plots,
    out_dir = out_dir
  )
}

# =========================
# Choose a peak file from a path (file or dir)
# =========================
resolve_peak_file <- function(path, cell_line = NULL) {
  p <- .clean_path(path)
  if (file.exists(p) && !dir.exists(p)) {
    return(tryCatch(normalizePath(p, winslash = "/", mustWork = TRUE), error = function(e) p))
  }
  if (!dir.exists(p)) stop(sprintf("Path does not exist (neither file nor dir): %s", p))
  cand <- list.files(p, pattern = "\\.(narrowPeak|bed)(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  if (!length(cand)) stop("No narrowPeak/bed files found in: ", p)
  if (!is.null(cell_line) && nzchar(cell_line)) {
    cl <- gsub("[^A-Za-z0-9]", ".*", cell_line)
    patt1 <- paste0("(?i)", cl, ".*_peaks\\.narrowPeak(\\.gz)?$")
    idx <- grep(patt1, basename(cand), perl = TRUE)
    if (!length(idx)) {
      patt2 <- paste0("(?i)", cl, ".*\\.(narrowPeak|bed)(\\.gz)?$")
      idx <- grep(patt2, basename(cand), perl = TRUE)
    }
    if (length(idx)) cand <- cand[idx]
  }
  if (length(cand) > 1) {
    info <- file.info(cand)
    cand <- cand[which.max(info$size)]
  }
  tryCatch(normalizePath(cand[1], winslash = "/", mustWork = TRUE), error = function(e) cand[1])
}

# =========================
# Batch runner with logging
# =========================
run_overlap_from_df <- function(
    samples_df,               # must have TF, Cell_Line, Path
    ref_tbl,
    proteoglycans,
    gag,
    parent_dir    = "Overlap",
    make_plots    = TRUE,
    save_outputs  = TRUE,
    write_run_log = TRUE
) {
  stopifnot(all(c("TF","Cell_Line","Path") %in% names(samples_df)))
  dir.create(parent_dir, showWarnings = FALSE, recursive = TRUE)
  
  n <- nrow(samples_df)
  results <- vector("list", n)
  run_log <- tibble::tibble(
    i = seq_len(n),
    TF = as.character(samples_df$TF),
    Cell_Line = as.character(samples_df$Cell_Line),
    Path = as.character(samples_df$Path),
    peaks_file = NA_character_,
    out_dir    = NA_character_,
    ok         = FALSE,
    error      = NA_character_
  )
  
  for (i in seq_len(n)) {
    tf  <- run_log$TF[i]
    cl  <- run_log$Cell_Line[i]
    pth <- run_log$Path[i]
    message(sprintf("[ %d/%d ] %s / %s", i, n, tf, cl))
    
    out_dir_i <- file.path(parent_dir, tf, cl)
    dir.create(out_dir_i, showWarnings = FALSE, recursive = TRUE)
    run_log$out_dir[i] <- out_dir_i
    
    peaks_file_i <- tryCatch(
      resolve_peak_file(pth, cell_line = cl),
      error = function(e) { run_log$error[i] <- conditionMessage(e); NA_character_ }
    )
    message("  b peaks: ", peaks_file_i)
    message("  b out:   ", out_dir_i)
    run_log$peaks_file[i] <- peaks_file_i
    
    res_i <- tryCatch({
      if (is.na(peaks_file_i) || !.file_exists(peaks_file_i))
        stop("Cannot resolve peaks file for ", tf, "/", cl)
      analyze_chip_pg_gag(
        peaks_file   = peaks_file_i,
        ref_tbl      = ref_tbl,
        proteoglycans= proteoglycans,
        gag          = gag,
        sample_id    = cl,
        tf_name      = tf,
        out_dir      = out_dir_i,
        make_plots   = make_plots,
        save_outputs = save_outputs
      )
    }, error = function(e) {
      run_log$error[i] <- conditionMessage(e)
      warning(sprintf("Failed %s/%s: %s", tf, cl, e$message), call. = FALSE)
      NULL
    })
    
    results[i] <- list(res_i)
    run_log$ok[i] <- !is.null(res_i)
  }
  
  names(results) <- paste(samples_df$TF, samples_df$Cell_Line, sep = "__")
  if (write_run_log) {
    readr::write_tsv(run_log, file.path(parent_dir, "batch_run_log.tsv"))
  }
  attr(results, "run_log") <- run_log
  results
}

# =========================
# Diagnostic (optional)
# =========================
check_one_overlap_pair <- function(peaks_path, ref_tbl, label = basename(peaks_path)) {
  peaks  <- parse_narrowpeak(peaks_path)
  ref_gr <- ref_to_granges_safe(ref_tbl)
  suppressWarnings({
    seqlevelsStyle(peaks)  <- "UCSC"
    seqlevelsStyle(ref_gr) <- "UCSC"
  })
  cat("\n=== ", label, " ===\n", sep = "")
  cat("peaks: used=", length(seqlevelsInUse(peaks)),  " total=", length(seqlevels(peaks)),  "\n", sep="")
  cat("ref  : used=", length(seqlevelsInUse(ref_gr)), " total=", length(seqlevels(ref_gr)), "\n", sep="")
  only_in_peaks <- setdiff(seqlevelsInUse(peaks),  seqlevels(ref_gr))
  only_in_ref   <- setdiff(seqlevelsInUse(ref_gr), seqlevels(peaks))
  cat("In peaks but not ref (in-use): ", ifelse(length(only_in_peaks), paste(only_in_peaks, collapse=", "), "<none>"), "\n", sep="")
  cat("In ref but not peaks (in-use): ", ifelse(length(only_in_ref),   paste(only_in_ref,   collapse=", "), "<none>"), "\n", sep="")
  if (length(only_in_peaks)) {
    cnt <- sort(table(as.character(seqnames(peaks))), decreasing = TRUE)
    cat("Counts in mismatched peaks levels:\n")
    print(cnt[names(cnt) %in% only_in_peaks])
  }
  invisible(list(only_in_peaks=only_in_peaks, only_in_ref=only_in_ref))
}





# =========================
# Batch QC summary for overlap results
# =========================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(ggplot2)
})

# Summarize a single result (one TF/Cell_Line)
summarize_one <- function(name, obj) {
  # name is like "SNAI1__HS578T"
  if (is.null(obj)) return(NULL)
  parts <- strsplit(name, "__", fixed = TRUE)[[1]]
  TF <- parts[1]; Cell_Line <- parts[2] %||% NA_character_
  
  # pull counts from summary table if present; else compute from per_gene
  if (!is.null(obj$summary)) {
    s <- obj$summary
    getN <- function(metric) {
      v <- s$N[match(metric, s$Metric)]
      ifelse(is.na(v), NA_integer_, as.integer(v))
    }
    called_peaks   <- getN("Called peaks")
    n_targets      <- getN("Target genes (PG+GAG)")
    n_promoter     <- getN("Promoter-bound targets")
    n_enhancer     <- getN("Enhancer-bound targets")
    n_both         <- getN("Both-labeled targets")
    n_any          <- getN("Any bound (promoter/enhancer/both)")
  } else {
    pg <- obj$per_gene %||% tibble()
    if (nrow(pg)) {
      n_promoter <- sum(pg$Promoter > 0, na.rm = TRUE)
      n_enhancer <- sum(pg$Enhancer > 0, na.rm = TRUE)
      n_both     <- if ("Both" %in% names(pg)) sum(pg$Both > 0, na.rm = TRUE) else 0L
      n_any      <- sum((pg$Promoter + pg$Enhancer + (pg$Both %||% 0L)) > 0, na.rm = TRUE)
    } else {
      n_promoter <- n_enhancer <- n_both <- n_any <- 0L
    }
    called_peaks <- length(obj$peaks %||% GRanges())
    n_targets    <- length(unique(obj$overlaps_targets$gene %||% character()))
  }
  
  tibble(
    TF, Cell_Line,
    called_peaks, n_targets,
    n_promoter, n_enhancer, n_both, n_any,
    frac_any = ifelse(n_targets > 0, n_any / n_targets, NA_real_)
  )
}

# Summarize entire batch result list returned by run_overlap_from_df()
summarize_batch <- function(res, write_csv = NULL) {
  stopifnot(is.list(res))
  nm <- names(res)
  out <- map_dfr(seq_along(res), ~ summarize_one(nm[.x], res[[.x]]))
  out <- out %>%
    arrange(desc(frac_any), desc(n_any)) %>%
    mutate(
      TF = factor(TF, levels = unique(TF)),
      Cell_Line = factor(Cell_Line, levels = unique(Cell_Line))
    )
  if (!is.null(write_csv)) readr::write_csv(out, write_csv)
  out
}

# -------- optional quick plots (bar + heatmap-ish tile) --------
plot_qc_bar <- function(qc_df, top_n = 30) {
  df <- qc_df %>%
    mutate(label = paste(TF, Cell_Line, sep="__")) %>%
    slice_head(n = min(top_n, n()))
  ggplot(df, aes(x = reorder(label, frac_any), y = frac_any)) +
    geom_col() +
    coord_flip() +
    labs(x = NULL, y = "Frac. targets with any binding", title = "Top samples by bound fraction")
}

plot_qc_tile <- function(qc_df) {
  ggplot(qc_df, aes(TF, Cell_Line, fill = frac_any)) +
    geom_tile() +
    scale_fill_viridis_c(na.value = "grey90") +
    labs(x = "TF", y = "Cell line", fill = "Frac any", title = "Binding fraction by TF / Cell line") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# -------- optional: top bound genes per sample --------
top_bound_genes <- function(res, sample_key, n = 20) {
  obj <- res[[sample_key]]
  stopifnot(!is.null(obj))
  pg <- obj$per_gene
  pg %>%
    mutate(total_hits = (Promoter %||% 0L) + (Enhancer %||% 0L) + (Both %||% 0L)) %>%
    arrange(desc(total_hits), desc(Promoter), gene) %>%
    slice_head(n = n)
}




