library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(janitor)
library(tibble)
library(writexl) 
library(purrr)


overlap_parent <- "Overlap"

files <- list.files(
  path = overlap_parent,
  pattern = "_Chipseq_Overlap\\.xlsx$",
  full.names = TRUE,
  recursive = TRUE
)

print(files)
length(files)


meta <- tibble(file = files) %>%
  mutate(
    base = basename(file),
    cell_line = str_match(base, "^(.*)_(SNAI1|SNAI2|ZEB1|ZEB2)_Chipseq_Overlap\\.xlsx$")[,2],
    TF        = str_match(base, "^(.*)_(SNAI1|SNAI2|ZEB1|ZEB2)_Chipseq_Overlap\\.xlsx$")[,3]
  )
print(meta)


out_list <- vector("list", nrow(meta))

for (i in seq_len(nrow(meta))) {
  f  <- meta$file[i]
  cl <- meta$cell_line[i]
  tf <- meta$TF[i]
  message("Reading: ", basename(f))
  
  raw <- read_xlsx(f)
  
  message("  Columns: ", paste(names(raw), collapse = ", "))
  
  # lower-case column names for matching
  ln <- tolower(names(raw))
  
  
  pick_first <- function(cands) {
    idx <- which(ln %in% cands)
    if (length(idx)) idx[1] else NA_integer_
  }
  
  i_gene       <- pick_first(c("gene","symbol","genes","gene_symbol"))
  i_feature    <- pick_first(c("feature","consensus_label","feat","type","label"))
  i_gh_type    <- pick_first(c("gh_type","type","class","gh_class"))
  i_region_uid <- pick_first(c("region_uid","uid","region_id"))
  i_chr        <- pick_first(c("chr","chrom","chromosome","seqnames"))
  i_start      <- pick_first(c("start","chromstart","start_coordinate","start_bp"))
  i_end        <- pick_first(c("end","chromend","end_coordinate","end_bp"))
  
  dat <- tibble(
    gene       = if (!is.na(i_gene))       raw[[i_gene]]       else NA_character_,
    feature    = if (!is.na(i_feature))    raw[[i_feature]]    else NA_character_,
    gh_type    = if (!is.na(i_gh_type))    raw[[i_gh_type]]    else NA_character_,
    region_uid = if (!is.na(i_region_uid)) raw[[i_region_uid]] else NA_character_,
    chr        = if (!is.na(i_chr))        raw[[i_chr]]        else NA_character_,
    start      = if (!is.na(i_start))      raw[[i_start]]      else NA_real_,
    end        = if (!is.na(i_end))        raw[[i_end]]        else NA_real_
  )
  
  # Normalize types
  dat <- dat %>%
    mutate(
      gene       = as.character(gene),
      gh_type    = as.character(gh_type),
      region_uid = as.character(region_uid),
      chr        = as.character(chr),
      start      = suppressWarnings(as.numeric(start)),
      end        = suppressWarnings(as.numeric(end))
    )
  
  # Fill missing region_uid from coords if available
  miss_uid <- is.na(dat$region_uid) | dat$region_uid == ""
  have_coords <- !is.na(dat$chr) & !is.na(dat$start) & !is.na(dat$end)
  if (any(miss_uid & have_coords)) {
    dat$region_uid[miss_uid & have_coords] <- paste0(dat$chr[miss_uid & have_coords], ":", dat$start[miss_uid & have_coords], "-", dat$end[miss_uid & have_coords])
  }
  still_miss <- is.na(dat$region_uid) | dat$region_uid == ""
  if (any(still_miss)) {
    dat$region_uid[still_miss] <- paste0("uid_", which(still_miss))
  }
  
 
  lx <- tolower(as.character(dat$feature))
  dat$feature <- dplyr::case_when(
    #str_detect(lx, "both|pe|promoter\\+enhancer") ~ "Both/PE",
    str_detect(lx, "promoter|tss")                ~ "Promoter",
    str_detect(lx, "enhancer")                    ~ "Enhancer",
    TRUE                                          ~ "Unclassified"
  )
  
  
  table(dat$feature, useNA = "ifany")
  
  
  
  dat$cell_line <- cl
  dat$TF        <- tf
  
  
  message("  Feature counts:"); print(table(dat$feature, useNA = "ifany"))
  
  out_list[[i]] <- dat
}

# Combine all files

hits_all <- bind_rows(out_list)


message("ALL features across all files:")
print(table(hits_all$feature, useNA = "ifany"))
message("By TF:")
print(hits_all %>% count(TF, feature) %>% arrange(TF, feature))

# Clean gene names
hits_all <- hits_all %>%
  mutate(gene = str_trim(gene), gene = na_if(gene, "")) %>%
  filter(!is.na(gene))

per_gene_cellline <- hits_all %>%
  group_by(gene, cell_line, TF) %>%
  summarise(
    Promoter       = n_distinct(region_uid[feature == "Promoter"]),
    Enhancer       = n_distinct(region_uid[feature == "Enhancer"]),
    n_total_hits   = n_distinct(region_uid[feature %in% c("Promoter","Enhancer","Both/PE","Unclassified")]),
    gh_type_unique = { u <- unique(na.omit(gh_type)); if (length(u)==0) NA_character_ else if (length(u)==1) u else "mixed" },
    .groups = "drop"
  ) %>%
  mutate(
    regulated_by_TF = n_total_hits > 0,
    Type = gh_type_unique
  ) %>%
  select(-gh_type_unique) %>%
  arrange(TF, gene, cell_line)

# checkpoint
per_gene_cellline %>% count(TF)  # rows per TF


# CHECKPOINT
per_gene_cellline %>%
  group_by(TF) %>%
  summarise(
    n_rows           = dplyr::n(),
    n_with_promoter  = sum(Promoter   > 0),
    n_with_enhancer  = sum(Enhancer   > 0),
    #n_with_both_pe   = sum(Both_PE    > 0),
    n_with_any       = sum(n_total_hits > 0),
    .groups = "drop"
  ) %>% print(n = 50)


presence_matrix <- per_gene_cellline %>%
  mutate(present = as.integer(n_total_hits > 0)) %>%
  select(TF, gene, cell_line, present) %>%
  pivot_wider(names_from = cell_line, values_from = present, values_fill = 0)


feature_across_lines <- per_gene_cellline %>%
  mutate(
    has_promoter     = Promoter   > 0,
    has_enhancer     = Enhancer   > 0,
    #has_both_pe      = Both_PE    > 0,
    #has_unclassified = Unclassified > 0,
    has_any          = n_total_hits > 0
  ) %>%
  group_by(TF, gene) %>%
  summarise(
    lines_with_any          = sum(has_any, na.rm = TRUE),
    lines_with_promoter     = sum(has_promoter, na.rm = TRUE),
    lines_with_enhancer     = sum(has_enhancer, na.rm = TRUE),
   # lines_with_both_pe      = sum(has_both_pe, na.rm = TRUE),
  #  lines_with_unclassified = sum(has_unclassified, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(TF, desc(lines_with_any),
         # desc(lines_with_both_pe),
          desc(lines_with_promoter + lines_with_enhancer))

# OPTIONAL CHECKS
feature_across_lines %>% filter(gene=="DCN") %>% print(n=50)
feature_across_lines %>% filter(grepl("MM001", TF) | TRUE) %>% head()


dir.create("TF_Summaries", showWarnings = FALSE)

for (tf in sort(unique(hits_all$TF))) {
  sub_hits   <- hits_all %>% filter(TF == tf) %>% arrange(gene, cell_line, feature)
  sub_pgcl   <- per_gene_cellline %>% filter(TF == tf)
  sub_pres   <- presence_matrix %>% filter(TF == tf) %>% select(-TF)
  sub_feat   <- feature_across_lines %>% filter(TF == tf) %>% select(-TF)
  
  out_file <- file.path("TF_Summaries", paste0("Summary_",tf, "_Chipseq_PG_GAG.xlsx"))
  writexl::write_xlsx(
    list(
      hits_long            = sub_hits,
      per_gene_cellline    = sub_pgcl,
      presence_matrix      = sub_pres,
      feature_across_lines = sub_feat
    ),
    path = out_file
  )
  message("[OK] Wrote: ", out_file)
}

table(hits_all$feature)

feature_across_lines <- feature_across_lines %>%
  mutate(category = dplyr::case_when(
    #lines_with_both_pe      > 0 ~ "Both/PE",
    lines_with_promoter     > 0 ~ "Promoter",
    lines_with_enhancer     > 0 ~ "Enhancer",
    TRUE                        ~ "None"
  ))

