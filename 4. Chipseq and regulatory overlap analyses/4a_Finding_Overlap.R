library(tidyverse)
source("Chipoverlap_Helpers.R")


ref_tbl <- "Chipseq_Annotation.xlsx"
ref_df <- read_xlsx(ref_tbl)

samples_df <- readxl::read_xlsx("All_Samples_Overlap.xlsx")
samples_df <- samples_df[complete.cases(samples_df$Path),]
samples_df <- samples_df %>%
  dplyr::mutate(
    Path = vapply(
      Path,
      function(p) normalizePath(.clean_path(p), winslash = "/", mustWork = TRUE),
      character(1)
    )
  )

# Optional: verify existence explicitly
print(dplyr::mutate(samples_df, Exists = file.exists(Path)), n = Inf)
stopifnot(all(file.exists(samples_df$Path)))


load("PG_GAG_Genes.RData")
stopifnot(file.exists(ref_tbl), exists("proteoglycans"), exists("gag"))

res <- run_overlap_from_df(
   samples_df   = samples_df,
   ref_tbl      = ref_df,
   proteoglycans= proteoglycans,
   gag          = gag,
   parent_dir   = "Overlap",
   make_plots   = TRUE,
   save_outputs = TRUE
 )


