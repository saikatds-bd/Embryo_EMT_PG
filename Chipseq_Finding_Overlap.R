source("Chipoverlap_Helpers.R")
ref_tbl <- "~/Chipseq_Annotation.xlsx"
ref_df <- read_xlsx(ref_tbl)


ref_df$Category <- ref_df$`Reannotated name`


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


print(dplyr::mutate(samples_df, Exists = file.exists(Path)), n = Inf)
stopifnot(all(file.exists(samples_df$Path)))

load("~/Interesting_GAG_Genes.RData")
#Making sure that the files exist
stopifnot(file.exists(ref_tbl), exists("proteoglycans"), exists("gag"))
ref_df$Category <- ref_df$`Reannotated name`

res <- run_overlap_from_df(
   samples_df   = samples_df,
   ref_tbl      = ref_df,
   proteoglycans= proteoglycans,
   gag          = gag,
   parent_dir   = "20251001_Overlap",
   make_plots   = TRUE,
   save_outputs = TRUE
 )
 
