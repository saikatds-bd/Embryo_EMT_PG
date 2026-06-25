library(dplyr)
library(purrr)
library(stringr)
library(writexl)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("intersect", "base")
conflict_prefer("reduce", "purrr")
# RData files containing DE objects are kept in separate folders according to their TFs
load_rdata_objects_local <- function(files) {
  
  obj_list <- list()
  
  for (f in files) {
    
    obj_name <- tools::file_path_sans_ext(basename(f))
    obj_name <- make.names(obj_name)
    
    tmp_env <- new.env()
    loaded_names <- load(f, envir = tmp_env)
    
    obj_list[[obj_name]] <- as.list(tmp_env)
    
    cat("Loaded:", basename(f), "as", obj_name, "\n")
    cat("  Contains:", paste(loaded_names, collapse = ", "), "\n")
  }
  
  return(obj_list)
}


extract_logfc_from_loaded_objects <- function(obj_list, pattern = "_DE_") {
  
  de_objects <- names(obj_list)[str_detect(names(obj_list), pattern)]
  
  cat("\nFound DE objects:\n")
  print(de_objects)
  
  extract_logfc <- function(obj_name) {
    
    obj <- obj_list[[obj_name]]
    
    for (el_name in names(obj)) {
      
      el <- obj[[el_name]]
      
      if (is.data.frame(el) && all(c("Symbol", "logFC") %in% colnames(el))) {
        
        cat("Found Symbol + logFC in:", obj_name, "->", el_name, "\n")
        
        return(
          el %>%
            dplyr::select(Symbol, logFC) %>%
            dplyr::rename(
              !!paste0(obj_name, "_logFC") := logFC
            )
        )
      }
    }
    
    warning(paste("Skipping", obj_name, "- no valid Symbol + logFC table found"))
    return(NULL)
  }
  
  logfc_list <- purrr::map(de_objects, extract_logfc)
  logfc_list <- logfc_list[!sapply(logfc_list, is.null)]
  
  return(logfc_list)
}


prepare_merged_logfc_simple <- function(logfc_list, genes_keep = NULL) {
  
  merged <- purrr::reduce(logfc_list, full_join, by = "Symbol")
  
  if (!is.null(genes_keep)) {
    merged <- merged %>%
      filter(Symbol %in% genes_keep)
  }
  
  return(merged)
}


# Column name cleaning


clean_dataset_colnames <- function(x) {
  
  x <- str_remove(x, "_logFC$")
  x <- str_replace_all(x, "_DE_", "_")
  x <- str_replace_all(x, "_", " ")
  
  # Fix common cell-line / dataset conventions
  x <- str_replace_all(x, regex("^Hacat", ignore_case = TRUE), "HaCaT")
  x <- str_replace_all(x, regex("^KLM1", ignore_case = TRUE), "KLM1")
  x <- str_replace_all(x, regex("^Keratinocytes", ignore_case = TRUE), "Keratinocytes")
  x <- str_replace_all(x, regex("^MCF10A", ignore_case = TRUE), "MCF10A")
  x <- str_replace_all(x, regex("^Pancreas", ignore_case = TRUE), "Pancreatic Cancer Spheroids")
  x <- str_replace_all(x, regex("^gastric", ignore_case = TRUE), "Gastric")
  x <- str_replace_all(x, regex("^hmec htert 56164", ignore_case = TRUE), "HMEC-hTERT 56164")
  x <- str_replace_all(x, regex("^hmec htert", ignore_case = TRUE), "HMEC-hTERT")
  x <- str_replace_all(x, regex("^huvec", ignore_case = TRUE), "HUVEC")
  x <- str_replace_all(x, regex("^neuroblastoma", ignore_case = TRUE), "Neuroblastoma")
  
  x <- str_replace_all(x, regex("^MDA MB 231", ignore_case = TRUE), "MDA-MB-231")
  x <- str_replace_all(x, regex("^MDA 231 D", ignore_case = TRUE), "MDA-231-D")
  x <- str_replace_all(x, regex("^SW480", ignore_case = TRUE), "SW480")
  x <- str_replace_all(x, regex("^MiaPaCa2", ignore_case = TRUE), "MiaPaCa2")
  x <- str_replace_all(x, regex("^RKO1", ignore_case = TRUE), "RKO1")
  x <- str_replace_all(x, regex("^LS147T", ignore_case = TRUE), "LS147T")
  x <- str_replace_all(x, regex("^SUM159", ignore_case = TRUE), "SUM159")
  x <- str_replace_all(x, regex("^HCC827", ignore_case = TRUE), "HCC827")
  x <- str_replace_all(x, regex("^HCEnC 21T", ignore_case = TRUE), "HCEnC-21T")
  x <- str_replace_all(x, regex("^C 09\\.10", ignore_case = TRUE), "C-09.10")
  x <- str_replace_all(x, regex("^ESC", ignore_case = TRUE), "ESC")
  x <- str_replace_all(x, regex("^H9", ignore_case = TRUE), "H9")
  x <- str_replace_all(x, regex("^CD4 Th1", ignore_case = TRUE), "CD4+ Th1")
  x <- str_replace_all(x, regex("^PEO1", ignore_case = TRUE), "PEO1")
  x <- str_replace_all(x, regex("^OVCA420", ignore_case = TRUE), "OVCA420")
  x <- str_replace_all(x, regex("^HGC27", ignore_case = TRUE), "HGC27")
  x <- str_replace_all(x, regex("^SGC7901", ignore_case = TRUE), "SGC7901")
  x <- str_replace_all(x, regex("^SH 5Y5Y", ignore_case = TRUE), "SH-SY5Y")
  x <- str_replace_all(x, regex("^U87", ignore_case = TRUE), "U87")
  x <- str_replace_all(x, regex("^PSC", ignore_case = TRUE), "PSC")
  x <- str_replace_all(x, regex("^GBM ECs", ignore_case = TRUE), "GBM ECs")
  x <- str_replace_all(x, regex("^SK N Be2C", ignore_case = TRUE), "SK-N-Be2C")
  x <- str_replace_all(x, regex("^MSC", ignore_case = TRUE), "MSC")
  x <- str_replace_all(x, regex("^MCF7", ignore_case = TRUE), "MCF7")
  x <- str_replace_all(x, regex("^Fibroblast", ignore_case = TRUE), "Fibroblast")
  x <- str_replace_all(x, regex("^THP1", ignore_case = TRUE), "THP1")
  x <- str_replace_all(x, regex("^SMS CTR", ignore_case = TRUE), "SMS-CTR")
  
  # Keep perturbation label readable
  x <- str_replace_all(x, " OE$", " OE")
  x <- str_replace_all(x, " KD$", " KD")
  x <- str_replace_all(x, " KO$", " KO")
  
  return(x)
}


make_converted_logfc <- function(df) {
  
  converted <- df
  
  kd_ko_cols <- colnames(converted)[
    str_detect(colnames(converted), "_KD_logFC$|_KO_logFC$")
  ]
  
  converted <- converted %>%
    mutate(
      across(
        all_of(kd_ko_cols),
        ~ .x * -1
      )
    )
  
  return(converted)
}


clean_logfc_colnames <- function(df) {
  
  old_cols <- colnames(df)
  new_cols <- old_cols
  
  non_symbol <- old_cols != "Symbol"
  new_cols[non_symbol] <- clean_dataset_colnames(old_cols[non_symbol])
  
  colnames(df) <- new_cols
  
  return(df)
}


process_one_tf <- function(tf_name, tf_path, genes_keep, out_dir = NULL) {
  
  cat("\n============================================\n")
  cat("Processing:", tf_name, "\n")
  cat("Path:", tf_path, "\n")
  cat("============================================\n")
  
  rdata_files <- list.files(
    path = tf_path,
    pattern = "\\.RData$",
    full.names = TRUE
  )
  
  if (length(rdata_files) == 0) {
    warning(paste("No .RData files found for", tf_name))
    return(NULL)
  }
  
  loaded_objects <- load_rdata_objects_local(rdata_files)
  
  logfc_list <- extract_logfc_from_loaded_objects(
    obj_list = loaded_objects,
    pattern = "_DE_"
  )
  
  if (length(logfc_list) == 0) {
    warning(paste("No valid logFC tables found for", tf_name))
    return(NULL)
  }
  
  regular_raw <- prepare_merged_logfc_simple(
    logfc_list = logfc_list,
    genes_keep = genes_keep
  )
  
  converted_raw <- make_converted_logfc(regular_raw)
  
  regular_clean <- clean_logfc_colnames(regular_raw)
  converted_clean <- clean_logfc_colnames(converted_raw)
  
  if (is.null(out_dir)) {
    out_dir <- tf_path
  }
  
  out_file <- file.path(
    out_dir,
    paste0(format(Sys.Date(), "%Y%m%d"), "_", tf_name, "_Interesting_Genes_logFC.xlsx")
  )
  
  sheet_list <- list(
    regular_clean,
    converted_clean
  )
  
  names(sheet_list) <- c(
    paste0(tf_name, "_DE_Regular"),
    paste0(tf_name, "_DE_Converted")
  )
  
  writexl::write_xlsx(
    x = sheet_list,
    path = out_file
  )
  
  cat("\nSaved:", out_file, "\n")
  
  return(
    list(
      regular = regular_clean,
      converted = converted_clean,
      output_file = out_file
    )
  )
}



load("C:/Users/ssa214/UiT Office 365/O365-Serglycin - General/Core protein asociation towards EMT/Saikat/Comments from Cathy/Pseudo_Files/Interesting_GAG_Genes.RData")

emt_tfs <- c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2")

proteoglycans <- readxl::read_xlsx("C:/Users/ssa214/UiT Office 365/O365-Serglycin - General/Core protein asociation towards EMT/Saikat/Comments from Cathy/Github_Scripts/Proteoglycans_GAGs_List.xlsx", sheet = "Proteoglycans")
gag <- readxl::read_xlsx("C:/Users/ssa214/UiT Office 365/O365-Serglycin - General/Core protein asociation towards EMT/Saikat/Comments from Cathy/Github_Scripts/Proteoglycans_GAGs_List.xlsx", sheet = "GAG")

genes2 <- unique(c(
  proteoglycans$Gene,
  gag$Gene,
  emt_tfs
))

tf_paths <- c(
  SNAI1  = "C:/Users/ssa214/UiT Office 365/O365-Serglycin - General/Core protein asociation towards EMT/Saikat/TF-RNA/SNAIL1/All_SNAI_DE",
  SNAI2  = "C:/Users/ssa214/UiT Office 365/O365-Serglycin - General/Core protein asociation towards EMT/Saikat/TF-RNA/SNAI2/SNAI2_DE/Subsetted",
  TWIST1 = "C:/Users/ssa214/UiT Office 365/O365-Serglycin - General/Core protein asociation towards EMT/Saikat/TF-RNA/TWIST1/All_DE",
  TWIST2 = "C:/Users/ssa214/UiT Office 365/O365-Serglycin - General/Core protein asociation towards EMT/Saikat/TF-RNA/TWIST2/All_De/Subsetted",
  ZEB1   = "C:/Users/ssa214/UiT Office 365/O365-Serglycin - General/Core protein asociation towards EMT/Saikat/TF-RNA/ZEB1/All_Zeb1_DE",
  ZEB2   = "C:/Users/ssa214/UiT Office 365/O365-Serglycin - General/Core protein asociation towards EMT/Saikat/TF-RNA/ZEB2/All_DE/Subsetted"
)

# ============================================================
# 6. Run all TFs
# ============================================================
tf_out_dir <- "C:/Users/ssa214/UiT Office 365/O365-Serglycin - General/Core protein asociation towards EMT/Saikat/Comments from Cathy/TF_RNA_Excel"

dir.create(tf_out_dir, recursive = TRUE, showWarnings = FALSE)
all_tf_results <- imap(
  tf_paths,
  ~ process_one_tf(
    tf_name = .y,
    tf_path = .x,
    genes_keep = genes2,
    out_dir = tf_out_dir
  )
)
