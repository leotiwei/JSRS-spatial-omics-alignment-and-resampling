#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(openxlsx)
})

project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
combined_dir <- file.path(project_dir, "resampled")
metadata_file <- file.path(project_dir, "metadata", "rna_bin_cluster_metadata.tsv")
outdir <- file.path(project_dir, "cluster_average")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sections <- strsplit(Sys.getenv("SECTIONS", unset = "cs12_3,cs12_4,cs12_5"), ",")[[1]]
methods <- list(
  conservative = "{section}_combined_multiomics_conservative.csv",
  nearest = "{section}_combined_multiomics_nearest.csv",
  IDW_k4_p2 = "{section}_combined_multiomics_IDW_k4_p2.csv"
)

make_name <- function(pattern, section) gsub("\\{section\\}", section, pattern)

detect_bin_col <- function(dt) {
  candidates <- c("bin_id", "_index", "barcode", "Barcode")
  hit <- intersect(candidates, colnames(dt))
  if (length(hit) == 0) stop("No bin identifier column found in resampled table.")
  hit[1]
}

get_feature_cols <- function(dt) {
  non_features <- c("bin_id", "_index", "barcode", "Barcode", "x", "y", "section_id", "cluster")
  setdiff(colnames(dt), non_features)
}

metadata <- fread(metadata_file) %>%
  mutate(bin_id = as.character(bin_id), section_id = as.character(section_id), cluster = as.character(cluster)) %>%
  filter(!is.na(cluster), cluster != "")

for (method_name in names(methods)) {
  message("Processing method: ", method_name)
  all_sections <- list()
  for (sec in sections) {
    file <- file.path(combined_dir, make_name(methods[[method_name]], sec))
    if (!file.exists(file)) {
      warning("Missing file: ", file)
      next
    }
    dt <- fread(file)
    bin_col <- detect_bin_col(dt)
    dt <- dt %>% mutate(bin_id = as.character(.data[[bin_col]]), section_id = sec)
    all_sections[[sec]] <- dt
  }
  if (length(all_sections) == 0) next

  combined <- bind_rows(all_sections) %>%
    left_join(metadata, by = c("bin_id", "section_id")) %>%
    filter(!is.na(cluster), cluster != "")
  feature_cols <- get_feature_cols(combined)

  cluster_average <- combined %>%
    group_by(cluster) %>%
    summarise(across(all_of(feature_cols), ~ mean(as.numeric(.x), na.rm = TRUE)), .groups = "drop")

  feature_cluster <- as.data.frame(t(as.matrix(cluster_average[, feature_cols, drop = FALSE])))
  colnames(feature_cluster) <- cluster_average$cluster
  feature_cluster <- data.frame(feature = rownames(feature_cluster), feature_cluster, check.names = FALSE)
  rownames(feature_cluster) <- NULL

  wb <- createWorkbook()
  addWorksheet(wb, "Metabolite_average")
  writeData(wb, "Metabolite_average", feature_cluster)
  addWorksheet(wb, "Cluster_summary")
  writeData(wb, "Cluster_summary", as.data.frame(table(combined$cluster)))
  saveWorkbook(wb, file.path(outdir, paste0("Metabolite_cluster_average_", method_name, ".xlsx")), overwrite = TRUE)
}
