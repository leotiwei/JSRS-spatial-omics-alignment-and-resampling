#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(openxlsx)
})

input_dir <- Sys.getenv("INPUT_DIR", unset = "input")
outdir <- Sys.getenv("OUTPUT_DIR", unset = "registration_qc_output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
sections <- strsplit(Sys.getenv("SECTIONS", unset = "cs12_3,cs12_4,cs12_5"), ",")[[1]]

read_coords <- function(file, section_id, source) {
  df <- fread(file)
  if (!all(c("x", "y") %in% colnames(df))) stop("Coordinate file must contain x and y: ", file)
  df %>%
    transmute(section_id = section_id, source = source, x = as.numeric(x), y = as.numeric(y)) %>%
    filter(is.finite(x), is.finite(y))
}

normalize_bbox <- function(df, ref_df) {
  df <- df %>% mutate(
    x_norm = (x - min(x)) / (max(x) - min(x) + 1e-12),
    y_norm = (y - min(y)) / (max(y) - min(y) + 1e-12)
  )
  ref_range <- ref_df %>% summarise(xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y))
  df %>% transmute(
    section_id, source,
    x = x_norm * (ref_range$xmax - ref_range$xmin) + ref_range$xmin,
    y = y_norm * (ref_range$ymax - ref_range$ymin) + ref_range$ymin
  )
}

rna_all <- list()
msi_before_all <- list()
msi_after_all <- list()

for (sec in sections) {
  rna_file <- file.path(input_dir, paste0("rna_coords_", sec, ".csv"))
  before_file <- file.path(input_dir, paste0("msi_coords_before_", sec, ".csv"))
  after_file <- file.path(input_dir, paste0("msi_coords_after_", sec, ".csv"))
  if (!file.exists(rna_file) || !file.exists(before_file) || !file.exists(after_file)) next
  rna <- read_coords(rna_file, sec, "RNA")
  before <- read_coords(before_file, sec, "MSI_before")
  after <- read_coords(after_file, sec, "MSI_after")
  before_norm <- normalize_bbox(before, rna)
  rna_all[[sec]] <- rna
  msi_before_all[[sec]] <- before_norm
  msi_after_all[[sec]] <- after
}

wb <- createWorkbook()
addWorksheet(wb, "RNA_mask")
writeData(wb, "RNA_mask", bind_rows(rna_all))
addWorksheet(wb, "MSI_before_mask")
writeData(wb, "MSI_before_mask", bind_rows(msi_before_all))
addWorksheet(wb, "MSI_after_mask")
writeData(wb, "MSI_after_mask", bind_rows(msi_after_all))
saveWorkbook(wb, file.path(outdir, "registration_masks.xlsx"), overwrite = TRUE)
