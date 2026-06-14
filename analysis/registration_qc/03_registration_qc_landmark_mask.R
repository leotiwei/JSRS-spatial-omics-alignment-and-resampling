#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

input_dir <- Sys.getenv("INPUT_DIR", unset = "input")
outdir <- Sys.getenv("OUTPUT_DIR", unset = "registration_qc_output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

landmark_file <- file.path(input_dir, "registration_landmarks.tsv")
mask_file <- file.path(outdir, "registration_masks.xlsx")
grid_size <- as.numeric(Sys.getenv("GRID_SIZE", unset = "100"))
dilation_radius <- as.integer(Sys.getenv("DILATION_RADIUS", unset = "1"))

landmarks <- read.delim(landmark_file, check.names = FALSE)
required <- c("section_id", "landmark_id", "rna_x", "rna_y", "met_before_x", "met_before_y", "met_after_x", "met_after_y")
if (!all(required %in% colnames(landmarks))) stop("Landmark table is missing required columns.")

landmark_errors <- landmarks %>%
  mutate(
    before_error = sqrt((rna_x - met_before_x)^2 + (rna_y - met_before_y)^2),
    after_error = sqrt((rna_x - met_after_x)^2 + (rna_y - met_after_y)^2)
  ) %>%
  pivot_longer(cols = c(before_error, after_error), names_to = "registration_state", values_to = "error") %>%
  mutate(registration_state = recode(registration_state, before_error = "before", after_error = "after"))

coords_to_keys <- function(df, grid_size = 100, dilation_radius = 1) {
  base <- df %>%
    mutate(x_grid = floor(x / grid_size), y_grid = floor(y / grid_size)) %>%
    distinct(section_id, x_grid, y_grid)
  offsets <- expand.grid(dx = -dilation_radius:dilation_radius, dy = -dilation_radius:dilation_radius)
  expanded <- tidyr::crossing(base, offsets) %>%
    transmute(section_id, key = paste(x_grid + dx, y_grid + dy, sep = "_")) %>%
    distinct(section_id, key)
  expanded
}

calc_overlap <- function(mask_a, mask_b, state) {
  sections <- intersect(unique(mask_a$section_id), unique(mask_b$section_id))
  bind_rows(lapply(sections, function(sec) {
    a <- mask_a %>% filter(section_id == sec) %>% pull(key) %>% unique()
    b <- mask_b %>% filter(section_id == sec) %>% pull(key) %>% unique()
    intersection <- length(intersect(a, b))
    union <- length(union(a, b))
    data.frame(
      section_id = sec,
      registration_state = state,
      dice = 2 * intersection / (length(a) + length(b)),
      jaccard = intersection / union,
      n_rna_cells = length(a),
      n_msi_cells = length(b)
    )
  }))
}

if (file.exists(mask_file)) {
  rna <- read.xlsx(mask_file, sheet = "RNA_mask")
  before <- read.xlsx(mask_file, sheet = "MSI_before_mask")
  after <- read.xlsx(mask_file, sheet = "MSI_after_mask")
  rna_keys <- coords_to_keys(rna, grid_size, dilation_radius)
  before_keys <- coords_to_keys(before, grid_size, dilation_radius)
  after_keys <- coords_to_keys(after, grid_size, dilation_radius)
  overlap <- bind_rows(calc_overlap(rna_keys, before_keys, "before"), calc_overlap(rna_keys, after_keys, "after"))
} else {
  overlap <- data.frame()
}

wb <- createWorkbook()
addWorksheet(wb, "landmark_errors")
writeData(wb, "landmark_errors", landmark_errors)
addWorksheet(wb, "mask_overlap")
writeData(wb, "mask_overlap", overlap)
saveWorkbook(wb, file.path(outdir, "registration_qc_results.xlsx"), overwrite = TRUE)

p1 <- ggplot(landmark_errors, aes(x = registration_state, y = error)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_point(position = position_jitter(width = 0.1), size = 0.8, alpha = 0.6) +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "Landmark error")
ggsave(file.path(outdir, "landmark_error_before_after.pdf"), p1, width = 4, height = 3)

if (nrow(overlap) > 0) {
  p2 <- ggplot(overlap, aes(x = registration_state, y = dice)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_point(position = position_jitter(width = 0.1), size = 0.8, alpha = 0.6) +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = "Dice overlap")
  ggsave(file.path(outdir, "mask_overlap_before_after.pdf"), p2, width = 4, height = 3)
}
