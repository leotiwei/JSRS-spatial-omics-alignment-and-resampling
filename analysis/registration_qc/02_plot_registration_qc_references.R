#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

mask_file <- Sys.getenv("MASK_FILE", unset = "registration_qc_output/registration_masks.xlsx")
outdir <- Sys.getenv("OUTPUT_DIR", unset = "registration_qc_output/reference_plots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

rna_mask <- read.xlsx(mask_file, sheet = "RNA_mask")
before_mask <- read.xlsx(mask_file, sheet = "MSI_before_mask")
after_mask <- read.xlsx(mask_file, sheet = "MSI_after_mask")

plot_section <- function(section_id) {
  rna <- rna_mask %>% filter(section_id == !!section_id)
  before <- before_mask %>% filter(section_id == !!section_id)
  after <- after_mask %>% filter(section_id == !!section_id)

  p1 <- ggplot(rna, aes(x, -y)) + geom_point(size = 0.1, alpha = 0.4) + coord_fixed() + theme_void() + ggtitle(paste(section_id, "RNA reference"))
  p2 <- ggplot() +
    geom_point(data = rna, aes(x, -y), size = 0.1, alpha = 0.25, color = "grey40") +
    geom_point(data = before, aes(x, -y), size = 0.1, alpha = 0.45, color = "red") +
    coord_fixed() + theme_void() + ggtitle("Before registration")
  p3 <- ggplot() +
    geom_point(data = rna, aes(x, -y), size = 0.1, alpha = 0.25, color = "grey40") +
    geom_point(data = after, aes(x, -y), size = 0.1, alpha = 0.45, color = "blue") +
    coord_fixed() + theme_void() + ggtitle("After registration")
  p1 + p2 + p3
}

for (sec in unique(rna_mask$section_id)) {
  p <- plot_section(sec)
  ggsave(file.path(outdir, paste0(sec, "_registration_reference.pdf")), p, width = 9, height = 3)
  ggsave(file.path(outdir, paste0(sec, "_registration_reference.png")), p, width = 9, height = 3, dpi = 300)
}
