#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
input_dir <- file.path(project_dir, "cluster_average")
outdir <- file.path(project_dir, "interpolation_comparison")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

method_files <- c(
  conservative = file.path(input_dir, "Metabolite_cluster_average_conservative.xlsx"),
  nearest = file.path(input_dir, "Metabolite_cluster_average_nearest.xlsx"),
  IDW_k4_p2 = file.path(input_dir, "Metabolite_cluster_average_IDW_k4_p2.xlsx")
)

cor_method <- Sys.getenv("COR_METHOD", unset = "spearman")
reference_method <- Sys.getenv("REFERENCE_METHOD", unset = "conservative")
cluster_a <- Sys.getenv("CLUSTER_A", unset = "CMC")
cluster_b <- Sys.getenv("CLUSTER_B", unset = "SCP")
top_n <- as.integer(Sys.getenv("TOP_N", unset = "20"))

read_feature_cluster_matrix <- function(file, sheet = "Metabolite_average") {
  if (!file.exists(file)) stop("Missing input file: ", file)
  df <- read.xlsx(file, sheet = sheet, check.names = FALSE)
  feature_col <- colnames(df)[1]
  features <- as.character(df[[feature_col]])
  keep <- !is.na(features) & features != ""
  df <- df[keep, , drop = FALSE]
  features <- make.unique(features[keep])
  mat <- as.matrix(df[, setdiff(colnames(df), feature_col), drop = FALSE])
  suppressWarnings(mode(mat) <- "numeric")
  rownames(mat) <- features
  mat <- mat[rowSums(!is.na(mat)) > 0, colSums(!is.na(mat)) > 0, drop = FALSE]
  mat
}

scale_by_feature <- function(mat) {
  out <- t(scale(t(mat)))
  out[is.na(out)] <- 0
  out
}

calc_log2fc <- function(mat, group1, group2, pseudocount = 1e-6) {
  if (!all(c(group1, group2) %in% colnames(mat))) return(NULL)
  log2((mat[, group1] + pseudocount) / (mat[, group2] + pseudocount))
}

top_abs_features <- function(x, n) {
  x <- x[is.finite(x) & !is.na(x)]
  names(sort(abs(x), decreasing = TRUE))[seq_len(min(n, length(x)))]
}

mats <- lapply(method_files, read_feature_cluster_matrix)
common_features <- Reduce(intersect, lapply(mats, rownames))
common_clusters <- Reduce(intersect, lapply(mats, colnames))
mats <- lapply(mats, function(m) m[common_features, common_clusters, drop = FALSE])

method_pairs <- expand.grid(method_1 = names(mats), method_2 = names(mats), stringsAsFactors = FALSE) %>%
  filter(method_1 < method_2)

pair_summary <- method_pairs %>%
  rowwise() %>%
  mutate(
    n_features = length(common_features),
    n_clusters = length(common_clusters),
    abundance_spearman = cor(
      as.vector(mats[[method_1]]),
      as.vector(mats[[method_2]]),
      method = cor_method,
      use = "pairwise.complete.obs"
    ),
    abundance_scaled_spearman = cor(
      as.vector(scale_by_feature(mats[[method_1]])),
      as.vector(scale_by_feature(mats[[method_2]])),
      method = cor_method,
      use = "pairwise.complete.obs"
    ),
    cluster_similarity_spearman = cor(
      as.vector(cor(mats[[method_1]], method = cor_method, use = "pairwise.complete.obs")),
      as.vector(cor(mats[[method_2]], method = cor_method, use = "pairwise.complete.obs")),
      method = cor_method,
      use = "pairwise.complete.obs"
    )
  ) %>%
  ungroup()

logfc_list <- lapply(mats, calc_log2fc, group1 = cluster_a, group2 = cluster_b)
if (!is.null(logfc_list[[reference_method]])) {
  ref_logfc <- logfc_list[[reference_method]]
  logfc_summary <- bind_rows(lapply(setdiff(names(logfc_list), reference_method), function(method_name) {
    comp <- logfc_list[[method_name]]
    if (is.null(comp)) return(NULL)
    keep <- is.finite(ref_logfc) & is.finite(comp) & !is.na(ref_logfc) & !is.na(comp)
    ref_top <- top_abs_features(ref_logfc[keep], top_n)
    comp_top <- top_abs_features(comp[keep], top_n)
    data.frame(
      comparison = paste0(reference_method, "_vs_", method_name),
      cluster_a = cluster_a,
      cluster_b = cluster_b,
      n_features = sum(keep),
      log2fc_spearman = cor(ref_logfc[keep], comp[keep], method = cor_method),
      direction_agreement = mean(sign(ref_logfc[keep]) == sign(comp[keep])),
      top_overlap_fraction = length(intersect(ref_top, comp_top)) / length(unique(c(ref_top, comp_top)))
    )
  }))
} else {
  logfc_summary <- data.frame()
}

wb <- createWorkbook()
addWorksheet(wb, "method_level_similarity")
writeData(wb, "method_level_similarity", pair_summary)
addWorksheet(wb, "selected_pair_log2FC")
writeData(wb, "selected_pair_log2FC", logfc_summary)
saveWorkbook(wb, file.path(outdir, "interpolation_method_comparison_summary.xlsx"), overwrite = TRUE)

plot_df <- pair_summary %>%
  pivot_longer(
    cols = c(abundance_spearman, abundance_scaled_spearman, cluster_similarity_spearman),
    names_to = "metric",
    values_to = "spearman"
  )

p1 <- ggplot(plot_df, aes(x = interaction(method_1, method_2, sep = " vs "), y = spearman)) +
  geom_col() +
  facet_wrap(~metric, nrow = 1) +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "Spearman correlation") +
  coord_cartesian(ylim = c(0, 1))

ggsave(file.path(outdir, "method_level_similarity_summary.pdf"), p1, width = 9, height = 3.5)
ggsave(file.path(outdir, "method_level_similarity_summary.png"), p1, width = 9, height = 3.5, dpi = 300)

if (nrow(logfc_summary) > 0) {
  p2 <- logfc_summary %>%
    pivot_longer(cols = c(log2fc_spearman, direction_agreement, top_overlap_fraction), names_to = "metric", values_to = "value") %>%
    ggplot(aes(x = comparison, y = value)) +
    geom_col() +
    facet_wrap(~metric, nrow = 1) +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = "Robustness metric") +
    coord_cartesian(ylim = c(0, 1))
  ggsave(file.path(outdir, "selected_cluster_pair_log2fc_robustness.pdf"), p2, width = 8, height = 3.5)
  ggsave(file.path(outdir, "selected_cluster_pair_log2fc_robustness.png"), p2, width = 8, height = 3.5, dpi = 300)
}
