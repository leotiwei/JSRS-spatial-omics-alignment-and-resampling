#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(grid)
})

project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
input_dir <- file.path(project_dir, "cluster_average")
outdir <- file.path(project_dir, "all_cluster_pair_robustness")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

method_files <- c(
  conservative = file.path(input_dir, "Metabolite_cluster_average_conservative.xlsx"),
  nearest = file.path(input_dir, "Metabolite_cluster_average_nearest.xlsx"),
  IDW_k4_p2 = file.path(input_dir, "Metabolite_cluster_average_IDW_k4_p2.xlsx")
)

reference_method <- Sys.getenv("REFERENCE_METHOD", unset = "conservative")
cor_method <- Sys.getenv("COR_METHOD", unset = "spearman")
min_features <- as.integer(Sys.getenv("MIN_FEATURES", unset = "50"))
top_n <- as.integer(Sys.getenv("TOP_N", unset = "60"))
drop_clusters <- strsplit(Sys.getenv("DROP_CLUSTERS", unset = "Otic,otic,OTIC"), ",")[[1]]

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

calc_log2fc <- function(mat, cluster_a, cluster_b, pseudocount = 1e-6) {
  log2((mat[, cluster_a] + pseudocount) / (mat[, cluster_b] + pseudocount))
}

top_abs_features <- function(x, n) {
  x <- x[is.finite(x) & !is.na(x)]
  names(sort(abs(x), decreasing = TRUE))[seq_len(min(n, length(x)))]
}

save_pheatmap <- function(ph, filename, width = 5, height = 4) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()
}

mats <- lapply(method_files, read_feature_cluster_matrix)
common_features <- Reduce(intersect, lapply(mats, rownames))
common_clusters <- setdiff(Reduce(intersect, lapply(mats, colnames)), drop_clusters)
mats <- lapply(mats, function(m) m[common_features, common_clusters, drop = FALSE])

method_names <- names(mats)
method_pairs <- expand.grid(method_1 = method_names, method_2 = method_names, stringsAsFactors = FALSE) %>%
  filter(method_1 < method_2)

method_level <- method_pairs %>%
  rowwise() %>%
  mutate(
    abundance_similarity = cor(as.vector(mats[[method_1]]), as.vector(mats[[method_2]]), method = cor_method, use = "pairwise.complete.obs"),
    abundance_scaled_similarity = cor(as.vector(scale_by_feature(mats[[method_1]])), as.vector(scale_by_feature(mats[[method_2]])), method = cor_method, use = "pairwise.complete.obs"),
    cluster_structure_similarity = cor(
      as.vector(cor(mats[[method_1]], method = cor_method, use = "pairwise.complete.obs")),
      as.vector(cor(mats[[method_2]], method = cor_method, use = "pairwise.complete.obs")),
      method = cor_method,
      use = "pairwise.complete.obs"
    )
  ) %>%
  ungroup()

cluster_pairs <- t(combn(common_clusters, 2))
all_pair_results <- list()
for (i in seq_len(nrow(cluster_pairs))) {
  cluster_a <- cluster_pairs[i, 1]
  cluster_b <- cluster_pairs[i, 2]
  lfc <- lapply(mats, calc_log2fc, cluster_a = cluster_a, cluster_b = cluster_b)
  keep <- Reduce("&", lapply(lfc, function(x) is.finite(x) & !is.na(x)))
  if (sum(keep) < min_features) next
  base <- lfc[[reference_method]][keep]
  for (method_name in setdiff(names(lfc), reference_method)) {
    comp <- lfc[[method_name]][keep]
    base_top <- top_abs_features(base, top_n)
    comp_top <- top_abs_features(comp, top_n)
    all_pair_results[[length(all_pair_results) + 1]] <- data.frame(
      cluster_a = cluster_a,
      cluster_b = cluster_b,
      comparison = paste0(reference_method, "_vs_", method_name),
      n_features = sum(keep),
      log2fc_spearman = cor(base, comp, method = cor_method),
      direction_agreement = mean(sign(base) == sign(comp)),
      top_overlap_fraction = length(intersect(base_top, comp_top)) / length(unique(c(base_top, comp_top)))
    )
  }
}
all_pair_summary <- bind_rows(all_pair_results)

all_lfc_vectors <- lapply(mats, function(mat) {
  out <- c()
  for (i in seq_len(nrow(cluster_pairs))) {
    cluster_a <- cluster_pairs[i, 1]
    cluster_b <- cluster_pairs[i, 2]
    out <- c(out, calc_log2fc(mat, cluster_a, cluster_b))
  }
  out[is.finite(out) & !is.na(out)]
})
common_len <- min(lengths(all_lfc_vectors))
all_lfc_vectors <- lapply(all_lfc_vectors, function(x) x[seq_len(common_len)])
all_pair_lfc_similarity <- outer(method_names, method_names, Vectorize(function(a, b) {
  cor(all_lfc_vectors[[a]], all_lfc_vectors[[b]], method = cor_method, use = "pairwise.complete.obs")
}))
dimnames(all_pair_lfc_similarity) <- list(method_names, method_names)

wb <- createWorkbook()
addWorksheet(wb, "method_level_similarity")
writeData(wb, "method_level_similarity", method_level)
addWorksheet(wb, "all_cluster_pair_summary")
writeData(wb, "all_cluster_pair_summary", all_pair_summary)
addWorksheet(wb, "all_pair_lfc_similarity")
writeData(wb, "all_pair_lfc_similarity", as.data.frame(all_pair_lfc_similarity), rowNames = TRUE)
saveWorkbook(wb, file.path(outdir, "all_cluster_pair_robustness_summary.xlsx"), overwrite = TRUE)

p <- all_pair_summary %>%
  pivot_longer(cols = c(log2fc_spearman, direction_agreement, top_overlap_fraction), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = comparison, y = value)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~metric, nrow = 1) +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "Robustness metric") +
  coord_cartesian(ylim = c(0, 1))

ggsave(file.path(outdir, "all_cluster_pair_robustness_summary.pdf"), p, width = 9, height = 3.5)
ggsave(file.path(outdir, "all_cluster_pair_robustness_summary.png"), p, width = 9, height = 3.5, dpi = 300)

ph <- pheatmap(all_pair_lfc_similarity, silent = TRUE, main = "All-pair log2FC similarity")
save_pheatmap(ph, file.path(outdir, "all_pair_log2fc_similarity_heatmap.pdf"), width = 4.5, height = 4)
