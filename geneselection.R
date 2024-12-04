tpm = master_harmonized_data[["training"]][["pbmc_gene_expression"]][["wide_tpm"]]
tpm_tr = tpm[, apply(tpm == 0, 2, mean) < 0.1]

library(matrixStats)
library(stats)
library(ggplot2)
library(dplyr)

expression_data <- tpm_tr[, !(colnames(tpm_tr) %in% 'specimen_id')]
specimin_id <- tpm_tr$specimen_id

gene_scores <- data.frame(
  gene = colnames(expression_data),
  variance = NA,
  mad = NA,
  mean_expression = NA,
  expression_range = NA,
  pca_loading = NA,
  correlation_score = NA,
  total_score = 0
)

gene_variances <- apply(expression_data, 2, var)
gene_scores$variance <- gene_variances
variance_threshold <- quantile(gene_variances, 0.75)
gene_scores$total_score <- gene_scores$total_score + as.integer(gene_scores$variance >= variance_threshold)

gene_mads <- apply(expression_data, 2, mad)
gene_scores$mad <- gene_mads
mad_threshold <- quantile(gene_mads, 0.75)
gene_scores$total_score <- gene_scores$total_score + as.integer(gene_scores$mad >= mad_threshold)

gene_means <- colMeans(expression_data)
gene_scores$mean_expression <- gene_means
mean_expression_threshold <- quantile(gene_means, 0.75)
gene_scores$total_score <- gene_scores$total_score + as.integer(gene_scores$mean_expression >= mean_expression_threshold)

gene_ranges <- apply(expression_data, 2, function(x) max(x) - min(x))
gene_scores$expression_range <- gene_ranges
range_threshold <- quantile(gene_ranges, 0.75)
gene_scores$total_score <- gene_scores$total_score + as.integer(gene_scores$expression_range >= range_threshold)

expression_data_scaled <- scale(expression_data)
pca_result <- prcomp(expression_data_scaled, center = TRUE, scale. = TRUE)
pc1_loadings <- abs(pca_result$rotation[, 1])
gene_scores$pca_loading <- pc1_loadings
pca_loading_threshold <- quantile(pc1_loadings, 0.75)
gene_scores$total_score <- gene_scores$total_score + as.integer(gene_scores$pca_loading >= pca_loading_threshold)

cor_matrix <- cor(expression_data, use = "pairwise.complete.obs")
average_correlations <- apply(cor_matrix, 1, function(x) mean(abs(x)))
gene_scores$correlation_score <- average_correlations
correlation_threshold <- quantile(average_correlations, 0.75)
gene_scores$total_score <- gene_scores$total_score + as.integer(gene_scores$correlation_score >= correlation_threshold)

max_score <- 6
score_threshold <- 4
selected_genes <- gene_scores$gene[gene_scores$total_score >= score_threshold]
selected_genes <- c(selected_genes, "specimen_id")
save(selected_genes, file = "selected_genes.RData")

selected_expression_data <- expression_data[, selected_genes]
selected_expression_data$specimin_id <- specimin_id

hist(gene_scores$total_score, breaks = max_score, main = "Distribution of Gene Scores",
     xlab = "Total Score", ylab = "Number of Genes", col = "lightblue")

library(pheatmap)
pheatmap(t(selected_expression_data[, selected_genes]), scale = "row",
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, show_colnames = FALSE,
         main = "Heatmap of Selected Genes")

pca_selected <- prcomp(selected_expression_data[, selected_genes], center = TRUE, scale. = TRUE)
plot(pca_selected$x[, 1], pca_selected$x[, 2], xlab = "PC1", ylab = "PC2",
     main = "PCA Plot of Samples Using Selected Genes", pch = 19)
