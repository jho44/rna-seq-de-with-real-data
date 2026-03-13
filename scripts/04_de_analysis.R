# Differential Expression Analysis with DESeq2
# Statistical testing and visualization of DE results

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)

cat(paste("=", paste(rep("-", 70), collapse = ""), "=\n", sep = ""))
cat("Differential Expression Analysis\n")
cat(paste("=", paste(rep("-", 70), collapse = ""), "=\n\n", sep = ""))

# Load data
cat("1. Loading data...\n")
dds <- readRDS("data/processed/dds.rds")
gene_info <- readRDS("data/processed/gene_info_filtered.rds")

cat(paste("   - Samples:", ncol(dds), "\n"))
cat(paste("   - Genes:", nrow(dds), "\n\n"))

# ========================================
# 1. Run DESeq2
# ========================================
cat("2. Running DESeq2 analysis...\n")

dds <- DESeq(dds)

cat("   ✓ Dispersion estimation complete\n")
cat("   ✓ Statistical testing complete\n\n")

# ========================================
# 2. Extract and process results
# ========================================
cat("3. Extracting results...\n")

# Get results with contrast: Tumor vs Normal
results_deseq <- results(dds, contrast = c("condition", "Tumor", "Normal"))

# Apply LFC shrinkage for better fold change estimates
# Use ashr method which works with contrast
results_deseq <- lfcShrink(
  dds,
  contrast = c("condition", "Tumor", "Normal"),
  type = "ashr"
)

# Convert to dataframe
de_results <- as.data.frame(results_deseq) %>%
  rownames_to_column("gene_id") %>%
  left_join(
    gene_info %>% select(gene_id, gene_name),
    by = "gene_id"
  ) %>%
  arrange(padj)

cat(paste("   - Total genes tested:", nrow(de_results), "\n"))
cat(paste("   - Significant genes (padj < 0.05):", sum(de_results$padj < 0.05, na.rm = TRUE), "\n\n"))

# ========================================
# 3. Save results tables
# ========================================
cat("4. Saving results tables...\n")

# All results
write.csv(de_results, "results/tables/de_results_all.csv", row.names = FALSE)

# Significant results
de_sig <- de_results %>%
  filter(padj < 0.05)

write.csv(de_sig, "results/tables/de_results_significant.csv", row.names = FALSE)

cat("   ✓ Saved: results/tables/de_results_all.csv\n")
cat("   ✓ Saved: results/tables/de_results_significant.csv\n\n")

# Print top DE genes
cat("Top 10 upregulated genes (Tumor > Normal):\n")
top_up <- de_sig %>% filter(log2FoldChange > 0) %>% head(10)
print(top_up %>% select(gene_id, gene_name, log2FoldChange, padj))

cat("\nTop 10 downregulated genes (Tumor < Normal):\n")
top_down <- de_sig %>% filter(log2FoldChange < 0) %>% head(10)
print(top_down %>% select(gene_id, gene_name, log2FoldChange, padj))

cat("\n")

# ========================================
# 4. MA plot
# ========================================
cat("5. Creating MA plot...\n")

ma_df <- de_results %>%
  mutate(
    sig = ifelse(padj < 0.05, "Significant", "Not Significant"),
    sig = factor(sig, levels = c("Not Significant", "Significant"))
  )

png("results/plots/ma_plot.png", width = 900, height = 600)
p <- ggplot(ma_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  theme_minimal() +
  labs(
    x = "log10(Mean Expression + 1)",
    y = "log2(Fold Change)",
    title = "MA Plot: Tumor vs Normal",
    color = "Status"
  ) +
  scale_color_manual(values = c("Not Significant" = "gray60", "Significant" = "#E41A1C"))
print(p)
dev.off()

cat("   ✓ Saved: results/plots/ma_plot.png\n\n")

# ========================================
# 5. Volcano plot
# ========================================
cat("6. Creating volcano plot...\n")

volcano_df <- de_results %>%
  mutate(
    sig = ifelse(padj < 0.05, "Significant", "Not Significant"),
    sig = factor(sig, levels = c("Not Significant", "Significant")),
    neg_log10_padj = -log10(padj)
  ) %>%
  replace_na(list(neg_log10_padj = 0))

png("results/plots/volcano_plot.png", width = 900, height = 600)
p <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj, color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(
    x = "log2(Fold Change)",
    y = "-log10(Adjusted P-value)",
    title = "Volcano Plot: Tumor vs Normal",
    color = "Status"
  ) +
  scale_color_manual(values = c("Not Significant" = "gray60", "Significant" = "#E41A1C"))
print(p)
dev.off()

cat("   ✓ Saved: results/plots/volcano_plot.png\n\n")

# ========================================
# 6. Heatmap of top DE genes
# ========================================
cat("7. Creating heatmap of top DE genes...\n")

# Get top 50 DE genes by adjusted p-value
top_genes <- de_results %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(50) %>%
  pull(gene_id)

# Load VST data
vst_data <- readRDS("data/processed/vst_data.rds")
sample_meta <- readRDS("data/processed/sample_metadata.rds")

# Extract expression for top genes
top_gene_expr <- vst_data[top_genes, ]

# Create sample annotation
sample_annotation <- data.frame(
  Condition = sample_meta$condition
)
rownames(sample_annotation) <- sample_meta$sample_id

# Create heatmap
png("results/plots/heatmap_top_de_genes.png", width = 1200, height = 800)
pheatmap(top_gene_expr,
  annotation_col = sample_annotation,
  annotation_colors = list(Condition = c("Tumor" = "#E41A1C", "Normal" = "#377EB8")),
  main = "Top 50 DE Genes: Tumor vs Normal",
  clustering_distance_cols = "correlation",
  clustering_distance_rows = "euclidean",
  fontsize = 8,
  fontsize_row = 7,
  show_colnames = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B35806"))(n = 100)
)
dev.off()

cat("   ✓ Saved: results/plots/heatmap_top_de_genes.png\n\n")

# ========================================
# Summary
# ========================================
cat(paste("=", paste(rep("=", 70), collapse = ""), "=\n", sep = ""))
cat("Differential Expression Analysis Complete!\n")
cat(paste("=", paste(rep("=", 70), collapse = ""), "=\n\n", sep = ""))

cat("Summary:\n")
cat(paste("  Genes tested:", nrow(de_results), "\n"))
cat(paste("  Significant genes (padj < 0.05):", nrow(de_sig), "\n"))
cat(paste("  Upregulated in Tumor:", sum(de_sig$log2FoldChange > 0, na.rm = TRUE), "\n"))
cat(paste("  Downregulated in Tumor:", sum(de_sig$log2FoldChange < 0, na.rm = TRUE), "\n\n"))

cat("Output tables:\n")
cat("  results/tables/de_results_all.csv\n")
cat("  results/tables/de_results_significant.csv\n\n")

cat("Output plots:\n")
cat("  results/plots/ma_plot.png\n")
cat("  results/plots/volcano_plot.png\n")
cat("  results/plots/heatmap_top_de_genes.png\n")
