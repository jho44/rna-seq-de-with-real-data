# Normalization and Exploration
# Size factor normalization, VST transformation, and exploratory visualizations

library(DESeq2)
library(tidyverse)
library(pheatmap)

cat(paste("=", paste(rep("-", 70), collapse = ""), "=\n", sep = ""))
cat("Normalization and Exploration\n")
cat(paste("=", paste(rep("-", 70), collapse = ""), "=\n\n", sep = ""))

# Load filtered data
cat("1. Loading filtered data...\n")
counts_filtered <- readRDS("data/processed/counts_filtered.rds")
sample_meta <- readRDS("data/processed/sample_metadata.rds")
gene_info <- readRDS("data/processed/gene_info_filtered.rds")

# ========================================
# 1. Create DESeqDataSet
# ========================================
cat("2. Creating DESeqDataSet object...\n")

# Ensure sample_meta rownames match counts colnames
rownames(sample_meta) <- sample_meta$sample_id

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = sample_meta,
  design = ~ condition
)

cat(paste("   - Samples:", ncol(dds), "\n"))
cat(paste("   - Genes:", nrow(dds), "\n\n"))

# ========================================
# 2. Size factor normalization
# ========================================
cat("3. Estimating size factors...\n")

dds <- estimateSizeFactors(dds)
size_factors <- sizeFactors(dds)

cat(paste("   - Size factors range: ", round(min(size_factors), 3), " to ",
          round(max(size_factors), 3), "\n\n", sep = ""))

# Get normalized counts
counts_normalized <- counts(dds, normalized = TRUE)

# ========================================
# 3. VST transformation
# ========================================
cat("4. Performing VST transformation...\n")

vst <- vst(dds, blind = TRUE)
vst_data <- assay(vst)

cat("   ✓ VST transformation complete\n\n")

# ========================================
# 4. PCA plot
# ========================================
cat("5. Creating PCA plot...\n")

pca_result <- prcomp(t(vst_data))
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  sample_id = rownames(pca_result$x),
  condition = sample_meta$condition[match(rownames(pca_result$x), sample_meta$sample_id)]
)

var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

png("results/plots/pca_vst.png", width = 800, height = 600)
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)"),
    title = "PCA Plot (VST-transformed counts)",
    color = "Condition"
  ) +
  scale_color_manual(values = c("Tumor" = "#E41A1C", "Normal" = "#377EB8"))
print(p)
dev.off()

cat(paste("   ✓ Saved: results/plots/pca_vst.png\n"))
cat(paste("   - PC1 explains ", var_explained[1], "% of variance\n"))
cat(paste("   - PC2 explains ", var_explained[2], "% of variance\n\n", sep = ""))

# ========================================
# 5. Sample correlation statistics
# ========================================
cat("6. Calculating sample correlation statistics...\n")

sample_cors <- cor(vst_data)
tumor_samples <- which(sample_meta$condition == "Tumor")
normal_samples <- which(sample_meta$condition == "Normal")

tumor_tumor_corr <- mean(sample_cors[tumor_samples, tumor_samples][upper.tri(sample_cors[tumor_samples, tumor_samples])])
normal_normal_corr <- mean(sample_cors[normal_samples, normal_samples][upper.tri(sample_cors[normal_samples, normal_samples])])
tumor_normal_corr <- mean(sample_cors[tumor_samples, normal_samples])

cat(paste("   - Tumor-Tumor (within-group): ",
          round(tumor_tumor_corr, 3), "\n", sep = ""))
cat(paste("   - Normal-Normal (within-group): ",
          round(normal_normal_corr, 3), "\n", sep = ""))
cat(paste("   - Tumor-Normal (between-group): ",
          round(tumor_normal_corr, 3), "\n\n", sep = ""))

# ========================================
# 6. Save normalized data for DE analysis
# ========================================
cat("7. Saving normalized data...\n")

saveRDS(dds, "data/processed/dds.rds")
saveRDS(vst_data, "data/processed/vst_data.rds")
saveRDS(counts_normalized, "data/processed/counts_normalized.rds")

cat("   ✓ Saved DESeqDataSet\n")
cat("   ✓ Saved VST-transformed counts\n")
cat("   ✓ Saved normalized counts\n\n")

# ========================================
# Summary
# ========================================
cat(paste("=", paste(rep("=", 70), collapse = ""), "=\n", sep = ""))
cat("Normalization and Exploration Complete!\n")
cat(paste("=", paste(rep("=", 70), collapse = ""), "=\n\n", sep = ""))

cat("Output plots:\n")
cat("  results/plots/pca_vst.png\n\n")
cat("Note: Sample correlation heatmap deferred to DE analysis\n")
cat("(will create heatmap of significant DE genes after DE testing)\n\n")

cat("Next step: Run 04_de_analysis.R\n")
