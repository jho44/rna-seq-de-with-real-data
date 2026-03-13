# QC and Low-Count Gene Filtering
# Performs quality control and filters out low-count genes

library(tidyverse)
library(edgeR)

cat(paste("=", paste(rep("-", 70), collapse = ""), "=\n", sep = ""))
cat("QC and Filtering\n")
cat(paste("=", paste(rep("-", 70), collapse = ""), "=\n\n", sep = ""))

# Load data
cat("Loading data...\n")
counts <- readRDS("data/processed/counts.rds")
sample_meta <- readRDS("data/processed/sample_metadata.rds")
gene_info <- readRDS("data/processed/gene_info.rds")

# ========================================
# 1. Library size QC
# ========================================
cat("\n1. Library Size Quality Control\n")
cat(paste(rep("-", 70), collapse = ""), "\n", sep = "")

lib_sizes <- colSums(counts)
lib_sizes_df <- data.frame(
  sample = names(lib_sizes),
  library_size = lib_sizes,
  condition = sample_meta$condition[match(names(lib_sizes), sample_meta$sample_id)]
)

cat(paste("Mean library size:", round(mean(lib_sizes)), "\n"))
cat(paste("Min library size:", round(min(lib_sizes)), "\n"))
cat(paste("Max library size:", round(max(lib_sizes)), "\n"))
cat(paste("SD:", round(sd(lib_sizes)), "\n\n"))

# Plot library sizes
png("results/plots/library_sizes.png", width = 800, height = 600)
p <- ggplot(lib_sizes_df, aes(x = condition, y = library_size, fill = condition)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 2) +
  theme_minimal() +
  theme(
    legend.position = "top"
  ) +
  labs(x = "Tissue Type", y = "Library Size (reads)", title = "Sequencing Depth Distribution", fill = "Condition") +
  scale_fill_manual(values = c("Tumor" = "#E41A1C", "Normal" = "#377EB8"))
print(p)
dev.off()
cat("✓ Saved: results/plots/library_sizes.png\n\n")

# ========================================
# 2. Count distribution
# ========================================
cat("2. Count Distribution\n")
cat(paste(rep("-", 70), collapse = ""), "\n", sep = "")

count_dist <- data.frame(
  counts_per_gene = rowSums(counts)
) %>%
  rownames_to_column("gene_id")

cat(paste("Mean counts per gene:", round(mean(count_dist$counts_per_gene)), "\n"))
cat(paste("Median counts per gene:", round(median(count_dist$counts_per_gene)), "\n\n"))

# ========================================
# 3. Low-count gene filtering
# ========================================
cat("3. Low-Count Gene Filtering\n")
cat(paste(rep("-", 70), collapse = ""), "\n", sep = "")

# Filter: keep genes with at least 10 total counts across all samples
min_count_threshold <- 10
keep <- rowSums(counts) >= min_count_threshold

cat(paste("Filtering threshold: >= ", min_count_threshold, " total counts\n", sep = ""))
cat(paste("Genes before filtering:", nrow(counts), "\n"))
cat(paste("Genes after filtering:", sum(keep), "\n"))
cat(paste("Genes removed:", nrow(counts) - sum(keep), "\n\n"))

# Filter counts and gene info
counts_filtered <- counts[keep, ]
gene_info_filtered <- gene_info[gene_info$gene_id %in% rownames(counts_filtered), ]

# ========================================
# 4. Visualization of filtering
# ========================================
png("results/plots/count_distribution.png", width = 800, height = 600)
p <- count_dist %>%
  ggplot(aes(x = log10(counts_per_gene + 1))) +
  geom_histogram(bins = 50, fill = "#377EB8", color = "black", alpha = 0.7) +
  geom_vline(xintercept = log10(min_count_threshold), linetype = "dashed", color = "#E41A1C", linewidth = 1.5) +
  theme_minimal() +
  labs(
    x = "log10(Total Counts + 1)",
    y = "Number of Genes",
    title = "Gene Count Distribution with Filtering Threshold"
  )
print(p)
dev.off()
cat("✓ Saved: results/plots/count_distribution.png\n\n")

# ========================================
# 5. Sample-level statistics
# ========================================
cat("4. Sample-Level Statistics After Filtering\n")
cat(paste(rep("-", 70), collapse = ""), "\n", sep = "")

sample_stats <- data.frame(
  sample = colnames(counts_filtered),
  total_counts = colSums(counts_filtered),
  genes_detected = colSums(counts_filtered > 0),
  condition = sample_meta$condition[match(colnames(counts_filtered), sample_meta$sample_id)]
)

print(sample_stats %>% arrange(condition))
cat("\n")

# ========================================
# 6. Save filtered data
# ========================================
cat("5. Saving Filtered Data\n")
cat(paste(rep("-", 70), collapse = ""), "\n\n", sep = "")

saveRDS(counts_filtered, "data/processed/counts_filtered.rds")
saveRDS(gene_info_filtered, "data/processed/gene_info_filtered.rds")
saveRDS(sample_meta, "data/processed/sample_metadata.rds")

cat("✓ Saved filtered count matrix\n")
cat("✓ Saved filtered gene info\n")
cat("✓ Saved sample metadata\n\n")

cat(paste("=", paste(rep("=", 70), collapse = ""), "=\n", sep = ""))
cat("QC and Filtering Complete!\n")
cat(paste("=", paste(rep("=", 70), collapse = ""), "=\n\n", sep = ""))

cat("Next step: Run 03_normalization.R\n")
