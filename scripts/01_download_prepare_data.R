# TCGA-LUAD Data Download and Preparation
# Downloads TCGA RNA-seq counts and filters for tumor vs normal tissue

library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)

# Create output directories
dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

cat(paste("=", paste(rep("-", 70), collapse = ""), "=\n", sep = ""))
cat("TCGA-LUAD Data Download & Preparation\n")
cat(paste("=", paste(rep("-", 70), collapse = ""), "=\n\n", sep = ""))

# ========================================
# 1. Download RNA-seq counts
# ========================================
cat("1. Downloading TCGA-LUAD RNA-seq counts...\n")

query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

# Download data
GDCdownload(query, method = "api")

# Prepare data
rnaseq_data <- GDCprepare(query)
cat(paste("✓ Downloaded", ncol(rnaseq_data), "samples\n\n"))

# ========================================
# 2. Extract sample metadata
# ========================================
cat("2. Processing sample metadata...\n")

sample_meta <- as.data.frame(colData(rnaseq_data))
sample_meta$sample_id <- rownames(sample_meta)
sample_meta$patient_id <- substr(sample_meta$sample_id, 1, 12)

# Create condition variable (Tumor vs Normal)
sample_meta$condition <- ifelse(sample_meta$sample_type == "Solid Tissue Normal", "Normal", "Tumor")

cat(paste("✓ Tumor samples:", sum(sample_meta$condition == "Tumor"), "\n"))
cat(paste("✓ Normal tissue samples:", sum(sample_meta$condition == "Normal"), "\n\n"))

# ========================================
# 3. Extract and save count matrix
# ========================================
cat("3. Extracting counts and saving data...\n")

# Get count matrix
counts <- assay(rnaseq_data)

# Ensure numeric
counts <- matrix(as.numeric(counts), nrow = nrow(counts))
rownames(counts) <- rownames(rnaseq_data)
colnames(counts) <- colnames(rnaseq_data)

# Get gene annotations
gene_info_raw <- as.data.frame(rowData(rnaseq_data))
cat("\nAvailable gene annotation columns:\n")
cat(paste(colnames(gene_info_raw), collapse = ", "), "\n\n")

# Select available columns for gene info
available_cols <- colnames(gene_info_raw)
if ("gene_name" %in% available_cols) {
  gene_info <- gene_info_raw %>%
    select(gene_name) %>%
    rownames_to_column("gene_id")
} else if ("external_gene_name" %in% available_cols) {
  gene_info <- gene_info_raw %>%
    select(external_gene_name) %>%
    rownames_to_column("gene_id")
} else {
  # If no gene name column, just use gene_id
  gene_info <- data.frame(
    gene_id = rownames(gene_info_raw)
  )
}

# Save raw counts
write.csv(counts, "data/raw/counts.csv", row.names = TRUE)
write.csv(gene_info, "data/raw/gene_info.csv", row.names = FALSE)

# Save and prepare sample metadata
sample_meta_save <- sample_meta %>%
  select(sample_id, patient_id, sample_type, condition)

write.csv(sample_meta_save, "data/raw/sample_metadata.csv", row.names = FALSE)

# ========================================
# 4. Save processed data for analysis
# ========================================
cat("4. Saving processed objects...\n\n")

# Create a clean count matrix with rownames
counts_clean <- counts
rownames(counts_clean) <- rownames(rnaseq_data)

# Save for next step
saveRDS(counts_clean, "data/processed/counts.rds")
saveRDS(sample_meta_save, "data/processed/sample_metadata.rds")
saveRDS(gene_info, "data/processed/gene_info.rds")

# ========================================
# Summary
# ========================================
cat(paste("=", paste(rep("=", 70), collapse = ""), "=\n", sep = ""))
cat("Data Preparation Complete!\n")
cat(paste("=", paste(rep("=", 70), collapse = ""), "=\n\n", sep = ""))

cat("Summary:\n")
cat(paste("  Total samples:", ncol(counts_clean), "\n"))
cat(paste("  Total genes:", nrow(counts_clean), "\n"))
cat(paste("  Tumor samples:", sum(sample_meta_save$condition == "Tumor"), "\n"))
cat(paste("  Normal tissue samples:", sum(sample_meta_save$condition == "Normal"), "\n\n"))

cat("Output files:\n")
cat("  data/raw/counts.csv\n")
cat("  data/raw/gene_info.csv\n")
cat("  data/raw/sample_metadata.csv\n")
cat("  data/processed/counts.rds\n")
cat("  data/processed/sample_metadata.rds\n")
cat("  data/processed/gene_info.rds\n\n")

cat("Next step: Run 02_qc_and_filtering.R\n")
