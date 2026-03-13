# RNA-seq Differential Expression Analysis: LUAD Tumor vs Normal Lung

## Project Overview

Differential expression analysis comparing lung adenocarcinoma (LUAD) tumors to normal lung tissue using TCGA-LUAD RNA-seq data.

**Research Goal**: Identify genes with significant expression differences between LUAD tumors and normal lung tissue, revealing tumor-specific expression signatures.

> _This analysis was developed with AI assistance (Claude) while applying concepts from Datacamp's [RNA-Seq with Bioconductor in R](https://app.datacamp.com/learn/courses/rna-seq-with-bioconductor-in-r) course and Liv Grant's [I Tested Claude Code on a Bioinformatics Project](https://www.youtube.com/watch?v=-hCyjQUF5RE) video_

## Data Source

**TCGA-LUAD (Lung Adenocarcinoma)**:

- RNA-seq counts from GDC (Genomic Data Commons)
- Samples: 540 tumor + 59 normal tissue samples
- Filter: Primary tumors + normal lung tissue only

## Analysis Workflow

```
1. Download & Prepare Data (01_download_prepare_data.R)
   ├── Fetch TCGA-LUAD RNA-seq counts
   ├── Filter for primary tumors + normal tissue
   └── Combine into analysis matrix
                ↓
2. QC + Filtering (02_qc_and_filtering.R)
   ├── Library sizes & count distributions
   └── Low-count gene filtering
                ↓
3. Normalization & Exploration (03_normalization.R)
   ├── Size factor normalization
   ├── VST transformation
   └── PCA & sample QC by tissue type
                ↓
4. Differential Expression Analysis (04_de_analysis.R)
   ├── DESeq2 testing (tumor vs normal)
   ├── LFC shrinkage
   └── Visualization (MA, volcano plots)
```

## Expected Outputs

### Tables

- `results/tables/de_results_all.csv` - All tested genes with DE statistics
- `results/tables/de_results_significant.csv` - Significant DE genes (FDR < 0.05)

### Plots

- `results/plots/library_sizes.png` - Library sizes by tissue type (boxplot)
- `results/plots/count_distribution.png` - Gene count distribution with filtering threshold
- `results/plots/pca_vst.png` - PCA plot (tumor vs normal separation)
- `results/plots/ma_plot.png` - MA plot (fold-change vs mean expression)
- `results/plots/volcano_plot.png` - Volcano plot (fold-change vs p-value)
- `results/plots/heatmap_top_de_genes.png` - Heatmap of top 50 DE genes

## Tools & Dependencies

**Language**: R (v4.0+)

**Key packages**:

- `TCGAbiolinks` - Download TCGA data
- `DESeq2` - Differential expression
- `ggplot2` - Visualization
- `pheatmap` - Heatmaps

## File Structure

```
rna-seq-de-with-real-data/
├── README.md                      (this file)
├── scripts/
│   ├── 01_download_prepare_data.R (fetch & prep TCGA data)
│   ├── 02_qc_and_filtering.R      (QC & low-count filtering)
│   ├── 03_normalization.R         (normalization & exploration)
│   └── 04_de_analysis.R           (DESeq2 & visualization)
├── data/
│   ├── raw/                       (downloaded TCGA data)
│   └── processed/                 (filtered, normalized data)
└── results/
    ├── tables/                    (DE results)
    └── plots/                     (visualizations)
```

## Running the Analysis

Run scripts sequentially:

```R
source("scripts/01_download_prepare_data.R")  # Download & prepare data
source("scripts/02_qc_and_filtering.R")       # QC & filtering
source("scripts/03_normalization.R")          # Normalization
source("scripts/04_de_analysis.R")            # DESeq2 analysis
```

## Results Summary

- **Genes tested**: 54,717
- **Significant DE genes (padj < 0.05)**: 27,845 (51%)
- **Upregulated in tumor**: 20,343
- **Downregulated in tumor**: 7,502

Strong separation between tumor and normal in PCA (PC1: 11.2%) indicates high biological signal.

## Next Steps & Future Analyses

### 1. Investigate Outlier Tumor Samples

- [ ] Identify tumor samples that cluster with normal samples in top 50 DE gene heatmap
- [ ] Compare expression profiles: Are these tumors biologically distinct subtypes?
- [ ] Check molecular annotations: Do they represent different LUAD subtypes (e.g., lepidic, acinar, solid)?
- [ ] Analyze separately: May require stratified DE analysis by subtype
- [ ] Investigate if outliers share specific mutation profiles

**Tools**: PCA/hierarchical clustering analysis, TCGA sample annotations

### 2. Pathway & Functional Enrichment Analysis

- [ ] Gene Ontology (GO) enrichment on upregulated genes (biological processes)
- [ ] KEGG pathway enrichment (cancer-related pathways)
- [ ] Hallmark gene sets (MSigDB - cancer signatures)
- [ ] Reactome pathway analysis (mechanistic pathways)
- [ ] Gene Set Enrichment Analysis (GSEA) - rank-based, directional

**Tools**: `clusterProfiler`, `enrichplot`, `fgsea`, `msigdbr`

### 3. Additional Downstream Analyses

- [ ] **Filter by biotype**: Separate protein-coding vs lncRNA genes
- [ ] **Tissue specificity**: Compare with other TCGA cancer types
- [ ] **Known cancer drivers**: Validate against published lung cancer signatures
- [ ] **Co-expression networks**: Gene module discovery with `WGCNA`
- [ ] **Immune infiltration**: Estimate immune cell abundance using `immunedeconv`
- [ ] **Survival correlation**: Link DE genes to patient outcomes (if available)
- [ ] **Cross-dataset validation**: Validate findings in independent lung cancer cohorts (GEO)

## Notes

- TCGA data download may take time (large files)
- Tumor samples: LUAD primary tumors
- Normal samples: normal lung tissue from TCGA
- Analysis uses standard best practices (Love et al. 2014)
- High proportion of significant genes (51%) suggests strong tumor vs normal signal
