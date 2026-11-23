# Single-Cell RNA-Seq Profiling of the Breast Invasive Ductal Carcinoma Immune Microenvironment (`GSE253587` / `GSM8023621`)

## Overview
End-to-end single-cell transcriptomics (`scRNA-seq`) analysis pipeline built in Python using **Scanpy** on a real human primary **Breast Invasive Ductal Carcinoma (`IDC`)** tumor-infiltrating immune (`CD45+`) dataset from the Dana-Farber Cancer Institute (**NCBI GEO `GSE253587`**, sample **`GSM8023621`**). This study resolves the cellular heterogeneity between anti-tumor cytotoxic effector lymphocytes and immunosuppressive regulatory T cells (`Tregs`) / myeloid populations driving tumor immune escape.

## Pipeline Workflow & Key Quantitative Results
1. **10x Genomics Matrix Ingestion & Clinical Quality Control (`scanpy.pp`):**
   - Loaded raw 10x Cell Ranger sparse count matrices (`1,362` raw barcodes $	imes$ `33,538` human `GRCh38` genes).
   - Filtered apoptotic/dying cells (`pct_counts_mt < 15%`), empty droplets (`min_genes >= 200`), and doublets (`n_genes_by_counts < 4,500`), retaining **956 high-confidence cells $	imes$ 16036 expressed genes**.
2. **Normalization, `adata.raw` Preservation & Feature Selection:**
   - Normalized UMI counts to `10,000` per cell, applied $\log(1+x)$ transformation, preserved the full `16,036`-gene transcriptome in `adata.raw`, and extracted **2,000 Highly Variable Genes (HVGs)**.
3. **Dimensionality Reduction & Graph Clustering (`scanpy.tl`):**
   - Performed Z-score scaling, PCA, k-Nearest Neighbor graph construction (`n_neighbors=12`, `n_pcs=15`), **UMAP** embedding, and **Leiden community detection** (`resolution=0.6`), resolving **8 distinct tumor-microenvironment cell populations**.
4. **Wilcoxon Marker Discovery & Biological Cell-Type Annotation:**
   - Executed Wilcoxon rank-sum differential expression (`rank_genes_groups`, Benjamini–Hochberg FDR < 0.05) across all `16,036` genes to annotate every cluster:
     - **B Cells (`CD79A+ MS4A1+`):** `24.6%` (`n = 235`)
     - **Activated Lymphocytes (`MALAT1-high`):** `19.4%` (`n = 185`)
     - **Cytotoxic CD8+ T / NK Cells (`NKG7+ CCL5+ CTSW+`):** `14.7%` (`n = 141`)
     - **Regulatory T Cells / Tregs (`TNFRSF18+ CD3D+ TRAC+`):** `13.1%` (`n = 125`) — *key immunosuppressive population mediating immune escape*
     - **CD4+ Helper T Cells (`IL7R+ LTB+`):** `11.6%` (`n = 111`)
     - **Inflammatory Stromal / CAFs (`NNMT+ SAA1+`):** `7.1%` (`n = 68`)
     - **Ductal Epithelial / Tumor Cells (`ELF3+ MGP+`):** `5.1%` (`n = 49`)
     - **Macrophages / Myeloid Cells (`SPI1+ TYROBP+ CST3+`):** `4.4%` (`n = 42`)

## Repository Structure
- `04_scRNAseq_Breast_Cancer_Immune_Microenvironment.ipynb`
- `README.md`
- `results/figures/` (Pre-QC violin plots, HVG dispersion, PCA/UMAP clusters, annotated composition bar chart, canonical marker DotPlot, 6-panel feature UMAP)
- `results/tables/idc_tme_marker_genes.csv`
