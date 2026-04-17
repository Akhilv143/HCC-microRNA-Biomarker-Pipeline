# Transcriptomic Meta-Analysis of miRNA-Gene Interactions in Hepatocellular Carcinoma (HCC)

## Project Overview
This repository contains a comprehensive bioinformatics workflow implemented in R for analyzing differential miRNA expression in Hepatocellular Carcinoma (HCC). By integrating multiple independent bulk RNA-seq datasets across different platforms, this project mitigates study-specific biases to uncover robust miRNA-gene regulatory networks driving HCC pathogenesis.

## Pipeline Implementation
The entire computational workflow—from raw data acquisition and batch correction to differential expression and target prediction—is implemented in a single, modular R script.

* **Main Analysis Script:** [`mirna_hcc.R`](mirna_hcc.R)

This script is designed to be reproducible, provided the required R dependencies and raw data accessions are available.

## Datasets Analyzed
Raw transcriptomic count data was sourced from the NCBI Gene Expression Omnibus (GEO). The meta-analysis integrates two independent cohorts to establish a robust combined miRNA dataset comparing normal tissue controls against HCC tumor samples.

* **GSE227378:** High-throughput miRNA profiling.
* **GSE76903:** Deep sequencing of small RNAs in paired normal and tumor liver tissues.

## Analytical Pipeline and Methodology

### 1. Data Acquisition and Preprocessing
* Raw count matrices and metadata were dynamically fetched and extracted using the `GEOquery` package.
* Expression counts corresponding to mature miRNAs (`hsa-miR-` and `hsa-let-`) were isolated and merged into a unified count matrix.
* Low-expressed miRNAs were filtered out to enhance statistical power.

### 2. Batch Effect Correction & Differential Expression
* `DESeq2` was utilized for robust differential expression analysis.
* To account for the technical variance of merging two independent GEO datasets, study origin was explicitly incorporated as a covariate in the DESeq2 design matrix (`~ batch + condition`).
* Significant differentially expressed miRNAs were isolated based on stringent adjusted p-value and log2 fold-change thresholds.

### 3. miRNA Target Prediction
* The `multiMiR` package was deployed to identify downstream target genes for the top up-regulated and down-regulated miRNAs.
* To ensure high biological relevance, the script intersects experimentally validated targets (from miRTarBase) with highly-scored predicted targets (from miRDB, score ≥ 80).

## Project Visualizations

### 1. Quality Control and PCA
Principal Component Analysis (PCA) was performed on Variance Stabilizing Transformed (VST) data to visualize sample clustering by biological condition and GEO batch.
<p align="center">
  <img src="results/plots/QC_PCA/PCA_by_Condition.png" width="48%" alt="PCA by Condition">
  <img src="results/plots/QC_PCA/PCA_by_Batch.png" width="48%" alt="PCA by Batch">
</p>

<br>

### 2. Differential Expression (Volcano Plot)
Visualizing the statistical significance and magnitude of change for differentially expressed miRNAs.
<p align="center">
  <img src="results/plots/DEG/Volcano_plot.png" width="600" alt="Volcano Plot">
</p>

<br>

### 3. Hierarchical Clustering (Heatmap)
Z-score scaled expression heatmap of the top 25 up-regulated and top 25 down-regulated miRNAs across all samples, grouped by condition.
<p align="center">
  <img src="results/plots/Heatmap/Heatmap_top25up_top25dn.png" width="800" alt="miRNA Heatmap">
</p>

## miRNA-Gene Interactions & Downstream Analysis
A critical component of this pipeline is the transition from miRNA profiling to functional gene networks. The script generates highly filtered target lists that serve as direct inputs for downstream network biology and clinical validation tools.

* **STRING Database Inputs:** High-confidence gene targets are exported as `STRING_input_UP.txt` and `STRING_input_DOWN.txt`. These files are pre-formatted for direct upload into the STRING-db web server for protein-protein interaction (PPI) network construction.
* **Network Node Matrices:** `targets_UP_degree_counts.csv` and `targets_DOWN_degree_counts.csv` provide quantitative mapping of how many distinct differentially expressed miRNAs target each specific gene.
* **Clinical Validation Preparation:** Generates `Top10_miRNAs_for_UALCAN.csv`, automatically stripping precursor suffixes to match the exact ID formatting required for external survival validation in the TCGA-LIHC cohort.

## Repository Structure

* **`mirna_hcc.R`**: The core analysis script executing the full transcriptomic and target prediction pipeline.
* **`/results/tables/`**: Contains the complete DESeq2 output, combined metadata, count matrices, and formatting for UALCAN.
* **`/results/tables/` (Interactions)**: Contains the highly filtered miRNA-gene target degrees and formatted text files for STRING-db.
* **`/results/plots/`**: Contains publication-quality outputs (PCA, Volcano, Heatmap) saved in PNG and high-resolution TIFF formats.
* **`/results/rds/`**: Serialized `multiMiR` database objects for rapid reproducible loading.
* **`.gitignore`**: Configured to exclude heavy raw data files (.tar, .gz) and massive TIFF images to maintain repository efficiency.

## R Dependencies
The pipeline requires the following core libraries:
`GEOquery`, `DESeq2`, `limma`, `ComplexHeatmap`, `multiMiR`, `tidyverse`, `ggplot2`, `ggrepel`.
