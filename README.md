# Integrative Transcriptomic Profiling of miRNA-Target Regulatory Networks in Hepatocellular Carcinoma (HCC)

## Project Overview
This repository contains a comprehensive bioinformatics R script designed to evaluate differential miRNA expression in Hepatocellular Carcinoma (HCC). By integrating multiple independent high-throughput RNA-seq datasets across different sequencing platforms, this project mitigates study-specific biases to uncover robust miRNA-gene regulatory networks driving HCC pathogenesis.

## Datasets Analyzed
Raw transcriptomic count data was sourced from the NCBI Gene Expression Omnibus (GEO). The integration of two independent cohorts establishes a robust combined miRNA dataset encompassing a total of 104 samples: 52 normal adjacent tissue controls and 52 HCC primary tumor samples. 

* **GSE227378:** High-throughput small RNA sequencing (BGISEQ-500) of paired tumor and adjacent normal tissues from 32 Chinese HCC patients (64 samples total). This study originally focused on identifying miRNAs like miR-3180 linked to overall survival and disease progression.
* **GSE76903:** Deep sequencing of small RNAs (Illumina HiSeq 2500) from HCC patients. While the original study included portal vein tumor thrombosis (PVTT) tissues, this specific analysis strictly filtered for the 20 primary tumor and 20 matched adjacent normal samples to maintain a clean control-versus-disease comparison.

## Analytical Workflow and Methodology

### 1. Data Acquisition and Preprocessing
* Raw count matrices and clinical metadata were dynamically fetched and extracted using the `GEOquery` package.
* Expression counts corresponding to mature miRNAs (`hsa-miR-` and `hsa-let-`) were isolated and merged into a unified count matrix.
* Low-expressed miRNAs were systematically filtered out to enhance statistical power across the combined dataset.

### 2. Batch Effect Correction and Differential Expression
* `DESeq2` was utilized for rigorous differential expression analysis.
* To account for the severe technical variance of merging distinct sequencing platforms (BGISEQ vs. Illumina), study origin was explicitly incorporated as a covariate in the DESeq2 design matrix (`~ batch + condition`).
* Significant differentially expressed miRNAs were isolated based on stringent adjusted p-value and log2 fold-change thresholds.

### 3. miRNA Target Prediction (miRTarBase and miRDB)
* The `multiMiR` package was deployed to identify downstream target genes for the top up-regulated and down-regulated miRNAs.
* To ensure the highest biological relevance, the script maps and intersects experimentally validated targets sourced directly from the **miRTarBase** database with highly-scored predicted targets sourced from the **miRDB** database (prediction score >= 80).

## Project Visualizations

### 1. Quality Control and Principal Component Analysis
Principal Component Analysis (PCA) was performed on Variance Stabilizing Transformed (VST) data to visualize sample clustering and validate the distribution of biological conditions against GEO batch effects.
<p align="center">
  <img src="results/plots/QC_PCA/PCA_by_Condition.png" width="48%" alt="PCA by Condition">
  <img src="results/plots/QC_PCA/PCA_by_Batch.png" width="48%" alt="PCA by Batch">
</p>

<br>

### 2. Differential Expression (Volcano Plot)
Visualizing the statistical significance and magnitude of change for differentially expressed miRNAs between the normal and tumor cohorts.
<p align="center">
  <img src="results/plots/DEG/Volcano_plot.png" width="600" alt="Volcano Plot">
</p>

<br>

### 3. Hierarchical Clustering (Heatmap)
Z-score scaled expression heatmap of the top 25 up-regulated and top 25 down-regulated miRNAs across all 104 samples, annotated by condition and dataset origin.
<p align="center">
  <img src="results/plots/Heatmap/Heatmap_top25up_top25dn.png" width="800" alt="miRNA Heatmap">
</p>

## miRNA-Gene Interactions and Downstream Outputs
A critical component of this project is the translation of miRNA profiling into functional gene networks. The R script generates highly filtered target lists that serve as direct inputs for network biology and validation tools.

* **STRING Database Inputs:** High-confidence gene targets are exported as `STRING_input_UP.txt` and `STRING_input_DOWN.txt`. These files are pre-formatted for direct upload into the STRING-db web server to construct protein-protein interaction (PPI) networks.
* **Network Node Matrices:** `targets_UP_degree_counts.csv` and `targets_DOWN_degree_counts.csv` provide quantitative mapping of the topological degree, detailing exactly how many distinct differentially expressed miRNAs target each specific gene.
* **Clinical Validation Preparation:** Generates `Top10_miRNAs_for_UALCAN.csv`, automatically stripping precursor suffixes to match the exact mature ID formatting required for external validation in the TCGA-LIHC cohort.

## Repository Structure

* **`mirna_hcc.R`**: The core R script executing the full transcriptomic and target prediction workflow.
* **`/results/tables/`**: Contains the complete DESeq2 output, combined metadata, count matrices, and formatting for external databases.
* **`/results/tables/` (Interactions)**: Contains the highly filtered miRNA-gene target degrees and formatted text files for STRING-db.
* **`/results/plots/`**: Contains publication-quality graphical outputs saved in PNG and high-resolution TIFF formats.
* **`/results/rds/`**: Serialized database objects for rapid, reproducible loading without requiring repeated web queries.
* **`.gitignore`**: Configured to exclude heavy raw data files (.tar, .gz) and massive TIFF images to maintain repository efficiency.

## R Dependencies
The script requires the following core libraries:
`GEOquery`, `DESeq2`, `limma`, `ComplexHeatmap`, `multiMiR`, `tidyverse`, `ggplot2`, `ggrepel`.
