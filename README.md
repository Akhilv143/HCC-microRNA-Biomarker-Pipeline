# Integrative Profiling of the miRNA-Gene Interactome in Hepatocellular Carcinoma (HCC)

## Project Overview
This repository hosts a high-throughput computational pipeline for investigating the differential expression of microRNAs (miRNAs) in Hepatocellular Carcinoma (HCC). By integrating transcriptomic data from independent cohorts, this workflow establishes a consistent biological signature of HCC while mitigating technical noise associated with different sequencing technologies.

## Biological Significance
Hepatocellular Carcinoma is characterized by complex regulatory disruptions. miRNAs serve as critical post-transcriptional regulators that can silence tumor suppressors or activate oncogenic pathways. This project identifies core miRNAs that are consistently dysregulated across multiple patient groups, providing a foundation for understanding the HCC regulatory landscape.

## Main Analysis Script
The entire analytical lifecycle—from raw data ingestion to high-confidence target prediction—is centralized in one modular script:

* **[`mirna_hcc.R`](mirna_hcc.R)**: This script automates library management, batch effect correction, differential expression statistics, and the integration of miRNA-gene databases.

## Datasets and Study Cohorts
Raw miRNA expression counts were retrieved from the NCBI Gene Expression Omnibus (GEO). To ensure biological consistency, the analysis utilizes a total of **104 samples**, consisting of 52 HCC tumors and 52 matched normal liver tissues.

| Dataset ID | Platform | Focus | Samples (N/T) |
|:---:|:---:|:---|:---:|
| **GSE227378** | BGISEQ-500 | Small RNA sequencing of Chinese HCC patients. | 32 Normal / 32 Tumor |
| **GSE76903** | Illumina HiSeq 2500 | Paired primary tumor and normal adjacent profiling. | 20 Normal / 20 Tumor |

## Top Differentially Expressed miRNAs
The following tables summarize the most significantly dysregulated miRNAs identified after cross-platform integration.

### Top 5 Up-regulated miRNAs
| miRNA | Log2 Fold Change | Adjusted p-value | Regulation |
|:---|:---:|:---:|:---:|
| hsa-miR-301a-3p | 1.79 | 3.08e-32 | Up |
| hsa-miR-421 | 1.76 | 2.47e-31 | Up |
| hsa-miR-3200-3p | 3.03 | 1.43e-29 | Up |
| hsa-miR-183-5p | 2.94 | 1.62e-29 | Up |
| hsa-miR-501-5p | 2.14 | 1.21e-27 | Up |

### Top 5 Down-regulated miRNAs
| miRNA | Log2 Fold Change | Adjusted p-value | Regulation |
|:---|:---:|:---:|:---:|
| hsa-miR-139-5p | -2.19 | 7.28e-30 | Down |
| hsa-miR-139-3p | -1.99 | 6.34e-25 | Down |
| hsa-miR-490-3p | -2.96 | 1.28e-21 | Down |
| hsa-miR-490-5p | -3.24 | 8.13e-19 | Down |
| hsa-miR-99a-3p | -1.54 | 1.08e-15 | Down |

## Computational Methodology

### 1. Multi-Cohort Integration
The pipeline dynamically fetches raw count data. Because the datasets were generated on different sequencing platforms (BGI vs. Illumina), the script utilizes **DESeq2** with a multi-factor design matrix. By including `batch` (Study ID) in the model, we isolate the biological signal (`condition`) from the technical variance.

### 2. Differential Expression Analysis
Statistical significance is calculated using the Wald test within DESeq2. We apply the Benjamini-Hochberg procedure to control the False Discovery Rate (FDR), ensuring that the identified miRNAs are robust across both studies.

### 3. Regulatory Network Discovery (miRTarBase & miRDB)
To determine the functional impact of these miRNAs, we perform a dual-database target search using the `multiMiR` package:
* **Validated Targets:** Sourced from **miRTarBase**, focusing on interactions with strong experimental evidence (Reporter assay, Western blot, etc.).
* **Predicted Targets:** Sourced from **miRDB**, utilizing only those with a target prediction score ≥ 80 to minimize false positives.

## Visualizations and QC
All plots are generated in both PNG and high-resolution TIFF formats for publication readiness.

* **PCA Analysis:** Validates that samples cluster primarily by biological status (Tumor vs. Normal) rather than their laboratory of origin.
* **Volcano Plot:** Provides a global view of the fold-change distribution and statistical significance.
* **Heatmap Clustering:** Demonstrates consistent miRNA expression patterns across the 104 samples using Z-score normalization.

## Repository Structure
* **`results/tables/`**: Final DESeq2 results, count matrices, and UALCAN-ready ID lists.
* **`results/plots/`**: Visual outputs for QC and Differential Expression.
* **`results/rds/`**: Pre-saved database query results for faster reproducibility.
* **`.gitignore`**: Optimized to keep the repository light by excluding large raw `.tar` files and intermediate data.

## R Environment Requirements
The pipeline depends on the following Bioconductor and CRAN packages:
`GEOquery`, `DESeq2`, `limma`, `ComplexHeatmap`, `multiMiR`, `tidyverse`, `ggplot2`, `ggrepel`, `BiocParallel`.
