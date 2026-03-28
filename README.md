# OA-Macrophage-Multiomics

**Code repository for the integrative multi-omics analysis of the osteoarthritis synovial microenvironment.**

This repository contains the custom computational scripts utilized in our study to seamlessly bridge macroscopic diagnostic signatures with microscopic single-cell pathogenesis in osteoarthritis (OA).

## 💻 Computational Environment
* **Primary Environment:** R version 4.5.2 (Native execution avoiding memory overflow via `dgCMatrix` sparse matrix architecture)
* **Key R Packages:** `glmnet` (v4.1-8), `randomForest` (v4.7-1.1), `Boruta` (v8.0.0), `Seurat` (v5.0.1), `CellChat` (v1.6.0)

## 🧬 Analytical Pipeline Overview
1. **Macroscopic Machine Learning Optimization:** Scripts for integrating LASSO, Random Forest, and Boruta algorithms to identify and validate the 3-gene diagnostic consensus core (*CD83*, *TNFAIP3*, *NFKBIA*).
2. **Microscopic Single-Cell Spatial Tracing:** Code for rigorous scRNA-seq quality control, non-linear dimensionality reduction (UMAP), and the precise mapping of macroscopic targets to pathogenic macrophage subpopulations.
3. **Mechanistic & Network Crosstalk Analysis:** Implementations for evaluating dual metabolic hyperactivation (Ferroptosis & SASP) and inferring intercellular communication networks via CellChat.

## 📂 Data Availability
The raw high-throughput transcriptomic datasets analyzed by these scripts are publicly available in the NCBI Gene Expression Omnibus (GEO) repository under accession numbers **GSE12021**, **GSE89408**, and **GSE216651**.
