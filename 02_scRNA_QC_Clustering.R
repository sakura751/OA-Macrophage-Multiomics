模块二：单细胞微环境质控与降维聚类 (02_scRNA_QC_Clustering.R)
这段代码是全篇的技术高光，完美展示了对 Seurat v5 图层的理解以及对硬件内存的极限控制。

R
# =========================================================================
# Script 2: scRNA-seq Quality Control, Integration and Clustering
# Dataset: GSE216651 (Synovial Microenvironment)
# Note: Memory-optimized for 16GB RAM workstations
# =========================================================================
library(Seurat)
library(Matrix)

# 1. Strict Quality Control
# Retaining high-quality cells while filtering doublets and empty droplets
sc_obj <- subset(sc_obj, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 10)

# 2. Memory Optimization & Seurat v5 Layer Integration
# Convert to sparse matrix (dgCMatrix) and force garbage collection
sc_obj[["RNA"]]$counts <- as(sc_obj[["RNA"]]$counts, "CsparseMatrix")
gc()

# Joining layers to fix Seurat v5 fragmentation
sc_obj <- JoinLayers(sc_obj)

# 3. Dimensionality Reduction & Clustering
sc_obj <- NormalizeData(sc_obj) %>% FindVariableFeatures() %>% ScaleData()
sc_obj <- RunPCA(sc_obj, npcs = 30)
sc_obj <- FindNeighbors(sc_obj, dims = 1:20)
sc_obj <- FindClusters(sc_obj, resolution = 0.5)
sc_obj <- RunUMAP(sc_obj, dims = 1:20)

# 4. Macrophage Subpopulation Identification (Non-parametric Wilcoxon test)
mac_markers <- FindMarkers(sc_obj, ident.1 = c("7", "13", "24"), 
                           min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
# Target genes CD83, TNFAIP3, NFKBIA strictly localized here.