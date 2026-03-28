# =========================================================================
# Script 3: Ferroptosis & SASP Pathway Scoring and Correlation
# =========================================================================

library(Seurat)
library(ggplot2)

# 1. Define Gene Signatures (FerrDb & MSigDB)
# Note: Full lists provided in Supplementary Materials
ferro_genes <- list(c("ACSL4", "PTGS2", "TFRC", "SLC7A11", "GPX4")) 
sasp_genes <- list(c("IL1B", "CXCL8", "CCL2", "MMP3", "TNF"))
panel_genes <- list(c("CD83", "TNFAIP3", "NFKBIA"))

# 2. Module Scoring
sc_obj <- AddModuleScore(sc_obj, features = ferro_genes, name = "Ferroptosis_Score")
sc_obj <- AddModuleScore(sc_obj, features = sasp_genes, name = "SASP_Score")
sc_obj <- AddModuleScore(sc_obj, features = panel_genes, name = "Panel_Score")

# 3. Pearson Correlation Analysis at Single-Cell Resolution
cor_ferro <- cor.test(sc_obj$Panel_Score1, sc_obj$Ferroptosis_Score1, method = "pearson")
cor_sasp <- cor.test(sc_obj$Panel_Score1, sc_obj$SASP_Score1, method = "pearson")

print(paste("Ferroptosis Correlation: R =", round(cor_ferro$estimate, 3), "p-value =", cor_ferro$p.value))
print(paste("SASP Correlation: R =", round(cor_sasp$estimate, 3), "p-value =", cor_sasp$p.value))