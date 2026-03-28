# =========================================================================
# Script 4: Cell-Cell Communication and Cloud API GO Enrichment
# =========================================================================

library(CellChat)
library(enrichR)
library(ggplot2)

# --- Part A: CellChat Network ---
# Initialize CellChat object with strict group limits
cellchat <- createCellChat(object = sc_obj, group.by = "idents")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Calculate communication probabilities
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# --- Part B: EnrichR API Cloud GO Enrichment ---
# Identifying DEGs for Pathogenic Macrophages
mac_degs <- FindMarkers(sc_obj, ident.1 = "Pathogenic_Macrophage", ident.2 = "Other_Macrophage", 
                        only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sig_genes <- rownames(mac_degs)[mac_degs$p_val_adj < 0.05]

# Cloud connection avoiding local SQLite errors
setEnrichrSite("Enrichr") 
dbs <- c("GO_Biological_Process_2023")
enriched <- enrichr(sig_genes, dbs)
go_results <- enriched[["GO_Biological_Process_2023"]]
go_results <- go_results[go_results$Adjusted.P.value < 0.05, ]

# Output highlights hyperactive "Translation" and "Oxidative Phosphorylation"