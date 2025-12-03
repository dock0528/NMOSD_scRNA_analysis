#devtools::install_github("jinworks/CellChat")
library(Seurat)
library(SeuratDisk)
library(CellChat)
library(NMF)
library(circlize)
library(ComplexHeatmap)

#----{Import RDS}
merge_data <- readRDS("C:/Users/Jane/Desktop/Wang實驗室/NMOSD研究計畫/scRNA/scRNA_DATA/My_merged_protein_coding_genes/My_merged_count_metadata(protein_coding).rds")

#----{Normalize:校正不同細胞測序深度的差異}
merge_data<-NormalizeData(merge_data, normalization.method = "LogNormalize", scale.factor = 10000)

#View(merge_data)

#----{scANVI metadat 加入 RDS}
# 讀scANVI後的metadata celltype
celltype_metadata <- read.csv("../scRNA_DATA/My_merged_protein_coding_genes/My_merged_cell_metadata_CellChat.csv")
head(celltype_metadata)

# 確保rownames為Cell ID
rownames(celltype_metadata)<- celltype_metadata$cell_id

# 加"scANVI分群" & "組別"加入merge_data
merge_data <- AddMetaData(
  object = merge_data,
  metadata = celltype_metadata$predicted_label,
  col.name = "Cell_label"
)
merge_data <- AddMetaData(
  object = merge_data,
  metadata = celltype_metadata$Condition,
  col.name = "Condition"
)

#----{Create Cellchat Object}
seurat_NMOSD <- subset(merge_data, subset = Condition == "NMOSD")
cellchat_NMOSD <- createCellChat(
  object  = seurat_NMOSD,      # Seurat object
  group.by = "Cell_label"    # cell type label
)

#----{L-R database}
CellChatDB<-CellChatDB.human
showDatabaseCategory(CellChatDB)
#library(dplyr)
#dplyr::glimpse(CellChatDB$interaction) #ligand–receptor資料庫

# 取database (可只選某部分)
#先全選
CellChatDB.use<-CellChatDB

#加入cellchat_NMOSD
cellchat_NMOSD@DB<-CellChatDB.use
cellchat_NMOSD <- subsetData(cellchat_NMOSD) # 資料過濾

# === 【Identify over-expressed "ligands" or "receptors" in each cell group】===
#devtools::install_github("immunogenomics/presto")
library(presto) #讓 Wilcoxon 在 dgCSparse 下跑
cellchat_NMOSD <- identifyOverExpressedGenes(cellchat_NMOSD)

# identify over-expressed "L–R interactions"
cellchat_NMOSD<- identifyOverExpressedInteractions(cellchat_NMOSD)

### Each cell types lognormalized average
computeAveExpr(cellchat_NMOSD, features = c("TNFRSF13C"),  
               type = "triMean") 

### 去除頭尾0.1%計算lognormalized average
computeAveExpr(
  cellchat_NMOSD,
  features = "TNFRSF13C",
  type = "truncatedMean",
  trim = 0.1
)

computeAveExpr(cellchat_NMOSD, features = c("CMPK2"),  
               type = "triMean")  # no result
computeAveExpr(cellchat_NMOSD, features = c("ISG15"),  
               type = "triMean")  # no result
computeAveExpr(cellchat_NMOSD, features = c("LY6E"),  
               type = "triMean")  # no result
computeAveExpr(cellchat_NMOSD, features = c("OASL"),  
               type = "triMean")  # no result

# ===【 Infer cell–cell communication at a L–R pair level】===
cellchat_NMOSD <- computeCommunProb(cellchat_NMOSD, type = "triMean", trim = NULL, 
                              raw.use = TRUE,nboot = 20) #約20mins
# raw.use = TRUE:用原始資料(lognormalized)做計算
# nboot:跑幾次permutation決定pvalue (default=100)

# ===【Filter fewer cells】===
cellchat_NMOSD <- filterCommunication(cellchat_NMOSD, min.cells = 10) #某celltype少於10cells->drop #default=10

# ===【 Infer cell–cell communication at a signaling pathway level】===
cellchat_NMOSD <- computeCommunProbPathway(cellchat_NMOSD,thresh = 0.05) #thresh:significant interaction P-value criteria
#saveRDS(cellchat_NMOSD, file = "../scRNA_DATA/My_merged_protein_coding_genes/cellchat_NMOSD.rds")

# ===【TNFRSF13C】===
lr.idx <- grep("TNFRSF13C", cellchat_NMOSD@LR$LRsig$interaction_name)
lr.use <- cellchat_NMOSD@LR$LRsig$interaction_name[lr.idx]
lr.use #"TNFSF13B_TNFRSF13C"
lr_TNFRSF13C_df<- data.frame(interaction_name = lr.use)
df_TNFRSF13C <- subsetCommunication(cellchat_NMOSD,  pairLR.use = lr_TNFRSF13C_df)

# ===【TNFRSF13C聚焦target = B naive】===
df_TNFRSF13C_Bnaive <- subset(df_TNFRSF13C, target == "B naive")
df_TNFRSF13C_Bnaive

sources.use <- unique(df_TNFRSF13C_Bnaive$source)
targets.use <- "B naive"
cellchat_TNFRSF13C_NMOSD <- aggregateNet(
  cellchat_NMOSD,
  sources.use = sources.use,
  targets.use = targets.use
)


par(mfrow = c(1, 1))
plot.new()
plot.window(xlim = c(-1, 1), ylim = c(-1, 1))
netVisual_aggregate(
  cellchat_TNFRSF13C_NMOSD,
  signaling = "BAFF",
  layout = "chord", 
  vertex.size.max = 6, 
  edge.width.max = 10,
  remove.isolate=TRUE,
  title.space = 4
)

