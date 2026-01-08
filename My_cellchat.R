#devtools::install_github("jinworks/CellChat")
library(Seurat)
library(SeuratDisk)
library(CellChat)
library(NMF)
library(circlize)
library(ComplexHeatmap)

#----{Import RDS}
merge_data <- readRDS("../scRNA_DATA/My_merged_protein_coding_genes/My_merged_count_metadata(protein_coding).rds")

#----{Normalize:校正不同細胞測序深度的差異}
merge_data<-NormalizeData(merge_data, normalization.method = "LogNormalize", scale.factor = 10000)

#View(merge_data)

#----{scANVI metadata 加入 RDS}
# 讀scANVI後的metadata celltype
celltype_metadata <- read.csv("../scRNA_DATA/My_merged_protein_coding_genes/My_merged_Azimuth_sub_celltype_CellChat.csv")
head(celltype_metadata)

# 確保rownames為Cell ID
rownames(celltype_metadata)<- celltype_metadata$cell_id

# 加"scANVI分群" & "組別"加入merge_data
merge_data <- AddMetaData(
  object = merge_data,
  metadata = celltype_metadata$sub_celltype,
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
#saveRDS(cellchat_NMOSD, file = "../scRNA_DATA/My_merged_protein_coding_genes/cellchat_NMOSD(Azimuth).rds")

# ===【Import CellChat_NMOSD RDS】===
cellchat_NMOSD <- readRDS("../scRNA_DATA/My_merged_protein_coding_genes/cellchat_NMOSD(Azimuth).rds")

# ===【TNFRSF13C】===
lr.idx <- grep("TNFRSF13C", cellchat_NMOSD@LR$LRsig$interaction_name)
lr.use <- cellchat_NMOSD@LR$LRsig$interaction_name[lr.idx]
lr.use #"TNFSF13B_TNFRSF13C"
lr_TNFRSF13C_df<- data.frame(interaction_name = lr.use)
df_TNFRSF13C <- subsetCommunication(cellchat_NMOSD,  pairLR.use = lr_TNFRSF13C_df)

# ===【TNFRSF13C聚焦target = B naive】===
df_TNFRSF13C_Bnaive <- subset(df_TNFRSF13C, target == "B naive")
df_TNFRSF13C_Bnaive

# ===【 Calculate the aggregated cell–cell communication network】===
# >>> for num of links & community probability
# {BAFF}
sources.use <- unique(df_TNFRSF13C$source)
targets.use <-  unique(df_TNFRSF13C$target )
cellchat_TNFRSF13C_BAFF_nmosd <- aggregateNet(
  cellchat_NMOSD,
  sources.use = sources.use,
  targets.use = targets.use
)

# {BAFF target B_naive}
sources.use <- unique(df_TNFRSF13C_Bnaive$source)
targets.use <- "B naive"
cellchat_TNFRSF13C_NMOSD <- aggregateNet(
  cellchat_NMOSD,
  sources.use = sources.use,
  targets.use = targets.use
)



# ===【Chord diagram - BAFF】===
par(mfrow = c(1, 1))
plot.new()
plot.window(xlim = c(-1, 1), ylim = c(-1, 1))
netVisual_aggregate(
  cellchat_TNFRSF13C_NMOSD,
  signaling = "BAFF", #signaling pathway
  signaling.name=NULL,
  sources.use =sources.use,
  layout = "chord", 
  vertex.size.max = 6, 
  edge.width.max = 10,
  remove.isolate=TRUE,
  title.space = 4
)

# ===【Chord diagram - BAFF in Bnaive】===
par(mfrow = c(1, 1))
plot.new()
plot.window(xlim = c(-1, 1), ylim = c(-1, 1))
netVisual_aggregate(
  cellchat_TNFRSF13C_NMOSD,
  signaling = "BAFF", #signaling pathway
  signaling.name='BAFF target B_naive',
  sources.use =sources.use,
  targets.use = targets.use,
  layout = "chord", 
  vertex.size.max = 6, 
  edge.width.max = 10,
  remove.isolate=TRUE,
  title.space = 4
)

# ===【Circle diagram - BAFF in Bnaive】===
par(mfrow = c(1, 1))
plot.new()
plot.window(xlim = c(-1, 1), ylim = c(-1, 1))
netVisual_aggregate(
  cellchat_TNFRSF13C_NMOSD,
  signaling = "BAFF", #signaling pathway
  signaling.name='BAFF target B_naive',
  sources.use =sources.use,
  targets.use = targets.use,
  layout = "circle", 
  vertex.size.max = 6, 
  edge.width.max = 10,
  remove.isolate=TRUE,
  title.space = 1
)

# ===【Hierarchy plot - BAFF in Bnaive】===
levels(cellchat_TNFRSF13C_NMOSD@idents) #[4] "B naive" 
par(mfrow = c(1, 1))
netVisual_aggregate(
  cellchat_TNFRSF13C_NMOSD, 
  signaling ="BAFF",
  layout = "hierarchy",
  vertex.receiver = seq(3,4), #設定target為"B memory",B naive"
  pt.title = 8, #整體大小
  edge.width.max = 12
)

# ===【Heatmap - BAFF 】===
par(mfrow=c(1,1)) 
ht1<-netVisual_heatmap(cellchat_TNFRSF13C_NMOSD, signaling = "BAFF",ylim.top=c(0,0.5),ylim.right=c(0,0.15),
                       color.heatmap = c('#ffecec','#ff4040'),
                       measure="weight",font.size = 9,col.show =c('B intermediate','B memory','B naive') )
draw(ht1, padding = unit(c(5, 25, 5, 15), "mm"))  #(上,左,下,右) 

# ===【Heatmap - BAFF in Bnaive】===
par(mfrow=c(1,1)) 
netVisual_heatmap(cellchat_TNFRSF13C_NMOSD, signaling = "BAFF",
                  measure="weight",font.size = 9,
                  targets.use = targets.use,
                  title.name='BAFF target B_naive')

#----------------------【BAFF NMOSD vs Control】----------------------

#----------------🔻{CellChat - Control}🔻------------------
seurat_control <- subset(merge_data, subset = Condition == "Control")
cellchat_control <- createCellChat(
  object  = seurat_control,      # Seurat object
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
cellchat_control@DB<-CellChatDB.use
cellchat_control <- subsetData(cellchat_control) # 資料過濾

# === 【Identify over-expressed "ligands" or "receptors" in each cell group】===
#devtools::install_github("immunogenomics/presto")
library(presto) #讓 Wilcoxon 在 dgCSparse 下跑
cellchat_control <- identifyOverExpressedGenes(cellchat_control)

# identify over-expressed "L–R interactions"
cellchat_control<- identifyOverExpressedInteractions(cellchat_control)

### Each cell types lognormalized average
computeAveExpr(cellchat_control, features = c("TNFRSF13C"),  
               type = "triMean") 

# ===【 Infer cell–cell communication at a L–R pair level】===
cellchat_control <- computeCommunProb(cellchat_control, type = "triMean", trim = NULL, 
                                    raw.use = TRUE,nboot = 20) #約20mins
# raw.use = TRUE:用原始資料(lognormalized)做計算
# nboot:跑幾次permutation決定pvalue (default=100)

# ===【Filter fewer cells】===
cellchat_control <- filterCommunication(cellchat_control, min.cells = 10) #某celltype少於10cells->drop #default=10

# ===【 Infer cell–cell communication at a signaling pathway level】===
cellchat_control<- computeCommunProbPathway(cellchat_control,thresh = 0.05) #thresh:significant interaction P-value criteria
#saveRDS(cellchat_control, file = "../scRNA_DATA/My_merged_protein_coding_genes/cellchat_control(Azimuth).rds")

# ===【Import CellChat_control RDS】===
cellchat_control <- readRDS("../scRNA_DATA/My_merged_protein_coding_genes/cellchat_control.rds")


# ===【TNFRSF13C】===
lr.idx_control <- grep("TNFRSF13C", cellchat_control@LR$LRsig$interaction_name)
lr.use_control <- cellchat_control@LR$LRsig$interaction_name[lr.idx_control]
lr.use_control #"TNFSF13B_TNFRSF13C"
lr_TNFRSF13C_df_control<- data.frame(interaction_name = lr.use_control)
df_TNFRSF13C_control <- subsetCommunication(cellchat_control,  pairLR.use = lr_TNFRSF13C_df_control)

# ===【 Calculate the aggregated cell–cell communication network】===
# >>> for num of links & community probability
# {BAFF}
sources.use <- unique(df_TNFRSF13C_control $source)
targets.use <-  unique(df_TNFRSF13C_control $target )
cellchat_TNFRSF13C_BAFF_control <- aggregateNet(
  cellchat_control,
  sources.use = sources.use,
  targets.use = targets.use
)

# ===【Heatmap - BAFF 】===
par(mfrow=c(1,1)) 
ht2<-netVisual_heatmap(cellchat_TNFRSF13C_BAFF_control, signaling = "BAFF",ylim.top=c(0,0.5),ylim.right=c(0,0.15),
                       color.heatmap = c('#ffecec','#ffa0a0'),
                       measure="weight",font.size = 9,col.show =c('B intermediate','B memory','B naive') )
draw(ht2, padding = unit(c(5, 25, 5, 15), "mm"))  #(上,左,下,右)

#----------------------------------------------------------------------------------------------

# ===【 Calculate the aggregated cell–cell communication network】===
# >>> for num of links & community probability
cellchat_NMOSD <- aggregateNet(cellchat_NMOSD)
cellchat_control <- aggregateNet(cellchat_control)

#===【Merge CellChat】===
object.list <- list(NMOSD = cellchat_NMOSD, Control = cellchat_control ) 
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))
#saveRDS(cellchat_merged  , file = "../scRNA_DATA/My_merged_protein_coding_genes/cellchat_merged.rds")

#===【Import CellChat_merge】===
cellchat_merged<-readRDS("../scRNA_DATA/My_merged_protein_coding_genes/cellchat_merged.rds")


#===【Compare num_interaction & total_strength】===
num_interaction<-compareInteractions(cellchat_merged , show.legend = F,  
                                     group = c(1,2))      
num_interaction

total_strength<- compareInteractions(cellchat_merged, show.legend = F,  
                                           group = c(1,2), measure = "weight") 
total_strength

#---🔻{Heatmap: NMOSD vs Control}🔻---
# ===【NMOSD: Heatmap - BAFF 】===
par(mfrow=c(1,1)) 
ht1<-netVisual_heatmap(cellchat_TNFRSF13C_NMOSD, signaling = "BAFF",ylim.top=c(0,0.5),ylim.right=c(0,0.15),
                      color.heatmap = c('#ffecec','#ff4040'),
                      measure="weight",font.size = 9,col.show =c('B intermediate','B memory','B naive') )
draw(ht1, padding = unit(c(5, 25, 5, 15), "mm"))  #(上,左,下,右) 


# ===【Control: Heatmap - BAFF 】===
par(mfrow=c(1,1)) 
ht2<-netVisual_heatmap(cellchat_TNFRSF13C_BAFF_control, signaling = "BAFF",ylim.top=c(0,0.5),ylim.right=c(0,0.15),
                  color.heatmap = c('#ffecec','#ffa0a0'),
                  measure="weight",font.size = 9,col.show =c('B intermediate','B memory','B naive') )
draw(ht2, padding = unit(c(5, 25, 5, 15), "mm"))  #(上,左,下,右)



