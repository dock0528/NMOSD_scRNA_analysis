library(SeuratObject)
library(Seurat)
nmosd_filtered_data<-readRDS("../scRNA_DATA/Filter_low_quality_cells_rds/NMOSD_count_metadata.rds")
control_filtered_data<-readRDS("../scRNA_DATA/Filter_low_quality_cells_rds/SRP349890_count_metadata.rds")

####################{檢查"NMOSD"資料內容}##################
View(nmosd_filtered_data)
head(nmosd_filtered_data@assays$RNA@layers$counts) #dgCMatrix
dim(nmosd_filtered_data@assays$RNA@layers$counts) #26135 x 193159
rownames(nmosd_filtered_data) #Features(HUGO)
colnames(nmosd_filtered_data) #Cells

colnames(nmosd_filtered_data@meta.data) #metadata
head(nmosd_filtered_data@meta.data)

####################{檢查"Control"資料內容}##################
View(control_filtered_data)
head(control_filtered_data@assays$RNA@layers$counts) #dgCMatrix
dim(control_filtered_data@assays$RNA@layers$counts) #38606 x 25053
rownames(control_filtered_data) #Features(ENSG)
colnames(control_filtered_data) #Cells

colnames(control_filtered_data@meta.data) #metadata
head(control_filtered_data@meta.data)

####################{篩選protein coding genes}##################
#----{protein coding gene list}
protein_coding_df<-read.csv("../scRNA_DATA/HUGO_with_ENSG_v32_v44(protein coding).csv")
colnames(protein_coding_df) <- c("HugoSymbol_merged", "ENSG") 
protein_coding_df<-protein_coding_df[!protein_coding_df$ENSG %in% c('ENSG00000269226','ENSG00000284024'),] #16489
protein_coding_ENSG<-protein_coding_df$ENSG #16489
protein_coding_HUGO<-protein_coding_df$HugoSymbol_merged #16489

#----{NMOSD}
hugo_map_ensg_df<-read.csv("../scRNA_DATA/nmosd_v32_HUGO_map_ENSG.csv")
hugo_map_ensg.protein_coding<-hugo_map_ensg_df[hugo_map_ensg_df$ENSG %in% protein_coding_ENSG,] #16506

#----重複的hugo symbol
dup_ensg<-hugo_map_ensg.protein_coding[
  duplicated(hugo_map_ensg.protein_coding$ENSG) | 
    duplicated(hugo_map_ensg.protein_coding$ENSG, fromLast = TRUE),  #顯示所有
]#17個重複

#去重複
dedup_map <- hugo_map_ensg.protein_coding[!duplicated(hugo_map_ensg.protein_coding$ENSG), ] #16489 去重複，若有重複保留第一個出現
merged_protein_coding_df <- merge(dedup_map, protein_coding_df, by = "ENSG", all = TRUE) #16489

nmosd_protein_data <- subset(nmosd_filtered_data, features = merged_protein_coding_df$HugoSymbol)
dim(nmosd_protein_data@assays$RNA@layers$counts) #16489 x 193159 怪怪的...
#merged_protein_coding_df[
#  duplicated(merged_protein_coding_df$HugoSymbol) | 
#    duplicated(merged_protein_coding_df$HugoSymbol, fromLast = TRUE),  #顯示所有] 
#出現4個不一樣的gene名(在原本nmosd是2個)

# create gene mappingsetNames(目標,原本)
gene_mapping_nmosd <- setNames(merged_protein_coding_df$HugoSymbol_merged,merged_protein_coding_df$HugoSymbol)

# 替換v32.HUGO 為 v44.HUGO
rownames(nmosd_protein_data) <- gene_mapping_nmosd[rownames(nmosd_protein_data )]

#----{Control}
control_protein_data <- subset(control_filtered_data, features = merged_protein_coding_df$ENSG)
dim(control_protein_data@assays$RNA@layers$counts) #16489 x 25053

# create gene mapping:setNames(目標,原本)
gene_mapping_control <- setNames(merged_protein_coding_df$HugoSymbol_merged,merged_protein_coding_df$ENSG)

# 替換為v44.HUGO
rownames(control_protein_data) <- gene_mapping_control[rownames(control_protein_data )]

#存rds
saveRDS(nmosd_protein_data , "../scRNA_DATA/My_merged_protein_coding_genes/My_NMOSD_count_metadata(protein_coding).rds")
saveRDS(control_protein_data, "../scRNA_DATA/Merged_protein_coding_genes/SRP349890_count_metadata(protein_coding).rds")

#################{合併nmosd & control rds}#################
nmosd_protein_data<-readRDS("../scRNA_DATA/My_merged_protein_coding_genes/My_NMOSD_count_metadata(protein_coding).rds")
control_protein_data<-readRDS("../scRNA_DATA/My_merged_protein_coding_genes/SRP349890_count_metadata(protein_coding).rds")

#刪除layer中多餘的data
nmosd_protein_data[["RNA"]]@layers$data <- NULL
control_protein_data[["RNA"]]@layers$data <- NULL
Layers(nmosd_protein_data[["RNA"]]) #counts
Layers(control_protein_data[["RNA"]]) #counts

#先補回原本dgCMatrix的gene和cell名
#{NMOSD}
nmosd_genes<-rownames(nmosd_protein_data) #16489
rownames(nmosd_protein_data[["RNA"]]@layers$counts) <- nmosd_genes
colnames(nmosd_protein_data[["RNA"]]@layers$counts) <-colnames(nmosd_protein_data)

#{Control}
idx <- match(nmosd_genes, rownames(control_protein_data))

control_protein_data[["RNA"]]@layers$counts <-
  control_protein_data[["RNA"]]@layers$counts[idx, , drop = FALSE]


#Seurat內層gene和cell名 (重要!!!)
rownames(control_protein_data[["RNA"]]@layers$counts) <- nmosd_genes
colnames(control_protein_data[["RNA"]]@layers$counts) <-colnames(control_protein_data)

# Seurat外層gene和cell名  (重要!!!)
rownames(control_protein_data) <- nmosd_genes

# 檢查
all(rownames(control_protein_data) == rownames(nmosd_protein_data)) #要為TRUE

#----{確認細胞名}
head(colnames(nmosd_protein_data))
head(colnames(control_protein_data))

#----{確認基因名}
head(rownames(nmosd_protein_data))
head(rownames(control_protein_data))

#----{只保留metadata為"orig.ident", "nCount_RNA","nFeature_RNA"}
nmosd_protein_data@meta.data <- nmosd_protein_data@meta.data[, (colnames(nmosd_protein_data@meta.data) %in% c("orig.ident", "nCount_RNA","nFeature_RNA"))]
control_protein_data@meta.data<- control_protein_data@meta.data[, (colnames(control_protein_data@meta.data) %in% c("orig.ident", "nCount_RNA","nFeature_RNA"))]


#合併
merged_protein_data <- merge(nmosd_protein_data, y = control_protein_data, merge.data = TRUE)

#確認layers
Layers(merged_protein_data[["RNA"]]) # "counts.NMOSD","counts.Control"
rownames(merged_protein_data) #Features(ENSG)
colnames(merged_protein_data) #Cells

#---{JoinLayers:合併多層layer}
merged_protein_data[["RNA"]] <- JoinLayers(merged_protein_data[["RNA"]])
Layers(merged_protein_data[["RNA"]])  #   "counts"

#----{補回Seurat內層gene和cell名}
rownames(merged_protein_data[["RNA"]]@layers$counts) <- rownames(merged_protein_data)
colnames(merged_protein_data[["RNA"]]@layers$counts)<-colnames(merged_protein_data)

####################{檢查"合併後"的資料內容}##################
View(merged_protein_data)
head(merged_protein_data@assays$RNA@layers) #dgCMatrix
dim(merged_protein_data@assays$RNA@layers$counts) #16489 x 218212
rownames(merged_protein_data) #Features(HUGO)
colnames(merged_protein_data) #Cells

colnames(merged_protein_data@meta.data) #metadata
head(merged_protein_data@meta.data)
tail(merged_protein_data@meta.data)

#----{畫圖確認Merge前後gene counts正確}
gene <- "NOC2L"  
plot(
  nmosd_protein_data[["RNA"]]@layers$counts[gene, ],
  merged_protein_data[["RNA"]]@layers$counts[gene, colnames(nmosd_protein_data)],
  xlab = "Original counts",
  ylab = "Merged counts",
  main = paste("Gene:", gene)
)
abline(0, 1, col = "red")

######################{存rds}##################
saveRDS(merged_protein_data , "../scRNA_DATA/My_merged_protein_coding_genes/My_merged_count_metadata(protein_coding).rds")

#####################{存成.mtx}###################
library(Matrix)

# 取 raw counts
m <- GetAssayData(merged_protein_data, assay = "RNA", layer = "counts")

# output folder
outdir <- "../scRNA_DATA/My_mtx_merge"
dir.create(outdir, showWarnings = FALSE)

#matrix.mtx
Matrix::writeMM(m, file.path(outdir, "My_merge_matrix(protein coding).mtx")) #回傳NULL代表有成功

#features.tsv (基因名)
write.table(
  data.frame(rownames(m)),
  file = file.path(outdir, "My_merge_features(protein coding).tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#barcodes.tsv (細胞ID)
write.table(
  data.frame(colnames(m)),
  file = file.path(outdir, "My_merge_barcodes(protein coding).tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#metadata
write.csv(
  merged_protein_data@meta.data,
  file = file.path(outdir, "My_merge_obs(protein coding).csv")
)
#-----------------------------------------------------