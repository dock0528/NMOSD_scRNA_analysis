library(SeuratObject)
library(Seurat)
nmosd_filtered_data<-readRDS("../scRNA_DATA/Filter_low_quality_cells_rds/NMOSD_count_metadata.rds")
control_filtered_data<-readRDS("../scRNA_DATA/Filter_low_quality_cells_rds/E-GEAD-551_count_metadata.rds")

####################{檢查"NMOSD"資料內容}##################
View(nmosd_filtered_data)
head(nmosd_filtered_data@assays$RNA@layers$counts) #dgCMatrix
dim(nmosd_filtered_data@assays$RNA@layers$counts) #26135 x 193159
rownames(nmosd_filtered_data) #Features(ENSG)
colnames(nmosd_filtered_data) #Cells

colnames(nmosd_filtered_data@meta.data) #metadata
head(nmosd_filtered_data@meta.data)

####################{檢查"Control"資料內容}##################
View(control_filtered_data)
head(control_filtered_data@assays$RNA@layers$counts) #dgCMatrix
dim(control_filtered_data@assays$RNA@layers$counts) #24541 x 39318
rownames(control_filtered_data) #Features(ENSG)
colnames(control_filtered_data) #Cells

colnames(control_filtered_data@meta.data) #metadata
head(control_filtered_data@meta.data)

####################{篩選protein coding genes}##################
#----{protein coding gene list}
protein_coding_gene_list<-read.csv("../scRNA_DATA/HUGO_with_ENSG_v32(protein coding).csv")$HugoSymbol #15906

#----{NMOSD}
nmosd_protein_data <- subset(nmosd_filtered_data, features = protein_coding_gene_list)
dim(nmosd_protein_data@assays$RNA@layers$counts) #15906 x 193159

#----{Control}
control_protein_data <- subset(control_filtered_data, features = protein_coding_gene_list)
dim(control_protein_data@assays$RNA@layers$counts) #15906 x 193159

#存rds
saveRDS(nmosd_protein_data , "../scRNA_DATA/Merged_protein_coding_genes/NMOSD_count_metadata(protein_coding).rds")
saveRDS(control_protein_data, "../scRNA_DATA/Merged_protein_coding_genes/E-GEAD-551_count_metadata(protein_coding).rds")

#################{合併nmosd & control rds}#################
nmosd_protein_data<-readRDS("../scRNA_DATA/Merged_protein_coding_genes/NMOSD_count_metadata(protein_coding).rds")
control_protein_data<-readRDS("../scRNA_DATA/Merged_protein_coding_genes/E-GEAD-551_count_metadata(protein_coding).rds")

#----{確認細胞名}
head(colnames(nmosd_protein_data))
head(colnames(control_protein_data))

#----{只保留metadata為"orig.ident", "nCount_RNA","nFeature_RNA"}
nmosd_protein_data@meta.data <- nmosd_protein_data@meta.data[, (colnames(nmosd_protein_data@meta.data) %in% c("orig.ident", "nCount_RNA","nFeature_RNA"))]
control_protein_data@meta.data<- control_protein_data@meta.data[, (colnames(control_protein_data@meta.data) %in% c("ID", "nCount_RNA","nFeature_RNA"))]
colnames(control_protein_data@meta.data)[colnames(control_protein_data@meta.data) == "ID"] <- "orig.ident"

#合併
merged_protein_data <- merge(nmosd_protein_data, y = control_protein_data, merge.data = FALSE)

#----{將原本分開的counts和data和}
mat_counts <- cbind(
  merged_protein_data@assays$RNA@layers$counts.NMOSD,
  merged_protein_data@assays$RNA@layers$counts.Control
)

mat_data <- cbind(
  merged_protein_data@assays$RNA@layers$data.NMOSD,
  merged_protein_data@assays$RNA@layers$data.Control
)

# 用 SetAssayData() 寫入 Assay5
merged_protein_data <- SetAssayData(
  merged_protein_data,
  assay = "RNA",
  slot = "counts",
  new.data = mat_counts
)

merged_protein_data <- SetAssayData(
  merged_protein_data,
  assay = "RNA",
  slot = "data",
  new.data = mat_data
)

#刪掉不要的layers
merged_protein_data@assays$RNA@layers <- list()

####################{檢查"合併後"的資料內容}##################
View(merged_protein_data)
head(merged_protein_data@assays$RNA@) #dgCMatrix
dim(control_filtered_data@assays$RNA@layers$counts) #24541 x 39318
rownames(control_filtered_data) #Features(ENSG)
colnames(control_filtered_data) #Cells

colnames(control_filtered_data@meta.data) #metadata
head(control_filtered_data@meta.data)####################{檢查"Control"資料內容}##################
View(control_filtered_data)
head(control_filtered_data@assays$RNA@layers$counts) #dgCMatrix
dim(control_filtered_data@assays$RNA@layers$counts) #24541 x 39318
rownames(control_filtered_data) #Features(ENSG)
colnames(control_filtered_data) #Cells

colnames(control_filtered_data@meta.data) #metadata
head(control_filtered_data@meta.data)