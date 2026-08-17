##################################【SoupX】##################################
#載入套件
install.packages("SoupX")
BiocManager::install("DropletUtils")
library(SoupX) #v1.6.2
library(DropletUtils)

#讀入所有樣本
ids <- list.dirs("../SoupX_cellranger_results_merged_v2", recursive = FALSE, full.names = FALSE) 
#recursive = FALSE:僅列出第一層資料夾 #full.names = FALSE:僅列出資料夾名稱(不要完整路徑)
base_dir <- "../SoupX_cellranger_results_merged_v2"   
outputdir<-"../SoupX_results"

for (id in ids) {
  message("Soupx正在處理:", id, " ...")
  
  #Read cellranger output
  sc <- load10X(file.path(base_dir, id, "outs"))
  
  #自動估算ambient RNA 
  soup <- autoEstCont(sc)
  
  #產出校正matrix
  adj.matrix <- adjustCounts(soup, roundToInt = TRUE)
  
  #設定feature name=ENSEMBL ID
  features <- read.table(
    file.path(base_dir, id, "outs","filtered_feature_bc_matrix","features.tsv.gz"),
    sep = "\t", header = FALSE, stringsAsFactors = FALSE
  )
  
  rownames(adj.matrix) <- features$V1
  
  #Output
  out_dir <- file.path(outputdir, id, "SoupX_corrected_outs")
  write10xCounts(out_dir, adj.matrix) 
  message("完成:", id)
}



##################################【所有sample合併】##################################
library(Seurat)

ids <- list.dirs("../SoupX_results", recursive = FALSE, full.names = FALSE) 
base_dir <- "../SoupX_results"

seurat_list <- list()

for (id in ids) {
  
  #Read SoupX ouput files
  data <- Read10X(file.path(base_dir, id, "SoupX_corrected_outs"), gene.column = 1)
  #gene.column = 1 把ENSG當作feature name
  
  
  # Create Seurat object
  obj <- CreateSeuratObject(data, project = id)
  
  # 加上 sample prefix 
  obj <- RenameCells(obj, add.cell.id = id)
  
  seurat_list[[id]] <- obj
}

# 合併
merged <- Reduce(function(x, y) merge(x, y), seurat_list) #Reduce()對所有seurat_list做merge

#合併多層layers(保留 cell names)
merged <- JoinLayers(merged, assay = "RNA")

#確認cell barcodes & features
colnames(merged) #cell barcodes
rownames(merged) #features(ENSG)
dim(merged) #38606 x 27909

# 存成 RDS
saveRDS(merged , "../scRNA_DATA_3control/Control_merged_SoupX_3samples.rds")


##################################【讀入Control_merged_soupx.rds】##################################
library(SeuratObject)
library(Seurat)
library(Matrix)
control_data=readRDS("../scRNA_DATA_3control/Control_merged_SoupX_3samples.rds")


# 取 raw counts
m <- GetAssayData(control_data, assay = "RNA", layer = "counts")

# OUTPUT folder
outdir <- "../scRNA_DATA_3control"
dir.create(outdir, showWarnings = FALSE)

#matrix.mtx
Matrix::writeMM(m, file.path(outdir, "control_matrix.mtx"))

#features.tsv (基因名)
write.table(
  data.frame(rownames(m)),
  file = file.path(outdir, "control_features.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#barcodes.tsv (細胞ID)
write.table(
  data.frame(colnames(m)),
  file = file.path(outdir, "control_barcodes.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#metadata
write.csv(
  control_data@meta.data,
  file = file.path(outdir, "control_obs.csv")
)
