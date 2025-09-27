##################################【SoupX】##################################
#載入套件
install.packages("SoupX")
BiocManager::install("DropletUtils")
library(SoupX)
library(DropletUtils)

#讀入所有樣本
ids <- list.dirs("../CellRanger_results", recursive = FALSE, full.names = FALSE) 
#recursive = FALSE:僅列出第一層資料夾 #full.names = FALSE:僅列出資料夾名稱(不要完整路徑)
base_dir <- "../CellRanger_results"   # 你的母資料夾

for (id in ids) {
  message("Soupx正在處理:", id, " ...")
  
  #Read cellranger output
  sc <- load10X(file.path(base_dir, id, "outs"))
  
  #自動估算ambient RNA 
  soup <- autoEstCont(sc)
  
  #產出校正matrix
  adj.matrix <- adjustCounts(soup, roundToInt = TRUE)
  
  #Output
  out_dir <- file.path(base_dir, id, "SoupX_corrected_outs")
  write10xCounts(out_dir, adj.matrix) 
  message("完成:", id)
}

##################################【所有sample合併】##################################
library(Seurat)

ids <- list.dirs("../SoupX_CellRanger_results", recursive = FALSE, full.names = FALSE) 
base_dir <- "../SoupX_CellRanger_results"

seurat_list <- list()

for (id in ids) {
  
  #Read SoupX ouput files
  data <- Read10X(file.path(base_dir, id, "SoupX_corrected_outs"))
  
  # Create Seurat object
  obj <- CreateSeuratObject(data, project = id)
  
  # 加上 sample prefix 
  obj <- RenameCells(obj, add.cell.id = id)
  
  seurat_list[[id]] <- obj
}

# 合併
merged <- Reduce(function(x, y) merge(x, y), seurat_list) #Reduce()對所有seurat_list做merge

# 取出所sample matrix(原本放在RNA@layers)
layers <- merged[["RNA"]]@layers  

# 把所有 counts 矩陣合併
all_counts <- do.call(cbind, layers)

# 建立一個新的 Seurat 物件，確保只有一個 counts matrix
merged_seurat <- CreateSeuratObject(all_counts, project = "control_merged")

# 存成 RDS
saveRDS(merged_seurat , "../scRNA_DATA_control_raw/Control_merged_soupx.rds")

##################################【讀入Control_merged_soupx.rds】##################################
library(SeuratObject)
library(Seurat)
library(Matrix)
control_data=readRDS("../scRNA_DATA_control_raw/Control_merged_soupx.rds")


# 取 raw counts
m <- GetAssayData(control_data, assay = "RNA", layer = "counts")

# OUTPUT folder
outdir <- "../scRNA_DATA_control_raw/mtx_control_raw"
dir.create(outdir, showWarnings = FALSE)

#matrix.mtx
Matrix::writeMM(m, file.path(outdir, "control_matrix_raw.mtx"))

#features.tsv (基因名)
write.table(
  data.frame(rownames(m)),
  file = file.path(outdir, "control_features_raw.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#barcodes.tsv (細胞ID)
write.table(
  data.frame(colnames(m)),
  file = file.path(outdir, "control_barcodes_raw.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#metadata
write.csv(
  control_data@meta.data,
  file = file.path(outdir, "control_obs_raw.csv")
)