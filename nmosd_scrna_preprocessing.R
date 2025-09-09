setwd("C:/Users/JANE/Desktop/Wang實驗室/NMOSD研究計畫/scRNA")

install.packages("SeuratObject")
install.packages("Seurat")
library(SeuratObject)
library(Seurat)
data=readRDS("./nmo25_pbmc.rds")
head(data@assays$RNA@counts)

colnames(data@meta.data)
head(data@meta.data)
########################【active.ident】################################
# 轉換為 DataFrame
cell_identities <- data.frame(
  Cell_ID = names(data@active.ident),  # 細胞名稱
  Identity = as.character(data@active.ident)  # 細胞分類標籤
)

write.csv(cell_identities, file = "cell_identities.csv", row.names = FALSE)
#-----------------------------------------------------------------------------

#E-GEAD-551
control_data<-readRDS("./E-GEAD-551_PBMC_scRNAseq.rds")
head(control_data@assays$RNA@counts)
control_cell_identities <- data.frame(
  Cell_ID = names(control_data@active.ident),  # 細胞名稱
  Identity = as.character(control_data@active.ident)  # 細胞分類標籤
)

write.csv(control_cell_identities , file = "cell_identities_control.csv", row.names = FALSE)

#取健康人的gene cell matrix
healthy_cell_identities<-control_cell_identities[grepl("^HC",control_cell_identities$Cell_ID),]


################【NMOSD scRNA data輸出】#############
# 0) 需要的套件
library(Seurat)
library(Matrix)

#---{取出 raw counts}
raw_count <- GetAssayData(data, assay = "RNA", slot = "counts")  # dgCMatrix 稀疏矩陣
print(dim(raw_count))     # 26135 x 195817 (gene x cell)

# 3) 10X 需要基因名稱唯一；若重複就自動去重
if (any(duplicated(rownames(raw_count)))) {
  message("Detected duplicated gene IDs; making them unique with suffixes")
  rownames(raw_count) <- make.unique(rownames(raw_count))
}

# （可選）如果你想去掉 Ensembl 版本號：ENSG000001.12 -> ENSG000001
# rownames(raw) <- sub("\\.\\d+$", "", rownames(raw))

# 4) 建立輸出資料夾
dir.create("10x_raw", showWarnings = FALSE)

# 5) 寫出 matrix.mtx.gz
Matrix::writeMM(raw, file = gzfile("10x_raw/matrix.mtx.gz"))

# 6) 寫出 features.tsv.gz（3 欄：gene_id, gene_name, feature_type）
features <- data.frame(
  gene_id     = rownames(raw),
  gene_name   = rownames(raw),         # 若你有 symbol 對照可放在這欄
  feature_type= "Gene Expression",
  check.names = FALSE
)
write.table(features, file = gzfile("10x_raw/features.tsv.gz"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 7) 寫出 barcodes.tsv.gz
write.table(colnames(raw), file = gzfile("10x_raw/barcodes.tsv.gz"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 8) （強烈建議）輸出每個細胞對應的樣本 ID，方便 Scrublet 逐樣本跑
meta_min <- data.frame(
  barcode = colnames(data),
  sample  = data$orig.ident,
  row.names = NULL, check.names = FALSE
)
write.csv(meta_min, "10x_raw/barcode_sample_map.csv", row.names = FALSE)

message("✅ Done. Files are in ./10x_raw/")

