#install.packages("SeuratObject")
#install.packages("Seurat")
library(SeuratObject)
library(Seurat)
data=readRDS("../scRNA_DATA/nmo25_pbmc.rds")
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


################【NMOSD scRNA data】#############
# 0) 需要的套件
library(Seurat)
library(Matrix)

#---{取出 raw counts}
raw_count <- GetAssayData(data, assay = "RNA", slot = "counts")  # dgCMatrix 稀疏矩陣
print(dim(raw_count))     # 26135 x 195817 (gene x cell)

library(Seurat)
library(dplyr)

# ----計算每個cell粒線體基因(以MT開頭)的reads占總reads比例 
if (!"percent.mt" %in% colnames(data@meta.data)) {
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
}
#----計算每個 cell 的血紅素基因（HBA1, HBA2, HBB）占總 reads 的比例
if (!"percent.hb" %in% colnames(data@meta.data)) {
  data[["percent.hb"]] <- PercentageFeatureSet(data, pattern = "^HB[AB][12]")
}

#----paper criteria
umi_min <- 1000
feature_min <- 200
mitochondrial_max  <- 12
hemoglobin_max  <- 10

# 99 百分位
p99_umi  <- quantile(data$nCount_RNA,   0.99, na.rm = TRUE)
p99_feature <- quantile(data$nFeature_RNA, 0.99, na.rm = TRUE)

cat(sprintf("99th percentile: UMI=%.0f, Features=%.0f\n", p99_umi, p99_feature))
#99th percentile: UMI=11798, Features=3240

#----計算不符合paper criteria的cell數
flag_lowUMI   <- data$nCount_RNA    < umi_min
flag_hiUMI    <- data$nCount_RNA    > p99_umi
flag_lowFeat  <- data$nFeature_RNA  < feature_min
flag_hiFeat   <- data$nFeature_RNA  > p99_feature
flag_hiMT     <- data$percent.mt    > mitochondrial_max
flag_hiHB     <- data$percent.hb    > hemoglobin_max 

flag_any <- flag_lowUMI | flag_hiUMI | flag_lowFeat | flag_hiFeat | flag_hiMT | flag_hiHB

#----計算結果
overall <- data.frame(
  criterion = c("UMI < 1000", "UMI > 99th pct", "Features < 200",
                "Features > 99th pct", "Mito % > 12", "Hemoglobin % > 10",
                "Any of above"),
  n_cells   = c(sum(flag_lowUMI), sum(flag_hiUMI), sum(flag_lowFeat),
                sum(flag_hiFeat), sum(flag_hiMT), sum(flag_hiHB),
                sum(flag_any))
)
overall$percent <- round(100 * overall$n_cells / ncol(data), 3)
print(overall, row.names = FALSE)
