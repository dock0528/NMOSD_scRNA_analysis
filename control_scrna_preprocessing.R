#install.packages("SeuratObject")
#install.packages("Seurat")
library(SeuratObject)
library(Seurat)
control_data=readRDS("../scRNA_DATA/E-GEAD-551_PBMC_scRNAseq.rds")
head(control_data@assays$RNA@counts)

colnames(control_data@meta.data)
head(control_data@meta.data)

################{原始data存成.mtx}##################
library(Matrix)

# 取 raw counts
m <- GetAssayData(control_data, assay = "RNA", layer = "counts")

# OUTPUT folder
outdir <- "../scRNA_DATA/mtx_control"
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

#----------------------------------------------------

#----{取 Healthy Control}
# 只保留 Condition == "Healthy" 的細胞 (Healthy Control)
healthy_data <- subset(control_data, subset = Condition == "Healthy")

# ----{存成RDS 格式}
#saveRDS(healthy_data  , "../scRNA_DATA/E-GEAD-551_PBMC_scRNAseq(Healthy Control).rds")

################{Healthy Control data存成.mtx}##################
library(Matrix)

# 取 raw counts
m_hc <- GetAssayData(healthy_data , assay = "RNA", layer = "counts")

# OUTPUT folder
outdir <- "../scRNA_DATA/mtx_control"
dir.create(outdir, showWarnings = FALSE)

#matrix.mtx
Matrix::writeMM(m_hc, file.path(outdir, "Healthy_Control_matrix.mtx"))

#features.tsv (基因名)
write.table(
  data.frame(rownames(m_hc)),
  file = file.path(outdir, "Healthy_Control_features.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#barcodes.tsv (細胞ID)
write.table(
  data.frame(colnames(m_hc)),
  file = file.path(outdir, "Healthy_Control_barcodes.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#metadata
write.csv(
  healthy_data@meta.data,
  file = file.path(outdir, "Healthy_Control_obs.csv")
)


################【Healthy Control scRNA data preprocessing】############
library(Seurat)
library(SeuratObject)
library(Matrix)
library(dplyr)
library(tidyr)
data =readRDS("../scRNA_DATA/E-GEAD-551_PBMC_scRNAseq(Healthy Control).rds")
doublet_cells<-read.csv("../Doublet_cells/E_GEAD_551_doublet_cells.txt",header=F)[[1]] #81
scrublet_data <- subset(data, cells = setdiff(Cells(data), doublet_cells))

#---{取出 raw counts}
raw_count <- GetAssayData(scrublet_data , assay = "RNA", slot = "counts")#從RNA提取counts
print(dim(raw_count))   # gene x cell: 24541 x 39804

# ----計算每個cell粒線體基因(以MT開頭)的reads占總reads比例 
#單位:%
if (!"percent.mt" %in% colnames(scrublet_data @meta.data)) {
  scrublet_data [["percent.mt"]] <- PercentageFeatureSet(scrublet_data , pattern = "^MT-")
}
#----計算每個 cell 的血紅素基因（HBA1, HBA2, HBB）占總 reads 的比例
#單位:%
if (!"percent.hb" %in% colnames(scrublet_data @meta.data)) {
  scrublet_data [["percent.hb"]] <- PercentageFeatureSet(scrublet_data , pattern = "^HB[AB][12]")
}    #^從字首匹配HB 、第3個字A or B、第4個字1 or 2

#----paper criteria
umi_min <- 1000
feature_min <- 200
mitochondrial_max  <- 12
hemoglobin_max  <- 10

#---{sample欄位}
head(scrublet_data @meta.data)
sample_col<-"ID"

#QC table
qc_df <- data.frame(
  sample       = scrublet_data @meta.data[[sample_col]],
  nCount_RNA   = scrublet_data $nCount_RNA,
  nFeature_RNA = scrublet_data $nFeature_RNA,
  percent.mt   = scrublet_data $percent.mt,
  percent.hb   = scrublet_data $percent.hb,
  row.names    = colnames(scrublet_data )
)

#Each sample calculate 99th percentile（UMI、Features）
p99_table <- qc_df %>%
  group_by(sample) %>%
  summarise(
    p99_umi     = quantile(nCount_RNA,   0.99, na.rm = TRUE), #99th 的cell UMIs
    p99_feature = quantile(nFeature_RNA, 0.99, na.rm = TRUE), #99th 的genes
    .groups = "drop" #不要保留分組形式
  )

print(n=25,p99_table)

#Each sample filter criteria 
qc_df <- qc_df %>%
  left_join(p99_table, by = "sample") %>% #依sample把p99_table和qc_df合併
  mutate( #加上 QC 的標記欄位(TRUE / FALSE)
    flag_lowUMI   = nCount_RNA   < umi_min,
    flag_hiUMI    = nCount_RNA   > p99_umi,
    flag_lowFeat  = nFeature_RNA < feature_min,
    flag_hiFeat   = nFeature_RNA > p99_feature,
    flag_hiMT     = percent.mt   > mitochondrial_max,
    flag_hiHB     = percent.hb   > hemoglobin_max,
    flag_any      = flag_lowUMI | flag_hiUMI | flag_lowFeat | flag_hiFeat | flag_hiMT | flag_hiHB 
    #if任何一個flag是TRUE，就把這個cell標記為TRUE
  )

# 把cell barcodes放在 rownames(確保之後subset data對得上)
rownames(qc_df) <- colnames(scrublet_data ) #每個sample_cell名

#[Each sample table] Row:criterion & Col:n_cells / percent
overall_by_sample_long <- qc_df %>%
  group_by(sample) %>%
  summarise(
    `UMI < 1000`          = sum(flag_lowUMI),
    `UMI > 99th pct`      = sum(flag_hiUMI),
    `Features < 200`      = sum(flag_lowFeat),
    `Features > 99th pct` = sum(flag_hiFeat),
    `Mito % > 12`         = sum(flag_hiMT),
    `Hemoglobin % > 10`   = sum(flag_hiHB),
    `Any of above`        = sum(flag_any),
    total                 = n(), #樣本中的 cell 總數
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(`UMI < 1000`:`Any of above`), #轉成長格式
               names_to = "criterion", values_to = "cell_counts") %>%
  mutate(percent = round(100 * cell_counts / total, 3)) %>%
  arrange(sample, match(criterion, c("UMI < 1000","UMI > 99th pct",
                                     "Features < 200","Features > 99th pct",
                                     "Mito % > 12","Hemoglobin % > 10",
                                     "Any of above"))) %>%
  #arrange排序
  select(sample, criterion, cell_counts, percent) #要保留的欄位

#Print sample table
invisible(lapply(split(overall_by_sample_long, overall_by_sample_long$sample), function(df) {
  #按照 sample分割成一個 list，每個 list 元素是一個樣本的 dataframe
  cat("\n=============== Sample:", unique(df$sample), "===============\n")
  print(df[, c("criterion","cell_counts","percent")], row.names = FALSE)
}))

#----{Filter cell}
# 保留通過 QC 的細胞
cells_to_keep <- rownames(qc_df)[!qc_df$flag_any]

# 建立新的 Seurat object(符合 paper criteria)
data_filtered <- subset(scrublet_data , cells = cells_to_keep)

# Filter counts matrix
qc_count <- GetAssayData(data_filtered, assay = "RNA", slot = "counts")

# Filter cell counts(Before vs After)
cat("Before filtering:", ncol(scrublet_data), "cells\n")
cat("After filtering:", ncol(data_filtered), "cells\n")

#----{存raw count matrix + metadata為.rds}
counts_filtered   <- GetAssayData(data_filtered, assay = "RNA", slot = "counts")
counts_filtered<- as(counts_filtered, "dgCMatrix")                      # 稀疏矩陣類型
#dim(counts_filtered)
metadata_filtered<-data_filtered@meta.data
#dim(metadata_filtered)

data_filtered_clean <- CreateSeuratObject(
  counts       = counts_filtered,
  meta.data    = metadata_filtered,
  assay        = "RNA",
  min.features = 0,  # 不因 features 太少丟 cell
  min.cells    = 0, # 不因出現細胞數太少丟基因
  project      = "Control"
)
DefaultAssay(data_filtered_clean) <- "RNA"

#驗證是否含counts
Layers(data_filtered_clean[["RNA"]])   #counts
dim(GetAssayData(data_filtered_clean, layer = "counts"))   # genes x cells:24541 x 39318

# 讓 data layer = counts（關鍵一步）
data_filtered_clean <-SetAssayData(
  object   = data_filtered_clean,
  assay    = "RNA",
  layer    = "data",
  new.data = GetAssayData(data_filtered_clean, assay = "RNA", layer = "counts")
)
Layers(data_filtered_clean[["RNA"]]) #"counts" "data" 
dim(GetAssayData(data_filtered_clean, layer="data"))  # genes x cells:24541 x 39318

# RDS 格式
saveRDS(data_filtered_clean , "../scRNA_DATA/Filter_low_quality_cells_rdsE-GEAD-551_count_metadata.rds")
