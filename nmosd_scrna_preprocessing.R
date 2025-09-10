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


################【NMOSD scRNA data preprocessing】############
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)

#---{取出 raw counts}
raw_count <- GetAssayData(data, assay = "RNA", slot = "counts")
print(dim(raw_count))   # gene x cell

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

#---{sample欄位}
sample_col<-"orig.ident"

#QC table
qc_df <- data.frame(
  sample       = data@meta.data[[sample_col]],
  nCount_RNA   = data$nCount_RNA,
  nFeature_RNA = data$nFeature_RNA,
  percent.mt   = data$percent.mt,
  percent.hb   = data$percent.hb,
  row.names    = colnames(data)
)

#Each ample calculate 99th percentile（UMI、Features）
p99_table <- qc_df %>%
  group_by(sample) %>%
  summarise(
    p99_umi     = quantile(nCount_RNA,   0.99, na.rm = TRUE),
    p99_feature = quantile(nFeature_RNA, 0.99, na.rm = TRUE),
    .groups = "drop"
  )

cat("Per-sample 99th percentiles:\n")
print(p99_table)

#Each sample filter criteria 
qc_df <- qc_df %>%
  left_join(p99_table, by = "sample") %>%
  mutate(
    flag_lowUMI   = nCount_RNA   < umi_min,
    flag_hiUMI    = nCount_RNA   > p99_umi,
    flag_lowFeat  = nFeature_RNA < feature_min,
    flag_hiFeat   = nFeature_RNA > p99_feature,
    flag_hiMT     = percent.mt   > mitochondrial_max,
    flag_hiHB     = percent.hb   > hemoglobin_max,
    flag_any      = flag_lowUMI | flag_hiUMI | flag_lowFeat | flag_hiFeat | flag_hiMT | flag_hiHB
  )

# 把cell barcodes放在 rownames(確保之後subset data對得上)
rownames(qc_df) <- colnames(data)

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
    total                 = n(),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(`UMI < 1000`:`Any of above`),
               names_to = "criterion", values_to = "n_cells") %>%
  mutate(percent = round(100 * n_cells / total, 3)) %>%
  arrange(sample, match(criterion, c("UMI < 1000","UMI > 99th pct",
                                     "Features < 200","Features > 99th pct",
                                     "Mito % > 12","Hemoglobin % > 10",
                                     "Any of above"))) %>%
  select(sample, criterion, n_cells, percent)

#Print sample table
invisible(lapply(split(overall_by_sample_long, overall_by_sample_long$sample), function(df) {
  cat("\n=== Sample:", unique(df$sample), "===\n")
  print(df[, c("criterion","n_cells","percent")], row.names = FALSE)
}))

#----{Filter cell}
# 保留通過 QC 的細胞
cells_to_keep <- rownames(qc_df)[!qc_df$flag_any]

# 建立新的 Seurat object(符合 paper criteria)
data_filtered <- subset(data, cells = cells_to_keep)

# Filter counts matrix
qc_count <- GetAssayData(data_filtered, assay = "RNA", slot = "counts")

# Filter cell counts(Before vs After)
cat("Before filtering:", ncol(data), "cells\n")
cat("After filtering:", ncol(data_filtered), "cells\n")
