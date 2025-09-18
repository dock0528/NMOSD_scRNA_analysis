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

#write.csv(cell_identities, file = "cell_identities.csv", row.names = FALSE)
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
raw_count <- GetAssayData(data, assay = "RNA", slot = "counts")#從RNA提取counts
print(dim(raw_count))   # gene x cell: 26135 x 195817

# ----計算每個cell粒線體基因(以MT開頭)的reads占總reads比例 
#單位:%
if (!"percent.mt" %in% colnames(data@meta.data)) {
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
}
#----計算每個 cell 的血紅素基因（HBA1, HBA2, HBB）占總 reads 的比例
#單位:%
if (!"percent.hb" %in% colnames(data@meta.data)) {
  data[["percent.hb"]] <- PercentageFeatureSet(data, pattern = "^HB[AB][12]")
}    #^從字首匹配HB 、第3個字A or B、第4個字1 or 2

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
    p99_umi     = quantile(nCount_RNA,   0.99, na.rm = TRUE), #99th 的cell UMIs
    p99_feature = quantile(nFeature_RNA, 0.99, na.rm = TRUE), #99th 的genes
    .groups = "drop" #不要保留分組形式
  )

cat("Per sample 99th percentiles:\n")
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
rownames(qc_df) <- colnames(data) #每個sample_cell名

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
data_filtered <- subset(data, cells = cells_to_keep)

# Filter counts matrix
qc_count <- GetAssayData(data_filtered, assay = "RNA", slot = "counts")

# Filter cell counts(Before vs After)
cat("Before filtering:", ncol(data), "cells\n")
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
  project      = "NMOSD"
)
DefaultAssay(data_filtered_clean) <- "RNA"

#驗證是否含counts
Layers(data_filtered_clean[["RNA"]])   #counts
dim(GetAssayData(data_filtered_clean, layer = "counts"))   # genes x cells:26135 x 193384

# 讓 data layer = counts（關鍵一步）
data_filtered_clean <-SetAssayData(
  object   = data_filtered_clean,
  assay    = "RNA",
  layer    = "data",
  new.data = GetAssayData(data_filtered_clean, assay = "RNA", layer = "counts")
)
Layers(data_filtered_clean[["RNA"]]) #"counts" "data" 
dim(GetAssayData(data_filtered_clean, layer="data"))  # genes x cells:26135 x 193384

# RDS 格式
saveRDS(data_filtered_clean , "../scRNA_DATA/NMOSD_count_metadata.rds")


################{存成.mtx}##################
library(Matrix)

# 取 raw counts
m <- GetAssayData(data_filtered_clean, assay = "RNA", layer = "counts")

# OUTPUT folder
outdir <- "../scRNA_DATA/mtx_nmosd"
dir.create(outdir, showWarnings = FALSE)

#matrix.mtx
Matrix::writeMM(m, file.path(outdir, "NMOSD_matrix.mtx"))

#features.tsv (基因名)
write.table(
  data.frame(rownames(m)),
  file = file.path(outdir, "NMOSD_features.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#barcodes.tsv (細胞ID)
write.table(
  data.frame(colnames(m)),
  file = file.path(outdir, "NMOSD_barcodes.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#metadata
write.csv(
  data_filtered_clean@meta.data,
  file = file.path(outdir, "NMOSD_obs.csv")
)

####################【"原始matrix"畫mitochondrial & Hemoglobin genes】##################
library(Seurat)
library(ggplot2)
library(reshape2)


#----{nFeature RNA}
p1 <- ggplot(qc_df, aes(x = sample, y = nFeature_RNA, fill = sample)) +
  geom_violin(trim = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
  ylab("nFeature_RNA") + xlab("")+
  ggtitle("nFeature RNA")

print(p1)


# nCount_RNA
p2 <- ggplot(qc_df, aes(x = sample, y = nCount_RNA, fill = sample)) +
  geom_violin(trim = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
  ylab("nCount_RNA") + xlab("")+
  ggtitle("nCount RNA")

print(p2)

# percent.mt
p3 <- ggplot(qc_df, aes(x = sample, y = percent.mt, fill = sample)) +
  geom_violin(trim = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
  ylab("percent.mt") + xlab("")+
  ggtitle("Percentage of Mitochondrial genes")

print(p3)

# percent.hb
p4 <- ggplot(qc_df, aes(x = sample, y = percent.hb, fill = sample)) +
  geom_violin(trim = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
  ylab("percent.hb") + xlab("")+
  ggtitle("Percentage of Hemoglobin genes")

print(p4)

####################【"過濾後matrix"畫mitochondrial & Hemoglobin genes】##################
library(Seurat)
library(ggplot2)
library(reshape2)
filtered_data=readRDS("../scRNA_DATA/NMOSD_count_metadata.rds")

# ----計算每個cell粒線體基因(以MT開頭)的reads占總reads比例 
#單位:%
if (!"percent.mt" %in% colnames(filtered_data@meta.data)) {
  filtered_data[["percent.mt"]] <- PercentageFeatureSet(filtered_data, pattern = "^MT-")
}
#----計算每個 cell 的血紅素基因（HBA1, HBA2, HBB）占總 reads 的比例
#單位:%
if (!"percent.hb" %in% colnames(filtered_data@meta.data)) {
  filtered_data[["percent.hb"]] <- PercentageFeatureSet(filtered_data, pattern = "^HB[AB][12]")
}    #^從字首匹配HB 、第3個字A or B、第4個字1 or 2


#---{sample欄位}
sample_col<-"orig.ident"

#QC table
qc_df_filtered <- data.frame(
  sample       = filtered_data@meta.data[[sample_col]],
  nCount_RNA   = filtered_data$nCount_RNA,
  nFeature_RNA = filtered_data$nFeature_RNA,
  percent.mt   = filtered_data$percent.mt,
  percent.hb   = filtered_data$percent.hb,
  row.names    = colnames(filtered_data)
)

#----{nFeature RNA}
p1 <- ggplot(qc_df_filtered, aes(x = sample, y = nFeature_RNA, fill = sample)) +
  geom_violin(trim = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
  ylab("nFeature_RNA") + xlab("")+
  ggtitle("nFeature RNA")

print(p1)


# nCount_RNA
p2 <- ggplot(qc_df_filtered, aes(x = sample, y = nCount_RNA, fill = sample)) +
  geom_violin(trim = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
  ylab("nCount_RNA") + xlab("")+
  ggtitle("nCount RNA")

print(p2)

# percent.mt
p3 <- ggplot(qc_df_filtered, aes(x = sample, y = percent.mt, fill = sample)) +
  geom_violin(trim = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
  ylab("percent.mt") + xlab("")+
  ggtitle("Percentage of Mitochondrial genes")

print(p3)

# percent.hb
p4 <- ggplot(qc_df_filtered, aes(x = sample, y = percent.hb, fill = sample)) +
  geom_violin(trim = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
  ylab("percent.hb") + xlab("")+
  ggtitle("Percentage of Hemoglobin genes")

print(p4)

#####################【處理gene一致的問題】#################
filtered_data=readRDS("../scRNA_DATA/NMOSD_count_metadata.rds")


# 取出 counts
counts <- filtered_data@assays$RNA@layers$counts

# 去掉版本號 (例如 AL627309.1 → AL627309)
gene_ids <- rownames(filtered_data@assays$RNA) #length(gene_ids):26135
gene_ids_clean <- gsub("\\..*", "", gene_ids)
gene_ids_clean <- as.factor(gene_ids_clean)

# 合併同一基因的 counts（保持稀疏格式）
#BiocManager::install("scran")
library(scran)
counts_clean <- sumCountsAcrossFeatures(counts, ids =gene_ids_clean )
class(counts_clean)
dim(counts_clean) #24092 X 193384

#加回gene & cell名
rownames(counts_clean) <- unique(gene_ids_clean)   
colnames(counts_clean) <- colnames(filtered_data)

counts_clean <- as(counts_clean, "dgCMatrix") #轉成"稀疏"矩陣
class(counts_clean)

# 取metadata
metadata_filtered<-filtered_data@meta.data
#dim(metadata_filtered)

data_filtered_clean <- CreateSeuratObject(
  counts       = counts_clean,
  meta.data    = metadata_filtered,
  assay        = "RNA",
  min.features = 0,  # 不因 features 太少丟 cell
  min.cells    = 0, # 不因出現細胞數太少丟基因
  project      = "NMOSD"
)
DefaultAssay(data_filtered_clean) <- "RNA"

#驗證是否含counts
Layers(data_filtered_clean[["RNA"]])   #counts
dim(GetAssayData(data_filtered_clean, layer = "counts"))   # genes x cells:24092 x 193384

# 讓 data layer = counts（重要!!!)
data_filtered_clean <-SetAssayData(
  object   = data_filtered_clean,
  assay    = "RNA",
  layer    = "data",
  new.data = GetAssayData(data_filtered_clean, assay = "RNA", layer = "counts")
)
Layers(data_filtered_clean[["RNA"]]) #"counts" "data" 
dim(GetAssayData(data_filtered_clean, layer="data"))  # genes x cells:24092 x 193384

# ----{存成RDS 格式}
saveRDS(data_filtered_clean , "../scRNA_DATA/NMOSD_count_metadata(gene unique).rds")


#----------------{存成.mtx}----------------
library(Matrix)

# 取 raw counts
m <- GetAssayData(data_filtered_clean, assay = "RNA", layer = "counts")

# output folder
outdir <- "../scRNA_DATA/mtx_nmosd"
dir.create(outdir, showWarnings = FALSE)

#matrix.mtx
Matrix::writeMM(m, file.path(outdir, "NMOSD_matrix(gene unique).mtx")) #回傳NULL代表有成功

#features.tsv (基因名)
write.table(
  data.frame(rownames(m)),
  file = file.path(outdir, "NMOSD_features(gene unique).tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#barcodes.tsv (細胞ID)
write.table(
  data.frame(colnames(m)),
  file = file.path(outdir, "NMOSD_barcodes(gene unique).tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#metadata
write.csv(
  data_filtered_clean@meta.data,
  file = file.path(outdir, "NMOSD_obs(gene unique).csv")
)
#-------------------------------------------
#----{確認是否有重複的Hugo Symbol}
nmosd_data<-readRDS("../scRNA_DATA/NMOSD_count_metadata(gene unique).rds")
#View(nmosd_data)
anyDuplicated(rownames(nmosd_data@assays$RNA@features)) #0

########################【篩protein coding gene】########################
Raw_count_merged<-read.csv('./RNA_DATA/Raw_count_merged_matrix.csv',row.names = 1,header=T) #dim(Raw_count_merged):78724 x 21

#Raw count 
count_df<-as.matrix(Raw_count_merged)
rownames(count_df)<-gsub("\\.\\d+$", "",rownames(count_df)) #去除小數點以後的值

#----連到Ensembl
library(biomaRt)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

#----My gene matrix
All_genes<-rownames(count_df)


#----抓註解
annot <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters    = "ensembl_gene_id",
  values     = All_genes,
  mart       = ensembl
)

#----篩 protein coding genes
protein_coding_genes <- annot[annot$gene_biotype == "protein_coding", ] #20091個
my_protein_coding_genes<-protein_coding_genes$ensembl_gene_id
#write.csv(my_protein_coding_genes,file='./RNA_DATA/My_protein_coding_genes.csv',row.names = FALSE)

#----proteing coding df
protein_coding_count_df<-count_df[rownames(count_df) %in% protein_coding_genes$ensembl_gene_id,]