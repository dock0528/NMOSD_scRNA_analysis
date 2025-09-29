#install.packages("SeuratObject")
#install.packages("Seurat")
library(SeuratObject)
library(Seurat)
data=readRDS("../scRNA_DATA/nmo25_pbmc.rds")
head(data@assays$RNA@counts)

colnames(data@meta.data)
head(data@meta.data)

##################【原始.rds存成.mtx】##############
# 取 raw counts
m <- GetAssayData(data, assay = "RNA", layer = "counts")

# OUTPUT folder
outdir <- "../scRNA_DATA/mtx_nmosd"
dir.create(outdir, showWarnings = FALSE)

#matrix.mtx
Matrix::writeMM(m, file.path(outdir, "nmosd_matrix_ori.mtx"))

#features.tsv (基因名)
write.table(
  data.frame(rownames(m)),
  file = file.path(outdir, "nmosd_features_ori.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#barcodes.tsv (細胞ID)
write.table(
  data.frame(colnames(m)),
  file = file.path(outdir, "nmosd_barcodes_ori.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#metadata
write.csv(
  control_data@meta.data,
  file = file.path(outdir, "nmosd_obs_ori.csv")
)

###################【取gene名，轉成Ensembl GENE ID】##################
hugo_genes <- rownames(data@assays$RNA@meta.features)
#write.table(hugo_genes, file = "../scRNA_DATA/HUGO_gene_list.csv", row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)



########################【active.ident】################################
# 轉換為 DataFrame
cell_identities <- data.frame(
  Cell_ID = names(data@active.ident),  # 細胞名稱
  Identity = as.character(data@active.ident)  # 細胞分類標籤
)

#write.csv(cell_identities, file = "cell_identities.csv", row.names = FALSE)
#-----------------------------------------------------------------------------

################【NMOSD scRNA data preprocessing】############
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)

# 建立新的 Seurat object(符合 paper criteria)
doublet_cells<-read.csv("../Doublet_cells/nmosd_doublet_cells.txt",header=F)[[1]]
scrublet_data <- subset(data, cells = setdiff(Cells(data), doublet_cells))

#---{取出 raw counts}
raw_count <- GetAssayData(scrublet_data , assay = "RNA", slot = "counts")#從RNA提取counts
print(dim(raw_count))   # gene x cell: 26135 x 195592

# ----計算每個cell粒線體基因(以MT開頭)的reads占總reads比例 
#單位:%
if (!"percent.mt" %in% colnames(scrublet_data @meta.data)) {
  scrublet_data [["percent.mt"]] <- PercentageFeatureSet(scrublet_data , pattern = "^MT-")
}
#----計算每個 cell 的血紅素基因（HBA1, HBA2, HBB）占總 reads 的比例
#單位:%
if (!"percent.hb" %in% colnames(scrublet_data@meta.data)) {
  scrublet_data[["percent.hb"]] <- PercentageFeatureSet(scrublet_data, pattern = "^HB[AB][12]")
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
  sample       = scrublet_data@meta.data[[sample_col]],
  nCount_RNA   = scrublet_data$nCount_RNA,
  nFeature_RNA = scrublet_data$nFeature_RNA,
  percent.mt   = scrublet_data$percent.mt,
  percent.hb   = scrublet_data$percent.hb,
  row.names    = colnames(scrublet_data)
)

#Each ample calculate 99th percentile（UMI、Features）
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
rownames(qc_df) <- colnames(scrublet_data) #每個sample_cell名

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
data_filtered <- subset(scrublet_data, cells = cells_to_keep)

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

###################【消失的ENSG篩protein coding gene】##############
library(clusterProfiler)
HUGO_gene_missed_list<-read.csv('../scRNA_DATA/HUGO_with_ENSG_v32(missing_ENSG).csv') #9個
HUGO_gene_missed_list$HugoSymbol<- sub("\\..*", "", HUGO_gene_missed_list$HugoSymbol) #去小數點

#----抓註解
annot <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters    = "external_gene_name",
  values     = HUGO_gene_missed_list$HugoSymbol,
  mart       = ensembl
)

#----篩 protein coding genes
protein_coding_genes <- annot[annot$gene_biotype == "protein_coding", ] #8個 ->"LINC01238"不是protein coding gene
protein_coding_genes

#補回ENSG
merged_df <- merge(
  HUGO_gene_missed_list, 
  protein_coding_genes, 
  by.x = "HugoSymbol", 
  by.y = "external_gene_name", 
  all.x = TRUE
)
merged_df <- subset(merged_df, select = -ENSG)
merged_df  <- subset(
  merged_df ,
  !(HugoSymbol %in% c("LINC01238"))
) #protein coding genes:16509個
merged_df$HugoSymbol <- paste0(merged_df$HugoSymbol, ".1") #加回原本.1
print(merged_df)
#write.csv(merged_df,'../scRNA_DATA/HUGO_with_ENSG_v32(missing_ENSG)(protein coding).csv',row.names =F)

########################【未消失ENSG的hugo symbol篩protein coding gene】########################
library(clusterProfiler)
HUGO_gene_list<-read.csv('../scRNA_DATA/HUGO_with_ENSG_v32(no_missing_ENSG).csv') #26179個
HUGO_gene_list$ENSG <- sub("\\..*", "", HUGO_gene_list$ENSG)

#----連到Ensembl
library(biomaRt)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

#----抓註解
annot <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters    = "ensembl_gene_id",
  values     = HUGO_gene_list$ENSG,
  mart       = ensembl
)

#----篩 protein coding genes
protein_coding_genes <- annot[annot$gene_biotype == "protein_coding", ] #16511個

#----重複的hugo symbol
dup_hugo <- protein_coding_genes[
  duplicated(protein_coding_genes$external_gene_name) | 
    duplicated(protein_coding_genes$external_gene_name, fromLast = TRUE),  #顯示所有
]
dup_hugo 
subset(dup_hugo, !external_gene_name=="") #2個hugo(4個ensg)

dedup_protein_coding_genes <- subset(
  protein_coding_genes,
  !(ensembl_gene_id %in% c("ENSG00000285437", "ENSG00000258724"))
) #protein coding genes:16509個


merged_df <- merge(
  dedup_protein_coding_genes, 
  HUGO_gene_list, 
  by.x = "ensembl_gene_id", 
  by.y = "ENSG", 
  all.x = TRUE
)
merged_df <- subset(merged_df, select = -external_gene_name) #16526(其中包含24個hugo重複)
dedup_merge_df <- merged_df[!duplicated(merged_df$HugoSymbol), ] #16503 去重複，若有重複保留第一個出現

#write.csv(dedup_merge_df,'../scRNA_DATA/HUGO_with_ENSG_v32(no_missing_ENSG)(protein coding).csv',row.names = F)

######################【2個protein coding df合併】################
protein_coding_df1<-read.csv("../scRNA_DATA/HUGO_with_ENSG_v32(missing_ENSG)(protein coding).csv") #8
protein_coding_df2<-read.csv("../scRNA_DATA/HUGO_with_ENSG_v32(no_missing_ENSG)(protein coding).csv") #16503

# 確保欄位順序相同
col_order <- c("HugoSymbol", "ensembl_gene_id", "gene_biotype")
protein_coding_df1 <- protein_coding_df1[, col_order]
protein_coding_df2 <- protein_coding_df2[, col_order]

# 合併
merged_df <- rbind(protein_coding_df1, protein_coding_df2)

# 移除 gene_biotype 欄位
merged_df <- merged_df[, c("HugoSymbol", "ensembl_gene_id")]

head(merged_df) 
#write.csv(merged_df,'../scRNA_DATA/HUGO_with_ENSG_v32(protein coding all).csv',row.names = F) #16511


#####################【protein coding matrix】#################
filtered_data=readRDS("../scRNA_DATA/NMOSD_count_metadata.rds")

# 取出 counts
counts <- filtered_data@assays$RNA@layers$counts

#加回gene & cell名
gene_ids <- rownames(filtered_data@assays$RNA)
rownames(counts) <- gene_ids   
colnames(counts) <- colnames(filtered_data)

#取protein coding gene
protein_coding_gene_list<-read.csv("../scRNA_DATA/HUGO_with_ENSG_v32(protein coding all).csv")$HugoSymbol #16511
keep_genes <- rownames(counts) %in% protein_coding_gene_list
counts_pc <- counts[keep_genes, ]
#dim(counts_pc) 16511x193159

# 取metadata
metadata_filtered<-filtered_data@meta.data
#dim(metadata_filtered)

data_pc <- CreateSeuratObject(
  counts       = counts_pc,
  meta.data    = metadata_filtered,
  assay        = "RNA",
  min.features = 0,  # 不因 features 太少丟 cell
  min.cells    = 0, # 不因出現細胞數太少丟基因
  project      = "NMOSD"
)
DefaultAssay(data_pc) <- "RNA"

#驗證是否含counts
Layers(data_pc[["RNA"]])   #counts
dim(GetAssayData(data_pc, layer = "counts"))   # genes x cells:16510 x 193384

# 讓 data layer = counts（重要!!!)
data_pc <-SetAssayData(
  object   = data_pc,
  assay    = "RNA",
  layer    = "data",
  new.data = GetAssayData(data_pc, assay = "RNA", layer = "counts")
)
Layers(data_pc[["RNA"]]) #"counts" "data" 
dim(GetAssayData(data_pc, layer="data"))  # genes x cells:16511 x 193384

# ----{存成RDS 格式}
saveRDS(data_pc , "../scRNA_DATA/NMOSD_count_metadata(protein coding).rds")


#----------------{存成.mtx}----------------
library(Matrix)

# 取 raw counts
m <- GetAssayData(data_pc, assay = "RNA", layer = "counts")

# output folder
outdir <- "../scRNA_DATA/mtx_nmosd"
dir.create(outdir, showWarnings = FALSE)

#matrix.mtx
Matrix::writeMM(m, file.path(outdir, "NMOSD_matrix(protein coding).mtx")) #回傳NULL代表有成功

#features.tsv (基因名)
write.table(
  data.frame(rownames(m)),
  file = file.path(outdir, "NMOSD_features(protein coding).tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#barcodes.tsv (細胞ID)
write.table(
  data.frame(colnames(m)),
  file = file.path(outdir, "NMOSD_barcodes(protein coding).tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#metadata
write.csv(
  data_pc@meta.data,
  file = file.path(outdir, "NMOSD_obs(protein coding).csv")
)
#-------------------------------------------



