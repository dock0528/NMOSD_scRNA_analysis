library(Seurat)
library(SeuratDisk)
#BiocManager::install("MAST")
library(MAST)

merge_data <- readRDS("C:/Users/Jane/Desktop/Wang實驗室/NMOSD研究計畫/scRNA/scRNA_DATA/My_merged_protein_coding_genes/My_merged_count_metadata(protein_coding).rds")

#----{Normalize:校正不同細胞測序深度的差異}
merge_data<-NormalizeData(merge_data, normalization.method = "LogNormalize", scale.factor = 10000)

#View(merge_data)

# 讀Azimuth後的metadata celltype
celltype_metadata <- read.csv("../scRNA_DATA/My_merged_protein_coding_genes/My_merged_Azimuth_sub_celltype_CellChat.csv")
head(celltype_metadata)

# 確保rownames為Cell ID
rownames(celltype_metadata)<- celltype_metadata$cell_id

# 加"分群" & "組別"加入merge_data
merge_data <- AddMetaData(
  object = merge_data,
  metadata = celltype_metadata$sub_celltype,
  col.name = "cluster_annotation"
)
merge_data <- AddMetaData(
  object = merge_data,
  metadata = celltype_metadata$Condition,
  col.name = "Condition"
)
#View(merge_data)

#Idents(merge_data) <- "Condition"

# 7個 IFN-I 基因
IFN_genes <- c("ISG15", "IFI6", "CMPK2", "LY6E", "OASL", "AKAP12", "TNFRSF13C")

# 建空dataframe
all_results <- data.frame()

library(dplyr)
#----{Differential expression analysis}
for (cell in unique(merge_data$cluster_annotation)) {
  
  cat("Processing:", cell, "\n")
  
  # celltype merge_data
  sub_obj <- subset(merge_data, subset = cluster_annotation == cell)
  
  #確認足夠的細胞數量
  cell_counts <- table(factor(sub_obj$Condition,
                              levels = c("NMOSD", "Control")))
  if (any(cell_counts < 3)) {
    cat(" 細胞數太少 ( < 3 ) >>>",
        "NMOSD :", cell_counts["NMOSD"],
        ", Control :", cell_counts["Control"], "\n")
    next
  }
  
  # FindMarkers()
  markers <- FindMarkers(
    sub_obj,
    slot= "data", 
    group.by = 'Condition',
    ident.1 = "NMOSD",
    ident.2 = "Control",
    features = IFN_genes,
    test.use = "MAST",
    logfc.threshold = 0,
    min.pct = 0.01,
    min.cells.group = 3
  )

  # 結果加上cell type
  markers <- markers %>%
    mutate(Celltype = cell, Gene = rownames(markers)) %>%  
    select(Celltype,Gene, everything())  
  # mutate:新增或修改欄位 #新增celltype欄位
  # Gene,Celltype放最前面，並列出所有統計結果值
  
  all_results <- bind_rows(all_results, markers)
}

#------ {檢查findmarkers p.adjust校正的方式}------
# 所有data上的基因一起做校正
n_genes_total <- nrow(merge_data) #16489 genes
all_results_test <- all_results %>%
  mutate(
    test_padj = p.adjust(p_val, method = "bonferroni", n = n_genes_total)
  )

# Findmarker results
head(all_results_test[order(all_results_test$test_padj), ])

# Upregulated genes
up_genes <- all_results_test[all_results_test$avg_log2FC > 1 & all_results_test$test_padj < 0.05, ]
up_genes <- up_genes %>%
  arrange(Gene, test_padj) 
up_genes

# Downregulated genes
down_genes <- all_results_test[all_results_test$avg_log2FC < -1 & all_results_test$test_padj < 0.05, ] #No
down_genes

### 輸出成csv給python
library(dplyr)
up_genes_MAST <- up_genes %>%
  dplyr::rename(
    `IFN-I gene`= Gene,
    log2FC = avg_log2FC,
    pval   = p_val,
    padj = test_padj
  ) %>%
  dplyr::select(Celltype, `IFN-I gene`, log2FC, pval, padj)

# 輸出 CSV
write.csv(up_genes_MAST, "../scRNA_DATA/My_upregulated_df_MAST(Azimuth).csv", row.names = FALSE)



# ------{各自Celltype做校正}------
all_results_new <- all_results %>%
  group_by(Celltype) %>%
  mutate(                                 
    my_padj = p.adjust(p_val, method = "bonferroni")
  ) 
# Findmarker results
head(all_results_new[order(all_results_new$my_padj), ])


# Upregulated genes
up_genes <- all_results_new[all_results_new$avg_log2FC > 1 & all_results_new$my_padj < 0.05, ]
up_genes <- up_genes %>%
  arrange(Gene, my_padj) 
up_genes

# Downregulated genes
down_genes <- all_results_new[all_results_new$avg_log2FC < -1 & all_results_new$my_padj < 0.05, ] #No
down_genes

### 輸出成csv給python
library(dplyr)
up_genes_renamed <- up_genes %>%
  dplyr::rename(
    `IFN-I gene`= Gene,
    log2FC = avg_log2FC,
    pval   = p_val,
    padj = my_padj
  ) %>%
  dplyr::select(Celltype, `IFN-I gene`, log2FC, pval, padj)

# 輸出 CSV
write.csv(up_genes_renamed, "../scRNA_DATA/My_upregulated_df_MAST_v2.csv", row.names = FALSE)

#------------------------------------
#檢查pytho有顯著，但R無顯著的值
# {CMPK2}
subset(all_results,
       Celltype %in% c("B memory", "Treg") & Gene == "CMPK2")

# {OASL}
subset(all_results,
       Celltype %in% c("CD8 Naive", "dnT") & Gene == "OASL")
# {TNFRSF13C}
subset(all_results,
       Celltype %in% c("CD4 Proliferating", "MAIT","CD8 Naive",'Treg') & Gene == "TNFRSF13C")

#############################【Mono/DC subcelltype marker genes】###################
library(Seurat)
library(SeuratDisk)
library(MAST)

merge_data <- readRDS("C:/Users/Jane/Desktop/Wang實驗室/NMOSD研究計畫/scRNA/scRNA_DATA/My_merged_protein_coding_genes/My_merged_count_metadata(protein_coding).rds")

#----{Normalize:校正不同細胞測序深度的差異}
merge_data<-NormalizeData(merge_data, normalization.method = "LogNormalize", scale.factor = 10000)

#View(merge_data)

# 讀scANVI後的metadata celltype
mono_dc_sub_df<- read.csv("../scRNA_DATA/Mono_DC_cluster_subtype_nn30.csv")
head(mono_dc_sub_df)

# Cell counts
sum(mono_dc_sub_df$cell_id %in% colnames(merge_data))  # 19463 cells

# 用 cell_id 當 rownames
rownames(mono_dc_sub_df) <- mono_dc_sub_df$cell_id

# 取 Mono_DC
common <- intersect(colnames(merge_data), rownames(mono_dc_sub_df))
merge_data <- AddMetaData(merge_data, mono_dc_sub_df[common, c("Mono_DC_cluster","Mono_DC_subtype")])
#View(merge_data)

mono_dc <- subset(merge_data, subset = !is.na(Mono_DC_cluster))
Idents(mono_dc) <- "Mono_DC_cluster" # 指定比較的欄位

Mono_DC_DEA<- FindMarkers(
  mono_dc,
  ident.1 = "8", # Cluster8
  ident.2 = NULL,        # vs other clusters
  test.use = "MAST",
  logfc.threshold = 0,
  min.pct = 0.01,
  min.cells.group = 3
)

# {Significant top30 genes}
head(Mono_DC_DEA[order(mk8$p_val_adj), ], 30)

# {看 pDC markers 是否顯著上升}
pDC<-c('ITM2C', 'PLD4', 'SERPINF1', 'LILRA4', 'IL3RA', 'TPM2', 'MZB1', 'SPIB', 'IRF4', 'SMPD3')
Mono_DC_DEA[pDC, c("avg_log2FC","pct.1","pct.2","p_val_adj")]

# {Mono/DC other sub celltype markers}
CD14Mono<-c('S100A9', 'CTSS', 'S100A8', 'LYZ', 'VCAN', 'S100A12', 'IL1B', 'CD14', 'G0S2', 'FCN1')
Mono_DC_DEA[CD14Mono, c("avg_log2FC","pct.1","pct.2","p_val_adj")]

CD16Mono<-c('CDKN1C', 'FCGR3A', 'PTPRC', 'LST1', 'IER5', 'MS4A7', 'RHOC', 'IFITM3', 'AIF1', 'HES4')
Mono_DC_DEA[CD16Mono, c("avg_log2FC","pct.1","pct.2","p_val_adj")]

cDC1<-c('CLEC9A', 'DNASE1L3', 'C1orf54', 'IDO1', 'CLNK', 'CADM1', 'FLT3', 'ENPP1', 'XCR1', 'NDRG2')
Mono_DC_DEA[cDC1, c("avg_log2FC","pct.1","pct.2","p_val_adj")]

cDC2<-c('FCER1A', 'HLA-DQA1', 'CLEC10A', 'CD1C', 'ENHO', 'PLD4', 'GSN', 'SLC38A1', 'NDRG2', 'AFF3')
Mono_DC_DEA[cDC2, c("avg_log2FC","pct.1","pct.2","p_val_adj")]
