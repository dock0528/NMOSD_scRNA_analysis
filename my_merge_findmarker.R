library(Seurat)
library(SeuratDisk)
#BiocManager::install("MAST")
library(MAST)

merge_data <- readRDS("C:/Users/Jane/Desktop/Wang實驗室/NMOSD研究計畫/scRNA/scRNA_DATA/My_merged_protein_coding_genes/My_merged_count_metadata(protein_coding).rds")

#----{Normalize:校正不同細胞測序深度的差異}
merge_data<-NormalizeData(merge_data, normalization.method = "LogNormalize", scale.factor = 10000)

#View(merge_data)

# 讀scANVI後的metadata celltype
celltype_metadata <- read.csv("../scRNA_DATA/My_merged_protein_coding_genes/My_merged_cell_metadata_HARMONY.csv")
head(celltype_metadata)

# 確保rownames為Cell ID
rownames(celltype_metadata)<- celltype_metadata$cell_id

# 加"分群" & "組別"加入merge_data
merge_data <- AddMetaData(
  object = merge_data,
  metadata = celltype_metadata$cluster_annotation,
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
write.csv(up_genes_MAST, "../scRNA_DATA/My_upregulated_df_MAST(padj_default).csv", row.names = FALSE)



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

