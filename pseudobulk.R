# Packages
#BiocManager::install(c("SingleCellExperiment", "scuttle"))
library(scuttle)
library(SingleCellExperiment)
library(Seurat)
library(SeuratDisk)
library(DESeq2)
library(edgeR)
library(ggplot2)

# My merge data (raw count)
merge_data <- readRDS("C:/Users/Jane/Desktop/Wang實驗室/NMOSD研究計畫/scRNA/scRNA_DATA/My_merged_protein_coding_genes/My_merged_count_metadata(protein_coding).rds")

# ---{確認sample名}
unique(merge_data@meta.data$orig.ident)

# ---{加入分群}
merge_data@meta.data$Condition <- ifelse(grepl("^nmo", merge_data@meta.data$orig.ident, ignore.case = TRUE),
                       "NMOSD",
                       "Control")

# 存回 Seurat
merge_data@meta.data <- merge_data@meta.data

pseudobulk_rna <- AggregateExpression(merge_data, 
                          group.by="orig.ident",
                          assays="RNA", 
                          slot="counts", 
                          fun="sum")
pseudobulk_matrix<- pseudobulk_rna$RNA
dim(pseudobulk_matrix) #gene x sample:16489 x 28

#--------------------------------🔻{DESeq2}🔻---------------------------------------------
########################{Sample metadata}########################
#Sample & Group
sample_df=data.frame(Sample=colnames(pseudobulk_matrix))
sample_df$Group <- ifelse(
  grepl("^nmo", sample_df$Sample),  # 如果 Sample 以 "SRR" 開頭
  "NMOSD",
  "Control"                   
)
# 把 sample_df 裡的 Group 欄轉成 factor，Control 當作 baseline
sample_df$Group <- factor(sample_df$Group,
                          levels = c("Control","NMOSD"))
rownames(sample_df) <- sample_df$Sample

########################{edgeR filterByExpr}########################
dge0 <- DGEList(counts = pseudobulk_matrix, group = sample_df$Group)

keep <- filterByExpr(dge0, group = sample_df$Group)   # 過濾低表達genes
dge  <- dge0[keep, , keep.lib.sizes = FALSE] #keep.lib.sizes = FALSE:重新計算每個樣本的 library size

counts_filter <- dge$counts    # 過濾後matrix -> gene x cell:13669 x 28

########################{DESeq2 分析}########################
# 確保樣本順序一致
counts_filter <- counts_filter[, rownames(sample_df)]

# ---- {建 DESeqDataSet} 
dds <- DESeqDataSetFromMatrix(countData = counts_filter,
                              colData   = sample_df,
                              design    = ~ Group)
# ----{執行DESeq }
dds <- DESeq(dds)


res<-results(dds)
head(results(dds,tidy=T)) #tidy=T整齊回傳
#表達量太低>>>padj=NA


#summary of differential gene expression
summary(res)

#sort summary list by padjust
res<-res[order(res$padj),]
head(res)



###########################{Volcano Plot}###########################

# res 欄位含 log2FoldChange, pvalue, padj
# 先把原本的 par 設定存起來，最後再還原
op <- par(no.readonly = TRUE) # no.readonly = TRUE:full list of parameters 

# 調整邊界，並允許在外部繪圖
par(mar = c(5, 4, 4, 6) , xpd = TRUE) #c(bottom, left, top, right)

# 基本火山圖：logFC vs. -log10(P.Value)
with(res, plot(
  log2FoldChange , -log10(pvalue),
  pch   = 20,
  main  = "NMOSD vs Healthy DEGs (pseudobulk)",
  xlim=c(-15,15),
  col   = "grey",
  xlab  = "log2 Fold Change",
  ylab  = "-log10(P Value)"
))


#----【Significants genes criteria: padj < 0.05 & |log2FoldChange|>1】-----
# 下調基因 (藍色)
with(subset(res,  padj < 0.05 & log2FoldChange   < -1),
     points(log2FoldChange, -log10(pvalue), pch=20, col='#8FA4FF')
)

# 上調基因 (紅色)
with(subset(res, padj < 0.05 & log2FoldChange   >1),
     points(log2FoldChange, -log10(pvalue), pch=20, col='#FF8080')
)

# 在圖外加圖例，往右外推 20%
legend("topright",
       inset   = c(-0.2, 0),
       legend  = c("Up", "Down"),
       title   = "Change",
       pch     = 20,
       col     = c("#FF8080", "#8FA4FF"),
       pt.cex  = 1.4, #點符號放大倍數
       bty     = "n" #the type of box 
)

# 還原原本的 par 設定
par(op)


#顯著基因
Up_DEGs<-row.names(subset(res,padj<0.05 & log2FoldChange>1)) #length(Up_DEGs):1326
Down_DEGs<-row.names(subset(res,padj<0.05 & log2FoldChange < -1)) #length(Down_DEGs):1576
DEGS<-row.names(subset(res,padj<0.05 & abs(log2FoldChange)>1)) #length(DEGS):2904
length(Up_DEGs);length(Down_DEGs);length(DEGS)
