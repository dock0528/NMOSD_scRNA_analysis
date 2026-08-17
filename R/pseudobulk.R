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


#===============================【NF-kB2 -> ICOSL、Mcl1、OTUD7B 】===========================
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

#---{加入sub celltype}
sub_celltype<-read.csv("../scRNA_DATA/My_merged_protein_coding_genes/My_merged_Azimuth_sub_celltype_CellChat.csv")
rownames(sub_celltype) <- sub_celltype$cell_id

common <- intersect(colnames(merge_data), rownames(sub_celltype))
length(common) # 218212(確認細胞皆有對到)

merge_data$sub_celltype <- NA
merge_data$sub_celltype[common] <-sub_celltype[common, "sub_celltype"]

#---{取B naive cells merge_data}
bnaive_merge_data<-subset(merge_data, subset = sub_celltype %in% c("B naive"))
#bnaive_merge_data:gene x cell:16489 x 14862

#---{pseudobulk RNA}
pseudobulk_rna <- AggregateExpression(bnaive_merge_data, 
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

res[rownames(res) %in% c('ICOSLG','OTUD7B'),]
res[rownames(res) %in% c('MCL1'),]

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
  main  = "NMOSD vs Healthy DEGs (pseudobulk-B naive)",
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
length(Up_DEGs);length(Down_DEGs);length(DEGS) #977;359;1336
ji5
c('ICOSLG','MCL1','OTUD7B')%in% rownames(bnaive_merge_data) #ICOSLG(即ICOSL) #TRUE TRUE TRUE
c('ICOSLG','MCL1','OTUD7B')%in% Up_DEGs #FALSE TRUE FALSE
c('ICOSLG','MCL1','OTUD7B')%in% Down_DEGs #FALSE FALSE FALSE

#====================【trascription factor stimulate genes】================
devtools::install_github("saezlab/dorothea")
library(dorothea)
library(dplyr)

data("dorothea_hs", package = "dorothea")  # human regulons

targets_NFKB2 <- dorothea_hs %>%
  filter(tf %in% c("NFKB2","RELB"),
         confidence %in% c("A","B")) %>%   # 先用高信心
  distinct(tf, target, confidence, mor)

targets_NFKB2 %>% count(tf, confidence)
head(targets_NFKB2)
length(unique(targets_NFKB2$target))
#---{在NFKB2 signal 刺激下是 up regulation}
targets_NFKB2_up<-targets_NFKB2 %>%
  filter(mor==1)
targets_NFKB2_up_genes<-unique(targets_NFKB2_up$target)
length(targets_NFKB2_up_genes) #18
targets_NFKB2_up_genes %in% rownames(bnaive_merge_data)

targets_NFKB2_up_genes[targets_NFKB2_up_genes %in% rownames(bnaive_merge_data)=='FALSE'] #CCL21

targets_NFKB2_up_genes[targets_NFKB2_up_genes %in% Up_DEGs=='TRUE'] #"TNFAIP3" "IL6" "IRF1"   

#---{在NFKB2 signal 刺激下是 down regulation}
targets_NFKB2_down<-targets_NFKB2 %>%
  filter(mor==-1)
targets_NFKB2_down_genes<-unique(targets_NFKB2_down$target)
length(targets_NFKB2_down_genes) #12
targets_NFKB2_down_genes %in% rownames(bnaive_merge_data)

targets_NFKB2_down_genes[targets_NFKB2_down_genes %in% rownames(bnaive_merge_data)=='FALSE'] # NO

targets_NFKB2_down_genes[targets_NFKB2_down_genes %in% Down_DEGs=='TRUE'] # NO   


#----[GSEA]

# res = as.data.frame(DESeq2結果)
ranks <- res$stat
names(ranks) <- rownames(res)

# 清理：去掉 NA/Inf，排序
ranks <- ranks[is.finite(ranks)]
ranks <- sort(ranks, decreasing = TRUE)

pathways <- list(NFKB2_RELB_targets = unique(targets_NFKB2$target))
library(fgsea)
fg1 <- fgsea(pathways = pathways, stats = ranks, nperm = 20000)
fg1
plotEnrichment(pathways$NFKB2_RELB_targets, ranks) +
  ggplot2::labs(title="GSEA: NFKB2/RELB targets")




