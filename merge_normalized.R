library(SeuratObject) # v5.2.0
library(Seurat) # v5.3.0
merge_data<-readRDS("../scRNA_DATA/Merged_protein_coding_genes/Merged_count_metadata(protein_coding).rds")

#----{Normalize:校正不同細胞測序深度的差異}
merge_data<-NormalizeData(merge_data, normalization.method = "LogNormalize", scale.factor = 10000)

#----{VST}
merge_data <- FindVariableFeatures(merge_data, selection.method = "vst", nfeatures = 2000)

# top2000 hvg
top2000_hvg<-VariableFeatures(merge_data)
#top10 hvg
top10_hvg<-head(top2000_hvg, 10)

#----{確認7個IFN-I都有在top2000 HVG}
IFN_I_DEGs<-c('ISG15','IFI6','CMPK2','LY6E','OASL','AKAP12','TNFRSF13C')
IFN_I_DEGs %in% top2000_hvg #皆為TRUE

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(merge_data)
plot2 <- LabelPoints(plot = plot1, points = top10_hvg, repel = TRUE)
x11()
plot1 + plot2
library(ggplot2)
ggsave("../scRNA_DATA/top10hvg_plot.png", plot = plot1 + plot2, width = 16, height = 8, dpi = 300)

#----{Scale:讓每個基因在細胞的平均為0，標準差為1}
all.genes <- rownames(merge_data)
merge_data <- ScaleData(merge_data, features = all.genes)

#----{使用top2000hvg畫pca}
merge_data <- RunPCA(merge_data, features = top2000_hvg) 
ElbowPlot(merge_data)

#----{儲存rds}
saveRDS(merge_data , "../scRNA_DATA/Merged_protein_coding_genes/Merged_PCA(protein_coding).rds")

#----{rds轉h5ad}
#if (!requireNamespace("remotes", quietly = TRUE)) {
#  install.packages("remotes")}
#remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)

merge_data<-readRDS("../scRNA_DATA/Merged_protein_coding_genes/Merged_PCA(protein_coding).rds")

# 保存Seurat物件為H5AD格式
SaveH5Seurat(merge_data, filename = "../scRNA_DATA/Merged_protein_coding_genes/Merge_PCA(protein_coding).h5Seurat")

# 使用 anndata 包將 .h5Seurat 轉換為 .h5ad 格式
Convert("../scRNA_DATA/Merged_protein_coding_genes/Merge_PCA(protein_coding).h5Seurat", dest = "h5ad")

