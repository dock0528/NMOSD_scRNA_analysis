library(SeuratObject) # v5.2.0
library(Seurat)
library(harmony)
# merge_data<-readRDS("/staging/biology/jane0528/NMOSD/scRNA/Dataset/My_merged_protein_coding_genes/RDS/My_merged_count_metadata(protein_coding).rds")

# #----{Normalize:校正不同細胞測序深度的差異}
# merge_data<-NormalizeData(merge_data, normalization.method = "LogNormalize", scale.factor = 10000)

# #----{VST}
# merge_data <- FindVariableFeatures(merge_data, selection.method = "vst", nfeatures = 2000)

# # top2000 hvg
# top2000_hvg<-VariableFeatures(merge_data)
# #top10 hvg
# top10_hvg<-head(top2000_hvg, 10)

# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(merge_data)
# plot2 <- LabelPoints(plot = plot1, points = top10_hvg, repel = TRUE)
# plot1 + plot2

# #----{Scale:讓每個基因在細胞的平均為0，標準差為1}
# all.genes <- rownames(merge_data)
# merge_data <- ScaleData(merge_data, features = all.genes)

# #----{使用top2000hvg畫pca}
# merge_data <- RunPCA(merge_data, features = top2000_hvg) 
# pca_elbow<-ElbowPlot(merge_data,ndims = 50)

# # ---- {Harmony 批次校正} ----
# #它會自動讀metadata
# merge_data <- RunHarmony(
#   object = merge_data,
#   group.by.vars = "orig.ident",   #要校正的批次(sample)  
#   reduction.use = "pca",
#   dims.use=1:15,
#   plot_convergence=TRUE,
#   early_stop=TRUE #允許迭代在cost不再下降時提前終止
# )
# # 檢查是否成功Harmony
# Reductions(merge_data)

# # Harmony cell embeddings
# harmony_embeddings <- Embeddings(merge_data, 'harmony')[, 1:2]
# head(harmony_embeddings)

# #----{存Harmony後的rds}
# library(SeuratDisk)
# saveRDS(
#   merge_data,
#   "/staging/biology/jane0528/NMOSD/scRNA/Dataset/My_merged_protein_coding_genes/RDS/My_merged_Harmony(protein_coding)_v2.rds"
# )

# #----{創建新Assay，改Assay5為Assay}
# merge_data_copy <- merge_data
# merge_data_copy[["RNA"]] <- as(object = merge_data_copy[["RNA"]], Class = "Assay")

# #----{保存Seurat物件為H5AD格式}
# SaveH5Seurat(
#   merge_data_copy,
#   filename = "/staging/biology/jane0528/NMOSD/scRNA/Dataset/My_merged_protein_coding_genes/Adata/My_Merge_Harmony(protein_coding).h5Seurat"
# )

# #----{將.h5Seurat 轉換為.h5ad 格式}
# SeuratDisk::Convert(
#   "/staging/biology/jane0528/NMOSD/scRNA/Dataset/My_merged_protein_coding_genes/Adata/My_Merge_Harmony(protein_coding).h5Seurat",
#   dest = "h5ad"
# )

#----{Plot harmony1 & 2}
# library(ggplot2)
# # PC15
# harmony_data_pc15<-readRDS("/staging/biology/jane0528/NMOSD/scRNA/Dataset/y_merged_protein_coding_genes/RDS/y_merged_Harmony(protein_coding)_v2.rds")
# harmony_plot_pc15 <- DimPlot(object = harmony_data_pc15, reduction = "harmony", pt.size = 1.0,raster=FALSE, group.by = "orig.ident") # raster=FALSE不會模糊化
# ggsave("/staging/biology/jane0528/NMOSD/scRNA/Dataset/Merged_protein_coding_genes/Harmony_plot_pc15.png", plot = harmony_plot_pc15 , width = 14, height = 8, dpi = 300,bg = "white")

# # PC50
# harmony_data_pc50<-readRDS("/staging/biology/jane0528/NMOSD/scRNA/Dataset/Merged_protein_coding_genes/RDS/Merged_Harmony(protein_coding).rds")
# harmony_plot_pc50 <- DimPlot(object = harmony_data_pc50, reduction = "harmony", pt.size = 1.0,raster=FALSE, group.by = "orig.ident") # raster=FALSE不會模糊化
# library(ggplot2)
# ggsave("/staging/biology/jane0528/NMOSD/scRNA/Dataset/Merged_protein_coding_genes/Harmony_plot_pc50.png", plot = harmony_plot_pc50 , width = 14, height = 8, dpi = 300,bg = "white")

# # PC50
# My_harmony_data_pc50<-readRDS("/staging/biology/jane0528/NMOSD/scRNA/Dataset/My_merged_protein_coding_genes/RDS/My_merged_Harmony(protein_coding).rds")
# My_harmony_plot_pc50 <- DimPlot(object = My_harmony_data_pc50, reduction = "harmony", pt.size = 1.0,raster=FALSE, group.by = "orig.ident") # raster=FALSE不會模糊化
# library(ggplot2)
# ggsave("/staging/biology/jane0528/NMOSD/scRNA/Dataset/My_merged_protein_coding_genes/My_Harmony_plot_pc50.png", plot = My_harmony_plot_pc50 , width = 14, height = 8, dpi = 300,bg = "white")

#----{Plot PCA vs UMAP vs Harmony_umap}
# PCA
library(ggplot2)
library(scCustomize)
colormap <- DiscretePalette_scCustomize(num_colors = 30, palette = "varibow")
ori_merged_data<-readRDS("/staging/biology/jane0528/NMOSD/scRNA/Dataset/My_merged_protein_coding_genes/RDS/My_merged_PCA(protein_coding).rds")
pca_plot <- DimPlot(
  ori_merged_data,
  reduction = "pca",
  dims = c(1, 2),
  group.by = "orig.ident",
  pt.size = 0.5,
  raster = FALSE,
  alpha=0.3,
  cols=colormap
) + 
  ggtitle("PCA (NMOSD vs SRP349890)") +
  theme(
  plot.title = element_text(hjust = 0.5, face = "bold"),
  # 去除背景與格線
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),

  # 加上外框
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

ggsave(
  "/staging/biology/jane0528/NMOSD/scRNA/Dataset/My_merged_protein_coding_genes/My_PCA_plot(pc50).png",
  plot = pca_plot,
  width = 10, height = 7, dpi = 300, bg = "white"
)
# UMAP
ori_merged_data <- RunUMAP(
  ori_merged_data,
  reduction = "pca",
  dims = 1:50,
  reduction.name = "umap_pca", #存入seurat物件的名字
  reduction.key = "UMAP_"
)

umap_plot <- DimPlot(
  ori_merged_data,
  reduction = "umap_pca",
  group.by = "orig.ident",
  pt.size = 0.5,
  raster = FALSE,
  alpha=0.3,
  cols=colormap
) +
  ggtitle("UMAP (NMOSD vs SRP349890)") +
  theme(
  plot.title = element_text(hjust = 0.5, face = "bold"),
  # 去除背景與格線
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),

  # 加上外框
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

ggsave(
  "/staging/biology/jane0528/NMOSD/scRNA/Dataset/My_merged_protein_coding_genes/My_UMAP_plot(pc50).png",
  plot = umap_plot,
  width = 10, height = 7, dpi = 300, bg = "white"
)

# Harmony UMAP
harmony_data_pc50<-readRDS("/staging/biology/jane0528/NMOSD/scRNA/Dataset/My_merged_protein_coding_genes/RDS/My_merged_Harmony(protein_coding).rds")
harmony_data_pc50 <- RunUMAP(
  harmony_data_pc50,
  reduction = "harmony",
  dims = 1:50,
  reduction.name = "umap_harmony",
  reduction.key = "Harmony_"
)

umap_harmony_plot <- DimPlot(
  harmony_data_pc50,
  reduction = "umap_harmony",
  group.by = "orig.ident",
  pt.size = 0.5,
  raster = FALSE,
  alpha=0.3,
  cols=colormap
) + 
  ggtitle("Harmony-corrected UMAP (NMOSD vs SRP349890)") +
  theme(
  plot.title = element_text(hjust = 0.5, face = "bold"),
  # 去除背景與格線
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),

  # 加上外框
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

ggsave(
  "/staging/biology/jane0528/NMOSD/scRNA/Dataset/My_merged_protein_coding_genes/My_Harmony_plot(pc50).png",
  plot = umap_harmony_plot,
  width = 10, height = 7, dpi = 300, bg = "white"
)
