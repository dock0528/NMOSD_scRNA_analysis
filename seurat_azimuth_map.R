#===【Install Packages】===
#devtools::install_github("satijalab/seurat-data", "seurat5")

#---{安裝相依檔(把原始gz檔手動載到電腦 -> 再install)}
#install.packages("../BSgenome.Hsapiens.UCSC.hg38_1.4.5.tar.gz",
                 repos = NULL, type = "source")

#remotes::install_github('satijalab/azimuth', ref = 'master')

#===【Load Packages】===
library(Seurat)
library(SeuratData)
library(Azimuth)
library(patchwork)

#===【Human PBMC reference】===
# 先手動載資料
human_pbmc_ref<- "../Azimuth_human_PBMC"   
dir.exists(human_pbmc_ref)
list.files(human_pbmc_ref)

#===【Harmony Data】===
#my_harmony_merged_data<-readRDS("../scRNA_DATA/My_merged_protein_coding_genes/My_merged_Harmony(protein_coding).rds")
my_PCA_merged_data<-readRDS("../scRNA_DATA/My_merged_protein_coding_genes/My_merged_PCA(protein_coding).rds")
#===【Azimuth annotation】===
my_harmony_merged_data <- RunAzimuth(my_harmony_merged_data, reference = human_pbmc_ref )
my_PCA_merged_data <- RunAzimuth(my_PCA_merged_data, reference = human_pbmc_ref )

#my_harmony_merged_data<-readRDS("../scRNA_DATA/My_merged_protein_coding_genes/My_merged_Azimuth(protein_coding).rds")

sub <- my_harmony_merged_data@meta.data[
  my_harmony_merged_data@meta.data$predicted.celltype.l1 == "B",
  c("predicted.celltype.l1", "predicted.celltype.l1.score",
    "predicted.celltype.l2", "predicted.celltype.l2.score")
]

sub[order(sub$predicted.celltype.l2.score), ][1:30, ]

#===【存成.h5ad】===
#----{存Harmony後的rds}
library(SeuratDisk)
saveRDS(
  my_harmony_merged_data,
  "../scRNA_DATA/My_merged_protein_coding_genes/My_merged_Azimuth(protein_coding).rds"
)

#----{創建新Assay，改Assay5為Assay}
my_harmony_merged_data_copy <- my_harmony_merged_data
my_harmony_merged_data_copy[["RNA"]] <- as(object = my_harmony_merged_data_copy[["RNA"]], Class = "Assay")

#----{保存Seurat物件為H5AD格式}
SaveH5Seurat(
  my_harmony_merged_data_copy,
  filename = "../scRNA_DATA/My_merged_protein_coding_genes/My_merged_Azimuth(protein_coding).h5Seurat"
)

#----{將.h5Seurat 轉換為.h5ad 格式}
SeuratDisk::Convert(
  "../scRNA_DATA/My_merged_protein_coding_genes/My_merged_Azimuth(protein_coding).h5Seurat",
  dest = "h5ad"
)
#-----------------------------------{未做Harmony校正，直Azimuth map PBMC}-----------------------------
#===【PCA Data】===
my_PCA_merged_data<-readRDS("../scRNA_DATA/My_merged_protein_coding_genes/My_merged_PCA(protein_coding).rds")

#===【Azimuth annotation】===
my_PCA_merged_data <- RunAzimuth(my_PCA_merged_data, reference = human_pbmc_ref )

sub <- my_PCA_merged_data@meta.data[
  my_PCA_merged_data@meta.data$predicted.celltype.l1 == "B",
  c("predicted.celltype.l1", "predicted.celltype.l1.score",
    "predicted.celltype.l2", "predicted.celltype.l2.score")
]

sub[order(sub$predicted.celltype.l2.score), ][1:30, ]

#===【存成.h5ad】===
#----{存Harmony後的rds}
library(SeuratDisk)
saveRDS(
  my_PCA_merged_data,
  "../scRNA_DATA/My_merged_protein_coding_genes/My_merged_Azimuth_no_correction(protein_coding).rds"
)

#----{創建新Assay，改Assay5為Assay}
my_PCA_merged_data_copy <- my_PCA_merged_data
my_PCA_merged_data_copy[["RNA"]] <- as(object = my_PCA_merged_data_copy[["RNA"]], Class = "Assay")

#----{保存Seurat物件為H5AD格式}
SaveH5Seurat(
  my_PCA_merged_data_copy,
  filename = "../scRNA_DATA/My_merged_protein_coding_genes/My_merged_Azimuth_no_correction(protein_coding).h5Seurat"
)

#----{將.h5Seurat 轉換為.h5ad 格式}
SeuratDisk::Convert(
  "../scRNA_DATA/My_merged_protein_coding_genes/My_merged_Azimuth_no_correction(protein_coding).h5Seurat",
  dest = "h5ad"
)
#-----------------------------------{不一定要跑，可以轉去python}-----------------------------------
#===【KNN graph】===
my_harmony_merged_data <- FindNeighbors(my_harmony_merged_data,k.param=20, reduction = "harmony", dims = 1:50)
# k.param=20(default)->n_neighbor

#===【Plot UMAP】===
my_harmony_merged_data <- RunUMAP(
  my_harmony_merged_data,
  reduction = "harmony",
  dims = 1:50,
  reduction.name = "umap_harmony",
  reduction.key = "Harmony_"
)
my_harmony_merged_data@meta.data$predicted.celltype.l2

library(ggplot2)
library(scCustomize)
colormap <- DiscretePalette_scCustomize(num_colors = 30, palette = "varibow")
DimPlot(
  my_harmony_merged_data,
  reduction = "umap_harmony",
  group.by = "predicted.celltype.l1",
  pt.size = 0.5,
  raster = FALSE,
  alpha=0.3,
  cols=colormap
) + 
  ggtitle("Harmony-corrected UMAP (group by celltype)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    # 去除背景與格線
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # 加上外框
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

