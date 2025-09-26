#install.packages("SeuratObject")
#install.packages("Seurat")
library(SeuratObject)
library(Seurat)
control_data=readRDS("../scRNA_DATA/E-GEAD-551_PBMC_scRNAseq.rds")
head(control_data@assays$RNA@counts)

colnames(control_data@meta.data)
head(control_data@meta.data)

#----{取 Healthy Control}
# 只保留 Condition == "Healthy" 的細胞 (Healthy Control)
healthy_data <- subset(control_data, subset = Condition == "Healthy")

# ----{存成RDS 格式}
#saveRDS(healthy_data  , "../scRNA_DATA/E-GEAD-551_PBMC_scRNAseq(Healthy Control).rds")
