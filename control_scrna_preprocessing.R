#install.packages("SeuratObject")
#install.packages("Seurat")
library(SeuratObject)
library(Seurat)
c=readRDS("../scRNA_DATA/E-GEAD-551_PBMC_scRNAseq.rds")
head(control_data@assays$RNA@counts)

colnames(control_data@meta.data)
head(control_data@meta.data)

################{原始data存成.mtx}##################
library(Matrix)

# 取 raw counts
m <- GetAssayData(control_data, assay = "RNA", layer = "counts")

# OUTPUT folder
outdir <- "../scRNA_DATA/mtx_control"
dir.create(outdir, showWarnings = FALSE)

#matrix.mtx
Matrix::writeMM(m, file.path(outdir, "control_matrix.mtx"))

#features.tsv (基因名)
write.table(
  data.frame(rownames(m)),
  file = file.path(outdir, "control_features.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#barcodes.tsv (細胞ID)
write.table(
  data.frame(colnames(m)),
  file = file.path(outdir, "control_barcodes.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#metadata
write.csv(
  control_data@meta.data,
  file = file.path(outdir, "control_obs.csv")
)

#----------------------------------------------------

#----{取 Healthy Control}
# 只保留 Condition == "Healthy" 的細胞 (Healthy Control)
healthy_data <- subset(control_data, subset = Condition == "Healthy")

# ----{存成RDS 格式}
#saveRDS(healthy_data  , "../scRNA_DATA/E-GEAD-551_PBMC_scRNAseq(Healthy Control).rds")

################{Healthy Control data存成.mtx}##################
library(Matrix)

# 取 raw counts
m_hc <- GetAssayData(healthy_data , assay = "RNA", layer = "counts")

# OUTPUT folder
outdir <- "../scRNA_DATA/mtx_control"
dir.create(outdir, showWarnings = FALSE)

#matrix.mtx
Matrix::writeMM(m_hc, file.path(outdir, "Healthy_Control_matrix.mtx"))

#features.tsv (基因名)
write.table(
  data.frame(rownames(m_hc)),
  file = file.path(outdir, "Healthy_Control_features.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#barcodes.tsv (細胞ID)
write.table(
  data.frame(colnames(m_hc)),
  file = file.path(outdir, "Healthy_Control_barcodes.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

#metadata
write.csv(
  healthy_data@meta.data,
  file = file.path(outdir, "Healthy_Control_obs.csv")
)
