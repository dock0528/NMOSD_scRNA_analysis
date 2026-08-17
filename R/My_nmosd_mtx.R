#####################{存成.mtx}###################
library(Matrix)
library(SeuratObject)
library(Seurat)
merged_protein_data<-readRDS("/staging/biology/jane0528/NMOSD/scRNA/Dataset/My_merged_protein_coding_genes/My_merged_count_metadata(protein_coding).rds")

# 取 raw counts
m <- GetAssayData(merged_protein_data, assay = "RNA", layer = "counts")

# output folder
outdir <- "/staging/biology/jane0528/NMOSD/scRNA/Dataset/My_merged_protein_coding_genes/My_mtx_merge_v2"
dir.create(outdir, showWarnings = FALSE)

#matrix.mtx
Matrix::writeMM(m, file.path(outdir, "My_merge_matrix(protein coding).mtx")) #回傳NULL代表有成功