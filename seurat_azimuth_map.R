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
human_pbmc_ref <- LoadReference('../Azimuth_human_PBMC/')

#===【Harmony Data】===
my_harmony_merged_data<-readRDS("../scRNA_DATA/My_merged_protein_coding_genes/My_merged_Harmony(protein_coding).rds")

mean(!grepl("^ENSG", rownames(my_harmony_merged_data[["RNA"]])))

my_harmony_merged_data <- RunAzimuth(my_harmony_merged_data, reference = human_pbmc_ref)