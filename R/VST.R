######################################【"VST" for each sample】############################################
#BiocManager::install("glmGamPoi")
#載入套件
library(Matrix)
library(matrixStats)
library(sctransform)
library(tools)
library(glmGamPoi)

#----{計算 row-wise mean 和 sd，支援稀疏矩陣 (CsparseMatrix)}
sparse_row_mean_sd <- function(mtx) {
  n_cells <- ncol(mtx)
  
  # row mean
  row_mean <- Matrix::rowMeans(mtx)
  
  # row variance = (E[X^2]) - (E[X])^2
  row_sq   <- Matrix::rowSums(mtx^2)
  row_var  <- (row_sq / n_cells) - (row_mean^2)
  row_var[row_var < 0] <- 0  # 避免浮點誤差
  
  row_sd   <- sqrt(row_var)
  
  return(data.frame(mean = row_mean, sd = row_sd))
}

#----{畫圖}
plot_vst_comparison <- function(mtx_file,barcode_file,gene_file, output_base_dir) {
  
  #create output folder
  if (!dir.exists(output_base_dir)) {
    dir.create(output_base_dir, recursive = TRUE)
  }
  
  # Read matrix.mtx
  data <- as(readMM(mtx_file), "CsparseMatrix")
  
  # Read barcodes (cells)
  barcodes <- read.delim(barcode_file, header = FALSE, stringsAsFactors = FALSE)[,1]
  colnames(data) <- barcodes
  
  # Read genes
  genes <- read.delim(gene_file, header = FALSE, stringsAsFactors = FALSE)[,1]
  rownames(data) <- genes
  
  #取sample名 
  sample_ids <- sub("(_[^_]+)$", "", barcodes)
  
  #每個sample分組計算VST
  for (sample in unique(sample_ids)) {
    cat("Processing sample:", sample, "\n")
    
    #取sample的data
    idx <- which(sample_ids == sample)
    sub_data <- as(data[, idx], "CsparseMatrix")
  
    #VST要確保gene在cell的raw count>0
    data_filter <- sub_data[rowSums(sub_data) > 0, ] 
  
    
    # raw counts
    df_raw <- sparse_row_mean_sd(data_filter)
    
    # log2(count+1)
    log_mat <- log1p(data_filter)   # log1p(x) = log(x+1)
    df_log  <- sparse_row_mean_sd(log_mat)
    
    # 輸出圖檔
    output_file <- file.path(output_base_dir, paste0(sample, "_VST_plot.png"))
    png(output_file, width = 1200, height = 600, res = 150)
    par(mfrow = c(1, 3), oma = c(2, 2, 5, 2))
    
    plot(df_raw$mean, df_raw$sd,
         pch = 20, col = "#00000033",
         main = "Raw counts",
         xlab = "Mean", ylab = "SD",
         xlim = c(0, 20), ylim = c(0, 100))
    
    plot(df_log$mean, df_log$sd,
         pch = 20, col = "#1f78b433",
         main = "log2(count + 1)",
         xlab = "Mean", ylab = "SD",
         xlim = c(0, 20), ylim = c(0, 100))
    
    vst_out <- sctransform::vst(data_filter, 
                                method = "glmGamPoi", #適合scrna
                                latent_var = c("log_umi"), 
                                return_gene_attr = TRUE, 
                                return_cell_attr = TRUE, 
                                verbosity = 1)
    
    #latent_var = c("log_umi"):每個 cell 的 UMI 總數取log ->模型考慮細胞的log UMI 數來校正定序深度差異
    #return_gene_attr = TRUE:輸出基因層級的統計資訊
    #verbosity = 1:控制訊息輸出，1 表示會顯示基本的運行過程
    
    gene_info <- vst_out$gene_attr #取vst後的資料
    
    gene_info$mean<-gene_info$amean #取每個gene的平均
    gene_info$sd<-sqrt(gene_info$residual_variance)
    
    
    #聚焦在0-20的Mean
    plot(gene_info$mean, gene_info$sd,
         pch = 20, col = "#e31a1c33",
         main = "VST",
         xlab = "Mean", ylab = "SD",
         xlim = c(0, 20),
         ylim=c(0,100))
    
    mtext(sample,
          outer = TRUE, line = 1.5, cex = 1.2, font = 2)
    
    dev.off()
  }
}

#----{執行輸出VST圖}
plot_vst_comparison(
  mtx_file = "/staging/biology/jane0528/NMOSD/scRNA/Dataset/My_merged_protein_coding_genes/My_mtx_merge/My_merge_matrix(protein coding).mtx",
  barcode_file = "/staging/biology/jane0528/NMOSD/scRNA/Dataset/My_merged_protein_coding_genes/My_mtx_merge/My_merge_barcodes(protein coding).tsv",
  gene_file = "/staging/biology/jane0528/NMOSD/scRNA/Dataset/My_merged_protein_coding_genes/My_mtx_merge/My_merge_features(protein coding).tsv",
  output_base_dir = "/staging/biology/jane0528/NMOSD/scRNA/Dataset/My_merged_protein_coding_genes/plot_VST_Merge"
)
cat("All VST plots finished!","\n")