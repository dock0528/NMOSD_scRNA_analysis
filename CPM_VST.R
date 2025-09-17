####################################【 "CPM" for each sample 】#######################################
library(Matrix)
library(tools)

#----{讀檔}
nmosd_matrix <- readMM("../scRNA_DATA/mtx_nmosd/NMOSD_matrix(gene unique).mtx")
genes <- read.delim("../scRNA_DATA/mtx_nmosd/NMOSD_features(gene unique).tsv", header = FALSE, stringsAsFactors = FALSE) #delim讀tab分隔
cells <- read.delim("../scRNA_DATA/mtx_nmosd/NMOSD_barcodes(gene unique).tsv", header = FALSE, stringsAsFactors = FALSE) #stringsAsFactors文字轉向量

rownames(nmosd_matrix) <- genes$V1   
colnames(nmosd_matrix) <- cells$V1  

dim(nmosd_matrix)   # gene × cell : 24092 x 193384

#----{畫圖}
plot_cpm_by_sample <- function(mtx_file, barcode_file, output_base_dir = "../scRNA_DATA/plot_CPM_NMOSD") {
  
  #create outout folder
  if (!dir.exists(output_base_dir)) {
    dir.create(output_base_dir, recursive = TRUE)
  }
  
  # Read matrix.mtx
  data <- readMM(mtx_file)
  
  # Read barcodes (cells)
  barcodes <- read.delim(barcode_file, header = FALSE, stringsAsFactors = FALSE)[,1]
  colnames(data) <- barcodes
  
  #取sample名 
  sample_ids <- sub("_.*", "", barcodes)
  
  #每個sample分組計算 CPM
  for (sample in unique(sample_ids)) {
    cat("Processing sample:", sample, "\n")
    
    #取sample的data
    idx <- which(sample_ids == sample)
    sub_data <- data[, idx]
    
    # 計算 CPM 和 log2(CPM+1)
    total_UMI <- colSums(sub_data)
    total_UMI[total_UMI == 0] <- 1
    cpm <- t(t(sub_data) / total_UMI) * 1e6
    cpm_values <- cpm@x
    log_cpm_values <- log2(cpm_values + 1)
    
    # 輸出圖檔
    output_file <- file.path(output_base_dir, paste0(sample, "_CPM_plot.png"))
    png(output_file, width = 1000, height = 500)
    par(mfrow = c(1, 2), oma = c(2, 2, 5, 2))
    hist(cpm_values, breaks = 100, main = "CPM Distribution", 
         xlab = "CPM", col = "skyblue2", border = "white")
    hist(log_cpm_values, breaks = 100, main = "log2(CPM + 1) Distribution", 
         xlab = "log2(CPM + 1)", col = "salmon", border = "white",xlim=c(0,20))
    mtext(sample, outer = TRUE, line = 1.5, cex = 1.5, font = 2)
    dev.off()
  }
}

#----{執行輸出CPM圖}
plot_cpm_by_sample(
  mtx_file = "../scRNA_DATA/mtx_nmosd/NMOSD_matrix(gene unique).mtx",
  barcode_file = "../scRNA_DATA/mtx_nmosd/NMOSD_barcodes(gene unique).tsv",
  output_base_dir = "../scRNA_DATA/plot_CPM_NMOSD"
)

######################################【"VST" for each sample】############################################

#載入套件
library(Matrix)
library(matrixStats)
library(sctransform)
library(tools)

#----{讀檔}
nmosd_matrix <- readMM("../scRNA_DATA/mtx_nmosd/NMOSD_matrix(gene unique).mtx")

#----{畫圖}
plot_vst_comparison <- function(mtx_file,barcode_file,gene_file, output_base_dir) {
  
  #create output folder
  if (!dir.exists(output_base_dir)) {
    dir.create(output_base_dir, recursive = TRUE)
  }
  
  # Read matrix.mtx
  data <- readMM(mtx_file)
  
  # Read barcodes (cells)
  barcodes <- read.delim(barcode_file, header = FALSE, stringsAsFactors = FALSE)[,1]
  colnames(data) <- barcodes
  
  # Read genes
  genes <- read.delim(gene_file, header = FALSE, stringsAsFactors = FALSE)[,1]
  rownames(data) <- genes
  
  #取sample名 
  sample_ids <- sub("_.*", "", barcodes)
  
  #每個sample分組計算VST
  for (sample in unique(sample_ids)) {
    cat("Processing sample:", sample, "\n")
    
    #取sample的data
    idx <- which(sample_ids == sample)
    sub_data <- data[, idx]
  
    #VST要確保gene在cell的raw count>0
    sub_data <- as(sub_data, "CsparseMatrix") #確保稀疏格式
    data_filter <- sub_data[rowSums(sub_data) > 0, ] 
  
    #row count & log2計算
    raw_mat <- as.matrix(data_filter)
    log2_mat <- log2(raw_mat + 1)
    
    get_mean_sd <- function(mat) {
      data.frame(mean = rowMeans(mat), sd = rowSds(mat))
    }
    df_raw <- get_mean_sd(raw_mat)
    df_log <- get_mean_sd(log2_mat)
    
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
    
    vst_out <- sctransform::vst(data_filter, latent_var = c("log_umi"), return_gene_attr = TRUE, 
                                return_cell_attr = TRUE, verbosity = 1)
    
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
  mtx_file = "../scRNA_DATA/mtx_nmosd/NMOSD_matrix(gene unique).mtx",
  barcode_file = "../scRNA_DATA/mtx_nmosd/NMOSD_barcodes(gene unique).tsv",
  gene_file = "../scRNA_DATA/mtx_nmosd/NMOSD_features(gene unique).tsv",
  output_base_dir = "../scRNA_DATA/plot_VST_NMOSD"
)
