####################################【 "CPM" for each sample 】#######################################
library(Matrix)
library(tools)

#----{讀檔}
nmosd_matrix <- readMM("../scRNA_DATA/mtx_nmosd/NMOSD_matrix.mtx")
genes <- read.delim("../scRNA_DATA/mtx_nmosd/NMOSD_features.tsv", header = FALSE, stringsAsFactors = FALSE) #delim讀tab分隔
cells <- read.delim("../scRNA_DATA/mtx_nmosd/NMOSD_barcodes.tsv", header = FALSE, stringsAsFactors = FALSE) #stringsAsFactors文字轉向量

rownames(nmosd_matrix) <- genes$V1   
colnames(nmosd_matrix) <- cells$V1  

dim(nmosd_matrix)   # gene × cell : 26135 x 193384

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
  mtx_file = "../scRNA_DATA/mtx_nmosd/NMOSD_matrix.mtx",
  barcode_file = "../scRNA_DATA/mtx_nmosd/NMOSD_barcodes.tsv",
  output_base_dir = "../scRNA_DATA/plot_CPM_NMOSD"
)
