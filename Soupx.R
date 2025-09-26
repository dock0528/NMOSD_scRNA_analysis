#載入套件
install.packages("SoupX")
BiocManager::install("DropletUtils")
library(SoupX)
library(DropletUtils)

#讀入所有樣本
ids <- list.dirs("../CellRanger_results", recursive = FALSE, full.names = FALSE) 
#recursive = FALSE:僅列出第一層資料夾 #full.names = FALSE:僅列出資料夾名稱(不要完整路徑)
base_dir <- "../CellRanger_results"   # 你的母資料夾

for (id in ids) {
  message("Soupx正在處理:", id, " ...")
  
  #Read cellranger output
  sc <- load10X(file.path(base_dir, id, "outs"))
  
  #自動估算ambient RNA 
  soup <- autoEstCont(sc)
  
  #產出校正matrix
  adj.matrix <- adjustCounts(soup, roundToInt = TRUE)
  
  #Output
  out_dir <- file.path(base_dir, id, "SoupX_corrected_outs")
  write10xCounts(out_dir, adj.matrix) 
  message("完成:", id)
  
}
