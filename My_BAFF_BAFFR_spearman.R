library(scuttle)
library(SingleCellExperiment)
library(Seurat)
library(SeuratDisk)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(dplyr)

# My merge data (raw count)
merge_data <- readRDS("C:/Users/Jane/Desktop/Wang實驗室/NMOSD研究計畫/scRNA/scRNA_DATA/My_merged_protein_coding_genes/My_merged_count_metadata(protein_coding).rds")

#----{Normalize:校正不同細胞測序深度的差異}
merge_data<-NormalizeData(merge_data, normalization.method = "LogNormalize", scale.factor = 10000)

# ---{確認sample名}
unique(merge_data@meta.data$orig.ident)

# ---{加入分群}
merge_data@meta.data$Condition <- ifelse(grepl("^nmo", merge_data@meta.data$orig.ident, ignore.case = TRUE),
                                         "NMOSD",
                                         "Control")

# 存回 Seurat
merge_data@meta.data <- merge_data@meta.data

#---{加入sub celltype}
sub_celltype<-read.csv("../scRNA_DATA/My_merged_protein_coding_genes/My_merged_Azimuth_sub_celltype_CellChat.csv")
rownames(sub_celltype) <- sub_celltype$cell_id

common <- intersect(colnames(merge_data), rownames(sub_celltype))
length(common) # 218212(確認細胞皆有對到)

merge_data$sub_celltype <- NA
merge_data$sub_celltype[common] <-sub_celltype[common, "sub_celltype"]

patient_col  <- "orig.ident"
group_col    <- "Condition"
celltype_col <- "sub_celltype"

x_gene <- "TNFSF13B"; x_celltype <- "cDC2"
y_gene <- "TNFRSF13C"; y_celltype <- "B naive"

#---{取cDC2 cells merge_data}
x_merge_data<-subset(merge_data, subset = sub_celltype %in% x_celltype)

x_n <- table(x_merge_data@meta.data$orig.ident)


#---{pseudobulk X RNA}
pseudobulk_rna_x <- AggregateExpression(x_merge_data, 
                                      group.by="orig.ident",
                                      assays="RNA", 
                                      slot="counts", 
                                      fun="sum",
                                      return.seurat=FALSE)

pseudobulk_matrix_x<- pseudobulk_rna_x$RNA
dim(pseudobulk_matrix_x)

pseudobulk_matrix_x <- as(pseudobulk_matrix_x, "dgCMatrix")
total_UMI_x <- Matrix::colSums(pseudobulk_matrix_x)
total_UMI_x[total_UMI_x == 0] <- 1

x_counts <- pseudobulk_matrix_x[x_gene, ]  # 1 x nPatients (仍是稀疏/向量形式)

# CPM per patient for this gene
x_vec <- as.numeric((x_counts / total_UMI_x) * 1e6)
x_vec <- log2(x_vec + 1)
names(x_vec) <- colnames(pseudobulk_matrix_x)

head(x_vec)


#---{取B naive cells}
y_merge_data <- subset(merge_data, subset = sub_celltype %in% y_celltype)

#---{pseudobulk Y}
pseudobulk_rna_y <- AggregateExpression(
  y_merge_data,
  group.by = "orig.ident",
  assays   = "RNA",
  fun      = "sum",
  return.seurat = FALSE
)
pseudobulk_matrix_y <- pseudobulk_rna_y$RNA




pseudobulk_matrix_y <- as(pseudobulk_matrix_y, "dgCMatrix")
total_UMI_y <- Matrix::colSums(pseudobulk_matrix_y)
total_UMI_y[total_UMI_y == 0] <- 1

y_counts <- pseudobulk_matrix_y[y_gene, ]  # 1 x nPatients (仍是稀疏/向量形式)

#---{log2CPM}
y_vec <- as.numeric((y_counts / total_UMI_y) * 1e6)
y_vec <- log2(y_vec + 1)
names(y_vec) <- colnames(pseudobulk_matrix_y)

head(y_vec)

#---{取Sample & group}
patient_group <- merge_data@meta.data %>%
  distinct(patient = orig.ident, group = Condition)

#---{交集皆有數值的samples}
common_patients <- intersect(names(x_vec), names(y_vec))

df_xy <- data.frame(
  patient = common_patients,
  X = x_vec[common_patients],
  Y = y_vec[common_patients],
  stringsAsFactors = FALSE
) %>%
  left_join(patient_group, by = "patient")

table(df_xy$group) #NMOSD:24 Control:3
head(df_xy)

#---{Spearman correlation}
result<- df_xy %>%
  group_by(group) %>%
  summarise(
    n = n(),
    r = cor(X, Y, method = "spearman", use = "complete.obs"),
    p = cor.test(X, Y, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("r = ", sprintf("%.3f", r),
                   "\np = ", format.pval(p, digits = 2, eps = 1e-3))
  )
result

#---{Plot:NMOSD & Control一起畫}

ggplot(df_xy, aes(X, Y, color = group)) +
  geom_point(size = 2.3, alpha = 0.85) +
  theme_classic(base_size = 14) +
  labs(
    title="Spearman's correlation between TNFSF13B and TNFRSF13C expression",
    x = paste0("TNFSF13B in ",x_celltype," (log2CPM)"),
    y = paste0("TNFRSF13C in ",y_celltype," (log2CPM)"),
    color = "Group"
  ) +
  theme(plot.title = element_text(size = 13, hjust = 0,face='bold'),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11)
        )+
  scale_x_continuous(limits = c(0, 9), breaks = seq(0, 9, by = 2)) +
  scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, by = 2))+
  scale_color_manual(values = c("NMOSD" = "#E3A19F", "Control" = "#67D6F0"))



#----------------------{Plot:NMOSD 和 Control 各自畫}------------------------
# 統一title & axis
base_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title  = element_text(size = 13, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  )

x_scale <- scale_x_continuous(limits = c(0, 9), breaks = seq(0, 9, by = 2))
y_scale <- scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, by = 2))

#---{Plot:NMOSD}
df_nmosd <- df_xy %>% filter(group == "NMOSD")
lab_nmosd <- result %>% filter(group == "NMOSD")

p_nmosd <- ggplot(df_nmosd, aes(X, Y)) +
  geom_point(size = 2.6, alpha = 0.8, color = "#E3A19F") +
  labs(
    title = "BAFF and BAFFR correlation - NMOSD",
    x = paste0("TNFSF13B in ",x_celltype," (log2CPM)"),
    y = paste0("TNFRSF13C in ",y_celltype," (log2CPM)")
  ) +
  x_scale + y_scale + base_theme +
  annotate("text",
           x = 7, y =2, hjust = 0, vjust = 1,
           label = lab_nmosd$label, size = 4)
p_nmosd

#---{Plot:Control}
df_ctrl <- df_xy %>% filter(group == "Control")
lab_ctrl <- result %>% filter(group == "Control")

p_ctrl <- ggplot(df_ctrl, aes(X, Y)) +
  geom_point(size = 2.8, alpha = 0.85, color = "#67D6F0") +
  labs(
    title = "BAFF and BAFFR correlation - Control",
    x = paste0("TNFSF13B in ",x_celltype," (log2CPM)"),
    y = paste0("TNFRSF13C in ",y_celltype," (log2CPM)")
  ) +
  x_scale + y_scale + base_theme +
  annotate("text",
           x = 7, y =2, hjust = 0, vjust = 1,
           label = lab_ctrl$label, size = 4)

p_ctrl


#----------------------{計算Sample細胞數}-----------------------

#---{patient;group;celltype;n_cells}
ct_counts <- merge_data@meta.data %>%
  count(
    patient  = .data[[patient_col]],
    group    = .data[[group_col]],
    celltype = .data[[celltype_col]],
    name = "n_cells"
  )

ct_counts %>% arrange(group, patient, desc(n_cells)) %>% head(20)

#---{patient;group;B naive cells;cDC2 cells}
ct_xy <- ct_counts %>%
  filter(celltype %in% c(x_celltype, y_celltype)) %>%
  pivot_wider(
    names_from  = celltype,
    values_from = n_cells,
    values_fill = 0
  )

ct_xy %>% arrange(group, patient) %>% head(20)

#---{Group cells distribution}
ct_xy %>%
  group_by(group) %>%
  summarise(
    n_patients = n(),
    cDC2_min = min(.data[[x_celltype]]),
    cDC2_q1  = quantile(.data[[x_celltype]], 0.25),
    cDC2_med = median(.data[[x_celltype]]),
    cDC2_q3  = quantile(.data[[x_celltype]], 0.75),
    cDC2_max = max(.data[[x_celltype]]),
    
    Bn_min = min(.data[[y_celltype]]),
    Bn_q1  = quantile(.data[[y_celltype]], 0.25),
    Bn_med = median(.data[[y_celltype]]),
    Bn_q3  = quantile(.data[[y_celltype]], 0.75),
    Bn_max = max(.data[[y_celltype]]),
    .groups="drop"
  )

#----------------------【Remove NMOSD outlier-nmo20】-----------------
# nmo020 的 cDC2 細胞數
x_n["nmo020"]

# nmo020 的 cDC2 pseudobulk 中，TNFSF13B 的 raw sum counts 到底是多少
pseudobulk_matrix_x["TNFSF13B", "nmo020"]
total_UMI_x["nmo020"]

df_xy_rm <- df_xy %>% 
  filter(patient != "nmo020")

# Spearman correlation
result_rm<- df_xy_rm %>%
  group_by(group) %>%
  summarise(
    n = n(),
    r = cor(X, Y, method = "spearman", use = "complete.obs"),
    p = cor.test(X, Y, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("r = ", sprintf("%.3f", r),
                   "\np = ", format.pval(p, digits = 2, eps = 1e-3))
  )

result_rm

# {Plot:NMOSD+Control}
ggplot(df_xy_rm, aes(X, Y, color = group)) +
  geom_point(size = 2.3, alpha = 0.85) +
  theme_classic(base_size = 14) +
  labs(
    title="Spearman's correlation between TNFSF13B and TNFRSF13C expression",
    x = "TNFSF13B in cDC2",
    y = "TNFRSF13C in B naive",
    color = "Group"
  ) +
  theme(plot.title = element_text(size = 13, hjust = 0,face='bold'),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11)
  )+
  scale_x_continuous(limits = c(0, 9), breaks = seq(0, 9, by = 2)) +
  scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, by = 2))+
  scale_color_manual(values = c("NMOSD" = "#E3A19F", "Control" = "#67D6F0"))

#---{Plot:NMOSD}
df_nmosd_rm <- df_xy_rm %>% filter(group == "NMOSD")
lab_nmosd_rm <- result_rm %>% filter(group == "NMOSD")

p_nmosd_rm <- ggplot(df_nmosd_rm, aes(X, Y)) +
  geom_point(size = 2.6, alpha = 0.8, color = "#E3A19F") +
  labs(
    title = "BAFF and BAFFR correlation - NMOSD",
    x = paste0("TNFSF13B in ",x_celltype," (log2CPM)"),
    y = paste0("TNFRSF13C in ",y_celltype," (log2CPM)")
  ) +
  scale_x_continuous(limits = c(6, 9), breaks = seq(0, 9, by = 1)) +
  scale_y_continuous(limits = c(6, 9), breaks = seq(0, 9, by = 1))+ 
  base_theme +
  annotate("text",
           x = 7, y =2, hjust = 0, vjust = 1,
           label = lab_nmosd_rm$label, size = 4)
p_nmosd_rm



#--------------------------------【建 function 一次跑所有celltypes】-----------------------------------
run_BAFF_BAFFR_one_Xcelltype <- function(
    merge_data,
    x_gene = "TNFSF13B",
    x_celltype = "cDC2",
    y_gene = "TNFRSF13C",
    y_celltype = "B naive",
    out_dir = "BAFF_BAFFR_corr_plots",
    x_lim = c(0, 9),
    y_lim = c(0, 9)
) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(dplyr)
    library(ggplot2)
    library(stringr)
  })
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------{你原本的 code：取 cDC2 cells}----------------------
  x_cells <- rownames(merge_data@meta.data)[merge_data@meta.data$sub_celltype %in% x_celltype]
  if (length(x_cells) == 0) {
    message("[SKIP] ", x_celltype, ": no cells")
    return(invisible(NULL))
  }
  x_merge_data <- subset(merge_data, cells = x_cells)
  
  x_n <- table(x_merge_data@meta.data$orig.ident)
  
  # ----------------------{pseudobulk X RNA}----------------------
  pseudobulk_rna_x <- AggregateExpression(
    x_merge_data,
    group.by="orig.ident",
    assays="RNA",
    slot="counts",
    fun="sum",
    return.seurat=FALSE
  )
  pseudobulk_matrix_x <- pseudobulk_rna_x$RNA
  pseudobulk_matrix_x <- as(pseudobulk_matrix_x, "dgCMatrix")
  
  total_UMI_x <- Matrix::colSums(pseudobulk_matrix_x)
  total_UMI_x[total_UMI_x == 0] <- 1
  
  if (!(x_gene %in% rownames(pseudobulk_matrix_x))) {
    message("[SKIP] ", x_celltype, ": x_gene not found in X pseudobulk: ", x_gene)
    return(invisible(NULL))
  }
  
  x_counts <- pseudobulk_matrix_x[x_gene, ]
  x_vec <- as.numeric((x_counts / total_UMI_x) * 1e6)
  x_vec <- log2(x_vec + 1)
  names(x_vec) <- colnames(pseudobulk_matrix_x)
  
  # ----------------------{你原本的 code：取 B naive cells}----------------------
  y_cells <- rownames(merge_data@meta.data)[merge_data@meta.data$sub_celltype %in% y_celltype]
  if (length(y_cells) == 0) {
    message("[SKIP] y_celltype no cells: ", y_celltype)
    return(invisible(NULL))
  }
  y_merge_data <- subset(merge_data, cells = y_cells)
  
  # ----------------------{pseudobulk Y}----------------------
  pseudobulk_rna_y <- AggregateExpression(
    y_merge_data,
    group.by = "orig.ident",
    assays   = "RNA",
    fun      = "sum",
    return.seurat = FALSE
  )
  pseudobulk_matrix_y <- pseudobulk_rna_y$RNA
  pseudobulk_matrix_y <- as(pseudobulk_matrix_y, "dgCMatrix")
  
  total_UMI_y <- Matrix::colSums(pseudobulk_matrix_y)
  total_UMI_y[total_UMI_y == 0] <- 1
  
  if (!(y_gene %in% rownames(pseudobulk_matrix_y))) {
    message("[SKIP] ", x_celltype, ": y_gene not found in Y pseudobulk: ", y_gene)
    return(invisible(NULL))
  }
  
  y_counts <- pseudobulk_matrix_y[y_gene, ]
  y_vec <- as.numeric((y_counts / total_UMI_y) * 1e6)
  y_vec <- log2(y_vec + 1)
  names(y_vec) <- colnames(pseudobulk_matrix_y)
  
  # ----------------------{取 Sample & group}----------------------
  patient_group <- merge_data@meta.data %>%
    distinct(patient = orig.ident, group = Condition)
  
  # ----------------------{交集皆有數值的samples}----------------------
  common_patients <- intersect(names(x_vec), names(y_vec))
  if (length(common_patients) < 3) {
    message("[SKIP] ", x_celltype, ": <3 common patients")
    return(invisible(NULL))
  }
  
  df_xy <- data.frame(
    patient = common_patients, #同時有 x_celltype 和 y_celltype 的 sample
    X = x_vec[common_patients],
    Y = y_vec[common_patients],
    stringsAsFactors = FALSE
  ) %>%
    left_join(patient_group, by = "patient")
  
  # ----------------------{Spearman correlation}----------------------
  result <- df_xy %>%
    group_by(group) %>%
    summarise(
      n = n(),
      r = cor(X, Y, method = "spearman", use = "complete.obs"),
      p = cor.test(X, Y, method = "spearman")$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0("n = ", n,
                     "\nr = ", sprintf("%.3f", r),
                     "\np = ", format.pval(p, digits = 2, eps = 1e-3))
    )
  
  # ----------------------{Plot:NMOSD 和 Control 各自畫}----------------------
  base_theme <- theme_classic(base_size = 14) +
    theme(
      plot.title  = element_text(size = 13, hjust = 0.5, face = "bold"),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10)
    )
  
  x_scale <- scale_x_continuous(limits = x_lim, breaks = seq(x_lim[1], x_lim[2], by = 2))
  y_scale <- scale_y_continuous(limits = y_lim, breaks = seq(y_lim[1], y_lim[2], by = 2))
  
  # NMOSD
  df_nmosd <- df_xy %>% filter(group == "NMOSD")
  lab_nmosd <- result %>% filter(group == "NMOSD")
  
  p_nmosd <- ggplot(df_nmosd, aes(X, Y)) +
    geom_point(size = 2.6, alpha = 0.8, color = "#E3A19F") +
    labs(
      title = "BAFF and BAFFR correlation - NMOSD",
      x = paste0(x_gene, " in ", x_celltype, " (log2CPM)"),
      y = paste0(y_gene, " in ", y_celltype, " (log2CPM)")
    ) +
    x_scale + y_scale + base_theme +
    annotate("text",
             x = x_lim[2] - 2, y = y_lim[1] + 2, hjust = 0, vjust = 1,
             label = ifelse(nrow(lab_nmosd) > 0, lab_nmosd$label, "NA"),
             size = 4)
  
  # Control
  df_ctrl <- df_xy %>% filter(group == "Control")
  lab_ctrl <- result %>% filter(group == "Control")
  
  p_ctrl <- ggplot(df_ctrl, aes(X, Y)) +
    geom_point(size = 2.8, alpha = 0.85, color = "#67D6F0") +
    labs(
      title = "BAFF and BAFFR correlation - Control",
      x = paste0(x_gene, " in ", x_celltype, " (log2CPM)"),
      y = paste0(y_gene, " in ", y_celltype, " (log2CPM)")
    ) +
    x_scale + y_scale + base_theme +
    annotate("text",
             x = x_lim[2] - 2, y = y_lim[1] + 2, hjust = 0, vjust = 1,
             label = ifelse(nrow(lab_ctrl) > 0, lab_ctrl$label, "NA"),
             size = 4)
  
  # ----------------------{存檔}----------------------
  safe_tag <- function(s) {
    s %>%
      str_replace_all("\\s+", "_") %>%
      str_replace_all("[/\\\\:*?\"<>|]+", "_")
  }
  
  tag <- paste0(
    safe_tag(x_gene), "_", safe_tag(x_celltype),
    "__", safe_tag(y_gene), "_", safe_tag(y_celltype)
  )
  
  ggsave(file.path(out_dir, paste0(tag, "__NMOSD.png")), p_nmosd, width = 5.5, height = 4.8, dpi = 300)
  ggsave(file.path(out_dir, paste0(tag, "__Control.png")), p_ctrl,  width = 5.5, height = 4.8, dpi = 300)
  
  message("[OK] ", x_celltype, " saved to ", out_dir)
  
  invisible(list(df_xy = df_xy, result = result, p_nmosd = p_nmosd, p_ctrl = p_ctrl))
}

# ----------------------------
# 跑所有 x_celltype 的 wrapper
# ----------------------------
run_all_x_celltypes <- function(
    merge_data,
    x_gene = "TNFSF13B",
    y_gene = "TNFRSF13C",
    y_celltype = "B naive",
    out_dir = "BAFF_BAFFR_corr_plots"
) {
  all_celltypes <- sort(unique(merge_data@meta.data$sub_celltype))
  all_celltypes <- all_celltypes[!is.na(all_celltypes)]
  
  res <- list()
  for (ct in all_celltypes) {
    res[[ct]] <- tryCatch(
      run_BAFF_BAFFR_one_Xcelltype (
        merge_data = merge_data,
        x_gene = x_gene,
        x_celltype = ct,
        y_gene = y_gene,
        y_celltype = y_celltype,
        out_dir = out_dir
      ),
      error = function(e) {
        message("[ERR] ", ct, " : ", e$message)
        NULL
      }
    )
  }
  invisible(res)
}
run_all_x_celltypes(
  merge_data = merge_data,
  x_gene = "TNFSF13B",
  y_gene = "TNFRSF13C",
  y_celltype = "B naive",
  out_dir = "BAFF_BAFFR_allX"
)

