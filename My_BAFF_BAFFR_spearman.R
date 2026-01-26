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
