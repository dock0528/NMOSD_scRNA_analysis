#################################################【scANVI Annotation】###########################################
#%%
#-----------------------------------------------【Install packages】----------------------------------
import scanpy as sc
import networkx as nx
import community   # python-louvain
import pandas as pd
import numpy as np
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import seaborn as sns
#-----------------------------------------------【將參考集和自己的資料合併】----------------------------------
#%%
adata = sc.read_h5ad(r"C:\Users\Jane\Desktop\Wang實驗室\NMOSD研究計畫\scRNA\NMOSD_scRNA_analysis\adata_cluster_v2.h5ad")#我的已分群後的adata

#%%
# 資料夾 GSE164378_RAW/底下，並改名為 barcodes.tsv、features.tsv、matrix.mtx
ref = sc.read_10x_mtx(
    "C:/Users/Jane/Desktop/Wang實驗室/NMOSD研究計畫/scRNA/NMOSD_scRNA_analysis/NMOSD_scRNA_DATA/Reference/GSE164378_RAW_renamed",                # 資料夾路徑
    var_names="gene_symbols",        # features 檔裡第二欄是 gene_symbols
    cache=True                       # 讀完會把轉好的 h5ad 存起來
)
# %%
import pandas as pd

#讀取已有的cell type reference
meta = pd.read_csv("C:/Users/Jane/Desktop/Wang實驗室/NMOSD研究計畫/scRNA/NMOSD_scRNA_analysis/NMOSD_scRNA_DATA/Reference/GSE164378_sc.meta.data_3P.csv", index_col=0)


# 把最細層的 celltype.l3 當作真標籤放進 ref.obs
# 確保 barcodes 一一對應
ref.obs["celltype"] = meta.loc[ref.obs_names, "celltype.l3"].astype(str)

# %%
#common genes
adata.var_names #adata gene name
common_genes = ref.var_names.intersection(adata.var_names)
ref = ref[:, common_genes].copy() #reference
query = adata[:, common_genes].copy() #自己的

# %%
#參考集 & 我的cell type做合併
adata_merged = ref.concatenate(
    query,
    batch_key="dataset",
    batch_categories=["ref", "query"],
    index_unique=None  # 保留原本的 barcode index
) #1213 geneid

#%%
#標記自己的cell type欄位為unknown
adata_merged.obs.loc[adata_merged.obs["dataset"] == "query", "celltype"]= "Unknown"

#%%
# 先把現存的這一欄轉成 categorical (重要!!!)
adata_merged.obs["celltype"] = adata_merged.obs["celltype"].astype("category")


# -----------------------------------------【scANVI preprocessing】-----------------------------------------
# %%
#pip install scvi-tools
import scvi
from scvi.model import SCANVI

scvi.model.SCANVI.setup_anndata(
    adata_merged,
    batch_key="dataset",  #concat欄位
    labels_key="celltype",  # reference 的真實標籤
    unlabeled_category="Unknown" # query 的標籤
)
#-----------------------------------------【Training scANVI】-----------------------------------------
#%%
#訓練scANVI
#【已跑過不用再重跑!!!】
scanvi = SCANVI(adata_merged)          
scanvi.train(
    max_epochs=400, #default=400
    accelerator='gpu'       
)

#%%
#存已訓練好的模型!!!
#scanvi.save("scanvi_model/")