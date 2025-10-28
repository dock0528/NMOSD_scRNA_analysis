#%%
#################################################【scANVI Annotation】###########################################
#-----------------------------------------------【Install packages】----------------------------------
import scanpy as sc
import networkx as nx
import community   # python-louvain
import pandas as pd
import numpy as np
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import seaborn as sns

#%%
#我的已分群後的adata
adata = sc.read_h5ad(r"C:\Users\JANE\Desktop\Wang實驗室\NMOSD研究計畫\scRNA\NMOSD_scANVI_DATA\Adata\Merge_PCA_2000HVG(Leiden)(raw_count).h5ad")
print(adata)

#%%
# 資料夾 GSE164378_RAW_renamed/底下，並改名為 barcodes.tsv、features.tsv、matrix.mtx
ref = sc.read_10x_mtx(
    r"C:\Users\JANE\Desktop\Wang實驗室\NMOSD研究計畫\scRNA\NMOSD_scANVI_DATA\Reference\GSE164378_RAW_renamed",                # 資料夾路徑
    var_names="gene_symbols",        # features 檔裡第二欄是 gene_symbols(已轉換成v32)
    cache=True                       # 讀完會把轉好的 h5ad 存起來
)
ref.layers["counts"] = ref.X.copy()

#%%
#讀取已有的cell type reference
meta = pd.read_csv(r"C:\Users\JANE\Desktop\Wang實驗室\NMOSD研究計畫\scRNA\NMOSD_scANVI_DATA\Reference\GSE164378_sc.meta.data_3P.csv", index_col=0)


# 把第二層的 celltype.l2 當作真標籤放進 ref.obs
# 確保 barcodes 一一對應
ref.obs["celltype"] = meta.loc[ref.obs_names, "celltype.l2"].astype(str)

#%%
#common genes
adata.var_names #adata gene name
common_genes = ref.var_names.intersection(adata.var_names)
ref_cp = ref[:, common_genes].copy() #reference
query = adata[:, common_genes].copy() #自己的

#%%
print(f'My adata gene counts: {len(adata.var_names)}')
print(f'Reference gene counts: {len(ref.var_names)}')
print(f'Intersection v32 gene counts: {len(common_genes)}')

#%%
#確認IFN-I皆有在內
IFN_I_list=['ISG15','IFI6','CMPK2','LY6E','OASL','AKAP12','TNFRSF13C']
pd.Series(IFN_I_list).isin(common_genes)

#%%
#參考集 & 我的cell type做合併
adata_merged = ref_cp.concatenate(
    query,
    batch_key="dataset",
    batch_categories=["ref", "query"],
    index_unique=None  # 保留原本的 barcode index
) #2000 geneid

#%%
#標記自己的cell type欄位為unknown
adata_merged.obs.loc[adata_merged.obs["dataset"] == "query", "celltype"]= "Unknown"

#%%
# 先把現存的這一欄轉成 categorical (重要!!!)
adata_merged.obs["celltype"] = adata_merged.obs["celltype"].astype("category")
print(adata_merged)

#%%
adata_merged.write_h5ad(r"C:\Users\JANE\Desktop\Wang實驗室\NMOSD研究計畫\scRNA\NMOSD_scANVI_DATA\adata_merged(Leiedn).h5ad")
#%%
#pip install scvi-tools==1.3.0
import scvi
from scvi.model import SCANVI

#%%
scvi.model.SCANVI.setup_anndata(
    adata_merged,
    layer="counts",
    batch_key="dataset",  #concat欄位
    labels_key="celltype",  # reference 的真實標籤
    unlabeled_category="Unknown" # query 的標籤
)

#%%
scanvi = SCANVI(adata_merged) 
scanvi.train(
    max_epochs=400, #default=400
    accelerator="gpu"
)
scanvi.save(r"C:\Users\JANE\Desktop\Wang實驗室\NMOSD研究計畫\scRNA\NMOSD_scANVI_DATA\scanvi_model")
# %%
