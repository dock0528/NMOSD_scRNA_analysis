#%%
#-----------------------------------【Install packages】----------------------------------

import scanpy as sc
import networkx as nx
import community   # pip install python-louvain
import pandas as pd
import numpy as np
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import seaborn as sns

#%%
#------------------------------【使用Networkx進行louvain分群】------------------------------
#louvain原理:讓群內相似度越大，群外相似度越小越好
def louvain_clustering(
    adata, 
    n_list=[10, 15, 20,25,30,35,40,45,50],#k個鄰居
    pca_dim=14,#使用14個PCs個PCs
    key_prefix="louvain_n",#分群結果寫入 adata.obs 的欄位名稱前綴
    group=None
):
    summary = []
    #取14PCs存在一個新欄位
    adata.obsm["X_louvain"] = adata.obsm["X_pca"][:, :pca_dim]

    for n in n_list:
        sc.pp.neighbors(adata, n_neighbors=n, use_rep="X_louvain") #return adata.obsp["connectivities"] & adata.obsp["distances"]
        key = f"{key_prefix}{n}"
        sc.tl.louvain(adata,
                    resolution=1.0,    # 固定用 1.0
                    key_added=key,
                    random_state=0)
        labels = adata.obs[key].astype(int)
        k_cls = labels.nunique()#分幾群
        modularity = community.modularity(
            dict(enumerate(labels)), nx.Graph(adata.obsp["connectivities"])
        ) # dict(enumerate(labels)->{cell:群}
        summary.append({"n": n, "k_clusters": k_cls, "modularity": modularity})
        print(f"clusters: {k_cls}, modularity: {modularity:.4f}")

    df = pd.DataFrame(summary)

    # 畫群數 & modularity
    fig, ax1 = plt.subplots(figsize=(8,5))
    ax1.plot(df["n"], df["k_clusters"], "o-", label="Clusters", color='#41C0F2')
    ax1.set_xlabel("n neighbors")
    ax1.set_ylabel("Number of Clusters")
    ax2 = ax1.twinx() #雙軸
    ax2.plot(df["n"], df["modularity"], "s--", label="Modularity", color='#C491F5')
    ax2.set_ylabel("Modularity Score")
    plt.title(f"{group} - Louvain Clustering vs n neighbors")
    fig.tight_layout()
    plt.show()

    
    return adata

#%%
#------------------------------【louvaing使用不同resolutions >>> silhouette score評估】------------------------------
def louvain_resolution_silhouette(
    adata,
    resolution_list=[0.6, 0.8, 1.0, 1.2, 1.5],
    pca_dim=14,
    n=25, #neighbors
    rep_key="louvain_kbest"
):
    # 建立 PCA 表示空間
    adata.obsm[rep_key] = adata.obsm["X_pca"][:, :pca_dim]
    sc.pp.neighbors(adata, n_neighbors=n, use_rep=rep_key)

    summary = []

    for res in resolution_list:
        print(f"\n Louvain resolution = {res}")
        key = f"louvain_res{res}"
        sc.tl.louvain(adata, resolution=res, key_added=key,random_state=0) #key_added->add the cluster labels
        #當前分群結果寫入 adata.obs["louvain_res{res}"]

        labels = adata.obs[key].astype(int)
        n_clusters = labels.nunique()

        # 當只有 1 群，無法計算 silhouette score
        if n_clusters <= 1 or n_clusters >= len(labels):
            sil_score = float("nan")
        else:
            sil_score = silhouette_score(adata.obsm[rep_key], labels)

        summary.append({
            "resolution": res,
            "n_clusters": n_clusters,
            "silhouette_score": sil_score
        })
        print(f"  -> Clusters: {n_clusters}, Silhouette Score: {sil_score:.4f}")

    df = pd.DataFrame(summary)

    # plot
    fig, ax1 = plt.subplots(figsize=(8, 5))

    ax1.plot(df["resolution"], df["n_clusters"], "o-", label="Clusters", color="#41C0F2")
    ax1.set_xlabel("Louvain Resolution")
    ax1.set_ylabel("Number of Clusters")
    ax1.tick_params(axis="y")

    ax2 = ax1.twinx()
    ax2.plot(df["resolution"], df["silhouette_score"], "s--", label="Silhouette Score", color="#F57A7A")
    ax2.set_ylabel("Silhouette Score")
    ax2.tick_params(axis="y")

    plt.title(f"Louvain Resolution  (n={n})")
    fig.tight_layout()
    plt.show()

    return df


#%%
# ---------------------【主程式：使用合併 pre + post 進行Louvain】----------------
adata = sc.read_h5ad("./adata_hvg_rawcount_pca50.h5ad")
adata= louvain_clustering(
    adata,
    n_list=[10, 15, 20,25,30,35,40,45,50],
    pca_dim=14,
    key_prefix="louvain_k",
    group="Pre/Post"
)

#%%
# ---------------------【主程式：使用合併 pre + post 調整Louvain-resolution】----------------
result_df = louvain_resolution_silhouette(
    adata,
    resolution_list=[0.6, 0.8, 1.0, 1.2, 1.5],
    pca_dim=14,
    n=25, #neighbors=25,cluster=16
    rep_key="louvain_kbest",
)
#%%
# ------------------------------【使用最佳 k 分群結果】-------------------------------
#由Moduality開始趨緩的k值畫UMAP >>> 由pre/post選k=15 
#Moduality範圍[-1,1]
best_n=25 #n取25，resolution:0.8，分群:14群
best_cluster_number = 14
best_res=0.8
adata.obs["cluster"] = adata.obs[f"louvain_res{best_res}"].astype(int)

#%%
#-----------------------------------【UMAP圖 標記pre/post】-----------------------------------
sc.tl.umap(adata)
sc.pl.umap(
    adata,
    color="dataset",
    title=f"Pre/Post UMAP",
    legend_loc="right margin",
    legend_fontsize=8,
    frameon=True,#圖加框
    legend_fontoutline=1, #圖例設置為圓點+標籤
    palette=["#A1D4E8","#EACDE4"]
)

#%%
#-----------------------------------【UMAP圖 標記cluster】-----------------------------------
colors_30 = sns.color_palette("tab20") + sns.color_palette("Paired")[:10]
adata.obs["cluster"] = adata.obs["cluster"].astype(str) #必須轉成字串在圖例時才不會呈現連續的bar
# 取得所有 cluster label（包含 pre+post）
unique_clusters = sorted(adata.obs["cluster"].unique()) 

# 指定固定顏色對應
colors_30 = sns.color_palette("tab20") + sns.color_palette("Paired")[:10]
cluster_color_map = dict(zip(unique_clusters, colors_30[:len(unique_clusters)]))
sc.tl.umap(adata)
sc.pl.umap(
    adata,
    color="cluster",
    title=f"Pre/Post Louvain (cluster={best_cluster_number})",
    legend_loc="right margin",
    legend_fontsize=6,
    frameon=True,#圖加框
    legend_fontoutline=1, #圖例設置為圓點+標籤
    palette=cluster_color_map 
)
#%%
#---------------------------【畫 pre/post各自UMAP】------------------------------
sc.tl.umap(adata)
sc.pl.umap(
    adata[adata.obs["dataset"] == "pre"],
    color="cluster",
    title=f"Pre Louvain (cluster={best_cluster_number})",
    legend_loc="right margin",
    legend_fontsize=6,
    frameon=True,#圖加框
    legend_fontoutline=1, #圖例設置為圓點+標籤
    palette=cluster_color_map 
)
sc.pl.umap(
    adata[adata.obs["dataset"] == "post"],
    color="cluster",
    title=f"Post Louvain (cluster={best_cluster_number})",
    legend_loc="right margin",
    legend_fontsize=6,
    frameon=True,
    legend_fontoutline=1, 
    palette=cluster_color_map
)

#%%
#---------------------------【 畫 pre/post 在各 cluster 的佔比】------------------------------
comb_df = adata.obs[["cluster", "dataset"]].rename(columns={"dataset": "group"})
ct = pd.crosstab(comb_df["cluster"].astype(int), comb_df["group"], normalize="index")
ct = ct.sort_index() #排序

#pre 藍色，post 粉色
colors = {"pre": "#A1D4E8", "post": "#EACDE4"}

# 畫圖
ax = ct.plot(
    kind="bar",
    stacked=True,
    color=[colors[col] for col in ct.columns],
    figsize=(12, 5)
)
plt.legend(
    title="Group",
    loc="center left",
    bbox_to_anchor=(1.0, 0.5),
    fontsize=10
)
plt.xlabel("Cluster", fontsize=12)
plt.ylabel("Proportion", fontsize=12)
plt.title(f"Pre/Post cluster composition (cluster={best_cluster_number})", fontsize=14)
plt.xticks(rotation=0)
plt.tight_layout()
plt.show()

#%%
#存adata
#adata.write("../scANVI_ref_dataset/adata_cluster_v2.h5ad")
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
#-----------------------------------------------【將參考集和自己的資料合併】----------------------------------
adata = sc.read_h5ad("../scANVI_ref_dataset/adata_cluster_v2.h5ad")#我的已分群後的adata

# 資料夾 GSE164378_RAW/底下，並改名為 barcodes.tsv、features.tsv、matrix.mtx
ref = sc.read_10x_mtx(
    "../scANVI_ref_dataset/GSE164378_RAW/",                # 資料夾路徑
    var_names="gene_symbols",        # features 檔裡第二欄是 gene_symbols
    cache=True                       # 讀完會把轉好的 h5ad 存起來
)
# %%
import pandas as pd

#讀取已有的cell type reference
meta = pd.read_csv("../scANVI_ref_dataset/GSE164378_sc.meta.data_3P.csv", index_col=0)


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
# %%
#訓練scANVI
#【已跑過不用再重跑!!!】
#scanvi = SCANVI(adata_merged)          
# scanvi.train(
#     max_epochs=400, #default=400
#     accelerator='gpu'       
# )

#%%
#存已訓練好的模型!!!
#scanvi.save("scanvi_model/")

#-----------------------------------------【Predicting scANVI】-----------------------------------------
# %%
#載入模型
scanvi_loaded = SCANVI.load("./scanvi_model/", adata=adata_merged) #!!!放入要用的adata

#模型預測
adata_merged.obs["predicted_label"] = scanvi_loaded.predict(adata_merged)
#_scvi_labels 原始label
#predicted_label 訓練後的label(用他!)

#-----------------------------------------【Mapping cell type】-----------------------------------------
# %%
cluster_key = "louvain_res0.8" #原始分群label
label_key = "predicted_label" #scVANI預測的結果

#計算每個 cluster裡最多的cell-type
cluster_annotation = (
    adata_merged.obs
        .groupby(cluster_key)[label_key] #以cluster_key做分群 #label_key做對照cell type
        .agg(lambda x: x.value_counts().idxmax()) #x.value_counts()：對同一個 cluster 裡的所有預測標籤做計數 #idxmax()取最多次的scANVI cell type
        .to_dict()
)

# %%
#將annotation後的cell type mapping 回原始資料
adata_merged.obs["cluster_annotation"] = adata_merged.obs[cluster_key].map(cluster_annotation)

# %%
#--------------------------------------------【UMAP圖 標記cluster】-----------------------------------------
adata_query = adata_merged[adata_merged.obs["dataset"] == "query"].copy()
adata_query.obsm["X_pca"]=adata.obsm["X_pca"] #shape=(24852, 50)
adata_query.obsm['scANVI_pca'] = adata_query.obsm["X_pca"][:, :14]

#%%
colors_30 = sns.color_palette("tab20") + sns.color_palette("Paired")[:10]
adata_query.obs["cluster_annotation"] = adata_query.obs["cluster_annotation"].astype(str) #必須轉成字串在圖例時才不會呈現連續的bar
# 取得所有 cluster label（包含 pre+post）
unique_clusters = sorted(adata_query.obs["cluster_annotation"].unique()) 

# 指定固定顏色對應
cluster_color_map = dict(zip(unique_clusters, colors_30[:len(unique_clusters)]))
sc.pp.neighbors(adata_query,n_neighbors=25, use_rep="scANVI_pca")
sc.tl.umap(adata_query)
sc.pl.umap(
    adata_query,
    color="cluster_annotation",
    title=f"Pre/Post Louvain after scANVI annotation",
    legend_loc="right margin",
    legend_fontsize=6,
    frameon=True,#圖加框
    legend_fontoutline=1, #圖例設置為圓點+標籤
    palette=cluster_color_map 
)
# %%
#-----------------------------------------【畫 pre/post各自UMAP】-----------------------------------------
adata_query.obs["dataset"]=adata.obs["dataset"]
sc.pl.umap(
    adata_query[adata_query.obs["dataset"] == "pre"],
    color="cluster_annotation",
    title=f"Pre Louvain after scANVI annotation",
    legend_loc="right margin",
    legend_fontsize=6,
    frameon=True,#圖加框
    legend_fontoutline=1, #圖例設置為圓點+標籤
    palette=cluster_color_map 
)
sc.pl.umap(
    adata_query[adata_query.obs["dataset"] == "post"],
    color="cluster_annotation",
    title=f"Post Louvain after scANVI annotation",
    legend_loc="right margin",
    legend_fontsize=6,
    frameon=True,
    legend_fontoutline=1, 
    palette=cluster_color_map
)

#%%
#存scANVI後的adata
#adata_query.write("../scANVI_ref_dataset/adata_scANVI_annotation.h5ad")
# %%
#-----------------------------------------【 畫 pre/post 在各 cluster 的佔比】-----------------------------------------
comb_df = adata_query.obs[["cluster_annotation", "dataset"]].rename(columns={"dataset": "group"})
ct = pd.crosstab(comb_df["cluster_annotation"], comb_df["group"], normalize="index")
ct = ct.sort_index() #排序

#pre 藍色，post 粉色
colors = {"pre": "#A1D4E8", "post": "#EACDE4"}

# 畫圖
ax = ct.plot(
    kind="bar",
    stacked=True,
    color=[colors[col] for col in ct.columns],
    figsize=(12, 5)
)
plt.legend(
    title="Group",
    loc="center left",
    bbox_to_anchor=(1.0, 0.5),
    fontsize=10
)
plt.xlabel("Cluster", fontsize=12)
plt.ylabel("Proportion", fontsize=12)
plt.title(f"Pre/Post cluster composition after scANVI annotation ", fontsize=14)
plt.xticks(rotation=0)
plt.tight_layout()
plt.show()

#%%
#存adata
#adata_query.write("../scANVI_ref_dataset/adata_annotation.h5ad")
#-----------------------------------------【 Adjusted Rand Index (ARI)】-----------------------------------------
#range:(-1,1)
# %%
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

#%%
adata_ann = sc.read_h5ad("../scANVI_ref_dataset/adata_annotation.h5ad")

#predict label
labels_true = adata_ann.obs["predicted_label"].values

#true label
labels_pred = adata_ann.obs["louvain_res0.8"].astype(str)
ari_score = adjusted_rand_score(labels_true, labels_pred)
print(f'Louvain ARI score:{ari_score:.4f}')

#-----------------------------------------【Normalized Mutual Information (NMI)】-----------------------------------------
#range:(0,1)
# %%
from sklearn.metrics.cluster import normalized_mutual_info_score
nmi_score=normalized_mutual_info_score(labels_true, labels_pred)
print(f'Louvain NMI score:{nmi_score:.4f}')
# %%
