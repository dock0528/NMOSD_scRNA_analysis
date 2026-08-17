#----{Packages}
import os
import scanpy as sc
import networkx as nx
import community   # pip install python-louvain
import pandas as pd
import numpy as np
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import seaborn as sns
import scvi
from scvi.model import SCANVI
from sklearn.metrics import adjusted_rand_score,v_measure_score
from sklearn.metrics.cluster import normalized_mutual_info_score

#----{Class}
class AZIMUTHvisualizer:
    def __init__(self, adata_celltype,use_rep='X_pca',ndim_pca=50, base_output="../Dataset/My_merged_protein_coding_genes"):
        """
        初始化
        adata_predicted: 含有 predicted_label 與 dataset='query' 的 AnnData
        adata_clustered: 含有 PCA 結果的 AnnData（通常是 Leiden clustering 後）
        base_output: 根目錄，底下會建立三個子資料夾：
                     Azimuth_plot / Azimuth_evaluation / Azimuth_celltype_composition
        """
        self.adata_celltype =  adata_celltype
        self.base_output = base_output

        self.adata_query = None
        self.ndim_pca=ndim_pca
        self.use_rep=use_rep

        # 建立輸出資料夾
        if self.use_rep=='X_pca': # X_pca
            self.plot_dir = os.path.join(base_output, "Azimuth_plot",f'PC{self.ndim_pca}')
            self.eval_dir = os.path.join(base_output, "Azimuth_evaluation",f'PC{self.ndim_pca}')
            self.comp_dir = os.path.join(base_output, "Azimuth_celltype_composition",f'PC{self.ndim_pca}')
              
        else: # X_harmony
            self.plot_dir = os.path.join(base_output, "Azimuth_plot",f'Harmony{self.ndim_pca}')
            #self.eval_dir = os.path.join(base_output, "Azimuth_evaluation",f'Harmony{self.ndim_pca}')
            #self.comp_dir = os.path.join(base_output, "Azimuth_celltype_composition",f'Harmony{self.ndim_pca}')
        os.makedirs(self.plot_dir, exist_ok=True)
        #os.makedirs(self.eval_dir, exist_ok=True)
        #os.makedirs(self.comp_dir, exist_ok=True)


    #----------------------------【UMAP 計算】----------------------------
    def prepare_umap(self, n_neighbors=20):
        adata_query = self.adata_celltype.copy()
        if self.use_rep=='X_pca':
            adata_query.obsm["scANVI_pca"] = adata_query.obsm["X_pca"][:, :self.ndim_pca]
            sc.pp.neighbors(adata_query, n_neighbors=n_neighbors, use_rep="scANVI_pca")
        else:
            adata_query.obsm["azimuth_harmony"] = adata_query.obsm["X_harmony"]
            sc.pp.neighbors(adata_query, n_neighbors=n_neighbors, use_rep="azimuth_harmony")

        
        sc.tl.umap(adata_query)
        self.adata_query = adata_query
        print(f"UMAP computed using {n_neighbors} neighbors and {self.ndim_pca} PCs.")
    #----------------------------【存高dpi圖】----------------------------
    def save_with_dpi(self,sc_path, dst_path, dpi=500):
        img = plt.imread(sc_path)
        fig, ax = plt.subplots(figsize=(10, 8), dpi=dpi)
        ax.imshow(img)
        ax.axis("off")
        plt.tight_layout()
        plt.savefig(dst_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        os.remove(sc_path)

    #----------------------------【UMAP繪圖與儲存】----------------------------
    def plot_umap_and_composition(self, cluster_key="predicted.celltype.l2"):
        if self.adata_query is None:
            raise ValueError("請先執行 prepare_umap()")

        sub_dir = os.path.join(self.plot_dir,cluster_key)
        os.makedirs(sub_dir, exist_ok=True)

        adata_query = self.adata_query
        colors_30 = sns.color_palette("tab20") + sns.color_palette("Paired")[:10]
        adata_query.obs["cluster_annotation"] = adata_query.obs[cluster_key]
        unique_clusters = sorted(adata_query.obs["cluster_annotation"].unique())
        cluster_color_map = dict(zip(unique_clusters, colors_30[:len(unique_clusters)]))

        # --- UMAP圖 標記cluster ---
        sc.pl.umap(
            adata_query,
            color="cluster_annotation",
            title="NMOSD/Control after Azimuth annotation",
            legend_loc="right margin",
            legend_fontsize=6,
            frameon=True,
            legend_fontoutline=1,
            palette=cluster_color_map,
            save=f"_{cluster_key}_All.png",
            show=False
            )
        sc_path = f"figures/umap_{cluster_key}_All.png"
        dst_path = f"{sub_dir}/UMAP_All.png"
        self.save_with_dpi(sc_path, dst_path, dpi=500)

        # --- UMAP圖 標記samples(batch) ---
        # 設定顏色(sample)
        adata_query.obs["orig.ident"] = adata_query.obs["orig.ident"].astype(str)
        unique_batches = sorted(adata_query.obs["orig.ident"].unique())
        batch_palette = sns.color_palette("husl", len(unique_batches))
        batch_color_map = dict(zip(unique_batches, batch_palette))
        sc.pl.umap(
            adata_query,
            color='orig.ident',
            title="UMAP colored by sample (after Harmony correction)",
            legend_loc="right margin",
            legend_fontsize=6,
            frameon=True,
            legend_fontoutline=1,
            palette=batch_color_map,
            save=f"_{cluster_key}_harmony_corrected_batch.png",
            show=False
            )
        sc_path = f"figures/umap_{cluster_key}_harmony_corrected_batch.png"
        dst_path = f"{sub_dir}/UMAP_harmony_corrected_batch.png"
        self.save_with_dpi(sc_path, dst_path, dpi=500)

        # --- 各組別 UMAP ---
        for cond in ["NMOSD", "Control"]:
            sc.pl.umap(
                adata_query[adata_query.obs["Condition"] == cond],
                color="cluster_annotation",
                title=f"{cond} after Azimuth annotation",
                legend_loc="right margin",
                legend_fontsize=6,
                frameon=True,
                legend_fontoutline=1,
                palette=cluster_color_map,
                save=f"_{cluster_key}_{cond}.png",
                show=False
            )
            sc_path = f"figures/umap_{cluster_key}_{cond}.png"
            dst_path = f"{sub_dir}/UMAP_{cond}.png"
            self.save_with_dpi(sc_path, dst_path, dpi=500)

        # --- 畫組別在各 cluster 的佔比 ---
        #comb_df = adata_query.obs[["cluster_annotation", "Condition"]].rename(columns={"Condition": "group"})
        #ct = pd.crosstab(comb_df["cluster_annotation"], comb_df["group"], normalize="index").sort_index()
        adata_query.obs["cluster_annotation"] = (
        adata_query.obs["cluster_annotation"]
        .astype("string") # 轉成字串，但保留 NaN
        .str.strip() # 去除前後空白
        )

        df = adata_query.obs[["orig.ident","cluster_annotation","Condition"]].copy()
        df = df.rename(columns={"Condition":"group"})

        # 每個 sample×cluster 的細胞數（缺的補 0）
        mat = (
            df.pivot_table(index=["orig.ident","group"],
                        columns="cluster_annotation",
                        aggfunc="size",
                        fill_value=0)
        )

        # 在同一 group 內：對 sample 取平均（每個 sample 等權重）
        mean_counts = mat.groupby("group").mean().T  # (cluster × group)

        # 每個 cluster 內兩組比例 (加總為1)
        ct = mean_counts.div(mean_counts.sum(axis=1), axis=0).fillna(0)
        ct = ct.reindex(sorted(ct.index))


        colors = {"Control": "#67D6F0", "NMOSD": "#E3A19F"}

        plt.figure(figsize=(28, 10))
        ct.plot(
            kind="bar",
            stacked=True,
            color=[colors[col] for col in ct.columns],
            ax=plt.gca()
        )
        plt.legend(
            title="Group",
            loc="center left",
            bbox_to_anchor=(1.0, 0.5),
            title_fontsize=16,
            fontsize=14
        )
        plt.xlabel("Celltype", fontsize=16)
        plt.ylabel("Proportion", fontsize=16)
        plt.title("Proportion of NMOSD/Control cells in each celltype after Azimuth annotation", fontsize=22)
        plt.xticks(rotation=90, fontsize=14)
        plt.tight_layout()
        plt.savefig(f"{sub_dir}/Composition.png", dpi=300)
        plt.close()
        print(f"All plots saved to: {sub_dir}")
        return self.adata_query
    
    