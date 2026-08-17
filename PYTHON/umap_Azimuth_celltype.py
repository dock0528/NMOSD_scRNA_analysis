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
    def __init__(self, adata_predicted, adata_clustered,use_rep='X_pca',ndim_pca=50, base_output="../Dataset/My_merged_protein_coding_genes"):
        """
        初始化
        adata_predicted: 含有 predicted_label 與 dataset='query' 的 AnnData
        adata_clustered: 含有 PCA 結果的 AnnData（通常是 Leiden clustering 後）
        base_output: 根目錄，底下會建立三個子資料夾：
                     Leiden_scANVI_plot / Leiden_scANVI_evaluation / Leiden_scANVI_celltype_composition
        """
        self.adata_predicted = adata_predicted
        self.adata_clustered = adata_clustered
        self.base_output = base_output
        self.adata_query = None
        self.ndim_pca=ndim_pca
        self.use_rep=use_rep

        # 建立輸出資料夾
        if self.use_rep=='X_pca': # X_pca
            self.plot_dir = os.path.join(base_output, "Leiden_scANVI_plot",f'PC{self.ndim_pca}')
            self.eval_dir = os.path.join(base_output, "Leiden_scANVI_evaluation",f'PC{self.ndim_pca}')
            self.comp_dir = os.path.join(base_output, "Leiden_scANVI_celltype_composition",f'PC{self.ndim_pca}')
            
        else: # X_harmony
            self.plot_dir = os.path.join(base_output, "Harmony_Azimuth_plot",f'Harmony{self.ndim_pca}')
            self.eval_dir = os.path.join(base_output, "Harmony_Azimuth_evaluation",f'Harmony{self.ndim_pca}')
            self.comp_dir = os.path.join(base_output, "Harmony_Azimuth_celltype_composition",f'Harmony{self.ndim_pca}')
        os.makedirs(self.plot_dir, exist_ok=True)
        os.makedirs(self.eval_dir, exist_ok=True)
        os.makedirs(self.comp_dir, exist_ok=True)


    #----------------------------【Cell-type mapping】----------------------------
    def map_celltype(self, cluster_key="leiden_nn20_res0.6", label_key="predicted_label"):
        cluster_annotation = (
            self.adata_predicted.obs
                .groupby(cluster_key)[label_key]
                .agg(lambda x: x.value_counts().idxmax())
                .to_dict()
        )
        self.adata_predicted.obs["cluster_annotation"] = self.adata_predicted.obs[cluster_key].map(cluster_annotation)
        print(f"Mapped {cluster_key} cell types done!")

    #----------------------------【UMAP 計算】----------------------------
    def prepare_umap(self, n_neighbors=20):
        if "dataset" in self.adata_predicted.obs.columns:
            adata_query = self.adata_predicted[self.adata_predicted.obs["dataset"] == "query"].copy()
        else:
            adata_query = self.adata_predicted.copy()
        if self.use_rep=='X_pca':
            adata_query.obsm['X_pca'] = self.adata_clustered.obsm['X_pca']
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
    def plot_umap_and_composition(self, cluster_key="leiden_nn20_res0.6",label_key="predicted_label"):
        if self.adata_query is None:
            raise ValueError("請先執行 prepare_umap()")

        sub_dir = os.path.join(self.plot_dir,cluster_key,label_key)
        os.makedirs(sub_dir, exist_ok=True)

        adata_query = self.adata_query
        colors_30 = sns.color_palette("tab20") + sns.color_palette("Paired")[:10]
        if "cluster_annotation" in self.adata_predicted.obs.columns:
            adata_query.obs["cluster_annotation"] = adata_query.obs["cluster_annotation"].astype(str)
        else:
            adata_query.obs["cluster_annotation"] = adata_query.obs[cluster_key]
        unique_clusters = sorted(adata_query.obs["cluster_annotation"].unique())
        cluster_color_map = dict(zip(unique_clusters, colors_30[:len(unique_clusters)]))

        # --- UMAP圖 標記cluster ---
        sc.pl.umap(
            adata_query,
            color="cluster_annotation",
            title="NMOSD/Control Leiden after Azimuth annotation",
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
                title=f"{cond} Leiden after Azimuth annotation",
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

        plt.figure(figsize=(22, 10))
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
        plt.xticks(rotation=45, fontsize=14)
        plt.tight_layout()
        plt.savefig(f"{sub_dir}/Composition.png", dpi=300)
        plt.close()
        print(f"All plots saved to: {sub_dir}")

    #----------------------------【ARI & NMI &Ｖ-measure評估指標】----------------------------
    def evaluate_clustering(self, cluster_key="leiden_nn20_res0.6",label_key="predicted_label"):
        if self.adata_query is None:
            raise ValueError("請先執行 prepare_umap()")

        labels_true = self.adata_query.obs[label_key].values
        labels_pred = self.adata_query.obs[cluster_key].astype(str).values

        ari_score = adjusted_rand_score(labels_true, labels_pred)
        nmi_score = normalized_mutual_info_score(labels_true, labels_pred)
        v_measure = v_measure_score(labels_true, labels_pred)

        metrics_df = pd.DataFrame({
            "Evaluation Index": ["ARI", "NMI", "V-measure"],
            "Score": [ari_score, nmi_score, v_measure]
        })
        out_path = f"{self.eval_dir}/{cluster_key}_{label_key}.csv"
        metrics_df.to_csv(out_path, index=False)

        print(f"ARI: {ari_score:.4f}, NMI: {nmi_score:.4f}, V-measure: {v_measure:.4f}")
        print(f"ARI & NMI & Ｖ-measure evaluation points saved to: {out_path}")
        return metrics_df

    #----------------------------【Cell-type composition + purity】----------------------------
    def export_celltype_composition(self, cluster_key="leiden_nn20_res0.6",label_key="predicted_label"):
        """
        計算每個 cluster 中各 cell type 的比例與 purity
        儲存至 Leiden_scANVI_celltype_composition/{cluster_key}.csv
        """
        if "cluster_annotation" not in self.adata_predicted.obs.columns:
            raise ValueError("請先執行 map_celltype()")

        comp_df = (
            self.adata_predicted.obs.groupby([cluster_key, label_key])
            .size()
            .reset_index(name="count")
        )

        # 計算比例與 purity
        comp_df["proportion"] = comp_df.groupby(cluster_key)["count"].transform(lambda x: x / x.sum())
        composition_df = (
            comp_df.groupby(cluster_key)
            .apply(lambda df: pd.Series({
                "Dominant_celltype": df.loc[df["proportion"].idxmax(), label_key],
                "Percentage": df["proportion"].max()
            }))
            .reset_index()
        )

        out_path = f"{self.comp_dir}/{cluster_key}_{label_key}.csv"
        composition_df.to_csv(out_path, index=False)
        print(composition_df)
        print(f"Cell-type composition saved to: {out_path}")
        return self.adata_query
    
    