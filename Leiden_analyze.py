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

#----{Class}
class Leidenanalyzer:
    def __init__(self, adata_hvg, ndim_pca=50,use_rep='X_pca',base_output="../Dataset/My_merged_protein_coding_genes"):
        
        self.adata_hvg = adata_hvg
        self.base_output = base_output
        self.ndim_pca = ndim_pca
        self.use_rep=use_rep
        # 建立輸出資料夾
        if self.use_rep=='X_pca': # X_pca
            self.plot_dir = os.path.join(base_output, "Leiden_plots",f'PC{self.ndim_pca}')
        else: #X_harmony
            self.plot_dir= os.path.join(base_output,"Harmony_plots",f'Harmony{self.ndim_pca}')
        os.makedirs(self.plot_dir, exist_ok=True)
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
    #----------------------------【Leiden 分群】----------------------------
    def perform_leiden_modularity(self,nn_list = [10, 15, 20, 25, 30, 35, 40]):
        summary = []

        for k in nn_list:
            print(f"\n處理 k = {k}")
                
            # 建立鄰近圖 (knn)
            #使用PC50
            if self.use_rep == 'X_pca':
                sc.pp.neighbors(self.adata_hvg, n_neighbors=k, n_pcs=self.ndim_pca, use_rep=self.use_rep)
            #使用Harmony 
            else:
                sc.pp.neighbors(self.adata_hvg, n_neighbors=k, use_rep=self.use_rep)
                
            # 執行 Leiden clustering
            col = f"leiden_nn{k}"
            sc.tl.leiden(self.adata_hvg, key_added=col, random_state=0)
                
            # 計算分群資訊
            labels = self.adata_hvg.obs[col]
            label_dict = dict(zip(range(len(labels)), labels.astype(str)))

            # 用 connectivities 建 NetworkX graph
            conn = self.adata_hvg.obsp["connectivities"]
            G = nx.Graph(conn)

            # 計算 modularity
            modularity = community.modularity(label_dict, G)

            n_cls = labels.nunique()
            summary.append({"n_nei": k, "n_clusters": n_cls, "modularity": modularity})
            print(f"clusters: {n_cls}, modularity: {modularity:.4f}")
            
        df = pd.DataFrame(summary) 

        import matplotlib.pyplot as plt

        # 畫群數 & modularity
        fig, ax1 = plt.subplots(figsize=(8, 5))
        ax1.plot(df['n_nei'], df["n_clusters"], "o-", label="Clusters", color='#41C0F2')
        ax1.set_xlabel("n neighbors")
        ax1.set_ylabel("Number of Clusters")
        ax2 = ax1.twinx()  # 共享x軸
        ax2.plot(df['n_nei'], df["modularity"], "s--", label="Modularity", color='#C491F5')
        ax2.set_ylabel("Modularity Score")
        plt.title(f"Leiden Clustering in different n_neighbors")
        fig.legend()
        plt.savefig(f"{self.plot_dir}/Leiden_vs_Modularity.png", dpi=300)
        plt.show()
        plt.close()
        print(f"Leiden_vs_Modularity plot saved done !")
        return self.adata_hvg
    # ---------------------------- 【尋找最佳 resolution】 ----------------------------
    def find_best_resolution(self, best_nn=20,resolutions = [0.6, 0.8, 1.0, 1.2, 1.4,1.6,1.8,2.0]):

        silhouette_scores = []
        summary_res = []
        #使用PC50
        if self.use_rep == 'X_pca':
            sc.pp.neighbors(self.adata_hvg, n_neighbors=best_nn, n_pcs=self.ndim_pca, use_rep=self.use_rep)
        #使用Harmony 
        else:
            sc.pp.neighbors(self.adata_hvg, n_neighbors=best_nn, use_rep=self.use_rep)

        for res in resolutions:
           
            print(f"\n處理 res = {res}")
            sc.tl.leiden(self.adata_hvg, resolution=res, key_added=f'leiden_res{res}')
            embed = self.adata_hvg.obsm[self.use_rep]  # 'X_harmony' 或 'X_pca'
            score = silhouette_score(embed, self.adata_hvg.obs[f'leiden_res{res}'])
            silhouette_scores.append((res, score))

            n_cls = self.adata_hvg.obs[f'leiden_res{res}'].nunique()
            summary_res.append({"resolution": res, "n_clusters": n_cls, "silhouette score": score})
            print(f"clusters: {n_cls}, silhouette score: {score:.4f}")

        best_res = sorted(silhouette_scores, key=lambda x: x[1], reverse=True)[0][0]
        print("Best resolution by silhouette score:", best_res)

        df = pd.DataFrame(summary_res)

        fig, ax1 = plt.subplots(figsize=(8, 5))
        ax1.plot(df['resolution'], df["n_clusters"], "o-", label="Clusters", color='#41C0F2')
        ax1.set_xlabel("resolution")
        ax1.set_ylabel("Number of Clusters")
        ax2 = ax1.twinx()  # 共享x軸
        ax2.plot(df['resolution'], df["silhouette score"], "s--", label="Silhouette Score", color="#EE4F33")
        ax2.set_ylabel("Silhouette Score")
        plt.title(f"Leiden Clustering in different resolutions(n_neighbors={best_nn})")
        fig.legend()
        plt.savefig(f"{self.plot_dir}/Leiden_vs_SilhoutteScore.png", dpi=300)
        plt.show()
        plt.close()
        print(f"Leiden_vs_SilhoutteScore plot saved done !")
        return self.adata_hvg
    # ---------------------------- 使用最佳參數組合畫UMAP ----------------------------
    def plot_leiden_umap(self, best_nn=20, best_res=0.8,best_cluster_number=24):
        
        sub_dir = os.path.join(self.plot_dir, f"leiden_nn{best_nn}_res{best_res}")
        os.makedirs(sub_dir, exist_ok=True)

        # 計算細胞之間的鄰近關係
        if self.use_rep == 'X_pca':
            sc.pp.neighbors(self.adata_hvg, n_neighbors=best_nn, n_pcs=self.ndim_pca, use_rep=self.use_rep) # method='umap'
        #使用Harmony 
        else:
            sc.pp.neighbors(self.adata_hvg, n_neighbors=best_nn, use_rep=self.use_rep) # method='umap'

        # 將圖中的節點（細胞）聚成群
        sc.tl.leiden(self.adata_hvg, key_added=f"leiden_nn{best_nn}_res{best_res}", resolution=best_res)  # 結果會存入 adata.obs[key_added]
        self.adata_hvg.obs["cluster"] =self.adata_hvg.obs[f"leiden_nn{best_nn}_res{best_res}"].astype(int)

        # 繪製聚類結果
        sc.tl.umap(self.adata_hvg)


        #----{UMAP 標記NMOSD/Control}
        sc.pl.umap(
            self.adata_hvg,
            color="Condition",
            title=f"NMOSD/Control UMAP (cluster={best_cluster_number})",
            legend_loc="right margin",
            legend_fontsize=8,
            frameon=True,#圖加框
            legend_fontoutline=1, #圖例設置為圓點+標籤
            palette={
                "Control":"#67D6F0", 
                "NMOSD": "#E3A19F"},
            save=f"_leiden_nn{best_nn}_res{best_res}_group.png"    
        )
        sc_path = f"figures/umap_leiden_nn{best_nn}_res{best_res}_group.png"
        dst_path = f"{sub_dir}/leiden_nn{best_nn}_res{best_res}_UMAP_group.png"
        self.save_with_dpi(sc_path, dst_path, dpi=500)

        #----{UMAP圖 標記cluster}
        colors_30 = sns.color_palette("tab20") + sns.color_palette("Paired")[:10]
        self.adata_hvg.obs["cluster"] = self.adata_hvg.obs["cluster"].astype(str) #必須轉成字串在圖例時才不會呈現連續的bar
        # 取得所有 cluster label（包含 NMOSD+Control）
        unique_clusters = sorted(self.adata_hvg.obs["cluster"].unique()) 

        # 指定固定顏色對應
        colors_30 = sns.color_palette("tab20") + sns.color_palette("Paired")[:10]
        cluster_color_map = dict(zip(unique_clusters, colors_30[:len(unique_clusters)]))
        sc.pl.umap(
            self.adata_hvg,
            color="cluster",
            title=f"NMOSD/Control Leiden (cluster={best_cluster_number})",
            legend_loc="right margin",
            legend_fontsize=6,
            frameon=True,#圖加框
            legend_fontoutline=1, #圖例設置為圓點+標籤
            save=f"_leiden_nn{best_nn}_res{best_res}_cluster.png",
            palette=cluster_color_map 
        )
        sc_path = f"figures/umap_leiden_nn{best_nn}_res{best_res}_cluster.png"
        dst_path = f"{sub_dir}/leiden_nn{best_nn}_res{best_res}_UMAP_cluster.png"
        self.save_with_dpi(sc_path, dst_path, dpi=500)

        #---------------------------【畫 NMOSD/Control各自UMAP】------------------------------
        #----{NMOSD}
        sc.pl.umap(
            self.adata_hvg[self.adata_hvg.obs["Condition"] == "NMOSD"],
            color="cluster",
            title=f"NMOSD Leiden (cluster={best_cluster_number})",
            legend_loc="right margin",
            legend_fontsize=6,
            frameon=True,#圖加框
            legend_fontoutline=1, #圖例設置為圓點+標籤
            save=f"_leiden_nn{best_nn}_res{best_res}_NMOSD.png",
            palette=cluster_color_map 
        )
        sc_path = f"figures/umap_leiden_nn{best_nn}_res{best_res}_NMOSD.png"
        dst_path = f"{sub_dir}/leiden_nn{best_nn}_res{best_res}_UMAP_NMOSD.png"
        self.save_with_dpi(sc_path, dst_path, dpi=500)

        #----{Control}
        sc.pl.umap(
            self.adata_hvg[self.adata_hvg.obs["Condition"] == "Control"],
            color="cluster",
            title=f"Control Leiden (cluster={best_cluster_number})",
            legend_loc="right margin",
            legend_fontsize=6,
            frameon=True,
            legend_fontoutline=1,
            save=f"_leiden_nn{best_nn}_res{best_res}_Control.png",
            palette=cluster_color_map
        )
        sc_path = f"figures/umap_leiden_nn{best_nn}_res{best_res}_Control.png"
        dst_path = f"{sub_dir}/leiden_nn{best_nn}_res{best_res}_UMAP_Control.png"
        self.save_with_dpi(sc_path, dst_path, dpi=500)
        print(f"Leiden_nn{best_nn}_res{best_res} UMAP all plots saved done !")
        return self.adata_hvg
    
    #---------------------------【 畫 NMOSD/Control 在各 cluster 的佔比】------------------------------
    def plot_cluster_composition(self, best_nn=20,best_res=0.8,best_cluster_number=24):

        sub_dir = os.path.join(self.plot_dir, f"leiden_nn{best_nn}_res{best_res}")
        os.makedirs(sub_dir, exist_ok=True)

        # 建立組合的 DataFrame
        #comb_df = self.adata_hvg.obs[["cluster", "Condition"]].rename(columns={"Condition": "group"})
        #ct = pd.crosstab(comb_df["cluster"].astype(int), comb_df["group"], normalize="index")
        #ct = ct.sort_index() #排序

        df = self.adata_hvg.obs[["orig.ident","cluster","Condition"]].copy()
        df = df.rename(columns={"Condition":"group"})
        df["cluster"] = df["cluster"].astype(int)

        # 每個 sample×cluster 的細胞數（缺的補 0）
        mat = (
            df.pivot_table(index=["orig.ident","group"],
                        columns="cluster",
                        aggfunc="size",
                        fill_value=0)
        )

        # 在同一 group 內：對 sample 取平均（每個 sample 等權重）
        mean_counts = mat.groupby("group").mean().T  # (cluster × group)

        # 每個 cluster 內兩組比例 (加總為1)
        ct = mean_counts.div(mean_counts.sum(axis=1), axis=0).fillna(0)
        ct = ct.sort_index() #排序

        #pre 藍色，post 粉色
        colors = {"Control": "#67D6F0", "NMOSD": "#E3A19F"}

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
        plt.title(f"Proportion of NMOSD/Control cells in each cluster (cluster={best_cluster_number})", fontsize=14)
        plt.xticks(rotation=0)
        plt.tight_layout()
        plt.savefig(f"{sub_dir}/leiden_nn{best_nn}_res{best_res}_Composition.png", dpi=300)
        plt.show()
        plt.close()
        print(f"Leiden_nn{best_nn}_res{best_res}_Composition plot saved done !")
        return self.adata_hvg