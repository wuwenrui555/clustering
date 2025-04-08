# %%
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import PyComplexHeatmap as pch

from subcluster import SubClusterSystem, ClusteringResult


# Create a sample dataset
def create_sample_dataset():
    data_f = "/mnt/nfs/home/wenruiwu/projects/bidmc-jiang-rcc/output/data/20250211_annotation/20250224_combined_dataScaleSize_rm=extreme_fix=CD68_normDAPI_arcsinh=0.01_quantile=1-99.csv"
    id = "RCC-TMA609(reg_4x5)-dst=reg019-src=reg005"
    data = pd.read_csv(data_f)

    # data for demo
    data_demo = data[data["id"].isin([id])].copy()
    data_demo["cell_id"] = (
        data_demo["id"].astype(str) + "_c" + data_demo["cellLabel"].astype(str)
    )
    data_demo = data_demo.set_index("cell_id")

    obs_columns = ["id", "cellLabel", "cellSize", "Y_cent", "X_cent", "tma"]
    data_demo_obs = data_demo[obs_columns]
    data_demo_X = data_demo.drop(columns=obs_columns)

    adata = ad.AnnData(data_demo_X)
    adata.obs = data_demo_obs
    adata.var_names = data_demo_X.columns.tolist()

    adata.write_h5ad("input/data_demo.h5ad", compression="gzip")


def simulate_sequential_clustering():
    """
    模拟多次连续聚类并生成summary，所有聚类结果都标记为clear
    """
    markers_all = [
        "CD45",
        "CD3e",
        "CD8",
        "CD4",
        "FoxP3",
        "CD20",
        "CD68",
        "CD163",
        "CD16",
        "CD11b",
        "MPO",
        "Cytokeratin",
        "CD31",
        "Podoplanin",
        "aSMA",
    ]

    # 创建输出目录
    output_dir = Path("output/sequential_clustering")
    output_dir.mkdir(parents=True, exist_ok=True)

    # 加载数据集
    print("Loading dataset...")
    adata = ad.read_h5ad("input/data_demo.h5ad")[0:1000]
    print(f"Dataset: {adata}")

    # 初始化子聚类系统
    print("Initializing subclustering system...")
    system = SubClusterSystem(adata, output_dir)

    # 第1轮聚类：使用所有细胞，K=3
    print("\n===== Round 1: Initial clustering with K=3 =====")
    unit_ids = adata.obs_names.tolist()

    clustering_1 = system.perform_clustering(
        unit_ids=unit_ids,
        features=markers_all,
        method="kmeans",
        method_params={"n_clusters": 15, "random_state": 42},
    )

    # 添加注释和标签（全部设为clear）
    cluster_names_1 = {
        "0": "Immune Cells",
        "1": "Epithelial Cells",
        "2": "Stromal Cells",
    }

    clustering_1 = system.add_annotations_to_clustering(clustering_1, cluster_names_1)

    # 所有标签设为clear
    cluster_tags_1 = {str(i): "clear" for i in range(3)}
    clustering_1 = system.add_tags_to_clustering(clustering_1, cluster_tags_1)

    # 保存聚类结果
    system.save_clustering_result(clustering_1)
    print("Saved round 1 clustering results")
    print(f"Clusters: {cluster_names_1}")

    # 创建每个聚类群体的信息
    cluster_counts = {}
    for cluster_id in clustering_1.cluster_df["cluster_id"].unique():
        count = sum(clustering_1.cluster_df["cluster_id"] == cluster_id)
        cluster_counts[cluster_id] = count

    print(f"Cluster counts: {cluster_counts}")

    # 第2轮聚类：只聚类免疫细胞，K=4
    print("\n===== Round 2: Subclustering Immune Cells with K=4 =====")
    # 选择免疫细胞亚群
    immune_cells = clustering_1.cluster_df[
        clustering_1.cluster_df["cluster_id"] == "0"
    ]["unit_id"].tolist()

    clustering_2 = system.perform_clustering(
        unit_ids=immune_cells,
        features=markers_all[:8],  # 使用免疫细胞标记物
        method="kmeans",
        method_params={"n_clusters": 4, "random_state": 42},
    )

    # 添加免疫细胞亚群注释和标签
    cluster_names_2 = {
        "0": "T Cells",
        "1": "B Cells",
        "2": "Macrophages",
        "3": "Other Immune",
    }

    clustering_2 = system.add_annotations_to_clustering(clustering_2, cluster_names_2)

    # 所有标签设为clear
    cluster_tags_2 = {str(i): "clear" for i in range(4)}
    clustering_2 = system.add_tags_to_clustering(clustering_2, cluster_tags_2)

    # 保存聚类结果
    system.save_clustering_result(clustering_2)
    print("Saved round 2 clustering results")
    print(f"Clusters: {cluster_names_2}")

    # 第3轮聚类：只聚类T细胞，K=3
    print("\n===== Round 3: Subclustering T Cells with K=3 =====")
    # 选择T细胞亚群
    t_cells = clustering_2.cluster_df[clustering_2.cluster_df["cluster_id"] == "0"][
        "unit_id"
    ].tolist()

    clustering_3 = system.perform_clustering(
        unit_ids=t_cells,
        features=["CD3e", "CD8", "CD4", "FoxP3"],  # T细胞相关标记物
        method="kmeans",
        method_params={"n_clusters": 3, "random_state": 42},
    )

    # 添加T细胞亚群注释和标签
    cluster_names_3 = {"0": "CD8+ T Cells", "1": "CD4+ T Cells", "2": "Tregs"}

    clustering_3 = system.add_annotations_to_clustering(clustering_3, cluster_names_3)

    # 所有标签设为clear
    cluster_tags_3 = {str(i): "clear" for i in range(3)}
    clustering_3 = system.add_tags_to_clustering(clustering_3, cluster_tags_3)

    # 保存聚类结果
    system.save_clustering_result(clustering_3)
    print("Saved round 3 clustering results")
    print(f"Clusters: {cluster_names_3}")

    # 第4轮聚类：只聚类上皮细胞，K=2
    print("\n===== Round 4: Subclustering Epithelial Cells with K=2 =====")
    # 选择上皮细胞亚群
    epithelial_cells = clustering_1.cluster_df[
        clustering_1.cluster_df["cluster_id"] == "1"
    ]["unit_id"].tolist()

    clustering_4 = system.perform_clustering(
        unit_ids=epithelial_cells,
        features=["Cytokeratin", "CD31"],  # 上皮细胞相关标记物
        method="kmeans",
        method_params={"n_clusters": 2, "random_state": 42},
    )

    # 添加上皮细胞亚群注释和标签
    cluster_names_4 = {"0": "Tumor Cells", "1": "Normal Epithelium"}

    clustering_4 = system.add_annotations_to_clustering(clustering_4, cluster_names_4)

    # 所有标签设为clear
    cluster_tags_4 = {str(i): "clear" for i in range(2)}
    clustering_4 = system.add_tags_to_clustering(clustering_4, cluster_tags_4)

    # 保存聚类结果
    system.save_clustering_result(clustering_4)
    print("Saved round 4 clustering results")
    print(f"Clusters: {cluster_names_4}")

    # 导出summary数据
    summary = system.get_summary()
    print("\n===== Final clustering summary =====")
    print(f"Summary shape: {summary.shape}")
    print(summary.head())

    # 整理列的顺序：unit_id, 聚类结果（按clustering_sequence），最后是annotation和tag
    desired_columns = ["unit_id", "latest_cluster_id", "annotation", "tag"]
    cluster_columns = [col for col in summary.columns if col not in desired_columns]
    ordered_columns = (
        ["unit_id"] + cluster_columns + ["latest_cluster_id", "annotation", "tag"]
    )

    # 按正确的顺序排序
    summary_ordered = summary[ordered_columns].sort_values("unit_id")

    # 导出最终结果到CSV
    summary_file = output_dir / "summary.csv"
    summary_ordered.to_csv(summary_file, index=False)
    print(f"Exported summary to {summary_file}")

    # 创建一个热图显示不同聚类层次的细胞数量分布
    print("\nCreating cluster distribution heatmap...")

    # 生成最终级别的聚类映射表
    summary_stats = pd.crosstab(summary["annotation"], columns="count").reset_index()
    summary_stats.columns = ["Cluster", "Cell Count"]
    summary_stats = summary_stats.sort_values("Cell Count", ascending=False)

    # 创建条形图
    plt.figure(figsize=(12, 6))
    sns.barplot(x="Cluster", y="Cell Count", data=summary_stats)
    plt.title("Cell Distribution Across Final Clusters")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_dir / "cluster_distribution.png")

    print(
        f"Saved cluster distribution plot to {output_dir / 'cluster_distribution.png'}"
    )

    return summary_ordered


# %%
if __name__ == "__main__":
    # 运行顺序聚类模拟
    simulate_sequential_clustering()


# %%
system.adata

unit_ids = clustering_1.unit_id
cluster_ids = clustering_1.cluster_id


# %%

features = markers_all
clustering_result = clustering_1


def plot_cluster_heatmap(
    adata: ad.AnnData,
    clustering_result: ClusteringResult,
    features: list[str],
    figsize=(10, 8),
    plot_value: str = "zscore",
    vmin: float = -2,
    vmax: float = 2,
):
    """
    Plot a heatmap of the cluster data

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object containing the cluster data.
    clustering_result : ClusteringResult
        ClusteringResult object containing the clustering result.
    features : list[str]
        List of features to plot.
    figsize : tuple, optional
        Figure size, by default (10, 8).
    plot_value : str, optional
        Value to plot, either "zscore" or "mean", by default "zscore"
    vmin : float, optional
        Minimum value for color scaling, by default -2
    vmax : float, optional
        Maximum value for color scaling, by default 2
    """
    metadata = adata.obs
    unit_ids = adata.obs_names
    # Fix: Convert cluster_ids to a pandas Series if it's not already
    cluster_ids = clustering_result.cluster_id
    clustering_data = adata.to_df().loc[unit_ids][features]

    # value to plot heatmap
    if plot_value == "zscore":
        cluster_mean = clustering_data.groupby(cluster_ids).mean()
        population_mean = clustering_data.mean(axis=0)
        population_std = clustering_data.std(axis=0)
        cluster_zscore = (cluster_mean - population_mean) / population_std
        heatmap_df = cluster_zscore.T
    elif plot_value == "mean":
        cluster_mean = clustering_data.groupby(cluster_ids).mean()
        heatmap_df = cluster_mean.T
        vmin = 0

    # value to plot subplot
    cluster_count = pd.Series(cluster_ids).value_counts().to_frame(name="count")
    cluster_mean_cellsize = metadata.groupby(cluster_ids)["cellSize"].mean().to_frame()

    # %%
    plt.figure(figsize=figsize)
    col_ha = pch.HeatmapAnnotation(
        cell_count=pch.anno_barplot(
            cluster_count, legend=False, colors="grey", height=20
        ),
        cell_size=pch.anno_barplot(
            cluster_mean_cellsize, legend=False, colors="grey", height=20
        ),
        verbose=0,
    )

    heatmap = pch.ClusterMapPlotter(
        data=heatmap_df,
        top_annotation=col_ha,
        col_cluster=True,
        row_cluster=False,
        label=plot_value,
        col_dendrogram=True,
        col_dendrogram_size=15,
        show_rownames=True,
        show_colnames=True,
        row_names_side="right",
        tree_kws={"colors": "blue"},
        verbose=0,
        legend_gap=5,
        cmap="RdBu_r",
        vmin=vmin,
        vmax=vmax,
        center=0,
        xticklabels_kws={"labelrotation": -90, "labelcolor": "blue"},
    )
    fig, ax = plt.subplots(figsize=figsize)
    heatmap.plot(ax=ax)
    fig.savefig("output/cluster_heatmap.png")
    plt.show()
    return heatmap


# Example usage:
# plot_cluster_heatmap(system.adata, clustering_1, markers_all)


# %%
