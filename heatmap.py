import matplotlib.pyplot as plt
import pandas as pd
import PyComplexHeatmap as pch
import seaborn as sns


def plot_cluster_heatmap(
    data: pd.DataFrame,
    metadata: pd.DataFrame,
    markers: list[str],
    figsize=(35, 15),
):
    """
    Plot a heatmap of the cluster data

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing the cluster data.
    metadata : pd.DataFrame
        DataFrame containing the metadata.
    markers : list[str]
        List of markers to plot.
    figsize : tuple, optional
        Figure size, by default (35, 15).
    """
    cmap_tab20 = plt.get_cmap("tab20")
    colors_tab20 = [cmap_tab20(i) for i in range(20)]

    clustering_data = system.raw_data.loc[clustering_result.unit_id][features_to_plot]

    # calculate the mean of each cluster
    cluster_mean = clustering_data.groupby(clustering_result.cluster_id).mean()

    # calculate the z-score of each cluster
    population_mean = clustering_data.mean(axis=0)
    population_std = clustering_data.std(axis=0)
    cluster_zscore = (cluster_mean - population_mean) / population_std

    # plot the heatmap
    heatmap_df = cluster_zscore.T
    cluster_n = pd.DataFrame(pd.Series(clustering_result.cluster_id).value_counts())

    plt.figure(figsize=figsize)
    cell_count_df = pd.DataFrame(data_cluster["cluster"].value_counts())
    cell_size_df = data.groupby("cluster")["cellSize"].mean().to_frame()
    tma_prop_df = (
        data_cluster.groupby("cluster")["tma"]
        .value_counts(normalize=True)
        .unstack()
        .fillna(0)
    )
    core_prop_df = (
        data_cluster.groupby("cluster")["id"]
        .value_counts(normalize=True)
        .unstack()
        .fillna(0)
    )
    tissue_type_df = (
        data.merge(metadata[["id", "tissue_type"]], on="id")
        .groupby("cluster")["tissue_type"]
        .value_counts(normalize=True)
        .unstack()
        .fillna(0)
    )

    col_ha = pch.HeatmapAnnotation(
        tma_prop=pch.anno_barplot(
            tma_prop_df,
            legend=True,
            colors=[colors_tab20[i % 20] for i in range(len(tma_prop_df.columns))],
            height=30,
            ylim=(0, 1),
        ),
        core_prop=pch.anno_barplot(
            core_prop_df,
            legend=False,
            colors=[colors_tab20[i % 20] for i in range(len(core_prop_df.columns))],
            height=30,
            ylim=(0, 1),
            linewidth=0,
        ),
        tissue_type=pch.anno_barplot(
            tissue_type_df,
            legend=True,
            colors=[colors_tab20[i % 20] for i in range(len(tissue_type_df.columns))],
            height=30,
            ylim=(0, 1),
        ),
        cell_count=pch.anno_barplot(
            cell_count_df, legend=False, colors="grey", height=20
        ),
        cell_size=pch.anno_barplot(
            cell_size_df, legend=False, colors="grey", height=20
        ),
    )
    pch.ClusterMapPlotter(
        data=df_heatmap,
        top_annotation=col_ha,
        col_cluster=True,
        row_cluster=False,
        label="Z-score",
        col_dendrogram=True,
        col_dendrogram_size=30,
        row_dendrogram=True,
        row_dendrogram_size=20,
        show_rownames=True,
        show_colnames=True,
        row_names_side="right",
        tree_kws={"row_cmap": "Set1", "colors": "blue"},
        verbose=0,
        legend_gap=5,
        cmap="RdYlBu_r",
        vmin=-2,
        vmax=2,
        center=0,
        xticklabels_kws={"labelrotation": -90, "labelcolor": "blue"},
    )


# %%

features_to_plot = system.raw_data.columns.tolist()


clustering_data = system.raw_data.loc[clustering_result.unit_id][features_to_plot]

# calculate the mean of each cluster
cluster_mean = clustering_data.groupby(clustering_result.cluster_id).mean()

# calculate the z-score of each cluster
population_mean = clustering_data.mean(axis=0)
population_std = clustering_data.std(axis=0)
cluster_zscore = (cluster_mean - population_mean) / population_std


cluster_n = pd.Series(clustering_result.cluster_id).value_counts()


# plot the heatmap
heatmap_df = cluster_zscore.T
sns.heatmap(heatmap_df, cmap="RdBu_r", vmin=-2, vmax=2)

bar_x = cluster_n.index
bar_y = heatmap_df_cluster_n[pheatmap_df.columns].values
sns.barplot(x=bar_x, y=bar_y)

sns.heatmap(heatmap_df_cluster_mean, cmap="Reds", vmin=0, vmax=1)
