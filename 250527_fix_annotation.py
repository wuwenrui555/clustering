# %%
import anndata as ad

from subcluster import (
    ClusteringResult,
    ClusteringResultManager,
    plot_clustering_heatmap_2,
)

# %%
adata = ad.read_h5ad(
    "/mnt/nfs/home/wenruiwu/projects/shuli_clustering/rcc/input/rcc_processed_data.h5ad"
)
output_dir = (
    "/mnt/nfs/home/wenruiwu/projects/shuli_clustering/rcc/output/clustering_rcc"
)
manager = ClusteringResultManager(output_dir=output_dir, unit_ids=adata.obs.index)


# %%
df_n = manager.summary_df.annotation.value_counts()
df_n.count
df_n.index

# %%
idx_non = manager.summary_df.annotation.isin(["Fibeo", "", "CD4T_CD8T"])
df_anno = manager.summary_df.loc[~idx_non, ["annotation"]]
for col in ["clustering_id", "method", "cluster_ids"]:
    df_anno[col] = "summary"
df_anno["unit_ids"] = df_anno.index
df_anno.to_csv(
    "/mnt/nfs/home/wenruiwu/projects/shuli_clustering/rcc/output/250527_rcc_new/clustering_results/summary.csv",
    index=False,
)

# %%
output_dir = (
    "/mnt/nfs/home/wenruiwu/projects/shuli_clustering/rcc/output/250527_rcc_new"
)
manager = ClusteringResultManager(output_dir=output_dir, unit_ids=adata.obs.index)

manager.summary_df.annotation.value_counts().count


manager.non_explicit_df
# %%
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
    "CD11c",
    "MPO",
    "Cytokeratin",
    "CD31",
    "Podoplanin",
    "aSMA",
]
result = ClusteringResult(
    clustering_id="whatever",
    method="whatever",
    unit_ids=manager.summary_df.index,
    cluster_ids=manager.summary_df.annotation,
)
plot_clustering_heatmap_2(
    adata,
    result,
    markers_all,
    figsize=(20, 8),
    kwargs_zscore={"vmin": -2, "center": 0, "vmax": 2},
    kwargs_mean={"vmin": 0, "vmax": 0.5},
)
