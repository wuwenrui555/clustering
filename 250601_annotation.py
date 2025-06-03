# %%
from pathlib import Path

import anndata as ad
import pandas as pd

from subcluster import (
    ClusteringResult,
    ClusteringResultManager,
    plot_clustering_heatmap,
)

TQDM_FORMAT = "{desc}: {percentage:3.0f}%|{bar:30}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]"


# %%
adata = ad.read_h5ad(
    "/mnt/nfs/home/wenruiwu/projects/shuli_clustering/rcc/input/rcc_processed_data.h5ad"
)
id_exclude = ["RCC-TMA609(reg_4x5)-dst=reg019-src=reg005"]
adata = adata[~adata.obs.id.isin(id_exclude), :]

output_dir = (
    "/mnt/nfs/home/wenruiwu/projects/shuli_clustering/rcc/output/250529_rcc_new_mono"
)
manager = ClusteringResultManager(output_dir=output_dir, unit_ids=adata.obs.index)

# %%

output_dir = Path(
    "/mnt/nfs/home/wenruiwu/projects/bidmc-jiang-rcc/output/data/20250601_annotation_2"
)
anno_mono = manager.summary_df.loc[manager.summary_df.annotation != ""].copy()
anno_mono["unit_ids"] = anno_mono.index
anno_mono[["unit_ids", "annotation"]].to_csv(
    output_dir / "annotation_mono.csv",
    index=False,
)

anno_old = pd.read_csv(
    "/mnt/nfs/home/wenruiwu/projects/bidmc-jiang-rcc/output/data/20250527_annotation/annotation.csv"
)
anno_old = anno_old.loc[~anno_old.unit_ids.isin(anno_mono.unit_ids)]

anno_all = pd.concat(
    [anno_mono[["unit_ids", "annotation"]], anno_old],
    ignore_index=True,
)
anno_all.annotation.value_counts().count

anno_mapping = {
    "Endo_Vas": "Endothelial",
    "Endo_Lym": "Lymphatics",
    "Epi": "Epithelial",
    "Fibro": "Fibroblast",
    "nothing": "Epithelial",
    "dnT": "CD4T",
    "M1": "M1-like",
    "M2": "M2-like",
}

annotation_fix = anno_all.annotation.map(anno_mapping).fillna(anno_all.annotation)
anno_all["annotation"] = annotation_fix
anno_all.to_csv(
    output_dir / "annotation.csv",
    index=False,
)
adata.obs["annotation"] = adata.obs.index.map(
    anno_all.set_index("unit_ids").annotation.to_dict()
)
adata.write_h5ad(output_dir / "adata_with_annotation.h5ad")

# %%
adata = ad.read_h5ad(output_dir / "adata_with_annotation.h5ad")
adata.obs.annotation.value_counts().count

anno_sequence = [
    "CD8T",
    "CD4T",
    "Treg",
    "B",
    "M1-like",
    "M2-like",
    "Monocyte",
    "Neutrophil",
    "DC",
    "Epithelial",
    "Endothelial",
    "Lymphatics",
    "Fibroblast",
]
anno_dict = {}
for i, anno in enumerate(anno_sequence):
    anno_dict[anno] = f"{i + 1:02d}  {anno}"


adata.obs["annotation"] = adata.obs["annotation"].map(anno_dict)
clustering_result = ClusteringResult(
    clustering_id="summary",
    method="summary",
    unit_ids=adata.obs.index,
    cluster_ids=adata.obs["annotation"],
)
markers_all = [
    # "CD45",
    "CD3e",
    "CD8",
    "CD4",
    "FoxP3",
    "CD20",
    "CD68",
    "CD163",
    # "CD16",
    "CD11b",
    "MPO",
    "CD11c",
    "Cytokeratin",
    "CD31",
    "Podoplanin",
    "aSMA",
]

# %%
fig = plot_clustering_heatmap(
    adata,
    clustering_result,
    markers_all,
    plot_value="zscore",
    figsize=(8, 8),
    vmin=-2,
    center=0,
    vmax=2,
)

# %%
