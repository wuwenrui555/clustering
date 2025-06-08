# %%
from pathlib import Path

import anndata as ad

from subcluster import (
    ClusteringResult,
    plot_clustering_heatmap,
    update_geojson_from_clustering_result,
)

TQDM_FORMAT = "{desc}: {percentage:3.0f}%|{bar:30}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]"

# %%
output_dir = Path(
    "/mnt/nfs/home/wenruiwu/projects/bidmc-jiang-rcc/output/data/20250601_annotation_2"
)
adata = ad.read_h5ad(output_dir / "adata_with_annotation.h5ad")

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
# generate geojson to visualize the clustering result
core_id = "RCC-TMA544-dst=reg014-src=reg008"
geojson_file = Path(
    f"/mnt/nfs/home/wenruiwu/projects/bidmc-jiang-rcc/output/data/20250324_update_geojson_R/geojson/{core_id}.geojson"
)
update_geojson_from_clustering_result(
    geojson_file,
    clustering_result,
    "/mnt/nfs/home/wenruiwu/projects/shuli_clustering/rcc/output/250605_final",
)

# %%
id_exclude = [
    "RCC-TMA001(reg4x3)-dst=reg001-src=reg001",
    "RCC-TMA001(reg4x5)-dst=reg003-src=reg004",
    "RCC-TMA001(reg4x5)-dst=reg006-src=reg001",
    "RCC-TMA001(reg4x5)-dst=reg011-src=reg001",
    "RCC-TMA002(reg4x6)-dst=reg001-src=reg002",
    "RCC-TMA002(reg4x6)-dst=reg002-src=reg001",
    "RCC-TMA002(reg4x6)-dst=reg003-src=reg001",
    "RCC-TMA002(reg4x6)-dst=reg007-src=reg004",
    "RCC-TMA003(reg4x4)-dst=reg001-src=reg001",
    "RCC-TMA003(reg4x5)-dst=reg006-src=reg003",
    "RCC-TMA003(reg5x6)-dst=reg001-src=reg002",
    "RCC-TMA003(reg5x6)-dst=reg002-src=reg001",
    "RCC-TMA004(reg4x4)-dst=reg005-src=reg001",
    "RCC-TMA004(reg4x5)-dst=reg002-src=reg004",
    "RCC-TMA004(reg4x5)-dst=reg003-src=reg001",
    "RCC-TMA004(reg4x6)-dst=reg001-src=reg005",
    "RCC-TMA542(reg4x5)-dst=reg015-src=reg012",
    "RCC-TMA542(reg4x6)-dst=reg001-src=reg001",
    "RCC-TMA542(reg4x6)-dst=reg002-src=reg002",
    "RCC-TMA542(reg4x6)-dst=reg003-src=reg005",
    "RCC-TMA542(reg4x6)-dst=reg004-src=reg001",
    "RCC-TMA542(reg4x6)-dst=reg006-src=reg011",
    "RCC-TMA542(reg4x6)-dst=reg008-src=reg014",
    "RCC-TMA543-dst=reg001-src=reg021",
    "RCC-TMA543-dst=reg006-src=reg022",
    "RCC-TMA543-dst=reg013-src=reg013",
    "RCC-TMA543-dst=reg018-src=reg014",
    "RCC-TMA543-dst=reg020-src=reg004",
    "RCC-TMA543-dst=reg021-src=reg025",
    "RCC-TMA543-dst=reg025-src=reg005",
    "RCC-TMA544-dst=reg001-src=reg021",
    "RCC-TMA544-dst=reg005-src=reg001",
    "RCC-TMA544-dst=reg006-src=reg022",
    "RCC-TMA544-dst=reg008-src=reg012",
    "RCC-TMA544-dst=reg013-src=reg013",
    "RCC-TMA544-dst=reg020-src=reg004",
    "RCC-TMA544-dst=reg025-src=reg005",
    "RCC-TMA609(reg_4x5)-dst=reg001-src=reg001",
    "RCC-TMA609(reg_4x5)-dst=reg007-src=reg013",
    "RCC-TMA609(reg_4x5)-dst=reg008-src=reg008",
    "RCC-TMA609(reg_4x5)-dst=reg014-src=reg004",
    "RCC-TMA609(reg_4x5)-dst=reg015-src=reg025",
    "RCC-TMA609(reg_4x5)-dst=reg016-src=reg020",
    "RCC-TMA609(reg_4x5)-dst=reg019-src=reg005",
]
adata_rcc = adata[~adata.obs.id.isin(id_exclude), :]
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
