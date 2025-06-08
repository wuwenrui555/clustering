# %%
from pathlib import Path

import anndata as ad
from tqdm import tqdm

from subcluster import (
    ClusteringResult,
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
geojson_dir = Path(
    "/mnt/nfs/home/wenruiwu/projects/bidmc-jiang-rcc/output/data/20250324_update_geojson_R/geojson/"
)
geojson_files = list(geojson_dir.glob("*.geojson"))

for geojson_file in tqdm(
    geojson_files,
    desc="Updating GeoJSON files",
    bar_format=TQDM_FORMAT,
):
    update_geojson_from_clustering_result(
        geojson_file,
        clustering_result,
        "/mnt/nfs/home/wenruiwu/projects/shuli_clustering/rcc/output/250605_final",
    )
