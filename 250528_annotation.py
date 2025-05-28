# %%
from pathlib import Path

import anndata as ad
from pyqupath.geojson import GeojsonProcessor
from tqdm import tqdm

from subcluster import (
    ClusteringResultManager,
    plot_clustering_heatmap,
    ClusteringResult,
)

TQDM_FORMAT = "{desc}: {percentage:3.0f}%|{bar:30}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]"


# %%
adata = ad.read_h5ad(
    "/mnt/nfs/home/wenruiwu/projects/shuli_clustering/rcc/input/rcc_processed_data.h5ad"
)
id_exclude = ["RCC-TMA609(reg_4x5)-dst=reg019-src=reg005"]
adata = adata[~adata.obs.id.isin(id_exclude), :]

output_dir = (
    "/mnt/nfs/home/wenruiwu/projects/shuli_clustering/rcc/output/250527_rcc_new"
)
manager = ClusteringResultManager(output_dir=output_dir, unit_ids=adata.obs.index)

# %%
annotation_dir = Path(
    "/mnt/nfs/home/wenruiwu/projects/bidmc-jiang-rcc/output/data/20250527_annotation"
)
annotation_dir.mkdir(parents=True, exist_ok=True)

manager.summary_df.annotation.value_counts().count
adata.obs["annotation"] = adata.obs.index.map(manager.summary_df.annotation.to_dict())


adata.write_h5ad(annotation_dir / "adata_with_annotation.h5ad")


# %%

geojson_dir = Path(
    "/mnt/nfs/home/wenruiwu/projects/bidmc-jiang-rcc/output/data/20250324_update_geojson_R/geojson/"
)
output_dir = Path(
    "/mnt/nfs/home/wenruiwu/projects/bidmc-jiang-rcc/output/data/20250527_annotation/"
)
geojson_files = list(geojson_dir.glob("*.geojson"))
for geojson_file in tqdm(
    geojson_files,
    desc="Processing GeoJSON files",
    bar_format=TQDM_FORMAT,
):
    geojson_file = Path(geojson_file)
    geojson = GeojsonProcessor.from_path(geojson_file)

    core_id = geojson_file.name.split(r".", 1)[0]
    adata_core = adata[adata.obs.id == core_id, :].copy()
    adata_core.obs["cellLabel"] = adata_core.obs["cellLabel"].astype(str)
    name_dict = adata_core.obs.set_index("cellLabel")["annotation"].to_dict()

    # select the cells involved in the clustering
    geojson.gdf = geojson.gdf[geojson.gdf["name"].isin(name_dict.keys())]
    geojson.update_classification(name_dict)

    # geojson.plot_classification()

    output_dir = Path(output_dir)
    output_file = output_dir / "geojson" / geojson_file.name
    output_file.parent.mkdir(parents=True, exist_ok=True)
    geojson.output_geojson(output_file)

# %%
adata = ad.read_h5ad(
    "/mnt/nfs/home/wenruiwu/projects/bidmc-jiang-rcc/output/data/20250527_annotation/adata_with_annotation.h5ad"
)
adata

# %%
anno_dict = {
    "CD8T": "01  CD8T",
    "CD4T": "02  CD4T",
    "Treg": "03  Treg",
    "dnT": "04  dnT",
    "B": "05  B",
    "M1": "06  M1",
    "M2": "07  M2",
    "Neutrophil": "08  Neutrophil",
    "Monocyte": "09  Monocyte",
    "DC": "10  DC",
    "Epi": "11  Epithelial",
    "Endo_Vas": "12  Endo_Vas",
    "Endo_Lym": "13  Endo_Lym",
    "Fibro": "14  Fibroblast",
    "nothing": "15  nothing",
}
adata.obs["annotation"] = adata.obs["annotation"].map(anno_dict)

clustering_result = ClusteringResult(
    clustering_id="summary",
    method="summary",
    unit_ids=adata.obs.index,
    cluster_ids=adata.obs["annotation"],
)

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
