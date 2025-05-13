# %%
import importlib
from pathlib import Path

import anndata as ad
from IPython.core.interactiveshell import InteractiveShell

import subcluster

importlib.reload(subcluster)
from subcluster import (  # noqa: E402
    ClusteringResultManager,
)

InteractiveShell.ast_node_interactivity = "all"


# %%
adata = ad.read_h5ad(
    "/mnt/nfs/home/wenruiwu/projects/shuli_clustering/rcc/input/rcc_processed_data.h5ad"
)
id_exclude = ["RCC-TMA609(reg_4x5)-dst=reg019-src=reg005"]
adata = adata[~adata.obs.id.isin(id_exclude), :]

output_dir = Path("./output/clustering_rcc")
manager = ClusteringResultManager(output_dir=output_dir, unit_ids=adata.obs.index)

# %%
manager.non_explicit_df.tag_sub_2.value_counts()
sum(
    manager.non_explicit_df[manager.non_explicit_df.tag_sub_2 == "CD4T_CD8T"].index
    == "RCC-TMA542(reg4x6)-dst=reg006-src=reg011_c654"
)
len(manager.non_explicit_df[manager.non_explicit_df.tag_sub_2 == "CD4T_CD8T"].index)
