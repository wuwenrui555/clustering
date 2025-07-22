# Clustering

## Installation

```sh
conda create -y -n clustering_3.11 python=3.11
conda activate clustering_3.1

pip install git+https://github.com/wuwenrui555/clustering.git@v0.1.0
```

## `ClusteringResultManager`

```python
# load and organize all saved clustering result in "clustering_sequence.txt"
manager = ClusteringResultManager(
    output_dir=output_dir,     # The directory to store clustering results.
    unit_ids=adata.obs.index,  # All unit ids to be clustered.
)

# units with no annotation need to be subclustered
manager.non_explicit_df
```

## `ClusteringResult`

### 01. Initiation

```python
# option 1: clustering
clustering_result = run_clustering(...)

# option 2: manually create
clustering_result = ClusteringResult(...)

# option 3: reload from stashed result
clustering_result = ClusteringResult.pop(...)

# option 4: reload from saved result
clustering_result = ClusteringResult.from_csv(...)
```

### 02. Reviewing

```python
# heatmap
plot_clustering_heatmap_2(...)

# update geojson
update_geojson_from_clustering_result(...)
```

### 03. Processing

```python
# add annotation for explicit clusters
clustering_result.add_annotation(...)

# add tag for non-explicit clusters (to make it easier to filter for sub-clustering)
clustering_result.add_tag(...)

# temperately stash and reload
clustering_result.stash(...)
clustering_result = ClusteringResult.pop(...)
```

### 04. Save

```python
clustering_result.save(...)
```
