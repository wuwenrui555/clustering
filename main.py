import pandas as pd
import numpy as np
from pathlib import Path
from subcluster import SubClusterSystem


# Create a sample dataset
def create_sample_dataset(n_units=100, n_features=20):
    """Create a sample dataset for demonstration."""
    np.random.seed(42)

    # Generate random data with some structure (3 clusters)
    cluster1 = np.random.normal(0, 1, size=(n_units // 3, n_features))
    cluster2 = np.random.normal(5, 1, size=(n_units // 3, n_features))
    cluster3 = np.random.normal(-5, 1, size=(n_units // 3, n_features))

    # Combine clusters
    data = np.vstack([cluster1, cluster2, cluster3])

    # Add some noise and variation
    data += np.random.normal(0, 0.5, size=data.shape)

    # Create unit IDs and feature names
    unit_ids = [f"unit_{i}" for i in range(data.shape[0])]
    feature_names = [f"feature_{i}" for i in range(data.shape[1])]

    # Create DataFrame
    df = pd.DataFrame(data, index=unit_ids, columns=feature_names)

    return df


def main():
    # Create output directory
    output_dir = Path("output")
    output_dir.mkdir(exist_ok=True)

    # Create sample dataset
    print("Creating sample dataset...")
    raw_data = create_sample_dataset()
    print(f"Dataset shape: {raw_data.shape}")

    # Initialize subclustering system
    print("Initializing subclustering system...")
    system = SubClusterSystem(raw_data, output_dir)

    # First clustering: cluster all units using all features
    print("\nPerforming initial clustering...")
    unit_ids = raw_data.index.tolist()
    features = raw_data.columns.tolist()

    # Use KMeans for simplicity (scanpy/phenograph requires more dependencies)
    system.perform_clustering(
        unit_ids=unit_ids,
        features=features,
        method="kmeans",
        method_params={"n_clusters": 3, "random_state": 42},
        annotations={"0": "Type A", "1": "Type B", "2": "Mixed"},
        tags={"0": "clear", "1": "clear", "2": "unclear"},
    )

    # Update metadata
    system.manager.export_metadata()

    # Print summary
    print("\nClustering summary:")
    summary = system.get_summary()
    print(summary.head())

    # Select units with "Mixed" annotation for subclustering
    mixed_units = system.select_units_by_annotation("Mixed")
    print(
        f"\nSelected {len(mixed_units)} units with 'Mixed' annotation for subclustering"
    )

    # Perform subclustering on selected units
    if mixed_units:
        print("\nPerforming subclustering on 'Mixed' units...")
        system.perform_clustering(
            unit_ids=mixed_units,
            features=features[:10],  # Use only first 10 features for subclustering
            method="kmeans",
            method_params={"n_clusters": 2, "random_state": 42},
            annotations={"0": "Subtype X", "1": "Subtype Y"},
            tags={"0": "clear", "1": "clear"},
        )

        # Update metadata
        system.manager.export_metadata()

        # Print updated summary
        print("\nUpdated clustering summary after subclustering:")
        summary = system.get_summary()
        print(summary.head())

        # Show latest cluster IDs
        print("\nLatest cluster IDs sample:")
        latest_ids = summary["latest_cluster_id"].unique()
        print(f"Unique latest cluster IDs: {latest_ids}")

    print("\nSubclustering workflow completed successfully!")
    print(f"Results are saved in: {output_dir.absolute()}")


if __name__ == "__main__":
    main()
