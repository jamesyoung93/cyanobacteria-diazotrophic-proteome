#!/usr/bin/env python3
"""Plot cluster size distribution.

Reads results/clusters/clusters_mapping.tsv and outputs a table of cluster sizes
along with a histogram figure. Run this after the clustering step has
completed to visualise cluster sizes.
"""
import pandas as pd
import matplotlib.pyplot as plt
import os

mapping_file = "results/clusters/clusters_mapping.tsv"
output_dir = "results/figures"

os.makedirs(output_dir, exist_ok=True)

# Load mapping
mapping = pd.read_csv(mapping_file, sep="\t", header=None, names=["clusterID", "member"])

# Compute sizes
sizes = mapping.groupby("clusterID").size()

sizes.to_csv(os.path.join(output_dir, "cluster_sizes.tsv"), sep="\t", header=False)

# Plot histogram
plt.figure(figsize=(6,4))
sizes.hist(bins=50)
plt.xlabel("Cluster size (# proteins)")
plt.ylabel("Count")
plt.title("Cluster size distribution")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "cluster_size_histogram.png"))
