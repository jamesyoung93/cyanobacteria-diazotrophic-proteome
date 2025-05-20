#!/usr/bin/env python3
import pandas as pd
import yaml
import os

# Load species → group mapping
with open("species.yaml") as f:
    cfg = yaml.safe_load(f)

sp_to_group = {}
for grp, entries in cfg.items():
    for e in entries:
        sp_to_group[e["name"]] = grp

# Read cluster → member mapping
mapping = pd.read_csv(
    "results/clusters/clusters_mapping.tsv",
    sep="\t",
    header=None,
    names=["clusterID","member"],
    dtype={"clusterID": int}
)

# Extract species from member (assumes names like "anabaena_ppc7120_protXYZ")
mapping["species"] = mapping["member"].str.split("_").str[0]

# Build presence/absence matrix
species_list = sorted(sp_to_group.keys())
clusters = sorted(mapping["clusterID"].unique())

mat = pd.DataFrame(0, index=clusters, columns=species_list)
for _, row in mapping.iterrows():
    sp = row["species"]
    cid = row["clusterID"]
    if sp in mat.columns:
        mat.at[cid, sp] = 1

# Save
os.makedirs("results/matrices", exist_ok=True)
mat.index.name = "clusterID"
mat.to_csv("results/matrices/pres_abs.csv")

# Also save species grouping
grp_df = pd.DataFrame.from_dict(
    sp_to_group, orient="index", columns=["group"]
)
grp_df.index.name = "species"
grp_df.to_csv("results/matrices/species_groups.csv")
