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

# Extract species from member IDs. Sequence headers were prefixed with the
# full species name followed by an underscore when building the BLAST
# database. Species names themselves contain underscores, so simply splitting
# on '_' would break. Instead, match against the known species list.
species_prefixes = sorted(sp_to_group.keys(), key=len, reverse=True)

def extract_species(seqid: str) -> str:
    for sp in species_prefixes:
        if seqid.startswith(sp + "_"):
            return sp
    return seqid.split("_")[0]

mapping["species"] = mapping["member"].apply(extract_species)

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
