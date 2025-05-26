#!/usr/bin/env python3
"""Identify diazotroph protein cohorts from filtered BLAST hits."""
import argparse
import pandas as pd
import yaml
import os

parser = argparse.ArgumentParser(description="Identify protein cohorts")
parser.add_argument("blast", help="Filtered BLAST tabular file")
parser.add_argument("output_dir", help="Directory for output cohort lists")
args = parser.parse_args()

if not os.path.isfile(args.blast):
    raise SystemExit(f"BLAST results not found: {args.blast}")

os.makedirs(args.output_dir, exist_ok=True)

# Load species grouping
with open("species.yaml") as fh:
    cfg = yaml.safe_load(fh)

sp_to_group = {}
for grp, entries in cfg.items():
    for e in entries:
        sp_to_group[e["name"]] = grp

unicell_diazo = [e["name"] for e in cfg.get("unicellular_diazotrophic", [])]
fil_diazo = [e["name"] for e in cfg.get("filamentous_diazotrophic", [])]
non_diazo = [e["name"] for g in cfg if "non_diazotrophic" in g for e in cfg[g]]
all_diazo = set(unicell_diazo + fil_diazo)

species_list = sorted(sp_to_group.keys(), key=len, reverse=True)

def extract_species(seqid: str) -> str:
    for sp in species_list:
        if seqid.startswith(sp + "_"):
            return sp
    return seqid.split("_")[0]

cols = ["qseqid","sseqid","pident","length","qlen","slen","evalue"]
hits = pd.read_csv(args.blast, sep="\t", header=None, names=cols)

hits["q_species"] = hits["qseqid"].apply(extract_species)
hits["s_species"] = hits["sseqid"].apply(extract_species)

cohort_uni = []
cohort_fil = []
cohort_all = []

for qseqid, grp in hits.groupby("qseqid"):
    q_sp = grp["q_species"].iloc[0]
    target_species = set(grp["s_species"]) - {q_sp}

    if q_sp in all_diazo:
        if (all_diazo - {q_sp}).issubset(target_species) \
                and target_species.isdisjoint(non_diazo):
            cohort_all.append(qseqid)

    if q_sp in unicell_diazo:
        required = set(unicell_diazo) - {q_sp}
        if required.issubset(target_species) \
                and target_species.isdisjoint(set(fil_diazo + non_diazo)):
            cohort_uni.append(qseqid)

    if q_sp in fil_diazo:
        required = set(fil_diazo) - {q_sp}
        if required.issubset(target_species) \
                and target_species.isdisjoint(set(unicell_diazo + non_diazo)):
            cohort_fil.append(qseqid)

pd.Series(sorted(cohort_uni)).to_csv(os.path.join(args.output_dir, "unicellular_specific.tsv"), index=False, header=False)
pd.Series(sorted(cohort_fil)).to_csv(os.path.join(args.output_dir, "filamentous_specific.tsv"), index=False, header=False)
pd.Series(sorted(cohort_all)).to_csv(os.path.join(args.output_dir, "diazotroph_common.tsv"), index=False, header=False)
