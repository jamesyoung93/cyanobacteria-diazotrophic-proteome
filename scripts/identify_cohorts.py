#!/usr/bin/env python3
"""Identify diazotroph protein cohorts from filtered BLAST hits.

The output tables now report, for every query sequence in each cohort,
how many hits were observed to non‑diazotrophic species as well as the
best hit identity and accession for every species configured in
``species.yaml``.  Tables are written both as TSV and CSV files.
"""
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

# Ordered list of all species for table columns
species_order = [e["name"] for grp in cfg.values() for e in grp]

cohort_uni = []
cohort_fil = []
cohort_all = []

for qseqid, grp in hits.groupby("qseqid"):
    q_sp = grp["q_species"].iloc[0]
    target_species = set(grp["s_species"]) - {q_sp}

    row = {"qseqid": qseqid}
    row["non_diazotroph_hits"] = int(grp[grp["s_species"].isin(non_diazo)].shape[0])

    for sp in species_order:
        sub = grp[grp["s_species"] == sp]
        if sub.empty:
            row[f"{sp}_pident"] = 0
            row[f"{sp}_acc"] = ""
        else:
            best = sub.loc[sub["pident"].idxmax()]
            row[f"{sp}_pident"] = best["pident"]
            row[f"{sp}_acc"] = best["sseqid"]

    if q_sp in all_diazo:
        if (all_diazo - {q_sp}).issubset(target_species):
            cohort_all.append(row)

    if q_sp in unicell_diazo:
        required = set(unicell_diazo) - {q_sp}
        if required.issubset(target_species) and target_species.isdisjoint(set(fil_diazo)):
            cohort_uni.append(row)

    if q_sp in fil_diazo:
        required = set(fil_diazo) - {q_sp}
        if required.issubset(target_species) and target_species.isdisjoint(set(unicell_diazo)):
            cohort_fil.append(row)

def write_tables(rows, basename):
    if not rows:
        df = pd.DataFrame(columns=["qseqid", "non_diazotroph_hits"] +
                               [f"{sp}_pident" for sp in species_order] +
                               [f"{sp}_acc" for sp in species_order])
    else:
        df = pd.DataFrame(rows)
        df = df.sort_values("qseqid")

    tsv_path = os.path.join(args.output_dir, basename + ".tsv")
    csv_path = os.path.join(args.output_dir, basename + ".csv")
    df.to_csv(tsv_path, sep="\t", index=False)
    df.to_csv(csv_path, index=False)

write_tables(cohort_uni, "unicellular_specific")
write_tables(cohort_fil, "filamentous_specific")
write_tables(cohort_all, "diazotroph_common")
