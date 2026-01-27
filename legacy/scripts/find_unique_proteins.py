#!/usr/bin/env python3
"""Identify proteins unique to diazotrophic species from BLAST hits."""
import pandas as pd
import yaml
import argparse
import os

parser = argparse.ArgumentParser(description="Find unique diazotroph proteins")
parser.add_argument("blast", help="Filtered BLAST tabular file")
parser.add_argument("output_dir", help="Directory for output tables")
args = parser.parse_args()

# Ensure the BLAST results exist before proceeding
if not os.path.isfile(args.blast):
    raise SystemExit(
        f"BLAST results not found: {args.blast}\n"
        "Run the `filter_blast_hits` step first (see README)."
    )

# Ensure the output directory exists
os.makedirs(args.output_dir, exist_ok=True)

# Load species groups
with open("species.yaml") as fh:
    cfg = yaml.safe_load(fh)

sp_to_group = {}
for grp, entries in cfg.items():
    for e in entries:
        sp_to_group[e["name"]] = grp

diazo_groups = [g for g in cfg if "diazotrophic" in g and "non_diazotrophic" not in g]
#diazo_groups = [g for g in cfg if "diazotrophic" in g]
non_diazo_groups = [g for g in cfg if "non_diazotrophic" in g]

diazo_species = [e["name"] for g in diazo_groups for e in cfg[g]]
non_diazo_species = [e["name"] for g in non_diazo_groups for e in cfg[g]]

cols = ["qseqid","sseqid","pident","length","qlen","slen","evalue"]
hits = pd.read_csv(args.blast, sep="\t", header=None, names=cols)

# Extract species names from sequence IDs. Headers were prefixed
# with the exact species name followed by an underscore when the BLAST
# database was built. Species names themselves contain underscores, so we
# can’t simply take the first field. Instead, match against the known
# species list.
species_list = sorted(sp_to_group.keys(), key=len, reverse=True)

# Build a mapping from sequence accessions (without the species prefix) back
# to their species by scanning the FASTA file used to create the BLAST
# database. This helps recover the species name even if the BLAST output
# lacks the explicit prefix for some reason.
id_to_species = {}
fasta_path = "data/db/all_cyano.faa"
if os.path.exists(fasta_path):
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]
                for sp in species_list:
                    pre = sp + "_"
                    if header.startswith(pre):
                        orig = header[len(pre):]
                        id_to_species[orig] = sp
                        id_to_species[header] = sp
                        # Also record just the accession if the header
                        # follows the typical "db|ACC|..." pattern so
                        # lookups work even when BLAST trims the prefix
                        if "|" in orig:
                            parts = orig.split("|")
                            if len(parts) > 1:
                                id_to_species[parts[1]] = sp
                        break

def extract_species(seqid: str) -> str:
    # First check explicit mapping from the FASTA database
    if seqid in id_to_species:
        return id_to_species[seqid]
    # BLAST may report IDs like "sp|ACC|..."; try the accession alone
    if "|" in seqid:
        acc = seqid.split("|")[1]
        if acc in id_to_species:
            return id_to_species[acc]
    for sp in species_list:
        if seqid.startswith(sp + "_"):
            return sp
    # If the sequence ID contains a known accession without the prefix,
    # try looking it up in the mapping as well
    acc = seqid.split("_")[0]
    if acc in id_to_species:
        return id_to_species[acc]
    # Fallback to the first token
    return seqid

hits["q_species"] = hits["qseqid"].apply(extract_species)
hits["s_species"] = hits["sseqid"].apply(extract_species)

# Hits from diazotroph queries to non-diazotroph subjects
hits_d2n = hits[hits["q_species"].isin(diazo_species) & hits["s_species"].isin(non_diazo_species)]
matched_to_non = set(hits_d2n["qseqid"])

all_diazo_proteins = hits["qseqid"][hits["q_species"].isin(diazo_species)].unique()
unique_diazo = sorted(set(all_diazo_proteins) - matched_to_non)

# Save list of proteins unique to diazotrophic species
pd.Series(unique_diazo).to_csv(f"{args.output_dir}/unique_diazotroph_proteins.tsv", index=False, header=False)

# Further classify by filamentous vs unicellular
filamentous_diazo = [e["name"] for e in cfg.get("filamentous_diazotrophic", [])]
unicellular_diazo = [e["name"] for e in cfg.get("unicellular_diazotrophic", [])]

hits_diazo = hits[hits["q_species"].isin(diazo_species) & hits["s_species"].isin(diazo_species)]

def unique_to_group(group_species):
    other_species = set(diazo_species) - set(group_species)
    grp_hits = hits_diazo[hits_diazo["q_species"].isin(group_species)]
    hits_to_other = grp_hits[grp_hits["s_species"].isin(other_species)]
    matched = set(hits_to_other["qseqid"]) | matched_to_non
    grp_all = grp_hits["qseqid"].unique()
    return sorted(set(grp_all) - matched)

unique_fil = unique_to_group(filamentous_diazo)
unique_uni = unique_to_group(unicellular_diazo)

pd.Series(unique_fil).to_csv(f"{args.output_dir}/unique_filamentous.tsv", index=False, header=False)
pd.Series(unique_uni).to_csv(f"{args.output_dir}/unique_unicellular.tsv", index=False, header=False)
