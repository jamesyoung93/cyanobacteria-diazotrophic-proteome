#!/usr/bin/env python3
import argparse
import pandas as pd
import yaml

COLS = ["qseqid","sseqid","pident","length","qlen","slen","evalue"]

def load_species_names(species_yaml: str):
    cfg = yaml.safe_load(open(species_yaml))
    names = []
    for _, entries in cfg.items():
        for e in entries:
            names.append(e["name"])
    names.sort(key=len, reverse=True)
    return names

def infer_species(id_series: pd.Series, species_names):
    s = id_series.astype(str)
    out = pd.Series([None]*len(s), index=s.index, dtype="object")
    for name in species_names:
        pref = name + "_"
        m = s.str.startswith(pref)
        out.loc[m] = name
    return out

def strip_species_prefix(pid: str, species_names):
    for name in species_names:
        pref = name + "_"
        if pid.startswith(pref):
            return pid[len(pref):]
    return pid

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True, help="BLAST filtered file (7 columns, no header)")
    ap.add_argument("--species_yaml", default="species.yaml")
    ap.add_argument("--query_species", default="anabaena_sp_pcc7120",
                    help="Species prefix for the query organism to make 1 row per query protein.")
    ap.add_argument("--outfile", required=True, help="Output CSV path")
    args = ap.parse_args()

    df = pd.read_csv(args.infile, sep="\t", header=None, names=COLS)

    # numeric types
    df["pident"] = df["pident"].astype(float)
    df["length"] = df["length"].astype(int)
    df["qlen"] = df["qlen"].astype(int)
    df["slen"] = df["slen"].astype(int)
    df["evalue"] = df["evalue"].astype(float)

    # coverages
    df["qcov"] = df["length"] / df["qlen"].replace(0, pd.NA)
    df["scov"] = df["length"] / df["slen"].replace(0, pd.NA)

    species_names = load_species_names(args.species_yaml)
    df["q_species"] = infer_species(df["qseqid"], species_names)
    df["s_species"] = infer_species(df["sseqid"], species_names)

    # keep only rows where query is the requested species
    df = df[df["q_species"] == args.query_species].copy()
    if df.empty:
        raise SystemExit(f"No rows with q_species == {args.query_species}. Check --query_species value.")

    # clean ids for readability
    df["q_id"] = df["qseqid"].apply(lambda x: strip_species_prefix(str(x), species_names))
    df["s_id"] = df["sseqid"].apply(lambda x: strip_species_prefix(str(x), species_names))

    # Enforce one hit per (q_id, s_species) by best ranking:
    # lowest evalue, then highest pident, then highest qcov
    df = df.sort_values(["q_id","s_species","evalue","pident","qcov"], ascending=[True,True,True,False,False])
    df = df.drop_duplicates(subset=["q_id","s_species"], keep="first")

    # Build wide columns per target species
    value_cols = ["s_id","pident","evalue","qcov","scov","length","slen"]
    wide = df.pivot(index="q_id", columns="s_species", values=value_cols)

    # Flatten MultiIndex columns: e.g., pident__nostoc_punctiforme_pcc73102
    wide.columns = [f"{metric}__{sp}" for metric, sp in wide.columns]
    wide = wide.reset_index().rename(columns={"q_id":"anabaena_sp_pcc7120_protein"})

    wide.to_csv(args.outfile, index=False)
    print(f"[DONE] wrote {args.outfile} with {len(wide)} query proteins")

if __name__ == "__main__":
    main()
