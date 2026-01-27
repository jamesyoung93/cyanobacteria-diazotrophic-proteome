#!/usr/bin/env python3
import argparse
import pandas as pd
import yaml

def load_species_groups(species_yaml: str):
    """
    Expects species.yaml to contain top-level groups like:
      filamentous_diazotrophic, unicellular_diazotrophic,
      filamentous_non_diazotrophic, unicellular_non_diazotrophic
    Each group is a list of dicts with key "name".
    """
    cfg = yaml.safe_load(open(species_yaml))

    diazo = set()
    nondiazo = set()

    for group_name, entries in cfg.items():
        names = {e["name"] for e in entries}
        g = group_name.lower()
        if "diazotrophic" in g and "non" not in g:
            diazo |= names
        elif "non_diazotrophic" in g or ("non" in g and "diazotrophic" in g):
            nondiazo |= names

    return diazo, nondiazo

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wide_csv", required=True, help="Wide table produced by blastp_to_wide_table.py")
    ap.add_argument("--species_yaml", default="species.yaml")
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--query_species", default="anabaena_sp_pcc7120")
    args = ap.parse_args()

    wide = pd.read_csv(args.wide_csv)

    diazo, nondiazo = load_species_groups(args.species_yaml)

    # pident columns are named like pident__<species>
    pident_cols = [c for c in wide.columns if c.startswith("pident__")]
    if not pident_cols:
        raise SystemExit("No pident__* columns found. Did you pass the correct wide CSV?")

    # Determine which species are represented in the wide file
    species_in_file = {c.split("__", 1)[1] for c in pident_cols}

    diazo_species = sorted(list(diazo & species_in_file))
    nondiazo_species = sorted(list(nondiazo & species_in_file))

    # Build lists of column names for each cohort
    diazo_cols = [f"pident__{sp}" for sp in diazo_species]
    nondiazo_cols = [f"pident__{sp}" for sp in nondiazo_species]

    # Compute mean and counts ignoring NaN (i.e., only where hit exists)
    out = pd.DataFrame()
    out["anabaena_sp_pcc7120_protein"] = wide["anabaena_sp_pcc7120_protein"]

    if diazo_cols:
        out["diazo_mean_pident"] = wide[diazo_cols].mean(axis=1, skipna=True)
        out["diazo_hit_count"] = wide[diazo_cols].notna().sum(axis=1)
    else:
        out["diazo_mean_pident"] = pd.NA
        out["diazo_hit_count"] = 0

    if nondiazo_cols:
        out["nondiazo_mean_pident"] = wide[nondiazo_cols].mean(axis=1, skipna=True)
        out["nondiazo_hit_count"] = wide[nondiazo_cols].notna().sum(axis=1)
    else:
        out["nondiazo_mean_pident"] = pd.NA
        out["nondiazo_hit_count"] = 0

    # Optional: also compute delta
    out["diazo_minus_nondiazo_mean_pident"] = out["diazo_mean_pident"] - out["nondiazo_mean_pident"]

    out.to_csv(args.out_csv, index=False)

    print("[DONE] wrote", args.out_csv)
    print("[INFO] diazo species used:", len(diazo_cols), " / nondiazo species used:", len(nondiazo_cols))
    if diazo_species:
        print("[INFO] diazo species:", ", ".join(diazo_species))
    if nondiazo_species:
        print("[INFO] nondiazo species:", ", ".join(nondiazo_species))

if __name__ == "__main__":
    main()
