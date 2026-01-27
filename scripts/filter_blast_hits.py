#!/usr/bin/env python3
"""
Filter BLAST tabular output by percent identity, e-value, and (optionally) coverage,
with an optional strict pairwise Reciprocal Best Hit (RBH) reduction.

Input format must match BLAST outfmt used by scripts/03_run_blastp.sh:
    outfmt "6 qseqid sseqid pident length qlen slen evalue"

Output format preserves the same 7 columns.

Species parsing:
- FASTA headers are prefixed as: >{species_name}_{original_header}
  (see scripts/02_build_blast_db.sh). Species names may contain underscores.
- We load species.yaml and match the *longest* known prefixes first.

RBH definition (strict pairwise RBH by species pair):
- For each unordered species pair {A,B}:
  - Consider only hits A->B and B->A.
  - For each query in A, keep its single best hit in B.
  - For each query in B, keep its single best hit in A.
  - Keep a hit (a_q -> b_s) iff b_s's best hit back points to a_q.
- Output is canonicalized to a single direction per pair to avoid duplicates.

Ranking for "best hit":
- lowest evalue, then highest pident, then longest alignment length.
"""

from __future__ import annotations

import argparse
from typing import Dict, List, Tuple

import pandas as pd
import yaml


BLAST_COLS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue"]


def load_species_names(species_yaml: str) -> List[str]:
    cfg = yaml.safe_load(open(species_yaml))
    names: List[str] = []
    for _, entries in cfg.items():
        for e in entries:
            names.append(e["name"])
    # Longest-first matching avoids prefix collisions when names share prefixes.
    names.sort(key=len, reverse=True)
    return names


def assign_species(ids: pd.Series, species_names: List[str]) -> pd.Series:
    """
    Assign species based on known name prefixes in qseqid/sseqid.

    For an ID like: anabaena_sp_pcc7120_WP_....
    We match species_name + "_" prefix.
    """
    s = ids.astype(str)
    out = pd.Series([None] * len(s), index=s.index, dtype="object")
    for name in species_names:
        prefix = name + "_"
        mask = s.str.startswith(prefix)
        out.loc[mask] = name
    return out


def sort_for_besthit(df: pd.DataFrame) -> pd.DataFrame:
    # Lower evalue first; then higher pident; then longer alignment.
    return df.sort_values(
        by=["evalue", "pident", "length"],
        ascending=[True, False, False],
        kind="mergesort",  # stable sort
    )


def strict_pairwise_rbh(df: pd.DataFrame, species_yaml: str) -> pd.DataFrame:
    """
    Strict RBH by species pair, returning canonical directional pairs.

    Canonical direction:
    - For each unordered pair {A,B}, output only A->B where A < B lexicographically.
    """
    species_names = load_species_names(species_yaml)

    df = df.copy()
    df["q_species"] = assign_species(df["qseqid"], species_names)
    df["s_species"] = assign_species(df["sseqid"], species_names)

    before_dropna = len(df)
    df = df.dropna(subset=["q_species", "s_species"])
    print(f"[RBH] species assigned for {len(df)}/{before_dropna} hits")

    # Remove within-species hits for RBH pairing
    df = df[df["q_species"] != df["s_species"]].copy()

    # Build unordered species-pair key
    def pair_key(a: str, b: str) -> Tuple[str, str]:
        return (a, b) if a < b else (b, a)

    df["pair"] = df.apply(lambda r: pair_key(r["q_species"], r["s_species"]), axis=1)

    out_blocks: List[pd.DataFrame] = []

    for (a, b), block in df.groupby("pair", sort=True):
        # Split into the two directions
        ab = block[(block["q_species"] == a) & (block["s_species"] == b)].copy()
        ba = block[(block["q_species"] == b) & (block["s_species"] == a)].copy()
        if ab.empty or ba.empty:
            continue

        # Best hit per query (A) into B
        ab = sort_for_besthit(ab)
        best_ab = ab.drop_duplicates(subset=["qseqid"], keep="first")[BLAST_COLS].copy()

        # Best hit per query (B) into A
        ba = sort_for_besthit(ba)
        best_ba = ba.drop_duplicates(subset=["qseqid"], keep="first")[BLAST_COLS].copy()

        # Map: B_query -> best A_subject
        back_map: Dict[str, str] = dict(zip(best_ba["qseqid"], best_ba["sseqid"]))

        # Keep only reciprocal pairs: A_query's best B_subject must map back to A_query
        mask = best_ab.apply(lambda r: back_map.get(r["sseqid"]) == r["qseqid"], axis=1)
        rbh_ab = best_ab.loc[mask].copy()

        # Canonical direction is a->b only (since (a,b) is ordered)
        out_blocks.append(rbh_ab)

    if not out_blocks:
        return df.iloc[0:0][BLAST_COLS].copy()

    out = pd.concat(out_blocks, ignore_index=True)

    # Ensure uniqueness and stable output
    out = out.drop_duplicates(subset=["qseqid", "sseqid"], keep="first")
    return out[BLAST_COLS]


def main() -> None:
    parser = argparse.ArgumentParser(description="Filter BLAST hits and optionally reduce to strict pairwise RBH.")
    parser.add_argument("input", help="BLAST tabular file (outfmt 6 with 7 columns: qseqid sseqid pident length qlen slen evalue)")
    parser.add_argument("output", help="Filtered output path")
    parser.add_argument("--pident", type=float, default=40.0, help="Minimum percent identity")
    parser.add_argument("--evalue", type=float, default=1e-10, help="Maximum e-value")
    parser.add_argument("--qcov", type=float, default=0.0, help="Minimum query coverage (length / qlen)")
    parser.add_argument("--scov", type=float, default=0.0, help="Minimum subject coverage (length / slen)")
    parser.add_argument("--use_rbh", action="store_true", help="Reduce to strict pairwise reciprocal best hits after filtering")
    parser.add_argument("--species_yaml", default="species.yaml", help="Path to species.yaml (required for --use_rbh)")
    parser.add_argument("--keep_self_hits", action="store_true", help="If set, do not remove qseqid==sseqid hits")

    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", header=None, names=BLAST_COLS)
    print(f"[FILTER] total hits read: {len(df)}")

    if not args.keep_self_hits:
        df = df[df["qseqid"] != df["sseqid"]]
        print(f"[FILTER] after removing self-hits: {len(df)}")

    # Coverage filters (avoid divide-by-zero)
    qlen = df["qlen"].replace(0, pd.NA)
    slen = df["slen"].replace(0, pd.NA)
    qcov = df["length"] / qlen
    scov = df["length"] / slen

    filt = df[
        (df["pident"] >= args.pident)
        & (df["evalue"] <= args.evalue)
        & (qcov >= args.qcov)
        & (scov >= args.scov)
    ].copy()

    print(
        "[FILTER] after thresholds: "
        f"{len(filt)} (pident>={args.pident}, evalue<={args.evalue}, qcov>={args.qcov}, scov>={args.scov})"
    )

    if args.use_rbh:
        before = len(filt)
        filt = strict_pairwise_rbh(filt, species_yaml=args.species_yaml)
        print(f"[FILTER] after RBH reduction: {len(filt)} (from {before})")

    filt.to_csv(args.output, sep="\t", header=False, index=False)
    print(f"[DONE] wrote: {args.output}")


if __name__ == "__main__":
    main()
