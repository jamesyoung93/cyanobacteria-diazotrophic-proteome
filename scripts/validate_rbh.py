#!/usr/bin/env python3
import argparse, re
from collections import defaultdict

def parse_args():
    ap = argparse.ArgumentParser(description="Validate RBH reduction behavior and reconstruct pairwise RBH from non-RBH hits.")
    ap.add_argument("--all_faa", required=True, help="FASTA used to build BLAST DB (e.g., data/db/all_cyano.faa)")
    ap.add_argument("--blast_all", required=True, help="Filtered hits without RBH (e.g., blastp_filtered.p40.rbh0.out)")
    ap.add_argument("--blast_rbh", required=True, help="Filtered hits with RBH (e.g., blastp_filtered.p40.rbh1.out)")
    ap.add_argument("--bitscore_col", type=int, default=12, help="1-indexed bitscore column in outfmt6 (default 12)")
    ap.add_argument("--skip_within_proteome", action="store_true", help="Skip hits where query and subject are from same proteome prefix.")
    return ap.parse_args()

def build_id2sp_from_all_faa(all_faa: str):
    id2sp = {}
    with open(all_faa) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            head = line[1:].strip()
            tok = head.split()[0]
            # Your headers are like: anabaena_sp_pcc7120_WP_...
            # species/proteome key can be inferred as prefix before last underscore chunk(s)
            # Here we use the first 3 underscore-delimited tokens as a stable "proteome key" heuristic.
            # If you want exact mapping, define it explicitly in species.yaml and parse that instead.
            parts = tok.split("_")
            if len(parts) >= 3:
                sp = "_".join(parts[:3])  # e.g., anabaena_sp_pcc7120
            else:
                sp = "UNKNOWN"
            id2sp[tok] = sp
    return id2sp

def main():
    args = parse_args()
    bitscore_idx = args.bitscore_col - 1

    id2sp = build_id2sp_from_all_faa(args.all_faa)

    # Load RBH pairs
    rbh_pairs = set()
    with open(args.blast_rbh) as f:
        for line in f:
            if not line.strip():
                continue
            q, s = line.split("\t", 2)[:2]
            rbh_pairs.add((q, s))

    # Reconstruct pairwise RBH from non-RBH file
    best_s_for_q_to_sp = {}
    best_q_for_s_to_sp = {}

    total_lines = 0
    used_lines = 0
    with open(args.blast_all) as f:
        for line in f:
            if not line.strip():
                continue
            total_lines += 1
            parts = line.rstrip("\n").split("\t")
            q, s = parts[0], parts[1]
            if q not in id2sp or s not in id2sp:
                continue
            used_lines += 1

            qsp = id2sp[q]
            ssp = id2sp[s]
            if args.skip_within_proteome and qsp != "UNKNOWN" and ssp != "UNKNOWN" and qsp == ssp:
                continue

            try:
                sc = float(parts[bitscore_idx])
            except Exception:
                sc = float(parts[-1])

            k1 = (q, ssp)
            cur = best_s_for_q_to_sp.get(k1)
            if cur is None or sc > cur[1]:
                best_s_for_q_to_sp[k1] = (s, sc)

            k2 = (s, qsp)
            cur2 = best_q_for_s_to_sp.get(k2)
            if cur2 is None or sc > cur2[1]:
                best_q_for_s_to_sp[k2] = (q, sc)

    pairwise_mutual = set()
    for (q, ssp), (s, sc) in best_s_for_q_to_sp.items():
        qsp = id2sp.get(q, "UNKNOWN")
        back = best_q_for_s_to_sp.get((s, qsp))
        if back and back[0] == q:
            pairwise_mutual.add((q, s))

    overlap = len(rbh_pairs & pairwise_mutual)
    print(f"[RECON] parsed {used_lines}/{total_lines} lines")
    print("RBH pairs (file):", len(rbh_pairs))
    print("Pairwise mutual-best (reconstructed):", len(pairwise_mutual))
    print("Overlap:", overlap)
    print("RBH covered by pairwise mutual-best:", overlap/len(rbh_pairs) if rbh_pairs else 0.0)
    print("Pairwise mutual-best covered by RBH:", overlap/len(pairwise_mutual) if pairwise_mutual else 0.0)

if __name__ == "__main__":
    main()
