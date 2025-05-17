#!/usr/bin/env bash
set -euo pipefail

mkdir -p results/annotations

echo "[InterProScan] annotating all_cyano.faa…"
interproscan.sh \
  -i data/db/all_cyano.faa \
  -f tsv \
  -o results/annotations/all_interpro.tsv \
  --goterms \
  --pathways \
  --iprlookup \
  --cpu ${SLURM_CPUS_ON_NODE:-16}

echo "[DONE] all_interpro.tsv"
