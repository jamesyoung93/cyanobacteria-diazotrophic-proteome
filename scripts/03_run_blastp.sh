#!/usr/bin/env bash
set -euo pipefail

mkdir -p results/blastp
CPU=${SLURM_CPUS_ON_NODE:-16}

echo "[BLASTP] all-vs-all → results/blastp/blastp_all.out"
blastp \
  -db data/db/cyano_db \
  -query data/db/all_cyano.faa \
  -evalue 1e-5 \
  -num_threads $CPU \
  -outfmt "6 qseqid sseqid pident length qlen slen" \
  -out results/blastp/blastp_all.out

echo "[DONE] blastp_all.out"
