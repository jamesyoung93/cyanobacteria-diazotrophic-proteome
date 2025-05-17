#!/usr/bin/env bash
set -euo pipefail
mkdir -p results/blastp
blastp \
  -db data/db/cyano_db \
  -query data/db/all_cyano.faa \
  -evalue 1e-5 \
  -num_threads ${SLURM_CPUS_ON_NODE:-16} \
  -outfmt 6 \
  -out results/blastp/blastp_all.out
