#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/db
cat data/raw/*/*.faa > data/db/all_cyano.faa

echo "[MAKEBLASTDB] building protein DB…"
makeblastdb \
  -in data/db/all_cyano.faa \
  -dbtype prot \
  -title cyano_db \
  -out data/db/cyano_db

echo "[DONE] BLAST DB in data/db/cyano_db.*"
