#!/usr/bin/env bash
set -euo pipefail
mkdir -p data/db
cat data/raw/*/*.faa > data/db/all_cyano.faa
makeblastdb -in data/db/all_cyano.faa -dbtype prot -out data/db/cyano_db
