#!/usr/bin/env bash
interproscan.sh -i data/db/all_cyano.faa \
  -f tsv -o results/annotations/all_interpro.tsv \
  --goterms --pathways
