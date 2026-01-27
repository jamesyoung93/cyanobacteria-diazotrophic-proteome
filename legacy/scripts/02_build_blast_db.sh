#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/db

# Build concatenated FASTA with species prefixes so downstream scripts can
# recover the species name from sequence IDs. Each sequence header is rewritten
# as ">speciesID_originalHeader" where "speciesID" comes from the FASTA file name.
> data/db/all_cyano.faa
for faa in data/raw/*/*.faa; do
    sp=$(basename "$faa" .faa)
    awk -v pre="${sp}_" '{
        if(substr($0,1,1)==">") {
            print ">" pre substr($0,2)
        } else {
            print
        }
    }' "$faa" >> data/db/all_cyano.faa
done

echo "[MAKEBLASTDB] building protein DB…"
makeblastdb \
  -in data/db/all_cyano.faa \
  -dbtype prot \
  -title cyano_db \
  -out data/db/cyano_db

echo "[DONE] BLAST DB in data/db/cyano_db.*"
