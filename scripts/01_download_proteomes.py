#!/usr/bin/env python3
import yaml
import os
import urllib.request
import gzip
import shutil

# Load species configuration
with open("config/species.yaml") as fh:
    species_cfg = yaml.safe_load(fh)

for group, entries in species_cfg.items():
    out_dir = f"data/raw/{group}"
    os.makedirs(out_dir, exist_ok=True)
    for e in entries:
        name = e["name"]
        url  = e["url"]
        gz_path = os.path.join(out_dir, f"{name}.faa.gz")
        faa_path = os.path.join(out_dir, f"{name}.faa")

        if os.path.exists(faa_path):
            print(f"[SKIP] {faa_path} already exists")
            continue

        print(f"Downloading {url} → {gz_path}")
        urllib.request.urlretrieve(url, gz_path)

        print(f"Decompressing to {faa_path}")
        with gzip.open(gz_path, "rb") as fin, open(faa_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)

        os.remove(gz_path)
        print(f"[DONE] {name}")
