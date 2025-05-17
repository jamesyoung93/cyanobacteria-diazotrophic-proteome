#!/usr/bin/env python3
import yaml, os, urllib.request, gzip, shutil

cfg = yaml.safe_load(open("config/species.yaml"))
for grp, entries in cfg.items():
    os.makedirs(f"data/raw/{grp}", exist_ok=True)
    for e in entries:
        out = f"data/raw/{grp}/{e['name']}.faa"
        gz = out + ".gz"
        print(f"Fetching {e['url']} → {gz}")
        urllib.request.urlretrieve(e['url'], gz)
        with gzip.open(gz, "rb") as f_in, open(out, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
