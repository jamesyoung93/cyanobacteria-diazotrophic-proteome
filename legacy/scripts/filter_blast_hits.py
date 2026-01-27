#!/usr/bin/env python3
"""Filter BLAST tabular output by percent identity and e-value."""
import argparse
import pandas as pd

# CLI
parser = argparse.ArgumentParser(description="Filter BLAST hits")
parser.add_argument("input", help="BLAST tabular file")
parser.add_argument("output", help="Filtered output path")
parser.add_argument("--pident", type=float, default=40.0, help="Minimum percent identity")
parser.add_argument("--evalue", type=float, default=1e-10, help="Maximum e-value")
args = parser.parse_args()

cols = ["qseqid","sseqid","pident","length","qlen","slen","evalue"]

df = pd.read_csv(args.input, sep="\t", header=None, names=cols)

filt = df[(df.pident >= args.pident) & (df.evalue <= args.evalue)]

filt.to_csv(args.output, sep="\t", header=False, index=False)
