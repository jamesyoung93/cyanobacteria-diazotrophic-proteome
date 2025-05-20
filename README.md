# cyanobacteria-diazotrophic-proteome

This repository contains a Snakemake workflow for collecting cyanobacterial proteomes,
building an all-vs-all BLAST database, clustering orthologs, annotating with
InterProScan, and running enrichment analysis.

## Usage

1. Create and activate the conda environment:

```bash
mamba env create -f environment.yaml -n cyano-diazotroph
conda activate cyano-diazotroph
```

2. Run the workflow (adjust `--cores` as appropriate):

```bash
snakemake --cores 16
```

To only fetch the raw FASTA files you can run:

```bash
snakemake data/raw/.download_complete
```

After clustering, you can inspect the cluster size distribution with:

```bash
python3 scripts/cluster_size_distribution.py
```

The Uniprot proteome links use HTTPS. If you encounter download issues,
verify the URLs in `species.yaml`.

### Adjusting BLAST filtering

The workflow now keeps the raw BLAST results and allows filtering by
percent identity and e-value prior to clustering. Default thresholds are
40% identity and `1e-10` for the e-value. You can tweak these by editing
the `filter_blast_hits` rule in the `Snakefile` or running the filtering
utility directly:

```bash
scripts/filter_blast_hits.py results/blastp/blastp_all.out \
    results/blastp/blastp_filtered.out --pident 50 --evalue 1e-20
```
