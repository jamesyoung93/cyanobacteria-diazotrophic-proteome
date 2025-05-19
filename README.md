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

The Uniprot proteome links use HTTPS. If you encounter download issues,
verify the URLs in `species.yaml`.
