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

The BLAST database is built from all downloaded FASTA files with each
sequence header prefixed by its species name (taken from `species.yaml`).
Downstream scripts rely on this convention to map proteins back to their
source species.

### Identifying diazotroph-specific proteins

If you only want the BLAST search and a list of proteins that lack hits in non-diazotrophic genomes, run the workflow up to the filtering step:

```bash
snakemake results/blastp/blastp_filtered.out --cores 16
```

Then run the helper script to summarise unique proteins:

```bash
python3 scripts/find_unique_proteins.py \
    results/blastp/blastp_filtered.out results/unique_lists
```

This will produce tables in `results/unique_lists` listing proteins found only in diazotrophic species, and subsets unique to filamentous or unicellular diazotrophs.
