## cyanobacteria-diazotrophic-proteome

This repository contains a minimal Snakemake workflow to fetch cyanobacterial proteomes,
build an all‑vs‑all BLAST database and identify proteins conserved within
different groups of diazotrophs.  The final output consists of three lists of
protein accessions describing the cohorts illustrated below:

1. `unicellular_specific.tsv` – proteins present in all unicellular diazotrophic
   cyanobacteria. Hits to filamentous diazotrophs may also be reported and can
   be filtered out downstream if desired.
2. `filamentous_specific.tsv` – proteins conserved across filamentous
   diazotrophic cyanobacteria. Hits to unicellular diazotrophs may likewise be
   included.
3. `diazotroph_common.tsv` – proteins found in both unicellular and filamentous
   diazotrophic species while absent from non‑diazotrophic cyanobacteria.

## Usage

1. Create and activate the conda environment:

```bash
mamba env create -f environment.yaml -n cyano-diazotroph
conda activate cyano-diazotroph
```

2. Run the workflow (adjust `--cores` as needed):

```bash
snakemake --cores 16
```

The workflow will download the configured proteomes, build the BLAST database,
run the all‑vs‑all search, filter the results and produce the three cohort files
in `results/cohorts/`.

### Adjusting BLAST filtering

By default BLAST hits are filtered at 80% identity with a maximum e‑value of
`1e-10`.  You can override these thresholds either by calling the filtering
script directly or by passing parameters to Snakemake:

```bash
snakemake --cores 16 --config pident=70 evalue=1e-20
```

This will propagate the values to the filtering step and subsequent cohort
calculation.
