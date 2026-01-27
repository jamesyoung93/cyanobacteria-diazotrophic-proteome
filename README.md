# cyanobacteria-diazotrophic-proteome

Comparative proteomics pipeline for generating conservation features used in FOX gene prediction. Identifies orthologous proteins across cyanobacterial species using strict pairwise reciprocal best hits (RBH).

**Associated Publication:** Young, J., Gu, L. & Zhou, R. "Predicting FOX gene candidates for oxic nitrogen fixation using multi-omic machine learning and comparative bioinformatics." *Nature Scientific Reports* (submitted).

## Setup

```bash
git clone https://github.com/jamesyoung93/cyanobacteria-diazotrophic-proteome.git
cd cyanobacteria-diazotrophic-proteome

conda create -n proteome python=3.10 -y
conda activate proteome
pip install -r requirements.txt
conda install -c bioconda blast
```

## Pipeline

```bash
# 1. Download proteomes
python scripts/01_download_proteomes.py --species_yaml species.yaml --outdir data/proteomes

# 2. Build BLAST database (prefixes headers with species names)
bash scripts/02_build_blast_db.sh

# 3. Run all-vs-all BLASTp (e-value 1e-5)
bash scripts/03_run_blastp.sh

# 4. Filter to RBH with manuscript parameters
python scripts/filter_blast_hits.py \
  results/blastp/blastp_all.out \
  results/blastp/filtered.p30.q70.s70.rbh1.tsv \
  --pident 30 --evalue 1e-5 --qcov 0.70 --scov 0.70 \
  --use_rbh --species_yaml species.yaml

# 5. Convert to wide format (one row per Anabaena protein)
python scripts/blastp_to_wide_table.py \
  --infile results/blastp/filtered.p30.q70.s70.rbh1.tsv \
  --species_yaml species.yaml \
  --query_species anabaena_sp_pcc7120 \
  --outfile results/blastp/anabaena7120_hits_wide.csv

# 6. Compute diazotroph vs non-diazotroph statistics
python scripts/wide_to_diazo_means.py \
  --wide_csv results/blastp/anabaena7120_hits_wide.csv \
  --species_yaml species.yaml \
  --out_csv results/blastp/anabaena7120_diazo_vs_nondiazo_means.csv
```

## RBH Methodology

For each species pair {A, B}, a hit A→B is kept as an RBH only if B's best hit back into A returns the original query. Best hit ranking uses lowest e-value, then highest percent identity, then longest alignment length.

**Filtering thresholds (from Methods):**
- E-value ≤ 1e-5
- Query coverage ≥ 70%
- Subject coverage ≥ 70%
- Percent identity ≥ 30%

## Species Panel

Defined in `species.yaml` (see Supplementary Table S3 for accessions):

| Group | Organisms |
|-------|-----------|
| Filamentous diazotrophs | *Anabaena* 7120, *A. variabilis* 29413, *Nostoc azollae* 0708, *N. punctiforme* 73102 |
| Unicellular diazotrophs | *Crocosphaera watsonii* 8501, *Cyanothece* 51142, plus taxonomic synonyms (see S3) |
| Non-diazotrophs | *Arthrospira maxima* CS-328, *Synechococcus elongatus* 7942 |

## Output Features

Conservation features computed for each *Anabaena* 7120 protein:

| Feature | Description |
|---------|-------------|
| `n_diazo_hits` | RBH count across diazotroph proteomes |
| `n_nondiazo_hits` | RBH count across non-diazotroph proteomes |
| `has_nondiazo_hit` | Binary indicator |
| `diazo_mean_pident` | Mean % identity to diazotroph RBHs |
| `diazo_minus_nondiazo_mean_pident` | Conservation bias toward diazotrophs |

These features feed into the ML pipeline at [FoxGenes_ML](https://github.com/jamesyoung93/FoxGenes_ML).

## Repository Structure

```
├── scripts/
│   ├── 01_download_proteomes.py
│   ├── 02_build_blast_db.sh
│   ├── 03_run_blastp.sh
│   ├── filter_blast_hits.py      # Core RBH implementation
│   ├── blastp_to_wide_table.py
│   ├── wide_to_diazo_means.py
│   └── validate_rbh.py
├── species.yaml
├── results/blastp/
└── legacy/                        # Previous implementation (archived)
```

## Related Repositories

- [FoxGenes_ML](https://github.com/jamesyoung93/FoxGenes_ML) — Machine learning pipeline
- [FoxGenesApp](https://github.com/jamesyoung93/FoxGenesApp) — Interactive web tool

## License

MIT
