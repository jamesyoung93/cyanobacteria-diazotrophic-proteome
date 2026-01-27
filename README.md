# cyanobacteria-diazotrophic-proteome

Comparative proteomics pipeline for generating conservation features used in FOX gene prediction. This pipeline identifies orthologous proteins across cyanobacterial species using strict pairwise reciprocal best hits (RBH) from all-vs-all BLASTp searches.

## Overview

This repository implements the comparative bioinformatics methodology described in the Methods section of our FOX gene prediction manuscript. The key output is a set of conservation features that indicate whether each *Anabaena* sp. PCC 7120 protein has homologs in diazotrophic versus non-diazotrophic cyanobacteria.

**Associated Publication:** Young, J.T. et al. "Predicting FOX gene candidates for oxic nitrogen fixation using multi-omic machine learning and comparative bioinformatics" (2026). *Nature Scientific Reports*.

## Methodology: Strict Pairwise Reciprocal Best Hits

The core of this pipeline is the strict pairwise RBH implementation in `scripts/filter_blast_hits.py`. This approach ensures robust ortholog identification across species.

### RBH Definition

For each unordered species pair {A, B}, we consider all BLAST hits in both directions (A→B and B→A). For each query protein in species A, we identify its single best hit in species B. Similarly, for each query in B, we identify its best hit in A. A hit from protein *a* in species A to protein *b* in species B is retained as an RBH if and only if *b*'s best hit back into species A is *a*. This bidirectional confirmation ensures that the relationship is symmetric and represents the most likely orthologous pair.

### Best Hit Ranking

When multiple hits exist from a query to a target species, we rank them using the following criteria in order: (1) lowest e-value, (2) highest percent identity, (3) longest alignment length. This ranking is applied consistently in both directions before checking reciprocity.

### Filtering Thresholds

The pipeline uses the following quality thresholds, which are specified in the Methods section of the manuscript:

| Parameter | Value | Description |
|-----------|-------|-------------|
| E-value | ≤ 1e-5 | Maximum e-value for BLAST hits |
| Query coverage | ≥ 70% | Minimum alignment length / query length |
| Subject coverage | ≥ 70% | Minimum alignment length / subject length |
| Percent identity | ≥ 30% | Minimum sequence identity |

### Output Canonicalization

To avoid duplicate pairs (A→B and B→A representing the same ortholog relationship), the output is canonicalized to a single direction per species pair. Specifically, for each unordered pair {A, B}, only A→B hits are output, where A < B lexicographically.

## Quick Start

### 1. Clone and Set Up Environment

```bash
git clone https://github.com/jamesyoung93/cyanobacteria-diazotrophic-proteome.git
cd cyanobacteria-diazotrophic-proteome

# Using conda (recommended)
conda create -n proteome python=3.10 -y
conda activate proteome
pip install -r requirements.txt

# Also need BLAST+ installed
conda install -c bioconda blast
```

### 2. Download Proteomes

```bash
python scripts/01_download_proteomes.py --species_yaml species.yaml --outdir data/proteomes
```

### 3. Build BLAST Database

```bash
bash scripts/02_build_blast_db.sh
```

### 4. Run All-vs-All BLASTp

```bash
bash scripts/03_run_blastp.sh
```

### 5. Filter to RBH with Quality Thresholds

```bash
python scripts/filter_blast_hits.py \
  results/blastp/blastp_all.out \
  results/blastp/filtered.p30.q70.s70.rbh1.tsv \
  --pident 30 \
  --evalue 1e-5 \
  --qcov 0.70 \
  --scov 0.70 \
  --use_rbh \
  --species_yaml species.yaml
```

### 6. Generate Wide Table and Conservation Features

```bash
# Convert to wide format (one row per Anabaena protein)
python scripts/blastp_to_wide_table.py \
  --infile results/blastp/filtered.p30.q70.s70.rbh1.tsv \
  --species_yaml species.yaml \
  --query_species anabaena_sp_pcc7120 \
  --outfile results/blastp/anabaena7120_hits_wide.csv

# Compute diazotroph vs non-diazotroph statistics
python scripts/wide_to_diazo_means.py \
  --wide_csv results/blastp/anabaena7120_hits_wide.csv \
  --species_yaml species.yaml \
  --out_csv results/blastp/anabaena7120_diazo_vs_nondiazo_means.csv
```

## Repository Structure

```
cyanobacteria-diazotrophic-proteome/
├── scripts/
│   ├── 01_download_proteomes.py    # Download proteomes from UniProt/NCBI
│   ├── 02_build_blast_db.sh        # Build combined BLAST database
│   ├── 03_run_blastp.sh            # Run all-vs-all BLASTp
│   ├── filter_blast_hits.py        # Core RBH filtering (main methodology)
│   ├── blastp_to_wide_table.py     # Convert to wide format
│   ├── wide_to_diazo_means.py      # Compute diazo vs non-diazo stats
│   ├── find_unique_proteins.py     # Identify lineage-specific proteins
│   ├── identify_cohorts.py         # Group proteins by conservation pattern
│   └── validate_rbh.py             # RBH validation utilities
├── results/
│   ├── blastp/                     # BLAST output and filtered results
│   └── cohorts/                    # Protein groupings by conservation
├── species.yaml                    # Proteome panel definition
├── environment.yaml                # Conda environment specification
├── requirements.txt                # Python dependencies
├── LICENSE
└── README.md
```

## Species Panel

The `species.yaml` file defines the proteome panel used for comparative analysis. The panel includes nine cyanobacterial proteomes organized into four groups based on morphology (filamentous vs unicellular) and nitrogen fixation capability (diazotrophic vs non-diazotrophic).

### Filamentous Diazotrophs (4 organisms)
These are heterocyst-forming cyanobacteria capable of nitrogen fixation:
- *Anabaena* sp. PCC 7120 (query organism)
- *Anabaena variabilis* ATCC 29413
- *Nostoc azollae* 0708
- *Nostoc punctiforme* PCC 73102

### Unicellular Diazotrophs (3 organisms, 5 taxonomy entries)
These are unicellular cyanobacteria capable of nitrogen fixation. Note that some have taxonomic synonyms that are counted separately for RBH feature generation:
- *Crocosphaera watsonii* WH 8501
- *Cyanothece* sp. ATCC 51142 / *Crocosphaera subtropica* ATCC 51142 (counted as 2)
- *Gloeothece* sp. 6803 / *Cyanosphaera ocelotica* WH 8902 (counted as 2)

### Non-Diazotrophs (2 organisms)
These serve as outgroup references lacking nitrogen fixation capability:
- *Arthrospira maxima* CS-328 (filamentous)
- *Synechococcus elongatus* PCC 7942 (unicellular)

## Conservation Features Generated

The pipeline produces several conservation features for each *Anabaena* sp. PCC 7120 protein:

| Feature | Description |
|---------|-------------|
| `n_diazo_hits` | Number of diazotroph species with an RBH hit |
| `n_nondiazo_hits` | Number of non-diazotroph species with an RBH hit |
| `has_nondiazo_hit` | Binary: protein has at least one non-diazotroph homolog |
| `diazo_mean_pident` | Mean percent identity across diazotroph RBH hits |
| `nondiazo_mean_pident` | Mean percent identity across non-diazotroph RBH hits |
| `diazo_minus_nondiazo_mean_pident` | Difference in mean identity (conservation bias) |

These features are incorporated into the machine learning feature matrix used by [FoxGenes_ML](https://github.com/jamesyoung93/FoxGenes_ML).

## Key Script: filter_blast_hits.py

The `filter_blast_hits.py` script is the core of the pipeline. It implements both quality filtering and strict pairwise RBH reduction.

### Usage

```bash
python scripts/filter_blast_hits.py INPUT OUTPUT [OPTIONS]

Required arguments:
  INPUT                 BLAST tabular file (outfmt 6 with 7 columns)
  OUTPUT                Filtered output path

Options:
  --pident FLOAT        Minimum percent identity (default: 40.0)
  --evalue FLOAT        Maximum e-value (default: 1e-10)
  --qcov FLOAT          Minimum query coverage (default: 0.0)
  --scov FLOAT          Minimum subject coverage (default: 0.0)
  --use_rbh             Reduce to strict pairwise RBH after filtering
  --species_yaml PATH   Path to species.yaml (required for --use_rbh)
  --keep_self_hits      Do not remove qseqid==sseqid hits
```

### Example with Manuscript Parameters

```bash
python scripts/filter_blast_hits.py \
  results/blastp/blastp_all.out \
  results/blastp/filtered.p30.q70.s70.rbh1.tsv \
  --pident 30 \
  --evalue 1e-5 \
  --qcov 0.70 \
  --scov 0.70 \
  --use_rbh \
  --species_yaml species.yaml
```

### Expected Output

The script prints filtering statistics at each step:

```
[FILTER] total hits read: 2847291
[FILTER] after removing self-hits: 2784562
[FILTER] after thresholds: 158432 (pident>=30, evalue<=1e-05, qcov>=0.7, scov>=0.7)
[RBH] species assigned for 158432/158432 hits
[FILTER] after RBH reduction: 24891 (from 158432)
[DONE] wrote: results/blastp/filtered.p30.q70.s70.rbh1.tsv
```

## Input Data Format

The pipeline expects BLAST output in tabular format with 7 columns, as produced by:

```bash
blastp -outfmt "6 qseqid sseqid pident length qlen slen evalue"
```

Columns are: query sequence ID, subject sequence ID, percent identity, alignment length, query length, subject length, and e-value.

### FASTA Header Convention

For species assignment to work correctly, FASTA headers in the combined database must be prefixed with the species name followed by an underscore. The `02_build_blast_db.sh` script handles this automatically:

```
>anabaena_sp_pcc7120_WP_010994653.1 hypothetical protein
>nostoc_punctiforme_pcc73102_WP_012409876.1 nitrogenase iron protein
```

## Related Repositories

The main machine learning pipeline is available at [FoxGenes_ML](https://github.com/jamesyoung93/FoxGenes_ML). An interactive web application for exploring predictions is available at [FoxGenesApp](https://github.com/jamesyoung93/FoxGenesApp).

## Citation

```
Young, J.T. et al. (2026). Predicting FOX gene candidates for oxic nitrogen fixation 
using multi-omic machine learning and comparative bioinformatics. 
Nature Scientific Reports. [DOI pending]
```

## License

MIT License - see [LICENSE](LICENSE) for details.
