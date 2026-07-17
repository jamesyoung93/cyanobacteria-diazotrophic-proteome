# Comparative proteomics for cyanobacterial diazotrophy

This repository builds reciprocal-best-hit conservation features used by the published FOX gene candidate-ranking workflow. It compares *Anabaena* sp. PCC 7120 proteins with a panel of filamentous diazotrophs, unicellular diazotrophs, and non-diazotrophs so that gene rankings can incorporate lineage-aware evidence while reducing ortholog and paralog ambiguity.

**Paper:** [Predicting FOX gene candidates for oxic nitrogen fixation using multi-omic machine learning and comparative bioinformatics](https://doi.org/10.1038/s41598-026-41873-w), *Scientific Reports* 16, 11412 (2026)

**Interactive app:** [FOX Gene Complement Explorer](https://foxgenesapp-h6mlszof5jdd6p2bpq22wr.streamlit.app/)

## Key result

This repository supplies the comparative-bioinformatics component of the associated study. The paper combines reciprocal-best-hit conservation with transcriptomics, proteomics, promoter architecture, Gene Ontology, and genomic context to rank FOX gene candidates. Full model results are reported in the paper and reproduced in [FoxGenes_ML](https://github.com/jamesyoung93/FoxGenes_ML).

## Methods and stack

- Bidirectional BLASTp with reciprocal-best-hit filtering
- Minimum 30% amino-acid identity
- Minimum 70% query and subject coverage
- Maximum e-value of `1e-5`
- Python, pandas, Biopython, Bash, Snakemake, and NCBI BLAST+ 2.14.1

## Install

```bash
git clone https://github.com/jamesyoung93/cyanobacteria-diazotrophic-proteome.git
cd cyanobacteria-diazotrophic-proteome
conda env create -f environment.yaml
conda activate proteome
```

If you install manually, use `requirements.txt` and install NCBI BLAST+ 2.14.1 from Bioconda.

## Run the pipeline

The repository includes a `Snakefile` and individual scripts. The explicit workflow is:

```bash
python scripts/01_download_proteomes.py --species_yaml species.yaml --outdir data/proteomes
bash scripts/02_build_blast_db.sh
bash scripts/03_run_blastp.sh

python scripts/filter_blast_hits.py \
  results/blastp/blastp_all.out \
  results/blastp/filtered.p30.q70.s70.rbh1.tsv \
  --pident 30 --evalue 1e-5 --qcov 0.70 --scov 0.70 \
  --use_rbh --species_yaml species.yaml

python scripts/blastp_to_wide_table.py \
  --infile results/blastp/filtered.p30.q70.s70.rbh1.tsv \
  --species_yaml species.yaml \
  --query_species anabaena_sp_pcc7120 \
  --outfile results/blastp/anabaena7120_hits_wide.csv

python scripts/wide_to_diazo_means.py \
  --wide_csv results/blastp/anabaena7120_hits_wide.csv \
  --species_yaml species.yaml \
  --out_csv results/blastp/anabaena7120_diazo_vs_nondiazo_means.csv
```

## Reciprocal-best-hit definition

For each species pair A and B, an A-to-B hit is retained only when B's best hit back into A returns the original query. Hits are ordered by lowest e-value, then highest percent identity, then longest alignment.

## Species panel

The exact accessions and taxonomic aliases are defined in `species.yaml` and Supplementary Table S3 of the paper.

| Group | Representative organisms |
|---|---|
| Filamentous diazotrophs | *Anabaena* 7120, *A. variabilis* 29413, *Nostoc azollae* 0708, *N. punctiforme* 73102 |
| Unicellular diazotrophs | *Crocosphaera watsonii* 8501 and *Cyanothece* 51142, including documented taxonomic synonyms |
| Non-diazotrophs | *Arthrospira maxima* CS-328 and *Synechococcus elongatus* 7942 |

## Outputs

For each *Anabaena* 7120 protein, the pipeline produces diazotroph and non-diazotroph RBH counts, mean percent identities, and lineage-bias features used by the machine-learning pipeline and app filters.

## Data availability

Reference proteome identifiers are recorded in `species.yaml`. Derived result tables needed to inspect the released features are included under `results/`. Proteomes can be downloaded from their public sources by the pipeline. Source database records retain their original terms.

## Related repositories

- [FoxGenes_ML](https://github.com/jamesyoung93/FoxGenes_ML): multi-omic candidate-ranking pipeline
- [FoxGenesApp](https://github.com/jamesyoung93/FoxGenesApp): interactive candidate and gene-complement explorer

## Citation

```bibtex
@article{Young2026FoxGenes,
  author  = {Young, James and Gu, Liping and Zhou, Ruanbao},
  title   = {Predicting FOX gene candidates for oxic nitrogen fixation using multi-omic machine learning and comparative bioinformatics},
  journal = {Scientific Reports},
  year    = {2026},
  volume  = {16},
  pages   = {11412},
  doi     = {10.1038/s41598-026-41873-w}
}
```

## License

Code is available under the [MIT License](LICENSE). Downloaded proteomes and other source data retain their original terms.
