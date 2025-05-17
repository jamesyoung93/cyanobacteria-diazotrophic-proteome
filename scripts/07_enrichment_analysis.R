#!/usr/bin/env Rscript
library(tidyverse)
mat <- read_csv("results/matrices/pres_abs.csv")
annot <- read_tsv("results/annotations/all_interpro.tsv")
# Fisher’s exact test per GO term / Pfam domain
# ggplot2 → results/figures/enrichment_bar.pdf
# heatmap of presence/absence → results/figures/presence_absence_heatmap.pdf
