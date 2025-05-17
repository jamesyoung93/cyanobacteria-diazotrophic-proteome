#!/usr/bin/env Rscript
library(data.table)
blast <- fread("results/blastp/blastp_all.out", header=FALSE)
# apply coverage & identity filters, then write MCL input
# ... then call system("mcl filtered.mci -I 1.5 -o results/clusters/clusters.txt")
