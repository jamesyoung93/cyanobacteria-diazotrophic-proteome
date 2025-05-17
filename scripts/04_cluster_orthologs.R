#!/usr/bin/env Rscript
library(data.table)

# Paths
blast_file   <- "results/blastp/blastp_all.out"
out_dir      <- "results/clusters"
mcl_input    <- file.path(out_dir, "blast.mci")
clusters_txt <- file.path(out_dir, "clusters.txt")
mapping_tsv  <- file.path(out_dir, "clusters_mapping.tsv")

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# 1) Read and filter BLAST hits
cols <- c("qseqid","sseqid","pident","length","qlen","slen")
blast <- fread(blast_file, header=FALSE, col.names=cols)

filtered <- blast[
  pident >= 30 &
  (length / qlen) >= 0.5 &
  (length / slen) >= 0.5
]

# 2) Prepare MCL input (tab-delimited edge list: node1 node2 weight)
edges <- filtered[, .(qseqid, sseqid, weight=pident)]
fwrite(edges, mcl_input, sep = "\t", col.names = FALSE)

# 3) Run MCL
message("[MCL] clustering…")
system(sprintf("mcl %s -I 1.5 -o %s", mcl_input, clusters_txt))

# 4) Parse clusters.txt → clusters_mapping.tsv
lines  <- readLines(clusters_txt)
mapping <- data.table(clusterID = integer(), member = character())

for (i in seq_along(lines)) {
  members <- strsplit(lines[i], "\t")[[1]]
  mapping <- rbind(mapping,
                   data.table(clusterID = i, member = members))
}

fwrite(mapping, mapping_tsv, sep="\t", col.names=FALSE)
message("[DONE] clusters and mapping written to ", out_dir)
