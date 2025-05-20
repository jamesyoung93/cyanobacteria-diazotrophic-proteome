#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(ggplot2)
library(yaml)

# 1) Load species groups
cfg <- yaml::read_yaml("species.yaml")
diazo      <- unlist(lapply(cfg[grep("diazo", names(cfg))], function(x) sapply(x, `[[`, "name")))
nondiazo   <- unlist(lapply(cfg[grep("non-?diazo", names(cfg), ignore.case=TRUE)], function(x) sapply(x, `[[`, "name")))

# 2) Read presence/absence matrix
mat <- fread("results/matrices/pres_abs.csv")
setDF(mat)
rownames(mat) <- mat$clusterID
mat$clusterID <- NULL

# 3) Identify diazo-specific clusters
is_diazo_sp <- apply(mat[, diazo] == 1, 1, all) &
               apply(mat[, nondiazo] == 0, 1, all)
diazo_clusters <- rownames(mat)[is_diazo_sp]
bg_clusters    <- rownames(mat)

# 4) Load InterProScan results & mapping
ipr_raw <- fread(
  "results/annotations/all_interpro.tsv",
  sep        = "\t",
  header     = FALSE,
  fill       = TRUE
)
col.names <- c("protein_acc","md5","seq_length","analysis",
               "signature_acc","signature_desc","start","end",
               "score","status","date","ipr_acc","ipr_desc",
               "GO","pathways")
ipr <- setnames(ipr_raw[, 1:length(col.names), with=FALSE], col.names)
                            
mapping <- fread("results/clusters/clusters_mapping.tsv",
                 header=FALSE, col.names=c("clusterID","member"))
mapping[, clusterID := as.character(clusterID)]

# 5) Map ipr → clusters
ipr_map <- merge(ipr, mapping, by.x="protein_acc", by.y="member")
ipr_map[, clusterID := as.character(clusterID)]

# 6) Build unique cluster–domain table
clus_dom <- unique(ipr_map[, .(clusterID, signature_acc, signature_desc)])
clus_dom[, signature := paste0(signature_acc,": ", signature_desc)]

# 7) Enrichment test per signature
fg <- diazo_clusters
bg <- bg_clusters
fg_size <- length(fg)
bg_size <- length(bg)

enrich <- clus_dom[, .(
  fg_with    = sum(clusterID %in% fg),
  total_with = .N
), by=signature]

enrich <- enrich %>%
  mutate(
    fg_without = fg_size - fg_with,
    bg_with    = total_with - fg_with,
    bg_without = bg_size - total_with + bg_with - fg_without,
    pvalue     = mapply(
                   function(a,b,c,d) {
                     fisher.test(matrix(c(a,b,c,d), 2))$p.value
                   },
                   fg_with, fg_without, bg_with, bg_without),
    adj_pvalue = p.adjust(pvalue, method="BH")
  ) %>%
  arrange(adj_pvalue)

# 8) Save results & plot
dir.create("results/figures", recursive=TRUE, showWarnings=FALSE)
fwrite(enrich, "results/figures/enrichment_results.tsv", sep="\t")

top20 <- head(enrich, 20)
ggplot(top20, aes(x = reorder(signature, -log10(adj_pvalue)),
                  y = -log10(adj_pvalue))) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(
    x = "InterPro Signature",
    y = "-log10(adj. p-value)",
    title = "Top 20 enriched domains in diazotroph-specific clusters"
  ) +
  theme_minimal() -> p

ggsave("results/figures/enrichment_bar.pdf", p, width=8, height=6)

# 9) Extract lead genes for Anabaena & Cyanothece
lead <- mapping[clusterID %in% diazo_clusters &
       grepl("anabaena_ppc7120|cyanothece_51142", member), ]
fwrite(lead, "results/figures/lead_genes.tsv", sep="\t")
