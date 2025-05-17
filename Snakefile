# Snakefile
import yaml
species = yaml.safe_load(open("config/species.yaml"))

rule all:
  input:
    expand("data/raw/{grp}/{name}.faa", grp=species, name=lambda w: w),
    "results/figures/presence_absence_heatmap.pdf"

rule download_proteome:
  output: "data/raw/{grp}/{name}.faa"
  params: url=lambda wildcards: species[wildcards.grp][wildcards.name]["url"]
  shell: """
    mkdir -p data/raw/{wildcards.grp}
    wget -O {output}.gz {params.url}
    gunzip -c {output}.gz > {output}
  """

rule make_blast_db:
  input: expand("data/raw/{grp}/{name}.faa", grp=species, name=lambda w:w)
  output: "data/db/cyano_db.phr"
  shell: """
    cat data/raw/*/*.faa > data/db/all_cyano.faa
    makeblastdb -in data/db/all_cyano.faa -dbtype prot -out data/db/cyano_db
  """

rule blastp_all:
  input:
    db="data/db/cyano_db",
    query="data/db/all_cyano.faa"
  output: "results/blastp/blastp_all.out"
  threads: 16
  shell: """
    mkdir -p results/blastp
    blastp -db {input.db} \
           -query {input.query} \
           -evalue 1e-5 \
           -num_threads {threads} \
           -outfmt 6 \
           -out {output}
  """


