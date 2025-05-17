import yaml

# Load species config (for reference if you need it in custom rules)
species_cfg = yaml.safe_load(open("config/species.yaml"))

# A little “sentinel” file to mark that all raw FASTAs have been fetched
DOWNLOAD_SENTINEL = "data/raw/.download_complete"

rule all:
    input:
        # final outputs
        "results/figures/enrichment_bar.pdf",
        "results/figures/enrichment_results.tsv",
        "results/figures/lead_genes.tsv"

rule download_proteomes:
    """
    Fetch all of the FASTA proteomes as configured in config/species.yaml
    and then touch a sentinel file so downstream rules know they can proceed.
    """
    output:
        DOWNLOAD_SENTINEL
    shell:
        """
        scripts/01_download_proteomes.py
        touch {output}
        """

rule build_blast_db:
    """
    Concatenate all downloaded FASTAs and build a BLAST protein DB.
    """
    input:
        DOWNLOAD_SENTINEL
    output:
        fasta="data/db/all_cyano.faa",
        phr="data/db/cyano_db.phr",
        pin="data/db/cyano_db.pin",
        psq="data/db/cyano_db.psq"
    shell:
        "scripts/02_build_blast_db.sh"

rule run_blastp:
    """
    Perform an all-vs-all BLASTP against the custom cyano_db.
    """
    input:
        fasta="data/db/all_cyano.faa",
        db="data/db/cyano_db.phr"
    output:
        "results/blastp/blastp_all.out"
    threads: 16
    shell:
        "scripts/03_run_blastp.sh"

rule cluster_orthologs:
    """
    Filter the BLAST hits and cluster with MCL into ortholog groups.
    """
    input:
        "results/blastp/blastp_all.out"
    output:
        clusters="results/clusters/clusters.txt",
        mapping="results/clusters/clusters_mapping.tsv"
    shell:
        "scripts/04_cluster_orthologs.R"

rule build_matrix:
    """
    Build presence/absence matrix of clusters × species.
    """
    input:
        mapping="results/clusters/clusters_mapping.tsv"
    output:
        pres_abs="results/matrices/pres_abs.csv",
        species_groups="results/matrices/species_groups.csv"
    shell:
        "scripts/05_build_matrix.py"

rule annotate_interpro:
    """
    Run InterProScan on the full proteome to pull in domain and GO annotations.
    """
    input:
        fasta="data/db/all_cyano.faa"
    output:
        interpro="results/annotations/all_interpro.tsv"
    threads: 16
    shell:
        "scripts/06_annotate_interpro.sh"

rule enrichment_analysis:
    """
    Identify diazotroph-specific clusters, run enrichment,
    and extract lead genes for Anabaena PCC 7120 & Cyanothece 51142.
    """
    input:
        pres_abs="results/matrices/pres_abs.csv",
        mapping="results/clusters/clusters_mapping.tsv",
        interpro="results/annotations/all_interpro.tsv",
        # config is read inside the R script, so we don’t need to list it as input
    output:
        pdf="results/figures/enrichment_bar.pdf",
        results="results/figures/enrichment_results.tsv",
        leads="results/figures/lead_genes.tsv"
    shell:
        "scripts/07_enrichment_analysis.R"
