import yaml

# Load species config (for reference if you need it in custom rules)
species_cfg = yaml.safe_load(open("species.yaml"))

# A little “sentinel” file to mark that all raw FASTAs have been fetched
DOWNLOAD_SENTINEL = "data/raw/.download_complete"

# Allow overriding BLAST filtering thresholds via `--config pident=XX evalue=YY`
PIDENT = config.get("pident", 80)
EVALUE = config.get("evalue", 1e-10)

rule all:
    input:
        "results/cohorts/unicellular_specific.tsv",
        "results/cohorts/filamentous_specific.tsv",
        "results/cohorts/diazotroph_common.tsv",
        "results/cohorts/unicellular_specific.csv",
        "results/cohorts/filamentous_specific.csv",
        "results/cohorts/diazotroph_common.csv"

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

rule filter_blast_hits:
    """
    Filter BLAST hits by percent identity and e-value.
    """
    input:
        "results/blastp/blastp_all.out"
    output:
        "results/blastp/blastp_filtered.out"
    params:
        pident=PIDENT,
        evalue=EVALUE
    shell:
        "scripts/filter_blast_hits.py {input} {output} --pident {params.pident} --evalue {params.evalue}"

rule identify_cohorts:
    """Summarise conserved protein cohorts across species groups."""
    input:
        "results/blastp/blastp_filtered.out"
    output:
        uni_tsv="results/cohorts/unicellular_specific.tsv",
        fil_tsv="results/cohorts/filamentous_specific.tsv",
        common_tsv="results/cohorts/diazotroph_common.tsv",
        uni_csv="results/cohorts/unicellular_specific.csv",
        fil_csv="results/cohorts/filamentous_specific.csv",
        common_csv="results/cohorts/diazotroph_common.csv"
    shell:
        "scripts/identify_cohorts.py {input} results/cohorts"

