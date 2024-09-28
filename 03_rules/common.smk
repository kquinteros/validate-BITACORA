#loading python libraries 
from snakemake.utils import validate
import pandas as pd
import yaml
from smart_open import open

##--- validate configuration file --##
validate(config, schema="../01_config/config.schema.yaml")

##--- load and validate protein database file ---##
proteins = pd.read_csv(config["protein_database"], sep="\t", dtype=str).set_index("Protein", drop=False)
proteins.index.names = ["Index"]
validate(proteins, schema= "../01_config/protein_database.schema.yaml")


##--- load and validate genome sequences and annotations ---##
genomes = pd.read_csv(config["target_genomes"], sep="\t", dtype = str).set_index("Sample", drop=False)
genomes.index.names = ["Index"]
validate(genomes, schema = "../01_config/target_genomes.schema.yaml")


##--- helper functions ---##
#This function retrieves protein PFAM accession ID and amino acid sequence length information
def protein_pfam_input(wildcards):
    return proteins.loc[wildcards.pro, "PFAM"]


#This function retrieves protein size selection. 
def protein_size_input(wildcards):
    return {
        "min": proteins.loc[wildcards.pro, "Min"],
        "max": proteins.loc[wildcards.pro, "Max"]
  }

#This function retrieves protein size selection. 
def protein_size_min(wildcards):
    return {
        "min": proteins.loc[wildcards.pro, "Min"],
        "max": proteins.loc[wildcards.pro, "Max"]
  }

#This function returns species abbreviation   
def genome_key(wildcards):
    return {
        genomes.loc[wildcards.sample, "Key"]
  }

##--- Output functions for bitacora pipeline ---##
# Rule: rule-run-interproscan
def pfam_output(wildcards):
    return expand("05_output/{sam}/{pro}/interproscan/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered.fasta.tsv", sam=genomes['Sample'], pro=proteins["Protein"])

# Rule: run-filter-protein-domain
def filter_protein_output(wildcards):
    return expand("05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_R1.fasta", sam=genomes['Sample'], pro=proteins["Protein"])

# Rule: run_TMbed
def run_tmbed(wildcards):
    return expand("05_output/{sam}/{pro}/TMbed/{pro}_transmembrane.pred", sam=genomes['Sample'], pro=proteins["Protein"])

# Rule: run_number_of_trans_domains
def run_awk_counts(wildcards):
    return expand("05_output/{sam}/{pro}/TMbed/{pro}_transmembrane_counts.csv", sam=genomes['Sample'], pro=proteins["Protein"])

# Rule: run_filter_trans_domains
def filter_trans_domains_output(wildcards):
    return expand("05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_R3.target.fasta", sam=genomes['Sample'], pro=proteins["Protein"])

# Rule: exclude_sequences
def exclude_sequences_output(wildcards):
    return expand("05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_excluded.gff3", sam=genomes['Sample'], pro=proteins["Protein"])

# Rule: ensure_feature_names
def ensure_feature_names_output(wildcards):
    return expand("05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_final.cds.fasta", sam=genomes['Sample'], pro=proteins["Protein"])

# Rule: run_filter_overlapping
def run_filter_overlapping_output(wildcards):
    return expand("05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_excluded_final.gff3", sam=genomes['Sample'], pro=proteins["Protein"])