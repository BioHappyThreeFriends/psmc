import pandas as pd
import re
from snakemake.utils import validate


# dataframe of info on samples #TODO deal not with only PE, but also SE
SAMPLES = pd.read_table(config["samples"], dtype=str)
# extract sample name from path of forward read if not filled
# SAMPLES['sample_id'].fillna(get_reads_names(SAMPLES['forward_read']), inplace=True)
# make column "sample" as identificator of dataframe
SAMPLES.set_index("sample_id", drop=False)
validate(SAMPLES, schema="../schemas/samples.schema.yaml")


def get_fasta(wildcards):
    return SAMPLES[SAMPLES["sample_id"]==wildcards.sample].fasta


def get_fasta_fai(wildcards):
    return SAMPLES[SAMPLES["sample_id"]==wildcards.sample].fasta_fai


def get_coverage(wildcards):
    return SAMPLES[SAMPLES["sample_id"]==wildcards.sample].coverage


def get_bam(wildcards):
    return SAMPLES[SAMPLES["sample_id"]==wildcards.sample].bam


def get_bai(wildcards):
    return SAMPLES[SAMPLES["sample_id"]==wildcards.sample].bai


rule create_sample_cluster_log_dirs:
    output:
        samples=directory(expand(cluster_log_dir_path / "{sample}", sample=SAMPLES.sample_id))
    shell:
        "mkdir -p {output.samples}; "


