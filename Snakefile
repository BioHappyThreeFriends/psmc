from snakemake.utils import min_version
from pathlib import Path
from os import walk

##### set minimum snakemake version #####
min_version("5.4.0")

##### setup config #####
configfile: "config/config.yaml"

# output dirs
alignment_dir_path = Path(config["alignment_dir"])
assembly_stats_dir_path = Path(config["assembly_stats_dir"])
variant_calling_dir_path = Path(config["variant_calling_dir"])
fastq_dir_path = Path(config["fastq_dir"])
psmc_dir_path = Path(config["psmc_dir"])
plots_dir_path = Path(config["plots_dir"])
# technical dirs
scripts_dir_path = str(config["scripts_dir"])
psmc_tool_dir_path = str(config["psmc_tool_dir"])
# log dirs
log_dir_path = Path(config["log_dir"])
cluster_log_dir_path = Path(config["cluster_log_dir"])
benchmark_dir_path = Path(config["benchmark_dir"])


#### load rules #####
include: "workflow/rules/common.smk"

include: "workflow/rules/Preparation/alignment.smk"
include: "workflow/rules/Preparation/variant_calling.smk"
include: "workflow/rules/PSMC/psmc.smk"
include: "workflow/rules/PSMC/drawing.smk"


##### target rules #####
localrules: all, create_sample_cluster_log_dirs

rule all:
    input:
        # create dirs for samples to avoid slurm error
        expand(cluster_log_dir_path / "{sample}", sample=SAMPLES.sample_id),

        # psmc plots
        expand(plots_dir_path / "{sample}.eps", sample=SAMPLES.sample_id)

