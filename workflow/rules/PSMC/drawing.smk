rule draw_psmc_one_sample: #TODO
    input:
        psmc=rules.psmc.output.psmc,
        psmc_tool=dir(psmc_tool_dir_path)
    output:
        psmc=plots_dir_path / "{sample}.eps"
    params:
        prefix=lambda w: (psmc_dir_path / "{sample}").format(assembly=w.assembly, sample=w.sample),
        generation_time=config["generation_time"],
        mutation_rate=config["mutation_rate"]
    log:
        std=log_dir_path / "{sample}.draw_psmc_one_sample.log",
        cluster_log=cluster_log_dir_path / "{sample}.draw_psmc_one_sample.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.draw_psmc_one_sample.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample}.draw_psmc_one_sample.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["draw_psmc_one_sample_threads"],
        time=config["draw_psmc_one_sample_time"],
        mem=config["draw_psmc_one_sample_mem_mb"],
    threads: config["draw_psmc_one_sample_threads"]
    shell:
        "{input.psmc_tool}/utils/psmc_plot.pl -g{params.generation_time} {params.prefix} {input.psmc}; "


# rule draw_psmc_all_samples: #TODO
#     input:
#         psmc=rules.psmc.output.psmc,
#         psmc_tool=dir(psmc_tool_dir_path)
#     output:
#         psmc=psmc_dir_path / "{sample}.psmc"
#     params:
#         prefix=$OUTDIR/$ASSEMBLY_NAME.$SAMPLE
#     log:
#         std=log_dir_path / "{sample}.psmc.log",
#         cluster_log=cluster_log_dir_path / "{sample}.psmc.cluster.log",
#         cluster_err=cluster_log_dir_path / "{sample}.psmc.cluster.err"
#     benchmark:
#          benchmark_dir_path / "{sample}.psmc.benchmark.txt"
#     conda:
#         "../../../%s" % config["conda_config"]
#     resources:
#         cpus=config["psmc_threads"],
#         time=config["psmc_time"],
#         mem=config["psmc_mem_mb"],
#     threads: config["psmc_threads"]
#     shell:
#         "{input.psmc_tool}/utils/psmc_plot.pl -g10 -M "female,male" $OUTDIR/$ASSEMBLY_NAME.pusa_sibirica $OUTDIR/$ASSEMBLY_NAME.female.diploid.psmc $OUTDIR/$ASSEMBLY_NAME.male.diploid.psmc; "


