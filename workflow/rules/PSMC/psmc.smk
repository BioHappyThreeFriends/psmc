rule fq_to_psmcfa: #TODO
    input:
        fastq_gz=rules.create_fastq.output.fastq_gz,
        psmc_tool=dir(psmc_tool_dir_path)
    output:
        psmcfa=psmc_dir_path / "{sample}.psmcfa"
    log:
        std=log_dir_path / "{sample}.fq_to_psmcfa.log",
        cluster_log=cluster_log_dir_path / "{sample}.fq_to_psmcfa.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.fq_to_psmcfa.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample}.fq_to_psmcfa.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["fq_to_psmcfa_threads"],
        time=config["fq_to_psmcfa_time"],
        mem=config["fq_to_psmcfa_mem_mb"],
    threads: config["fq_to_psmcfa_threads"]
    shell:
        "{input.psmc_tool}/utils/fq2psmcfa -q20 {input.fastq_gz} > {output.psmcfa}; "


# rule splitfa: #TODO
#     input:
#         psmcfa=rules.fq_to_psmcfa.output.psmcfa,
#         psmc_tool=dir(psmc_tool_dir_path)
#     output:
#         split_psmcfa=psmc_dir_path / "{sample}.split.psmcfa"
#     log:
#         std=log_dir_path / "{sample}.splitfa.log",
#         cluster_log=cluster_log_dir_path / "{sample}.splitfa.cluster.log",
#         cluster_err=cluster_log_dir_path / "{sample}.splitfa.cluster.err"
#     benchmark:
#          benchmark_dir_path / "{sample}.splitfa.benchmark.txt"
#     conda:
#         "../../../%s" % config["conda_config"]
#     resources:
#         cpus=config["splitfa_threads"],
#         time=config["splitfa_time"],
#         mem=config["splitfa_mem_mb"],
#     threads: config["splitfa_threads"]
#     shell:
#         "{input.psmc_tool}/utils/splitfa {input.psmcfa} > {output.split_psmcfa}; "


rule psmc: #TODO
    input:
        psmcfa=rules.fq_to_psmcfa.output.psmcfa,
        psmc_tool=dir(psmc_tool_dir_path)
    output:
        psmc=psmc_dir_path / "{sample}.psmc"
    log:
        std=log_dir_path / "{sample}.psmc.log",
        cluster_log=cluster_log_dir_path / "{sample}.psmc.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.psmc.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample}.psmc.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["psmc_threads"],
        time=config["psmc_time"],
        mem=config["psmc_mem_mb"],
    threads: config["psmc_threads"]
    shell:
        '{input.psmc_tool}/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o {output.psmc} {input.psmcfa}; '


# rule psmc_seq: #TODO
#     input:
#         split_psmcfa=rules.splitfa.output.split_psmcfa,
#         psmc_tool=dir(psmc_tool_dir_path)
#     output:
#         psmc=psmc_dir_path / "{sample}.psmc"
#     log:
#         std=log_dir_path / "{sample}.psmc_seq.log",
#         cluster_log=cluster_log_dir_path / "{sample}.psmc_seq.cluster.log",
#         cluster_err=cluster_log_dir_path / "{sample}.psmc_seq.cluster.err"
#     benchmark:
#          benchmark_dir_path / "{sample}.psmc_seq.benchmark.txt"
#     conda:
#         "../../../%s" % config["conda_config"]
#     resources:
#         cpus=config["psmc_seq_threads"],
#         time=config["psmc_seq_time"],
#         mem=config["psmc_seq_mem_mb"],
#     threads: config["psmc_seq_threads"]
#     shell:
#         "seq 100 | xargs -i echo {input.psmc_tool}/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o round-{{}}.psmc {input.split_psmcfa} | sh; "


