rule mosdepth: #TODO
    input:
        bam=get_bam,
        bai=get_bai
    output:
        out=alignment_dir_path / "{sample}.coverage.per-base.bed.gz"
    params:
        min_mapping_quality=config["min_mapping_quality"],
        prefix=lambda w: alignment_dir_path / ("{sample}.coverage".format(assembly = w.assembly, sample=w.sample))
    log:
        std=log_dir_path / "{sample}.mosdepth.log",
        cluster_log=cluster_log_dir_path / "{sample}.mosdepth.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.mosdepth.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample}.mosdepth.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["mosdepth_threads"],
        time=config["mosdepth_time"],
        mem=config["mosdepth_mem_mb"],
    threads: config["mosdepth_threads"]
    shell:
        "mosdepth -t {threads} --mapq {params.min_mapping_quality} {params.prefix} {input.bam} > {log.std} 2>&1"


rule samtools_faidx:
    input:
        assembly=get_fasta
    output:
        assembly_stats_dir_path / "{sample}.fai"
    log:
        std=log_dir_path / "{sample}/samtools_faidx.log",
        cluster_log=cluster_log_dir_path / "{sample}/samtools_faidx.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/samtools_faidx.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/samtools_faidx.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["samtools_faidx_threads"],
        mem=config["samtools_faidx_mem_mb"],
        time=config["samtools_faidx_time"]
    threads:
        config["samtools_faidx_threads"]
    shell:
        "samtools faidx {input.assembly}; "
        "mv {input.assembly}.fai {output} "


