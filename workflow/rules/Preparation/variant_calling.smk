rule variant_calling: #TODO
    input:
        fasta=get_fasta,
        fai=get_fasta_fai,
        bam=get_bam
    output:
        bcf_gz=variant_calling_dir_path / "{sample}.bcf.gz",
        vcf_gz=variant_calling_dir_path / "{sample}.vcf.gz"
    params:
        prefix=lambda w: variant_calling_dir_path / ("split/bcf/{sample}/".format(assembly = w.assembly, sample=w.sample)),
        bcf=lambda wildcards, output: output["bcf_gz"][:-3],
        vcf=lambda wildcards, output: output["vcf_gz"][:-3],
        pigz_threads=config["pigz_threads"] 
    log:
        std=log_dir_path / "{sample}.variant_calling.log",
        cluster_log=cluster_log_dir_path / "{sample}.variant_calling.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.variant_calling.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample}.variant_calling.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["variant_calling_threads"],
        time=config["variant_calling_time"],
        mem=config["variant_calling_mem_mb"],
    threads: config["variant_calling_threads"]
    shell:
        "mkdir -p {params.prefix}; "
        "prepare_region_list.py -r {input.fai} -s -m 1500000 -n 1 -g samtools -x 1000 2>/dev/null | "
        "parallel -j 10 'samtools mpileup -C 50 -uf '{input.fasta}' -r {} '{input.bam}' | "
        "bcftools view -b -c - > {params.prefix}/tmp.{{#}}.bcf' > {log.std} 2>&1; "
        "bcftools cat `ls {params.prefix}/tmp.*.bcf | sort -V` >> {output.bcf}; "
        "bcftools view {input.bcf} >> {output.vcf}; "
        "pigz -p {params.pigz_threads} {output.bcf}; "
        "pigz -p {params.pigz_threads} {output.vcf}; "


rule create_fastq: #TODO
    input:
        vcf_gz=rules.variant_calling.output.vcf_gz
    output:
        fastq_gz=fastq_dir_path / "{sample}.fq.gz"
    params:
        d=config["d"],
        D=config["D"]
    log:
        std=log_dir_path / "{sample}.create_fastq.log",
        cluster_log=cluster_log_dir_path / "{sample}.create_fastq.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.create_fastq.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample}.create_fastq.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["create_fastq_threads"],
        time=config["create_fastq_time"],
        mem=config["create_fastq_mem_mb"],
    threads: config["create_fastq_threads"]
    shell:
        "zcat {input.vcf_gz} | vcfutils.pl vcf2fq -d {params.d} -D {params.D} | "
        "gzip > {output.fastq_gz}; "


