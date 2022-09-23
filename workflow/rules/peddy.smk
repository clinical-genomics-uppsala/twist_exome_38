rule glnexus:
    input:
        gvcf=expand("parabricks/pbrun_deepvariant/{sample}.g.vcf", sample=get_samples(samples)),
    output:
        bcf=temp("qc/peddy/all.bcf"),
        glnexus=temp(directory("GLnexus.DB"))
    params:
        bed=config["glnexus"]["bedfile"],
        in_gvcf=get_in_gvcf
    log:
        "qc/peddy/all.bcf.log",
    benchmark:
        repeat(
            "qc/peddy/all.bcf.benchmark.tsv",
            config.get("glnexus", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("glnexus", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("glnexus", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("glnexus", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("glnexus", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("glnexus", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("glnexus", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("glnexus", {}).get("container", config["default_container"])
    message:
        "{rule}: Run GLNexus for joint genotyping of Deepvariant gVCFs"
    shell:
        "glnexus_cli --config DeepVariant --bed {params.bed} -i {params.in_gvcf} > {output.bcf}"


rule bcftools_view:
    input:
        "qc/peddy/all.bcf",
    output:
        "qc/peddy/all.vcf.gz",
    log:
        "qc/peddy/all.vcf.gz.log"
    benchmark:
        repeat(
            "qc/peddy/all.bcf.benchmark.tsv",
            config.get("bcftools_view", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=config.get("bcftools_view", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_view", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcftools_view", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcftools_view", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcftools_view", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools_view", {}).get("container", config["default_container"])
    message:
        "{rule}: Run bcftools view to convert glNexus bcf to vcf and then bgzip and tabix"
    shell:
        """
        bcftools view {input} | bgzip -c > {output}
        tabix {output}
        """


rule create_ped:
    input:
        config["peddy"]["samples"],
    output:
        'qc/peddy/all.ped'
    log:
        'qc/peddy/all.ped'
    resources:
        mem_mb=config.get("create_ped", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_ped", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("create_ped", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("create_ped", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("create_ped", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("create_ped", {}).get("container", config["default_container"])
    message:
        "{rule}: Create a peddy ped/FAM file from the SampleSheet.csv file"
    script:
        "../scripts/create_peddy_fam.py"


rule peddy:
    input:
        vcf="qc/peddy/all.vcf.gz",
        ped="qc/peddy/all.ped",
    output:
        temp("qc/peddy/peddy_mqc.peddy.ped"),
        temp("qc/peddy/peddy_mqc.ped_check.csv"),
        temp("qc/peddy/peddy_mqc.sex_check.csv"),
        temp("qc/peddy/peddy_mqc.het_check.csv"),
        temp("qc/peddy/peddy_mqc.html")
    params:
        pre="qc/peddy/peddy_mqc",
    log:
        "qc/peddy/peddy.log",
    benchmark:
        repeat(
            "qc/peddy/all.peddy.ped.benchmark.tsv",
            config.get("peddy", {}).get("benchmark_repeats", 1),)
    resources:
        mem_mb=config.get("peddy", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("peddy", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("peddy", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("peddy", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("peddy", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("peddy", {}).get("container", config["default_container"])
    message:
        "{rule}: Run peddy analysis to check relatedeness in trios and check the sex of samples"
    shell:
        "peddy --procs {resources.threads} --plot --loglevel 'WARNING' --prefix {params.pre} {input.vcf} {input.ped}"


rule create_peddy_mqc_tsv:
    input:
        peddy_rel_check = "qc/peddy/peddy_mqc.ped_check.csv",
        peddy_sex_check = "qc/peddy/peddy_mqc.sex_check.csv",
        ped="qc/peddy/all.ped"
    output:
        rel_check_mqc="qc/peddy/peddy_rel_check_mqc.tsv",
        sex_check_mqc="qc/peddy/peddy_sex_check_mqc.tsv"
    params:
        pre="qc/peddy/peddy_mqc",
    log:
        "qc/peddy/peddy.log",
    benchmark:
        repeat(
            "qc/peddy/create_peddy_mqc_tsv.benchmark.tsv",
            config.get("create_peddy_mqc_tsv", {}).get("benchmark_repeats", 1),)
    resources:
        mem_mb=config.get("create_peddy_mqc_tsv", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_peddy_mqc_tsv", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("create_peddy_mqc_tsv", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("create_peddy_mqc_tsv", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("create_peddy_mqc_tsv", {}).get("time", config["default_resources"]["time"]),
    message:
        "{rule}: Create multiqc custom content embedded config tsv files from peddy sex_check and ped_check files"
    script:
        "../scripts/create_peddy_mqc_config.py"
