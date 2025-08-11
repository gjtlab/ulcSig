rule bwa_mem2_mem:
    input:
        reads = get_trimmed_reads,
        idx = multiext(config["resources"]["genome"], 
                      ".0123", ".amb", ".bwt.2bit.64", ".ann", ".pac")
    output:
        bam = temp("results/mapped/{sample}{unit}.sorted.bam"),
        index = temp("results/mapped/{sample}{unit}.sorted.bam.bai")
    log:
        "logs/mapping/{sample}{unit}.bwa_mem2.log"
    params:
        extra = lambda wildcards: config["params"]["bwaMem2"]["mem"] + " " + get_read_group(wildcards),
        sort_extra = "--tmpdir {}".format(config["params"]["tmpdir"])
    threads: min(workflow.cores, 16)
    wrapper:
        "file:wrappers/bwa-mem2/mem-sambamba"


rule sambamba_markdup:
    input:
        "results/mapped/{sample}{unit}.sorted.bam"
    output:
        bam = temp("results/mapped/{sample}{unit}.mkdup.bam"),
        index = temp("results/mapped/{sample}{unit}.mkdup.bam.bai")
    params:
        extra = config["params"]["sambamba"]["markdup_extra"] + " --tmpdir {}".format(config["params"]["tmpdir"])
    log:
        "logs/mapping/{sample}{unit}.markdup.log"
    threads: min(workflow.cores, 4)
    wrapper:
        "file:wrappers/sambamba/markdup"


rule sambamba_merge:
    input:
        get_sample_bams
    output:
        bam = temp("results/mapped/{sample}.mkdup.bam")
    params:
        extra = ""
    log:
        "logs/mapping/{sample}.merge.log"
    threads: min(workflow.cores, 8)
    wrapper:
        "file:wrappers/sambamba/merge"


rule recalibrate_base_qualities:
    input:
        bam = get_recal_input(),
        bai = get_recal_input(bai=True),
        ref = config["resources"]["genome"],
        known = config["resources"]["known_vcfs"]
    output:
        recal_table = temp("results/mapped/{sample}.grp")
    log:
        "logs/mapping/{sample}.bqsr.recal.log"
    params:
        extra = get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
        java_opts = "",
        spark_runner = "LOCAL" if config["params"]["gatk"]["use_spark"] else None,
        spark_master = "local[{threads}]" if config["params"]["gatk"]["use_spark"] else None
    resources:
        mem_mb = config["params"]["gatk"]["mem"]
    wrapper:
        "file:wrappers/gatk/baserecalibratorspark" if config["params"]["gatk"]["use_spark"] else "file:wrappers/gatk/baserecalibrator"


rule apply_base_quality_recalibration:
    input:
        bam = get_recal_input(),
        bai = get_recal_input(bai=True),
        ref = config["resources"]["genome"],
        recal_table = "results/mapped/{sample}.grp"
    output:
        bam = protected("results/mapped/{sample}.bqsr.bam")
    log:
        "logs/mapping/{sample}.bqsr.apply.log"
    params:
        extra = get_regions_param(),
        java_opts = "",
        spark_runner = "LOCAL" if config["params"]["gatk"]["use_spark"] else None,
        spark_master = "local[{threads}]" if config["params"]["gatk"]["use_spark"] else None
    resources:
        mem_mb = config["params"]["gatk"]["mem"]
    wrapper:
        "file:wrappers/gatk/applybqsrspark" if config.get("use_spark", False) else "file:wrappers/gatk/applybqsr"


rule sambamba_index:
    input:
        "{prefix}.bam"
    output:
        temp("{prefix}.bam.bai")
    params:
        extra = ""
    log:
        "logs/mapping/{prefix}.index.log"
    threads: 8
    wrapper:
        "file:wrappers/sambamba/index"


rule clean:
    output:
        temp(directory("results/mapped/tmps"))
    log:
        "logs/mapping/clean.log"
    shell:
        "rm -rf {output}"

