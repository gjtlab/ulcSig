rule bwa_mem2_mem:
    input:
        reads = get_trimmed_reads,
        # Index can be a list of (all) files created by bwa-mem2, or one of them
        # idx = rules.bwa_index.output,
        idx = multiext(config["resources"]["genome"], ".0123", ".amb", ".bwt.2bit.64", ".ann",".pac"),
    output:
        bam = temp("results/fq2bams/{sample}{unit}.sorted.bam"),
        index = temp("results/fq2bams/{sample}{unit}.sorted.bam.bai")
    log:
        "logs/fq2bams/{sample}{unit}.bwa_mem2.log",
    params:
        extra = get_read_group,
        sort_extra = "--tmpdir results/fq2bams/tmps",  # Extra args for sambamba.
    threads: 32
    wrapper:
        "file:wrappers/bwa-mem2/mem-sambamba"


rule sambamba_markdup:
    input:
        "results/fq2bams/{sample}{unit}.sorted.bam",
    output:
        bam = temp("results/fq2bams/{sample}{unit}.mkdup.bam"),
        #bam = protected("results/fq2bams/{sample}{unit}.mkdup.bam"),
        index = temp("results/fq2bams/{sample}{unit}.mkdup.bam.bai")
    params:
        extra = "--hash-table-size 2621440 --sort-buffer-size 102400 --overflow-list-size 600000" # optional parameters
    log:
        "logs/fq2bams/{sample}{unit}.mkdup.log"
    threads: 8
    wrapper:
        "file:wrappers/sambamba/markdup"


rule sambamba_merge:
    input:
        get_sample_bams
    output:
        bam = temp("results/fq2bams/{sample}.mkdup.bam"),
        #bam = protected("results/fq2bams/{sample}.mkdup.bam")
    params:
        extra = ""  # optional parameters
    log:
        "logs/sambamba-merge/{sample}.log"
    threads: 8
    wrapper:
        "file:wrappers/sambamba/merge"


rule recalibrate_base_qualities:
    input:
        bam = get_recal_input(),
        bai = get_recal_input(bai=True),
        ref = config["resources"]["genome"],
        known = config["resources"]["known_vcfs"],
    output:
        recal_table = temp("results/fq2bams/{sample}.grp"),
    log:
        "logs/fq2bams/{sample}.bqsr1.log",
    params:
        extra = get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
        java_opts = ""
    resources:
        mem_mb = config["params"]["gatk"]["mem"],
    wrapper:
        "file:wrappers/gatk/baserecalibrator"


rule apply_base_quality_recalibration:
    input:
        bam = get_recal_input(),
        bai = get_recal_input(bai=True),
        ref = config["resources"]["genome"],
        recal_table = "results/fq2bams/{sample}.grp",
    output:
        bam = protected("results/fq2bams/{sample}.bqsr.bam"),
    log:
        "logs/fq2bams/{sample}.bqsr2.log",
    params:
        extra = get_regions_param(),
    resources:
        mem_mb = config["params"]["gatk"]["mem"],
    wrapper:
        "file:wrappers/gatk/applybqsr"


rule sambamba_index:
    input:
        "{prefix}.bam"
    output:
        temp("{prefix}.bam.bai")
    params:
        extra = ""  # optional parameters
    log:
        "logs/sambamba-index/{prefix}.log"
    threads: 8
    wrapper:
        "file:wrappers/sambamba/index"


# rule clean:
#     output:
#         [x for x in glob.glob("results/fq2bams/*.bai") if "bqsr" not in x],
#         glob.glob("results/fq2bams/*.grp"),
#         "results/fq2bams/tmps",
#     log:
#         "logs/clean.log"
#     shell:
#         "echo {output}"
#         "rm -rf {output}"

