rule gatk_mutect2:
    input:
        map = get_mutect_input,
        fasta = config["resources"]["genome"],
        germline = config["resources"]["snpdb"],
        pon = config["resources"]["pon"],
        intervals = config["resources"]["intervals"]
    output:
        vcf = "results/{sample}.mutect.vcf.gz",
    threads: 10
    log:
        "logs/{sample}.mutect.log",
    params:
        extra = get_mutect_sample,
        java_opts = ""
    resources:
        mem_mb = config["params"]["gatk"]["mem"],
    wrapper:
        "file:wrappers/gatk/mutect"


rule gatk_filtermutectcalls:
    input:
        vcf = "results/{sample}.mutect.vcf.gz",
        ref = config["resources"]["genome"],
    output:
        vcf = "results/{sample}.mutect.filtered.vcf.gz",
    log:
        "logs/{sample}.mutect.filtered.log",
    params:
        extra="",  # optional arguments, see GATK docs
        java_opts="",  # optional
    resources:
        mem_mb = config["params"]["gatk"]["mem"],
    wrapper:
        "file:wrappers/gatk/filtermutectcalls"


rule gatk_leftalignandtrimvariants:
    input:
        vcf = "results/{sample}.mutect.filtered.vcf.gz",
        ref = config["resources"]["genome"],
    output:
        vcf = "results/{sample}.mutect.filtered.split.norm.vcf.gz",
    log:
        "logs/{sample}.mutect.filtered.split.norm.log",
    params:
        extra="--split-multi-allelics", # optional
        java_opts="",  # optional
    resources:
        mem_mb = config["params"]["gatk"]["mem"],
    wrapper:
        "file:wrappers/gatk/leftalignandtrimvariants"