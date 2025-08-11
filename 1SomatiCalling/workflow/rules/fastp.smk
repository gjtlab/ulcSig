rule fastp_se:
    input:
        sample=get_fastq
    output:
        trimmed="results/fastps/{sample}{unit}.cleaned.fq.gz",
        html="results/fastps/{sample}{unit}_se.fastp.html",
        json="results/fastps/{sample}{unit}_se.fastp.json"
    log:
        "logs/fastps/{sample}{unit}.log"
    params:
        adapters="",
        extra=""
    threads: 16
    wrapper:
        "file:wrappers/fastp"


rule fastp_pe:
    input:
        sample=get_fastq,
    output:
        trimmed=["results/fastps/{sample}{unit}_R1.cleaned.fq.gz",
                "results/fastps/{sample}{unit}_R2.cleaned.fq.gz"],
        html="results/fastps/{sample}{unit}_pe.fastp.html",
        json="results/fastps/{sample}{unit}_pe.fastp.json"
    log:
        "logs/fastps/{sample}{unit}.log"
    params:
        adapters="",
        extra=""
    threads: 16
    wrapper:
        "file:wrappers/fastp"

# rule fastp:
#     input:
#         sample=get_fastq,
#     output:
#         trimmed = fastp_output_dict[("{sample}", "{unit}")]['fastqs'],
#         html="results/fastps/{sample}{unit}.fastp.html",
#         json="results/fastps/{sample}{unit}.fastp.json"
#         # trimmed=list(filter(lambda x: "clean" in x and "{sample}{unit}" in x, fastp_output)),
#         # html=list(filter(lambda x: "html" in x and "{sample}{unit}" in x, fastp_output)),
#         # json=list(filter(lambda x: "json" in x and "{sample}{unit}" in x, fastp_output)),
#     log:
#         "logs/fastps/{sample}{unit}.log"
#     params:
#         adapters="",
#         extra=""
#     threads: 16
#     wrapper:
#         "file:wrappers/fastp"

rule multiqc:
    input:
        # list(itertools.chain.from_iterable([x['json'] for x in fastp_output_dict.values()]))
        [x for x in fastp_output if "json" in x]
    output:
        report(
            "results/fastps_multiqc.html",
            caption="../report/multiqc.rst",
            category="Quality control",
        ),
    log:
        "logs/multiqc.log",
    wrapper:
        "file:wrappers/multiqc"
