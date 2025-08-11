if not config["params"]["fastp"]["skip"]:
    # 执行fastp的代码逻辑
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
            adapters=config["params"]["fastp"].get("adapters", ""),
            extra=get_fastp_params
        threads: min(workflow.cores, 16)
        wrapper:
            "file:wrappers/fastp"
    
    rule fastp_pe:
        input:
            sample=get_fastq
        output:
            trimmed=["results/fastps/{sample}{unit}_R1.cleaned.fq.gz",
                    "results/fastps/{sample}{unit}_R2.cleaned.fq.gz"],
            html="results/fastps/{sample}{unit}_pe.fastp.html",
            json="results/fastps/{sample}{unit}_pe.fastp.json"
        log:
            "logs/fastps/{sample}{unit}.log"
        params:
            adapters=config["params"]["fastp"].get("adapters", ""),
            extra=get_fastp_params
        threads: min(workflow.cores, 16)
        wrapper:
            "file:wrappers/fastp"
else:
    rule fastp_se:
        input:
            sample=get_fastq
        output:
            trimmed="results/fastps/{sample}{unit}.cleaned.fq.gz",
            html="results/fastps/{sample}{unit}_se.fastp.html",
            json="results/fastps/{sample}{unit}_se.fastp.json"
        log:
            "logs/fastps/{sample}{unit}.log"
        shell:
            """
            # 跳过fastp处理，直接复制原始文件
            cp {input.sample} {output.trimmed}
            touch {output.html} {output.json}
            """
    
    rule fastp_pe:
        input:
            sample=get_fastq
        output:
            trimmed=["results/fastps/{sample}{unit}_R1.cleaned.fq.gz",
                    "results/fastps/{sample}{unit}_R2.cleaned.fq.gz"],
            html="results/fastps/{sample}{unit}_pe.fastp.html",
            json="results/fastps/{sample}{unit}_pe.fastp.json"
        log:
            "logs/fastps/{sample}{unit}.log"
        shell:
            """
            cp {input.sample[0]} {output.trimmed[0]}
            cp {input.sample[1]} {output.trimmed[1]}
            touch {output.html} {output.json}
            """


rule multiqc:
    input:
        expand(
            "results/fastps/{sample}{unit}_{layout}.fastp.json",
            zip,
            sample=input_data["sample"],
            unit=input_data["unit"],
            layout=["pe" if x == "Paired" else "se" for x in input_data["layout"]]
        )
    output:
        report(
            "results/fastps_multiqc.html",
            caption="../report/multiqc.rst",
            category="Quality control",
        )
    log:
        "logs/fastps_multiqc.log"
    wrapper:
        "file:wrappers/multiqc"