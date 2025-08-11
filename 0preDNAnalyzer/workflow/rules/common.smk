import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version
import itertools
import glob

min_version("7.0")

###### Config file and sample sheets #####
configfile: "configs/config.yaml"


# validate(config, schema="../schemas/config.schema.yaml")

input_data = pd.read_table(config["input_data"], sep="\t", dtype=str).drop_duplicates()

# 确保 'unit' 列存在，如果不存在或全为空则创建一个默认值
if 'unit' not in input_data.columns or input_data['unit'].isnull().all():
    input_data['unit'] = [f"_u{x}" for x in input_data.groupby('sample').cumcount().add(1).astype(str)]
else:
    input_data['unit'] = "_" + input_data['unit'].astype(str)

input_data.loc[:, "count"] = input_data.groupby('sample')["sample"].transform('count')
input_data.loc[:, "fastp_ID"] = input_data['sample'] + input_data['unit']

samples = input_data["sample"].unique()
units = input_data["unit"].unique()
fastp_IDs = input_data["fastp_ID"].unique()

input_data = input_data.set_index(["sample", "unit"], drop=False)

def get_fastp_output(x):
    if x["layout"] == "Paired":
        ext = ["_R1.cleaned.fq.gz", "_R2.cleaned.fq.gz", "_pe.fastp.html", "_pe.fastp.json"]
    else:
        ext = [".cleaned.fq.gz", "_se.fastp.html", "_se.fastp.json"]
    return(expand("results/fastps/{sample}{unit}{ext}",
            sample = x["sample"],
            unit = x["unit"],
            ext = ext))

fastp_output = list(itertools.chain.from_iterable(input_data.apply(get_fastp_output, axis=1)))

##### Wildcard constraints #####
wildcard_constraints:
    sample = "|".join(samples),
    unit = "|".join(units),
    fastp_ID = "|".join(fastp_IDs),


##### Helper functions #####

# contigs in reference genome
def get_contigs():
    with checkpoints.genome_faidx.get().output[0].open() as fai:
        return pd.read_table(fai, header=None, usecols=[0], squeeze=True, dtype=str)


def get_fastq(wildcards):
    """Get fastq files of given sample_unit."""
    fastqs = input_data.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return [fastqs.fq1, fastqs.fq2]
    return [fastqs.fq1]


def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(input_data.loc[(sample, unit), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}\tLB:{library}'".format(
        sample=wildcards.sample,
        platform=input_data.loc[(wildcards.sample, wildcards.unit), "platform"],
        library=input_data.loc[(wildcards.sample, wildcards.unit), "library"],
    )


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if is_single_end(wildcards.sample, wildcards.unit):
        return "results/fastps/{sample}{unit}.cleaned.fq.gz".format(
            sample=wildcards.sample, unit=wildcards.unit
        )
    else:
        return [
            "results/fastps/{sample}{unit}_R1.cleaned.fq.gz".format(
                sample=wildcards.sample, unit=wildcards.unit
            ),
            "results/fastps/{sample}{unit}_R2.cleaned.fq.gz".format(
                sample=wildcards.sample, unit=wildcards.unit
            ),
        ]


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    units = input_data.loc[wildcards.sample]["unit"].unique()
    return expand(
        "results/mapped/{sample}{unit}.mkdup.bam",
        sample=wildcards.sample,
        unit=units,
    )


def get_regions_param(regions=config["resources"].get("interval_region"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["resources"].get("interval_padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


def get_recal_input(bai=False):
    """Get the input file for BQSR"""
    if bai:
        return "results/mapped/{sample}.mkdup.bam.bai"
    return "results/mapped/{sample}.mkdup.bam"


def get_fastp_params(wildcards):
    """Get fastp parameters."""
    extra = config["params"]["fastp"].get("extra", "")
    length_required = config["params"]["fastp"].get("length_required")
    if length_required:
        extra += f" --length_required {length_required}"
    return extra

