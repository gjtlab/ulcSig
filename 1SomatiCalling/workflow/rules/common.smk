import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version
import itertools
import glob

min_version("5.18.0")


#report: "../report/workflow.rst"


# container: "continuumio/miniconda3:4.8.2"

###### Config file and sample sheets #####
configfile: "configs/config.yaml"


# validate(config, schema="../schemas/config.schema.yaml")

input_data = pd.read_table(config["input_data"], sep="\t", dtype=str).drop_duplicates()

samples = input_data["tumor_sample"].unique()

input_data = input_data.set_index(["tumor_sample"], drop=False)

##### Wildcard constraints #####
wildcard_constraints:
    sample = "|".join(samples),

##### Helper functions #####

def get_mutect_input(wildcards):
    """Get bams files of given sample."""
    bams = input_data.loc[(wildcards.sample), ["tumor_bam", "normal_bam"]].dropna()
    if len(bams) == 2:
        return [bams.tumor_bam, bams.normal_bam]
    return [bams.tumor_bam]

def get_mutect_sample(wildcards):
    """Get bams files of given sample."""
    bams = input_data.loc[(wildcards.sample), ["tumor_bam", "normal_bam"]].dropna()
    if len(bams) == 2:
        return "--tumor-sample " + bams.tumor_bam + " --normal-sample " + bams.normal_bam
    return ""
