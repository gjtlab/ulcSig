__author__ = "Christopher Schröder"
__copyright__ = "Copyright 2020, Christopher Schröder"
__email__ = "christopher.schroeder@tu-dortmund.de"
__license__ = "MIT"


from os import path

from snakemake.shell import shell


# Extract arguments.
extra = snakemake.params.get("extra", "")
sort_extra = snakemake.params.get("sort_extra", "")

index = snakemake.input.get("index", "")
if isinstance(index, str):
    index = path.splitext(snakemake.input.idx)[0]
else:
    index = path.splitext(snakemake.input.idx[0])[0]

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Check inputs/arguments.
if not isinstance(snakemake.input.reads, str) and len(snakemake.input.reads) not in {
    1,
    2,
}:
    raise ValueError("input must have 1 (single-end) or 2 (paired-end) elements")

# 根据总线程数分配各步骤的线程数
bwa_threads = snakemake.threads
view_threads = max(2, snakemake.threads // 4)  # view操作分配较少线程
sort_threads = max(4, snakemake.threads // 2)  # sort操作分配中等线程

shell(
    "(bwa-mem2 mem"
    " -t {bwa_threads}"  # bwa-mem2使用全部线程
    " {extra}"
    " {index}"
    " {snakemake.input.reads}"
    " | sambamba view -S -f bam /dev/stdin"
    " -t {view_threads}"  # view使用较少线程
    " | sambamba sort /dev/stdin"
    " -t {sort_threads}"  # sort使用中等线程
    " -o {snakemake.output.bam}"
    " {sort_extra}"
    ") {log}"
)
