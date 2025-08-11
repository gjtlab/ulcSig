__author__ = "Jan Forster"
__copyright__ = "Copyright 2021, Jan Forster"
__email__ = "j.forster@dkfz.de"
__license__ = "MIT"


import os
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if len(snakemake.input) > 1:
    shell(
        "sambamba merge {snakemake.params.extra} -t {snakemake.threads} "
        "{snakemake.output[0]} {snakemake.input} "
        "{log}"
    )
else:
    shell(
        """
        mv -v {snakemake.input[0]} {snakemake.output[0]} {log};
        if [ -f {snakemake.input[0]}.bai ]; then
            mv -v {snakemake.input[0]}.bai {snakemake.output[0]}.bai >> {snakemake.log} 2>&1;
        else
            echo "{snakemake.input[0]}.bai does not exist, skipping..." >> {snakemake.log}
        fi
        """
    )
