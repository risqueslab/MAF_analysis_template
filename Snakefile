import sys

rule all:
    input:
        in_config = "processing_config.txt"
    output:
        temp(touch(".done"))
    conda:
        "envs/R_env.yaml"
    shell:
        """
        set -e
        set -o pipefail
        set -x

        Rscript {sys.path[0]}/scripts/MAF_processing.R \
        processing_config.txt
        """