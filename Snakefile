import sys

rule all:
    input:
        ".maf.done",
        ".tp53.done"
    output:
        temp(touch(".all.done"))

rule MAF_processing:
    input:
        in_config = "processing_config.txt"
    output:
        temp(touch(".maf.done"))
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

rule TP53_processing:
    input:
        in_config = "processing_config.txt"
    output:
        temp(touch(".tp53.done"))
    conda:
        "envs/R_env.yaml"
    shell:
        """
        set -e
        set -o pipefail
        set -x

        Rscript {sys.path[0]}/scripts/TP53_processing.R \
        processing_config.txt
        """