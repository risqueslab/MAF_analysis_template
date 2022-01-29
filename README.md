# MAF_analysis_template
A MAF analysis template to be used with multi-sample MAF files in the Risques lab.  This is designed for use with output from the [Duplex-Seq-Pipeline](https://github.com/Kennedy-Lab-UW/Duplex-Seq-Pipeline).

## Basic use:

This template mostly uses Snakemake for environment control.  If you are using Snakemake, to run the associated R script:
* `cd` to this directory on the command line
* run `bash run_processing.sh`

If you don't want to use Snakemake, then set up your own environment with (at minimum):
* R=4.1* 
* r-tidyverse
* r-gridextra
* r-gridbase
* r-patchwork
* r-ggpubr
* r-cowplot

While most of these packages aren't used in the template, they are imported by the R script so that they are there if you want to use them.  
Then you can run the R script by:
* `cd` to this directory on the command line
* run `Rscript scripts/MAF_processing.R`

## Configuration
The text file processing_config.txt is used for configuration.  All rows currently in the file are required; more can be added as key-value pairs:
`key=value`, with no spaces around the `=`.  These values can be referenced in the R script using `inputs$key`.  

The required inputs are: 
* A multi-sample MAF file containing output data from the Duplex-Seq-Pipeline
* A tab-deliniated file containing extra data for all the samples in question.  This file MUST be present, but if there is no extra data, could just be a list of the sample names with the header **Sample** (capitalization matters)
* A minimum depth to use for analysis; this will override the *low_depth* filter from the Duplex-Seq-Pipeline.  
* A maximum clonality level to use; this can override the *SNP* filter from the Duplex-Seq-Pipeline
* A prefix for any output files generated; suggested that this start with **outputs/** to place outputs in the outputs directory.  
* A comma-delimited list of filters to apply (no spaces)
* A long-form [Seshat](https://p53.fr/tp53-database/seshat) file (for TP53 analysis only)

It is suggested that all input files go in the **inputs** folder, and that their definition in the processing_config.txt file reflect that.  
