library("tidyverse")
library("patchwork")
library("ggpubr")
library("cowplot")

# Start by setting your working directory
#setwd()
#==================================================================================================#
# Set up values and files to use in processing
#==================================================================================================#
inputs <- read_delim("processing_config.txt", delim="=", col_names = FALSE, comment = "#", skip_empty_rows = TRUE) %>%
  pivot_wider(names_from = X1, values_from = X2) %>% 
  type_convert()

FiltersToApply = str_split(inputs$filters_to_apply, pattern = ",")[[1]]

# input multisample MAF file
in_maf_file = inputs$concat_maf_file
# File translating sample names, with (at minimum) column "Tumor_Sample_Barcode"
# reflecting the names in the MAF file
sample_data_file = inputs$extra_data_file

# minimum depth to filter to (not implemented yet)
min_depth = inputs$min_depth

# Prefix for naming output files
out_prefix = inputs$out_prefix

# Set MAF limit
MAF_lim = inputs$MAF_lim

# Set up non-standard inputs

#==================================================================================================#
# GENERAL PROCESSING CODE:
#==================================================================================================#



# set up filter string
filter_string = do.call(paste, as.list(c(FiltersToApply, sep="|")))

# set up variant classification table
variant_clasification_table = read_delim(inputs$varClassTrans, delim="\t") 
# Filter

sample_data <- read_delim(sample_data_file, delim="\t")

# colnames(maf_table)
# Load MAF file
maf_table <- read_delim(in_maf_file, delim="\t", skip=1) %>% 
  select(Hugo_Symbol,
         NCBI_Build,
         Chromosome,
         Start_Position,
         Variant_Classification,
         Variant_Type,
         Reference_Allele,
         Tumor_Seq_Allele2,
         Tumor_Sample_Barcode,
         HGVSc,
         HGVSp,
         HGVSp_Short,
         Exon_Number,
         t_depth,
         t_alt_count,
         Consequence,
         Existing_variation,
         IMPACT,
         FILTER 
         ) %>% 
  mutate(MAF=t_alt_count/t_depth) %>% 
  left_join(variant_clasification_table, by=c("Variant_Classification")) %>%
  left_join(sample_data, by=c("Tumor_Sample_Barcode"="Sample")) %>% 
  mutate(coding2 = if_else(
    Variant_Classification == "Splice_Region", 
    if_else(is.na(Exon_Number), 
            "non-coding", 
            "coding"),
    coding)) %>% 
  mutate(coding = coding2) %>% 
  select(!coding2)
  

# Save intermediate
maf_table  %>%
  write_delim(paste(out_prefix,".pre_anal.tsv", sep = ""), 
              delim = "\t")
maf_table %>% 
  filter(!grepl(filter_string, FILTER)) %>% 
  filter(MAF <= MAF_lim) %>% 
  filter(t_depth >= inputs$min_depth) %>% 
  write_delim(paste(out_prefix,".pre_anal.filt.tsv", sep = ""), 
              delim = "\t")

#==================================================================================================#
# After this point, things get less general.  Put any other analyses here. 
#==================================================================================================#


