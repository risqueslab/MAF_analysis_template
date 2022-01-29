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
variant_clasification_table = read_delim(inputs$varClassTrans) 

seshat_col_trans = read_delim(inputs$seshat_com1_translate)
seshat_path_trans <- read_delim(inputs$seshat_path_translate)

# Filter

sample_data <- read_delim(sample_data_file, delim="\t")

# Load seshat file
seshat_long <- read_delim(inputs$seshat_long_file, delim="\t") %>% 
  # remove unnecessary columns
  select(cDNA_Variant, Comment_1_Frequency, Pathogenicity) %>% 
  # Pull out first non-space letter of Comment_1_Frequency
  mutate(C1F = str_sub(str_replace_all(Comment_1_Frequency, " ",""), 1,1)) %>%
  # Join on the classification by C1F
  left_join(seshat_col_trans, by="C1F") %>% 
  # Join on the classification by Pathogenicity
  left_join(seshat_path_trans, by="Pathogenicity") 

# Load MAF file
TP53_maf_table <- read_delim(in_maf_file, delim="\t", skip=1) %>% 
  # Remove unnecessary columns
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
  # Filter to just TP53 mutants
  filter(Hugo_Symbol == "TP53") %>% 
  # Join with Seshat output
  bind_cols(seshat_long) %>% 
  # Add MAF
  mutate(MAF=t_alt_count/t_depth) %>% 
  # Join on Variant Classification
  left_join(variant_clasification_table, by=c("Variant_Classification")) %>%
  # Join on extra sample information
  left_join(sample_data, by=c("Tumor_Sample_Barcode"="Sample")) %>% 
  # Fix coding classification for Splice Regions
  mutate(coding2 = if_else(
    Variant_Classification == "Splice_Region", 
    if_else(is.na(Exon_Number), 
            "non-coding", 
            "coding"),
    coding)) %>% 
  mutate(coding = coding2) %>% 
  select(!coding2)
  

# Save intermediate
TP53_maf_table  %>%
  write_delim(paste(out_prefix,".TP53.pre_anal.tsv", sep = ""), 
              delim = "\t")
TP53_coding_maf <- TP53_maf_table %>% 
  filter(!grepl(filter_string, FILTER)) %>% 
  filter(MAF <= MAF_lim) %>% 
  filter(t_depth >= inputs$min_depth) %>% 
  filter(coding=="coding") %>%
  write_delim(paste(out_prefix,".TP53.pre_anal.filtCoding.tsv", sep = ""), 
              delim = "\t")

TP53_noncoding_maf <- TP53_maf_table %>% 
  filter(!grepl(filter_string, FILTER)) %>% 
  filter(MAF <= MAF_lim) %>% 
  filter(t_depth >= inputs$min_depth) %>% 
  filter(!coding=="coding") %>%
  write_delim(paste(out_prefix,".TP53.pre_anal.filtNonCoding.tsv", sep = ""), 
              delim = "\t")

#==================================================================================================#
# After this point, things get less general.  Put any other analyses here. 
#==================================================================================================#

