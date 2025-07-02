# Load necessary packages
library(tidyverse)
library(DESeq2)
library(readxl)
library(ashr)
library(ggrepel)

# Clean up the environment
rm(list=ls())

# Define script directories
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
home_dir <- dirname(script_dir)
output_dir <- file.path(home_dir, "output")
rna_dir <- dirname(home_dir)
raw_dir <- file.path(home_dir, "input", "DE_per_condition")
plot_dir <- file.path(output_dir, "deseq", "plots")
table_dir <- file.path(output_dir, "deseq", "tables")
res_dir <- file.path(home_dir, "resource")

# Read necessary files
count_table <- read.csv(file.path(res_dir, "2024_05_22_STRATEGIC_combined_exon_featureCounts.raw.txt"), sep = "\t")
metadata_full <- read_xlsx(file.path(res_dir, "240613_STRATEGIC_PERSCO_Overview_Hartwig_WGS.xlsx"))
metadata_rna <- read_xlsx(file.path(res_dir, "STRATEGIC_count_metadata.xlsx"))
RNA_top5 <- readRDS(file.path(raw_dir, "PDO_RNA_screen_top5.rds"))
RNA_top5_chemonaive <- readRDS(file.path(raw_dir, "PDO_RNA_screen_top5_chemonaive.rds"))
RNA_top5_pretreated <- readRDS(file.path(raw_dir, "PDO_RNA_screen_top5_pretreated.rds"))
RNA <- readRDS(file.path(raw_dir, "PDO_RNA_screen.rds"))
WGS <- readRDS(file.path(raw_dir, "PDO_WGS_screen.rds"))

RNA <- RNA %>%
  mutate(weekFrozen = case_when(
    org_name == "OPT0032" ~ 117,
    org_name == "RAS34" ~ 117, 
    .default = `Week frozen`
  )
  )

WGS <- WGS %>%
  mutate(weekFrozen = case_when(
    org_name == "OPT0032" ~ 117,
    org_name == "RAS34" ~ 117, 
    .default = `Week frozen`
  )
  )

# Assume feature_counts is a subset of count_table for DESeq analysis
feature_counts <- count_table %>%
  dplyr::select(-gene_name, -Chr, -Start, -End, -Strand, -Length)

gene_names <- count_table %>%
  rownames_to_column(var = 'Ensembl_ID') %>%
  dplyr::select(Ensembl_ID, gene_name)

# write_xlsx(gene_names, file.path(res_dir, "gene_names_df.xlsx"))
# RNA <- RNA %>% 
#   mutate(ox_resp_all = case_when(
#     Palliative_oxaliplatin_response_all == "PD" ~ "Refractory",
#     Palliative_oxaliplatin_response_all == "SD" ~ "Not refractory",
#     Palliative_oxaliplatin_response_all == "PR" ~ "Not refractory",
#     Palliative_oxaliplatin_response_all == "SD/PR" ~ "Not refractory",
#     .default = Palliative_oxaliplatin_response_all
#   ), 
#   ox_resp_all = as.factor(ox_resp_all))

# contrast_var = "response_oxaliplatin_chemo_naive" # contrast var to test clinically responsive vs unresponsive
