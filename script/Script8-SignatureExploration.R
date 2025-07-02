library(tidyverse)
library(pheatmap)
library(readxl)

# Clear the workspace
rm(list=ls())

# Define directories relative to the script location
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)


home_dir <- dirname(dirname(script_dir))

WGS_dir <- file.path(home_dir, "WGS")
Arne_signature_dir <- file.path(WGS_dir, "analysis_files_arne2", "analysis", "MutSig")

rna_dir <- file.path(home_dir, "RNA sequencing", "RNA seq analyses Maarten Lidwien")
raw_dir <- file.path(rna_dir, "input", "DE_signature")
plot_dir <- file.path(rna_dir, "plots")
res_dir <- file.path(rna_dir, "resource")

metadata <- read_xlsx(file.path(res_dir, "240515_STRATEGIC_PERSCO_Overview_Hartwig_WGS.xlsx"))
count_table <- read_tsv(file.path(raw_dir, "2024_05_22_STRATEGIC_counts_normalized.tsv"))
sig_PDO <- read.csv(file.path(Arne_signature_dir, "SBS_signature_contribution.txt"), sep= "\t")

apo_sig_PDO <- sig_PDO %>% 
  select(SampleID, SBS2, SBS13)

# Replace missing gene names with Ensembl_ID
count_table <- count_table %>%
  mutate(gene_name = ifelse(is.na(gene_name), Ensembl_ID, gene_name))

# Check for duplicates and handle them
make_unique <- function(names) {
  dup <- duplicated(names)
  while (any(dup)) {
    names[dup] <- paste0(names[dup], "_", seq_len(sum(dup)))
    dup <- duplicated(names)
  }
  names
}
count_table$gene_name <- make_unique(count_table$gene_name)

metadata_sig <- metadata %>% 
  mutate(SampleID = str_replace(ProjectID, "MetNaive", "MetNa"), 
         SampleID = str_replace(SampleID, "MetPretreated", "MetPret"), 
         SampleID = paste0(SampleID, "T")) %>%
  right_join(apo_sig_PDO, by="SampleID") %>%
  select(SampleID, ProjectID, `Sample type`, Palliative_oxaliplatin_response, SBS2, SBS13) %>%
  filter(`Sample type` == "pretreated")

# Create the boxplot
sbs2_ox <- ggplot(metadata_sig, aes(x = Palliative_oxaliplatin_response, y = SBS2)) +
  geom_boxplot() +
  labs(title = "Boxplot of SBS2 by Palliative Oxaliplatin Response",
       x = "Palliative Oxaliplatin Response",
       y = "SBS2") +
  theme_minimal()

ggsave(file.path(plot_dir, "sbs2vsoxresponse.png"), plot=sbs2_ox)

# Create the boxplot
sbs13_ox <- ggplot(metadata_sig, aes(x = Palliative_oxaliplatin_response, y = SBS13)) +
  geom_boxplot() +
  labs(title = "Boxplot of SBS13 by Palliative Oxaliplatin Response",
       x = "Palliative Oxaliplatin Response",
       y = "SBS2") +
  theme_minimal()

ggsave(file.path(plot_dir, "sbs13vsoxresponse.png"), plot=sbs13_ox)

count_table_no_names <- count_table %>%
  select(-Ensembl_ID) %>%
  rename(gene = gene_name) %>%
  column_to_rownames(var="gene") %>%
  mutate(across(everything(), as.numeric))



# Get the list of genes of interest from DE results
genes_of_interest <- DE_results %>%
  filter(abs(log2FoldChange) > 1 & padj < 0.05) %>%
  pull(gene_name)

# Filter count table for genes of interest
count_table_filtered <- count_table_no_names[rownames(count_table_no_names) %in% genes_of_interest, ]

# Prepare sample names to merge with metadata
sample_names <- tibble(
  sample = colnames(count_table_no_names)
)

metadata_select <- metadata %>%
  mutate(
    sample = gsub("-", ".", `Organoid line`),
    Sample_type = factor(`Sample type`), 
    MAPK_status = factor(`MAPK status`), 
    Adjuvant_Ox = factor(Adjuvant_Ox) 
    # Folfox_AUClow = `5-FU_DS_AUClow` + Oxaliplatin_DS_AUClow,
    # Folfox_sensitive = factor(case_when(
    #   is.na(`5-FU_DS_AUClow`) | is.na(Oxaliplatin_DS_AUClow) ~ NA_character_,
    #   (`5-FU_DS_AUClow` + Oxaliplatin_DS_AUClow) < 1 ~ "Sensitive",
    #   TRUE ~ "Resistant"))
  ) %>%
  filter(STRATEGIC_analysis == "yes") %>%
  select(sample, Sample_type, MAPK_status, Adjuvant_Ox, Palliative_oxaliplatin_response
         #, Folfox_sensitive
  ) %>%
  inner_join(sample_names, by = "sample")


