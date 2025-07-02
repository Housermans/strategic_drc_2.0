
### version control
# v2 (12-06-2024):
# toegevoegd: rm statement (line 11)
# adjusted file.path callings in loading instead of paste0

#---- load libraries and data ----
# Load dplyr and readxl packages for data manipulation and reading Excel files
library(tidyverse)
library(readxl)
library(data.table)

rm(list=ls())

# Define necessary directories for accessing and saving files
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
analysis_dir <- dirname(script_dir)
home_dir <- dirname(dirname(script_dir))
WGS_dir <- file.path(home_dir, "WGS")
WGS_data_dir <- file.path(WGS_dir, "Data_HMF")
WGS_plot_dir <- file.path(analysis_dir, "7_WGS")
Arne_driver_dir <- file.path(WGS_dir, "analysis_files_arne", "Drivers") 
resources_dir <- file.path(analysis_dir, "resources")
metrics_dir <- file.path(analysis_dir, "6_metrics", "240607_QC_includes_RAS05")
RNA_seq <- file.path(home_dir, "RNA sequencing")
metadata_RNA <- file.path(analysis_dir, "STRATEGIC_table1")

#-----Import the WGS data PDOs, screening data and metadata----- 

mat_PDO <- read.csv(file.path(Arne_driver_dir, "Driveroverview_rastric.txt"), sep = "\t") #PDO data drivers
metrics <- read_xlsx(file.path(metrics_dir, "240607_QC_includes_RAS05_metrics_normalized.xlsx")) #PDO data drug screen metrics
metadata <- read_xlsx(file.path(resources_dir, "STRATEGIC_PERSCO_Overview_Hartwig_WGS.xlsx")) #PDO metadata
purple_PDO <- read.csv(file.path(WGS_dir, "organoid_purple_purity_data_prefilter.tsv"), sep = "\t") #PDO data ploidy, purity, TMB/TML
ssig_PDO <- read.csv(file.path(paste0(WGS_dir, "/analysis_files_arne/MutSig"), "SBS_signature_contribution.txt"), sep = "\t")
dsig_PDO <- read.csv(file.path(paste0(WGS_dir, "/analysis_files_arne/MutSig"), "DBS_signature_contribution.txt"), sep = "\t")
CN_PDO <- read.csv(file.path(paste0(WGS_dir, "/analysis_files_arne/karyotype"), "Results_CopyNumber.tsv"), sep = "\t")
RNA_seq_counts <- read.csv(file.path(RNA_seq, "2024_06_10_STRATEGIC_counts_normalized_weekFrozen_batch_corrected.tsv"), sep="\t")
RNA_metadata <- read_xlsx(file.path(metadata_RNA, "STRATEGIC_count_metadata.xlsx"))

#-----Clean and merge WGS data-----

# Clean copy number (karyotype) data, change columns to rows and vice versa
CN_PDO <- as.data.frame(CN_PDO)  # Ensure CN_PDO is a dataframe
rownames(CN_PDO) <- CN_PDO$chromArms
CN_PDO$chromArms <- NULL  # Remove the first column
CN_PDO_transposed <- t(CN_PDO)
CN_PDO_transposed <- as.data.frame(CN_PDO_transposed)
CN_PDO_transposed$sampleId <- rownames(CN_PDO_transposed)
rownames(CN_PDO_transposed) <- NULL
CN_PDO <- CN_PDO_transposed %>%
  dplyr::select(sampleId, everything())

# Clean signature data
names(dsig_PDO)[names(dsig_PDO) == "sample_ID"] <- "sampleId"
names(ssig_PDO)[names(ssig_PDO) == "SampleID"] <- "sampleId"
sig_PDO <- merge(ssig_PDO, dsig_PDO, by="sampleId")
sig_PDO$SBS17 <- sig_PDO$SBS17a + sig_PDO$SBS17b
sig_PDO$SBSplatinum <- sig_PDO$SBS35 + sig_PDO$DBS5

# Function to categorize columns as "present" or "absent"
categorize_column <- function(col) {
  categories <- rep("present", length(col))
  categories[col == 0] <- "absent"
  return(categories)
}

# Columns to categorize
cols_to_categorize <- c("SBS17a", "SBS17b", "SBS25", "SBS35", "DBS5", "SBS2", "SBS13")

# Apply the function to each specified column and create new categorized columns
for (col_name in cols_to_categorize) {
  sig_PDO[[paste0(col_name, "_cat")]] <- categorize_column(sig_PDO[[col_name]])
  
  #sig_PDO <- sig_PDO %>%
  #  mutate(across(-sampleId, ~ifelse(grepl("Na", sampleId), NA, .)))
}

sig_PDO_categorized <- sig_PDO
sig_PDO_categorized$SBS17_threshold <- ifelse(sig_PDO_categorized$SBS17>1000, "present", "absent")

# Rename sample IDs for merging with WGS
metadata <- metadata %>% rename(sampleId = ProjectID) #Rename sample IDs metadata for merging
metadata$sampleId <- gsub("MetNaive", "MetNa", metadata$sampleId)
metadata$sampleId <- gsub("MetPretreated", "MetPret", metadata$sampleId)

sig_CN_PDO <- merge(sig_PDO_categorized, CN_PDO, by = "sampleId") #merge signatures, and copy numbers
mat_sig_PDO <- merge(sig_CN_PDO, mat_PDO, by = "sampleId") #merge signatures, copy number and PDO data drivers
mat_sig_PDO$sampleId <- gsub("T$", "", mat_sig_PDO$sampleId) # Rename sample IDs for merging
mat_HMF_treatment <- merge(mat_sig_PDO, metadata, by = "sampleId") #merge WGS PDO data (drivers + signatures + copy nr) and metadata

# Rename sample IDs for merging
purple_PDO$sampleId <- gsub("MetNaive", "MetNa", purple_PDO$sampleId)
purple_PDO$sampleId <- gsub("MetPretreated", "MetPret", purple_PDO$sampleId)

mat_HMF_treatment <- merge(mat_HMF_treatment, purple_PDO, by = "sampleId") #merge WGS PDO data ((drivers + signatures + copy nr + metadata) and PDO data ploidy, purity, TMB/TML
mat_HMF_treatment$`Organoid line` <- gsub("OPT379-(\\d{2})(\\d+)(\\d{2})", "OPT\\1\\3", mat_HMF_treatment$`Organoid line`, perl = TRUE) # Rename Organoid line OPTIC
mat_HMF_treatment <- mat_HMF_treatment[!grepl("tumor biopsy", mat_HMF_treatment$`Organoid line`), ] #Remove tumor biopsies from the dataset
names(mat_HMF_treatment)[names(mat_HMF_treatment) == "Organoid line"] <- "org_name" #Rename Organoid line column to org_name

metadata$`Organoid line` <- gsub("OPT379-(\\d{2})(\\d+)(\\d{2})", "OPT\\1\\3", metadata$`Organoid line`, perl = TRUE) # Rename Organoid line OPTIC
metadata <- metadata[!grepl("tumor biopsy", metadata$`Organoid line`), ] #Remove tumor biopsies from the dataset
metadata <- metadata[!grepl("Tumor biopsy", metadata$`Organoid line`), ] #Remove tumor biopsies from the dataset
metadata <- metadata[!grepl("Percorso", metadata$`Organoid line`), ] #Remove tumor biopsies from the dataset
names(metadata)[names(metadata) == "Organoid line"] <- "org_name" #Rename Organoid line column to org_name



# Select samples without contamination
mat_HMF_treatment <- subset(
  mat_HMF_treatment,
  grepl("yes", STRATEGIC_analysis) | org_name == "RAS05"
)

# Select samples without contamination
metadata <- subset(
  metadata,
  grepl("yes", STRATEGIC_analysis) | org_name == "HUB-02-C2-089"
)

# ---- clinical data categories ----
# Make new categories for clinical response data
mat_HMF_treatment$response_5FU <- ifelse(mat_HMF_treatment$Palliative_5FU_response != "PD" & mat_HMF_treatment$`Sample type` == "pretreated", "no PD",
                                         ifelse(mat_HMF_treatment$Palliative_5FU_response == "No" & mat_HMF_treatment$`Sample type` == "chemonaive", "resistant naive", NA))
mat_HMF_treatment$response_SN38 <- ifelse(mat_HMF_treatment$Palliative_irinotecan_response != "PD" & mat_HMF_treatment$`Sample type` == "pretreated", "no PD",
                                          ifelse(mat_HMF_treatment$Palliative_irinotecan_response == "No" & mat_HMF_treatment$`Sample type` == "chemonaive", "resistant naive", NA))
mat_HMF_treatment$response_oxaliplatin<- ifelse(mat_HMF_treatment$Palliative_oxaliplatin_response != "PD" & mat_HMF_treatment$`Sample type` == "pretreated", "no PD",
                                                ifelse(mat_HMF_treatment$Palliative_oxaliplatin_response == "No" & mat_HMF_treatment$`Sample type` == "chemonaive", "resistant naive", NA))
mat_HMF_treatment$Adjuvant_received_chemonaive <- ifelse(mat_HMF_treatment$Adjuvant_received == "Yes" & mat_HMF_treatment$`Sample type` == "chemonaive", "Yes", mat_HMF_treatment$Adjuvant_received)
mat_HMF_treatment$Adjuvant_chemo_naive <- ifelse(mat_HMF_treatment$Adjuvant_received_chemonaive == "Yes","adjuvant", mat_HMF_treatment$`Sample type`)
mat_HMF_treatment$response_5FU_chemo_naive<- ifelse(!is.na(mat_HMF_treatment$response_5FU),mat_HMF_treatment$response_5FU, mat_HMF_treatment$`Sample type`)
mat_HMF_treatment$response_SN38_chemo_naive <- ifelse(!is.na(mat_HMF_treatment$response_SN38),mat_HMF_treatment$response_SN38, mat_HMF_treatment$`Sample type`)
mat_HMF_treatment$response_oxaliplatin_chemo_naive <- ifelse(!is.na(mat_HMF_treatment$response_oxaliplatin),mat_HMF_treatment$response_oxaliplatin, mat_HMF_treatment$`Sample type`)
mat_HMF_treatment$response_5FU_adjuvant_selection <- ifelse(mat_HMF_treatment$Adjuvant_chemo_naive == "chemonaive" | mat_HMF_treatment$Adjuvant_chemo_naive == "pretreated", 
                                                            mat_HMF_treatment$Adjuvant_chemo_naive, "NA")
mat_HMF_treatment$response_5FU_chemo_naive_selection <- ifelse(mat_HMF_treatment$response_5FU_chemo_naive == "chemonaive" | mat_HMF_treatment$response_5FU_chemo_naive == "pretreated", 
                                                               mat_HMF_treatment$response_5FU_chemo_naive, "NA")
mat_HMF_treatment$response_5FU_chemo_naive_adjuvant_selection <- ifelse((mat_HMF_treatment$response_5FU_adjuvant_selection == "chemonaive" & mat_HMF_treatment$response_5FU_chemo_naive_selection == "chemonaive") | 
                                                                          (mat_HMF_treatment$response_5FU_adjuvant_selection == "pretreated" & mat_HMF_treatment$response_5FU_chemo_naive_selection == "pretreated"),
                                                                        mat_HMF_treatment$response_5FU_chemo_naive, "NA")

# ---- RNA seq data ----
RNA_metadata <- RNA_metadata %>%
  mutate(
    sampleId = str_replace(ProjectID, "MetNaive", "MetNa"), 
    sampleId = str_replace(sampleId, "MetPretreated", "MetPret"),
    RNA_Id = organoid_count_name
  ) %>%
  # filter(potentially_contaminated == FALSE) 
  # %>%
  dplyr::select(RNA_Id, organoid_name, sampleId)

RNA_combined <- mat_HMF_treatment %>%
  left_join(RNA_metadata, by="sampleId")

#----drugscreen data import----

# Rename metrics sample IDs for merging
metrics$org_name <- gsub("_Demi|_2|_3", "", metrics$org_name)
metrics$org_name <- gsub("HUB-02-C2-89", "HUB-02-C2-089", metrics$org_name)

#check which PDOs do not have WGS or drugscreen data
no_WGS <- anti_join(metrics, mat_HMF_treatment, by = "org_name") 
unique(no_WGS$org_name) #"HUB-02-C2-89" "OPT0426"  "OPT0502"  "RAS34" ("RAS06" wordt opnieuw gedaan! iets raars in pipeline zie orange file)   
no_drugscreen <- anti_join(mat_HMF_treatment, metrics, by = "org_name")
unique(no_drugscreen$org_name) #"OPT0044" "OPT0404"

PDO_WGS_screen <- merge(metrics, mat_HMF_treatment, by = "org_name") #merge WGS & drugscreen data
PDO_metadata_screen <- merge(metrics, metadata, by = "org_name") #merge WGS & drugscreen data

#new dataframe with only pretreated samples
PDO_WGS_screen_pretreatedfiltered <- PDO_WGS_screen %>%
  filter(chemo_naive != "yes")

#deep deletions
character_or_factor_columns <- sapply(PDO_WGS_screen_pretreatedfiltered, function(column) is.character(column))
subset_df <- PDO_WGS_screen_pretreatedfiltered[, character_or_factor_columns]
columns_with_deep_deletion <- colnames(subset_df)[sapply(subset_df, function(column) any(column == "deep deletion"))]
unique(columns_with_deep_deletion) #check in which columns a deep deletion is present

#make 1 variable for presence of deep deletions (in CFS)
columns_to_check <- c("FHIT","NAALADL2","CCSER1","PRKN","MACROD2","WWC3")
PDO_WGS_screen$has_deep_deletion <- apply(PDO_WGS_screen[, columns_to_check], 1, function(row) {
  if ('deep deletion' %in% row) {
    return('present')
  } else {
    return('absent')
  }
})

# Filter for pretreated samples
PDO_WGS_screen_pretreatedonly <- PDO_WGS_screen %>%
  mutate(across(contains(c("SBS", "DBS", "has_deep_deletion")), ~ ifelse(chemo_naive == "yes", NA, .)))

# Create AUC categories
median(PDO_WGS_screen$AUC)
median(PDO_metadata_screen$AUC)
PDO_WGS_screen$AUC_category <- ifelse(PDO_WGS_screen$AUC < 0.5, "sensitive", "resistant")
PDO_metadata_screen$AUC_category <- ifelse(PDO_metadata_screen$AUC < 0.5, "sensitive", "resistant")

# Make new df with 5 highest and lowest AUC for all treatments, regardless of treatment history
PDO_WGS_screen_top5 <- PDO_WGS_screen %>%
  group_by(condition) %>%
  arrange(condition, AUC) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 5 | rank > n() - 5) %>%
  ungroup() %>%
  dplyr::select(-rank) #gedaan voor de samples waarvan we ook WGS data hebben

# Save dataframes (rds)
file_path <- file.path(WGS_plot_dir, "PDO_metadata_screen.rds")
saveRDS(PDO_metadata_screen, file_path)

file_path <- file.path(WGS_plot_dir, "sig_PDO_categorized.rds")
saveRDS(sig_PDO_categorized, file_path)

file_path <- file.path(WGS_plot_dir, "mat_HMF_treatment.rds")
saveRDS(mat_HMF_treatment, file_path)

file_path <- file.path(WGS_plot_dir, "PDO_WGS_screen.rds")
saveRDS(PDO_WGS_screen, file_path)

file_path <- file.path(WGS_plot_dir, "PDO_WGS_screen_pretreatedonly.rds")
saveRDS(PDO_WGS_screen_pretreatedonly, file_path)

file_path <- file.path(WGS_plot_dir, "PDO_WGS_screen_pretreatedfiltered.rds")
saveRDS(PDO_WGS_screen_pretreatedfiltered, file_path)

file_path <- file.path(WGS_plot_dir, "PDO_WGS_screen_top5.rds")
saveRDS(PDO_WGS_screen_top5, file_path)
