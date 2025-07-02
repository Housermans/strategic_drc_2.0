### version control
# v6 (01-08-2024)
# change sample names

# v5 (26-06-2024):
# Import all purple.purity.tsv files in the WGS folder toegevoegd met nieuwe organoid data
# Import new sig dataframe already merged SBS and DBS in Script12.MutSig
# Import new CNV dataframe made in Script13

# v4 (17-06-2024):
# een 'maak top 5' functie die standaard aanpassingen bevat die elke
# df moet ondergaan en die gelinkt aan een filter stap, zodat we top5 pre-
# treated en chemonaive kunnen maken voor RNA en WGS databases

# v3 (13-06-2024):
# - Switched to loading the tidyverse library for broader functionality.
# - Integrated RNA ID within all metadata files to make comparisons
# - Created clinical response data categories within metadata instead of within 
# WGS data to make sure that this can also apply to datasets that do nothave WGS
# data
# - Changed order where different WGS mutations were done
# - Calculated median AUC per condition instead of overall and categorized AUC 
# values accordingly.
# - Applied the new mutations to the metrics dataframe before merging 
# so these changes are applied to all daughter dataframes
# - Merged drug screen data with WGS and RNA data.

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
WGS_analysis_files_arne2_dir <- file.path(WGS_dir, "analysis_files_arne2/GCRuns_organoidWGSdata")

# Import all purple.purity.tsv files in the WGS folder
dirs<-Sys.glob(file.path(WGS_analysis_files_arne2_dir, "/"))
sampleslist=list.files(dirs, pattern = "purple.purity.tsv$",full.names=TRUE, recursive = TRUE)
sampleslist_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
datalist = lapply(sampleslist, function(x){read.table(file=x,sep="\t", header=T)})
names(datalist) <- sampleslist_names

datalist1 <- do.call(rbind, datalist)
datalist1 <- tibble::rownames_to_column(datalist1, "sampleId")

write.table(datalist1, file = file.path(WGS_dir, "organoid_purple_purity_data_prefilter.tsv"), sep="\t", col.names = TRUE, quote = FALSE, row.names = FALSE)

#-----Import the WGS data PDOs, screening data and metadata----- 
# 54 samples in total, 37 non contaminated etc

mat_PDO <- read.csv(file.path(WGS_plot_dir, "Driveroverview.txt"), sep = "\t") #PDO data drivers
metrics <- read_xlsx(file.path(metrics_dir, "240607_QC_includes_RAS05_metrics_normalized.xlsx")) #PDO data drug screen metrics
metadata <- read_xlsx(file.path(resources_dir, "STRATEGIC_PERSCO_Overview_Hartwig_WGS.xlsx")) #PDO metadata
purple_PDO <- read.csv(file.path(WGS_dir, "organoid_purple_purity_data_prefilter.tsv"), sep = "\t") #PDO data ploidy, purity, TMB/TML
sig_PDO <- read.csv(file.path(WGS_plot_dir, "SBS_DBS_signature_contribution.txt"), sep = "\t")
CN_PDO <- read.csv(file.path(WGS_plot_dir, "Results_CopyNumber.tsv"), sep = "\t")
# RNA_seq_counts <- read.csv(file.path(RNA_seq, "2024_06_10_STRATEGIC_counts_normalized_weekFrozen_batch_corrected.tsv"), sep="\t")
RNA_metadata <- read_xlsx(file.path(metadata_RNA, "STRATEGIC_count_metadata.xlsx"))

#-----Clean and merge WGS data-----

# Rename sample IDs for merging with WGS
metadata <- metadata %>% mutate(sampleId = ProjectID) #Rename sample IDs metadata for merging
metadata$sampleId <- gsub("MetNaive", "MetNa", metadata$sampleId)
metadata$sampleId <- gsub("MetPretreated", "MetPret", metadata$sampleId)

metadata$`Organoid line` <- gsub("OPT379-(\\d{2})(\\d+)(\\d{2})", "OPT\\1\\3", metadata$`Organoid line`, perl = TRUE) # Rename Organoid line OPTIC
metadata <- metadata[!grepl("tumor biopsy", metadata$`Organoid line`), ] #Remove tumor biopsies from the dataset
metadata <- metadata[!grepl("Tumor biopsy", metadata$`Organoid line`), ] #Remove tumor biopsies from the dataset
metadata <- metadata[!grepl("Percorso", metadata$`Organoid line`), ] #Remove tumor biopsies from the dataset
names(metadata)[names(metadata) == "Organoid line"] <- "org_name" #Rename Organoid line column to org_name

# # Select samples without contamination
metadata <- subset(
  metadata,
  grepl("yes", STRATEGIC_analysis) | org_name == "HUB-02-C2-089"
)

# ---- clinical data categories ----
# Make new categories for clinical response data
metadata$response_5FU <- ifelse(metadata$Palliative_5FU_response != "PD" & metadata$`Sample type` == "pretreated", "no PD",
                                         ifelse(metadata$Palliative_5FU_response == "No" & metadata$`Sample type` == "chemonaive", "resistant naive", NA))
metadata$response_SN38 <- ifelse(metadata$Palliative_irinotecan_response != "PD" & metadata$`Sample type` == "pretreated", "no PD",
                                          ifelse(metadata$Palliative_irinotecan_response == "No" & metadata$`Sample type` == "chemonaive", "resistant naive", NA))
metadata$response_oxaliplatin<- ifelse(metadata$Palliative_oxaliplatin_response != "PD" & metadata$`Sample type` == "pretreated", "no PD",
                                                ifelse(metadata$Palliative_oxaliplatin_response == "No" & metadata$`Sample type` == "chemonaive", "resistant naive", NA))
metadata$Adjuvant_received_chemonaive <- ifelse(metadata$Adjuvant_received == "Yes" & metadata$`Sample type` == "chemonaive", "Yes", metadata$Adjuvant_received)
metadata$Adjuvant_chemo_naive <- ifelse(metadata$Adjuvant_received_chemonaive == "Yes","adjuvant", metadata$`Sample type`)
metadata$response_5FU_chemo_naive<- ifelse(!is.na(metadata$response_5FU),metadata$response_5FU, metadata$`Sample type`)
metadata$response_SN38_chemo_naive <- ifelse(!is.na(metadata$response_SN38),metadata$response_SN38, metadata$`Sample type`)
metadata$response_oxaliplatin_chemo_naive <- ifelse(!is.na(metadata$response_oxaliplatin),metadata$response_oxaliplatin, metadata$`Sample type`)
metadata$response_5FU_adjuvant_selection <- ifelse(metadata$Adjuvant_chemo_naive == "chemonaive" | metadata$Adjuvant_chemo_naive == "pretreated", 
                                                            metadata$Adjuvant_chemo_naive, "NA")
metadata$response_5FU_chemo_naive_selection <- ifelse(metadata$response_5FU_chemo_naive == "chemonaive" | metadata$response_5FU_chemo_naive == "pretreated", 
                                                               metadata$response_5FU_chemo_naive, "NA")
metadata$response_5FU_chemo_naive_adjuvant_selection <- ifelse((metadata$response_5FU_adjuvant_selection == "chemonaive" & metadata$response_5FU_chemo_naive_selection == "chemonaive") | 
                                                                          (metadata$response_5FU_adjuvant_selection == "pretreated" & metadata$response_5FU_chemo_naive_selection == "pretreated"),
                                                                        metadata$response_5FU_chemo_naive, "NA")

# ---- RNA seq data ----
# deze naar voren gehaald om aan metadata RNA count id toe te voegen, zodat we WGS ook aan RNA seq kunnen koppelen
RNA_metadata <- RNA_metadata %>%
  mutate(
    sampleId = str_replace(ProjectID, "MetNaive", "MetNa"), 
    sampleId = str_replace(sampleId, "MetPretreated", "MetPret"),
    RNA_Id = organoid_count_name
  ) %>%
  # filter(potentially_contaminated == FALSE) 
  # %>%
  dplyr::select(RNA_Id, organoid_name, sampleId)

RNA_combined <- metadata %>%
  left_join(RNA_metadata, by="sampleId") %>%
  dplyr::filter(RNA_Id !="OPT379.0000032", 
                RNA_Id != "RAS34")

metadata <- RNA_combined 

RNA_combined <- RNA_combined %>%
  dplyr::filter(org_name != "HUB-02-C2-089")

# ---- WGS data ----

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
names(sig_PDO)[names(sig_PDO) == "SampleID"] <- "sampleId"
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

sig_CN_PDO <- merge(sig_PDO_categorized, CN_PDO, by = "sampleId") #merge signatures, and copy numbers
mat_sig_PDO <- merge(sig_CN_PDO, mat_PDO, by = "sampleId") #merge signatures, copy number and PDO data drivers
mat_sig_PDO$sampleId <- gsub("T$", "", mat_sig_PDO$sampleId) # Rename sample IDs for merging
mat_HMF_treatment <- merge(mat_sig_PDO, metadata, by = "sampleId") #merge WGS PDO data (drivers + signatures + copy nr) and metadata

# Rename sample IDs for merging
purple_PDO$sampleId <- gsub("T", "", purple_PDO$sampleId)

mat_HMF_treatment <- merge(mat_HMF_treatment, purple_PDO, by = "sampleId") #merge WGS PDO data ((drivers + signatures + copy nr + metadata) and PDO data ploidy, purity, TMB/TML
# mat_HMF_treatment$`Organoid line` <- gsub("OPT379-(\\d{2})(\\d+)(\\d{2})", "OPT\\1\\3", mat_HMF_treatment$`Organoid line`, perl = TRUE) # Rename Organoid line OPTIC
# mat_HMF_treatment <- mat_HMF_treatment[!grepl("tumor biopsy", mat_HMF_treatment$`Organoid line`), ] #Remove tumor biopsies from the dataset
# names(mat_HMF_treatment)[names(mat_HMF_treatment) == "Organoid line"] <- "org_name" #Rename Organoid line column to org_name

# Select samples without contamination
mat_HMF_treatment <- subset(
  mat_HMF_treatment,
  grepl("yes", STRATEGIC_analysis) | org_name == "RAS05"
)


#----drugscreen data import----


# Rename metrics sample IDs for merging
metrics$org_name <- gsub("_Demi|_2|_3", "", metrics$org_name)
metrics$org_name <- gsub("HUB-02-C2-89", "HUB-02-C2-089", metrics$org_name)

#check which PDOs do not have WGS or drugscreen data
no_WGS <- anti_join(metrics, mat_HMF_treatment, by = "org_name") 
unique(no_WGS$org_name) #"HUB-02-C2-89" "OPT0426"  "OPT0502"  "RAS34" ("RAS06" wordt opnieuw gedaan! iets raars in pipeline zie orange file)   
no_drugscreen <- anti_join(mat_HMF_treatment, metrics, by = "org_name")
unique(no_drugscreen$org_name) #"OPT0044" "OPT0404"

# --- calculate median per condition and apply them to drugscreen data ----
# Calculate the median AUC for each condition
condition_medians <- metrics %>%
  group_by(condition) %>%
  summarize(median_AUC = median(AUC, na.rm = TRUE))

# Join the median AUC values back to the original data
metrics <- metrics %>%
  left_join(condition_medians, by = "condition")

# Use the condition-specific median to categorize AUC values
metrics <- metrics %>%
  mutate(AUC_category = ifelse(AUC < median_AUC, "sensitive", "resistant"))

# Print a summary to verify the categorization
summary(metrics$AUC_category)

# ---- merge PDO and drugscreen data ----

PDO_WGS_screen <- merge(metrics, mat_HMF_treatment, by = "org_name") #merge WGS & drugscreen data
PDO_metadata_screen <- merge(metrics, metadata, by = "org_name") #merge WGS & drugscreen data
PDO_RNA_screen <- merge(metrics, RNA_combined, by="org_name")

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

# Change sample names
PDO_metadata_screen$sampleId <- gsub("Met|a|ret", "", PDO_metadata_screen$sampleId)
sig_PDO_categorized$sampleId <- gsub("Met|a|ret", "", sig_PDO_categorized$sampleId)
PDO_WGS_screen$sampleId <- gsub("Met|a|ret", "", PDO_WGS_screen$sampleId)
PDO_RNA_screen$sampleId <- gsub("Met|a|ret", "", PDO_RNA_screen$sampleId)
mat_HMF_treatment$sampleId <- gsub("Met|a|ret", "", mat_HMF_treatment$sampleId)

# Filter for pretreated samples
PDO_WGS_screen_pretreatedonly <- PDO_WGS_screen %>%
  mutate(across(contains(c("SBS", "DBS", "has_deep_deletion")), ~ ifelse(chemo_naive == "yes", NA, .)))

#new dataframe with only pretreated samples, now including deep deletion column
PDO_WGS_screen_pretreatedfiltered <- PDO_WGS_screen %>%
  filter(chemo_naive != "yes")

# # Create AUC categories
# median(PDO_WGS_screen$AUC)
# median(PDO_metadata_screen$AUC)
# PDO_WGS_screen$AUC_category <- ifelse(PDO_WGS_screen$AUC < 0.5, "sensitive", "resistant")
# PDO_metadata_screen$AUC_category <- ifelse(PDO_metadata_screen$AUC < 0.5, "sensitive", "resistant")

get_top5 <- function(df) {
  df %>%
    group_by(condition) %>%
    arrange(condition, AUC) %>%
    mutate(rank = row_number()) %>%
    filter(rank <= 5 | rank > n() - 5) %>%
    ungroup() %>%
    dplyr::select(-rank)
}

# Make new df with 5 highest and lowest AUC for all treatments, regardless of treatment history
PDO_WGS_screen_top5 <- PDO_WGS_screen %>%
  get_top5()

# Make new df with 5 highest and lowest AUC for all treatments, regardless of treatment history
PDO_WGS_screen_top5_pretreated <- PDO_WGS_screen %>%
  filter(`Sample type` == "pretreated") %>%
  get_top5()

# Make new df with 5 highest and lowest AUC for all treatments, regardless of treatment history
PDO_WGS_screen_top5_chemonaive <- PDO_WGS_screen %>%
  filter(`Sample type` == "chemonaive") %>%
  get_top5()

# Make new df with 5 highest and lowest AUC for all treatments, regardless of treatment history
PDO_RNA_screen_top5 <- PDO_RNA_screen %>%
  get_top5()

# Make new df with 5 highest and lowest AUC for all treatments, regardless of treatment history
PDO_RNA_screen_top5_pretreated <- PDO_RNA_screen %>%
  filter(`Sample type` == "pretreated") %>%
  get_top5()

# Make new df with 5 highest and lowest AUC for all treatments, regardless of treatment history
PDO_RNA_screen_top5_chemonaive <- PDO_RNA_screen %>%
  filter(`Sample type` == "chemonaive") %>%
  get_top5()


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

file_path <- file.path(WGS_plot_dir, "PDO_WGS_screen_top5_pretreated.rds")
saveRDS(PDO_WGS_screen_top5_pretreated, file_path)

file_path <- file.path(WGS_plot_dir, "PDO_WGS_screen_top5_chemonaive.rds")
saveRDS(PDO_WGS_screen_top5_chemonaive, file_path)

file_path <- file.path(WGS_plot_dir, "PDO_RNA_screen.rds")
saveRDS(PDO_RNA_screen, file_path)

file_path <- file.path(WGS_plot_dir, "PDO_RNA_screen_top5.rds")
saveRDS(PDO_RNA_screen_top5, file_path)

file_path <- file.path(WGS_plot_dir, "PDO_RNA_screen_top5_pretreated.rds")
saveRDS(PDO_RNA_screen_top5_pretreated, file_path)

file_path <- file.path(WGS_plot_dir, "PDO_RNA_screen_top5_chemonaive.rds")
saveRDS(PDO_RNA_screen_top5_chemonaive, file_path)


# # Function to compare top selections between two data frames for a given condition
# compare_top_selections <- function(condition, df1, df2, key_column) {
#   top_df1 <- df1 %>% filter(condition == !!condition) %>% pull(!!key_column)
#   top_df2 <- df2 %>% filter(condition == !!condition) %>% pull(!!key_column)
#   
#   list(
#     in_df1_not_in_df2 = setdiff(top_df1, top_df2),
#     in_df2_not_in_df1 = setdiff(top_df2, top_df1)
#   )
# }
# 
# # Get unique conditions
# conditions <- unique(c(PDO_WGS_screen_top5$condition, PDO_RNA_screen_top5$condition))
# 
# # Compare top selections for each condition
# comparison_results <- lapply(conditions, function(cond) {
#   compare_top_selections(cond, PDO_WGS_screen_top5, PDO_RNA_screen_top5, key_column = "org_name")
# })
# 
# names(comparison_results) <- conditions
# 
# # Function to summarize differences
# summarize_differences <- function(results) {
#   summary <- lapply(names(results), function(cond) {
#     res <- results[[cond]]
#     list(
#       condition = cond,
#       in_WGS_not_in_RNA = length(res$in_df1_not_in_df2),
#       in_RNA_not_in_WGS = length(res$in_df2_not_in_df1),
#       WGS_only = res$in_df1_not_in_df2,
#       RNA_only = res$in_df2_not_in_df1
#     )
#   })
#   return(summary)
# }
# 
# # Summarize the differences
# difference_summary <- summarize_differences(comparison_results)
# 
# # Print the summary
# print(difference_summary)
