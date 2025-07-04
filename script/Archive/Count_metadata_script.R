#load packages
library(tidyverse)
library(boot) 
library(table1, quietly=TRUE)
library(rmarkdown)
library(writexl)

rm(list=ls())

#TODO: OPTIC RASTRIC cohort


#Load files
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# db_clinical <- readRDS(file = file.path(script_dir,"db_clinic_only_table1.Rds"))
db_counts <- read.csv(file.path(script_dir, "STRATEGIC_raw_counts.txt"), sep= "\t")
db_optic_original <- readRDS(file= file.path(script_dir, "OPTIC_clinical_data.rds"))
db_rastric_original <- read_delim(file.path(script_dir, "Gegevens_RASTRIC.csv"), delim=",")
db_rastric_biopsyloc <- read_delim(file.path(script_dir, "RASTRIC_biopsy_locations.csv"), delim=";")
db_strategic_original <- read_delim(file.path(script_dir, "STRATEGIC_overview2.csv"), delim=";")

length(colnames(db_counts))
sample_names <- colnames(db_counts)[7:49]


db_strategic <- db_strategic_original %>% 
    mutate(organoid_name = str_trim(`Organoid line`), 
         MAPK_status = str_trim(`MAPK status`),
         Palliative_duration = as.numeric(Palliative_duration)) %>%
  mutate(
    Palliative_duration = ifelse(is.na(Palliative_duration), 0, Palliative_duration),
    Palliative_duration = case_when(
      Palliative_duration > 100 ~ floor(Palliative_duration / 30.5),
      TRUE ~ Palliative_duration
      ),
  ) %>%
  dplyr::select(organoid_name, `Sample type`, MAPK_status, Adjuvant_received, 
                Adjuvant_5FU, Adjuvant_Ox, Palliative_duration, STRATEGIC_analysis, Palliative_5FU_cycles, Palliative_5FU_response, 
                Palliative_irinotecan_cycles, Palliative_irinotecan_response, Palliative_oxaliplatin_cycles, Palliative_oxaliplatin_response)
  
db_optic <- db_optic_original %>% mutate(
  creation_date = dmy_hms(Participant.Creation.Date),
  creation_year = year(creation_date),
  Age = creation_year - Birth_year, 
  Sex = Gender_Pt, 
  organoid_name = ID_OPTIC,
  origin_location = Biopsy_tissue
) %>%
  dplyr::select(organoid_name, Age, Sex, Sidedness, origin_location)

db_rastric <- db_rastric_original %>%
  mutate(
    organoid_name = str_replace(SUBJID, "RAS0", "RAS"),
    regisdat = dmy(regisdat),
    birthday= dmy(dob), 
    Age = as.integer(floor(interval(birthday, regisdat) / years(1))), 
    # PS = str_sub(whops_bl, -1, -1), 
    Sidedness = case_when(
      initmhloc %in% c("sigmoid", "descending colon") ~ "Left-sided (splenic flexure-sigmoid)",
      initmhloc == "cecum" ~ "Right-sided (coecum-transverse colon)",
      initmhloc %in% c("rectum", "Rectosigmoid") ~ "Rectum (rectosigmoid/rectal)",
      TRUE ~ NA_character_
    ),
    # RASMut = case_when(
    #   str_detect(krasmut, "G12") ~ "KRAS codon 12 mutation", 
    #   str_detect(krasmut, "Gly12") ~ "KRAS codon 12 mutation",
    #   str_detect(krasmut, "c.34G") ~ "KRAS codon 12 mutation", 
    #   str_detect(krasmut, "c.35G") ~ "KRAS codon 12 mutation", 
    #   str_detect(krasmut, "exon 2") ~ "KRAS codon 12 mutation", 
    #   str_detect(krasmut, "G13") ~ "KRAS codon 13 mutation", 
    #   str_detect(krasmut, "c.38G") ~ "KRAS codon 13 mutation",
    #   str_detect(krasmut, "436G") ~ "KRAS non-exon 2 mutation",
    #   str_detect(krasmut, "183A") ~ "KRAS non-exon 2 mutation",
    #   str_detect(krasmut, "873G") ~ "KRAS wildtype amplification",
    #   str_detect(nrasmut, "61") ~ "NRAS mutation"
    # ), 
    # nmeta_sites = case_when(
    #   nmeta > 3 ~ "≥ 4",
    #   .default = as.character(nmeta)
    # ),
    # nmeta_sites = factor(nmeta_sites, levels=c("1", "2", "3", "≥ 4")),
    Sex = sex
  ) %>% 
  left_join(db_rastric_biopsyloc, c("organoid_name" = "Study Number")) %>%
  mutate(
    origin_location = str_replace(Location, "metastasis", ""),
    origin_location = str_replace(origin_location, "l biopsy", ""),
    origin_location = case_when(
      str_detect(origin_location,"(skin only)") ~ NA_character_, 
      .default = origin_location
    )) %>%
  mutate(
    origin_location = str_replace(origin_location, "\\(.*?\\)", "")
  ) %>%
  dplyr::select(organoid_name, Age, Sex, Sidedness, origin_location, organoid_count_name)

# Combine all data frames
db_combined <- db_strategic %>%
  left_join(db_optic, by = "organoid_name") %>%
  left_join(db_rastric, by = "organoid_name") %>%
  mutate(
    Age = coalesce(Age.x, Age.y),
    Sex = coalesce(Sex.x, Sex.y),
    Sidedness = coalesce(Sidedness.x, Sidedness.y),
    origin_location = coalesce(origin_location.x, origin_location.y),
    `Sample type` = str_to_title(`Sample type`)
  ) %>%
  dplyr::select(-Age.x, -Age.y, -Sex.x, -Sex.y, -Sidedness.x, -Sidedness.y, -origin_location.x, -origin_location.y) %>%
  mutate(
    Sidedness = ifelse(Sidedness == "", "Unknown", Sidedness), 
    MAPK_status = case_when(
      MAPK_status == "KRAS" ~ "RAS-mutant", 
      MAPK_status == "NRAS" ~ "RAS-mutant", 
      MAPK_status == "WT" ~ "RAS/BRAF-wildtype",
      MAPK_status == "WT KRAS amplification" ~ "KRAS wildtype amplification"),
    origin_location = str_trim(origin_location), 
    organoid_count_name = str_replace(organoid_name, "-", ".")
  ) %>%
  mutate(
    organoid_count_name = case_when(
      organoid_count_name == "OPT379.0400002" ~ "OPT379.0400002.110014", 
      organoid_count_name == "OPT379.0400004" ~  "OPT379.0400004.110008", 
      organoid_count_name == "OPT379.0400008" ~ "OPT379.0400008.110007",
      .default = organoid_count_name
    )
  )

print(sample_names %in% db_combined$organoid_count_name)

# Assuming sample_names is your vector and db_combined is your data frame
missing_names <- sample_names[!(sample_names %in% db_combined$organoid_count_name)]

# Print the missing names
print(missing_names)

# Names to filter and their corresponding new organoid_count_names
names_to_filter <- c("OPT379-0000032", "OPT379-0000067", "RAS34")
new_count_names <- c("OPT00.32.PDO", "OPT00.67.PDO", "RAS34.PDO")

# Filter the rows
rows_to_duplicate <- db_combined[db_combined$organoid_name %in% names_to_filter, ]

# Change the organoid_count_name for the filtered rows
rows_to_duplicate$organoid_count_name <- new_count_names[match(rows_to_duplicate$organoid_name, names_to_filter)]

# Append the modified rows back to the original data frame
db_combined <- rbind(db_combined, rows_to_duplicate)

print(sample_names %in% db_combined$organoid_count_name)

# Assuming sample_names is your vector and db_combined is your data frame
missing_names <- sample_names[!(sample_names %in% db_combined$organoid_count_name)]

# Print the missing names
print(missing_names)

# Add the new column 'potentially_contaminated' and set to FALSE by default
db_combined$potentially_contaminated <- FALSE

# List of organoid_count_name values to set 'potentially_contaminated' to TRUE
contaminated_samples <- c("OPT379.0000024", "RAS05", "RAS25")  # Replace this with your actual list

# Update the 'potentially_contaminated' column
db_combined$potentially_contaminated[db_combined$organoid_count_name %in% contaminated_samples] <- TRUE

# Display the updated data frame
print("Updated Data Frame:")
print(db_combined)

# Filter `db_combined` and rearrange columns
db_combined_filtered <- db_combined %>%
  filter(organoid_count_name %in% sample_names) %>%
  select(organoid_count_name, everything())

# Write the data frame to an Excel file
write_xlsx(db_combined_filtered, file.path(script_dir, "STRATEGIC_count_metadata.xlsx"))

