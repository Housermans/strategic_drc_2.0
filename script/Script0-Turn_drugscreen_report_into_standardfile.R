# Drug Screening Data Processing Script
#
# This script processes drug screening data from a 384-well plate in R. The primary goals of this script are:

# To use this script, ensure that you have the dplyr and readxl packages installed. You may also need to uncomment
# the appropriate block of code depending on the plate type being used (full 384-well plate, half 384-well plate, or
# inverted half 384-well plate). Make sure to adjust the file paths as needed for your system.
#
# For any questions or issues, please consult the documentation or reach out to the script author.

# Load dplyr and readxl packages for data manipulation and reading Excel files
library(dplyr)
library(readxl)

# Clear existing objects in the workspace
rm(list=ls())

# Define necessary directories for accessing and saving files
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
home_dir <- dirname(script_dir)
resource_dir <- file.path(home_dir, "resources")


# Read the tabular detail tab from the drug screen report (.xml file) from a tecan 300D drug printer) as a data frame. 
D <- read_excel(file.path(resource_dir,"G12Cscreen_top.xlsx"))   


# Rename columns for better readability
# Add a Timepoint column with value "D5"
# Add conc_condition, Value, value_corr, and GR columns with NA values
# Adjust conc_ to appropriate names
D <- D %>%
  select(-Plate, -`Non-random\r\nwell`, -`Non-random\r\nrow`, -`Non-random\r\ncol`, -starts_with("Volume"), -starts_with("Total"))    %>%
  rename(
    Organoid = "Plate ID",
    condition = "Fluid name",
    Dispensed_well = `Dispensed\r\nwell`,
    Dispensed_row = `Dispensed\r\nrow`,
    Dispensed_col = `Dispensed\r\ncol`,
    conc_Sotorasib = "Conc. (µM)\r\nSotorasib",
    conc_Adagrasib = "Conc. (µM)\r\nAdagrasib",
    conc_Divarasib = "Conc. (µM)\r\nDivarasib",
    conc_Lapatinib = "Conc. (µM)\r\nLapatinib",
    conc_Vinorelbine = "Conc. (µM)\r\nVinorelbine",
    DMSO_pct = "DMSO %"
  ) %>%
  mutate(Timepoint = "D5", .after=Organoid) %>%
  mutate(conc_condition = NA, .after=Concentration) %>%
  mutate(Value = NA, value_corr = NA, GR = NA, .after=Value)
  
# Modify condition column based on the given rules for different plate types
# (uncomment the appropriate block of code depending on the plate type being used)

# USE WITH FULL 384 WELL PLATE
# D <- D %>%
#   mutate(condition = case_when(
#     Dispensed_col == 1 & Dispensed_row %in% 1:8 ~ "Fluorescence",
#     Dispensed_col == 2 & Dispensed_row %in% 1:8 ~ "Tween",
#     Dispensed_col == 3 & Dispensed_row %in% 1:8 ~ "DMSO_1",
#     TRUE ~ condition # keep the original value if none of the above match
#   ))

# USE WITH HALF 384 WELL PLATE
# D <- D %>%
#   mutate(condition = case_when(
#     Dispensed_col == 1 & Dispensed_row %in% 1:8 ~ "Fluorescence",
#     Dispensed_col == 2 & Dispensed_row %in% 1:8 ~ "Tween",
#     Dispensed_col == 3 & Dispensed_row %in% 1:8 ~ "DMSO_1",
#     TRUE ~ condition # keep the original value if none of the above match
#   ))

D <- D %>%
  mutate(condition = case_when(
    Dispensed_col == 1 & Dispensed_row %in% 1:8 ~ "Fluorescence",
    Dispensed_col == 2 & Dispensed_row %in% 1:8 ~ "DMSO_1",
    TRUE ~ condition # keep the original value if none of the above match
  ))


# USE WITH BOTTOM HALF 384 WELL PLATE
# D <- D %>%
#   mutate(condition = case_when(
#     Dispensed_col == 24 & Dispensed_row %in% 9:16 ~ "Fluorescence",
#     Dispensed_col == 1 & Dispensed_row %in% 9:16 ~ "Tween",
#     Dispensed_col == 2 & Dispensed_row %in% 9:16 ~ "DMSO_1",
#     TRUE ~ condition # keep the original value if none of the above match
#   ))

# D <- D %>%
#   mutate(condition = case_when(
#     Dispensed_col == 23 & Dispensed_row %in% 9:16 ~ "Fluorescence",
#     Dispensed_col == 24 & Dispensed_row %in% 9:16 ~ "DMSO_1",
#     TRUE ~ condition # keep the original value if none of the above match
#   ))

# Update the condition column based on the combinations of fluid concentrations
D <- D %>%
  mutate(condition = case_when(
    condition == "2 Fluids" & !is.na(conc_Sotorasib) & !is.na(conc_Lapatinib) ~ "Sotorasib_Lapatinib",
    condition == "2 Fluids" & !is.na(conc_Adagrasib) & !is.na(conc_Lapatinib) ~ "Adagrasib_Lapatinib",
    condition == "2 Fluids" & !is.na(conc_Divarasib) & !is.na(conc_Lapatinib) ~ "Divarasib_Lapatinib",
    condition == "3 Fluids" & !is.na(conc_Sotorasib) & !is.na(conc_Lapatinib) & !is.na(conc_Vinorelbine) ~ "Sotorasib_Lapatinib_Vinorelbine",
    condition == "3 Fluids" & !is.na(conc_Adagrasib) & !is.na(conc_Lapatinib) & !is.na(conc_Vinorelbine) ~ "Adagrasib_Lapatinib_Vinorelbine",
    condition == "3 Fluids" & !is.na(conc_Divarasib) & !is.na(conc_Lapatinib) & !is.na(conc_Vinorelbine) ~ "Divarasib_Lapatinib_Vinorelbine",
    TRUE ~ condition
  ))

# Create a new column "conc_condition" based on the max of columns starting with "conc_" (excluding "conc_condition")
conc_cols <- grep("^conc_", names(D), value = TRUE)[-1]
D <- D %>%
  mutate(conc_condition = pmax(!!!syms(conc_cols), na.rm = TRUE))

# Update the conc_condition column based on a specific condition (combination of conc_Lapatinib, conc_Binimetinib, and conc_Vinorelbine)
# D <- D %>%
#   mutate(conc_condition = ifelse(condition == "vi_bi_la",
#                                  conc_Vinorelbine,
#                                  conc_condition))

# save the changed file back to the original location
write.xlsx(D, file.path(resource_dir,"G12Cscreen_top.xlsx"))