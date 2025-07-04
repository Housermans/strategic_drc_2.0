# Load required libraries
library(dplyr)

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

# Read the OPTIC clinical data
OPTIC <- readRDS(file.path(script_dir, "OPTIC_clinical_data.rds"))

# Read the Project-OPTICid.csv file
project_ids <- read.csv(file.path(script_dir, "Project-OPTICid.csv"), sep = ";", stringsAsFactors = FALSE)

# Filter OPTIC data for participants that appear in Project-OPTICid.csv
# The Project-OPTICid.csv has OPTIC_ID column that should match Participant.Id
filtered_OPTIC <- OPTIC %>%
  filter(Participant.Id %in% project_ids$OPTIC_ID)

selected_OPTIC <- filtered_OPTIC %>%
  select(
    Participant.Id,  # OPTIC ID
    DA_Prim,         # Date of primary tumor
    DA_Met           # Date of metastasis
  )

# Convert date variables to Date format
selected_OPTIC$DA_Prim <- as.Date(selected_OPTIC$DA_Prim, format = "%d-%m-%Y")
selected_OPTIC$DA_Met <- as.Date(selected_OPTIC$DA_Met, format = "%d-%m-%Y")

# Calculate the interval between primary and metastasis
selected_OPTIC$interval_prim_met <- selected_OPTIC$DA_Met - selected_OPTIC$DA_Prim

# Categorize the interval
selected_OPTIC$interval <- ifelse(selected_OPTIC$interval_prim_met > 182, "Meta", "Sync")

# Display the results
print("Filtered and processed OPTIC data:")
print(selected_OPTIC)

# Save the filtered data
write.csv(selected_OPTIC, "filtered_OPTIC_data.csv", row.names = FALSE)
cat("Filtered data saved to 'filtered_OPTIC_data.csv'\n") 
