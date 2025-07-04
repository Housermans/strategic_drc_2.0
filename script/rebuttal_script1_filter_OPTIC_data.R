# Load required libraries
library(dplyr)

# Set paths relative to the Analyse 2.0 directory
resources <- "resources"
table1_data <- "STRATEGIC_table1"

# Read the OPTIC clinical data
OPTIC <- readRDS(file.path(resources, "OPTIC_clinical_data.rds"))

# Read the Project-OPTICid.csv file
project_ids <- read.csv(file.path(resources, "Project-OPTICid.csv"), sep = ";", stringsAsFactors = FALSE)

# Clean the OPTIC_ID column by removing trailing/leading whitespace
project_ids$OPTIC_ID <- trimws(project_ids$OPTIC_ID)

# Filter OPTIC data for participants that appear in Project-OPTICid.csv
# The Project-OPTICid.csv has OPTIC_ID column that should match Participant.Id
filtered_OPTIC <- OPTIC %>%
  filter(Participant.Id %in% project_ids$OPTIC_ID)

selected_OPTIC <- project_ids %>%
  left_join(
    filtered_OPTIC %>%
      select(
        Participant.Id,  # OPTIC ID
        DA_Prim,         # Date of primary tumor
        DA_Met,          # Date of metastasis
        interval
      ) %>%
      rename(
        OPTIC_ID = Participant.Id
      ),
    by="OPTIC_ID"
  ) %>%
  arrange(ProjectID)


# Display the results
print("Filtered and processed OPTIC data:")
print(selected_OPTIC)

# Save the filtered data
write.csv(selected_OPTIC, file.path(table1_data, "filtered_OPTIC_data.csv"), row.names = FALSE)
cat("Filtered data saved to 'filtered_OPTIC_data.csv'\n") 
