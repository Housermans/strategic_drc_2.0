# Load the dplyr package
library(dplyr)

# Read the .xml file as a data frame
drug_screen <- read_xml("drug_screen.xml")

# Select the relevant columns from the tabular tab
drug_screen <- drug_screen %>%
  select(Dispensedwell, Dispensedrow, Dispensedcol, Fluidname, starts_with("Conc."), ends_with("DMSO %"))

# Add new columns for Organoid and Timepoint
drug_screen <- drug_screen %>%
  mutate(Organoid = paste0("X", Dispensedwell),
         Timepoint = "D5")

# Add new columns for Value, value_corr and GR
drug_screen <- drug_screen %>%
  mutate(Value = 0.01,
         value_corr = 0.01,
         GR = 0.01)

# Add new column for conc_condition
drug_screen <- drug_screen %>%
  mutate(conc_condition = paste0(Fluidname, "_", Conc))

# # Rename some columns
# drug_screen <- drug_screen %>%
#   rename(condition = Fluidname,
#          conc_<drug> = Conc.<Um>.<drug>,
#          DMSO_pct = DMSO.%)