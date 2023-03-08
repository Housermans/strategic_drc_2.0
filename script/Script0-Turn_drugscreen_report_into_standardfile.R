# Load the dplyr package
library(dplyr)
library(readxl)

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
home_dir <- dirname(script_dir)
resource_dir <- file.path(home_dir, "resources")

# Read the .xml file as a data frame
D <- read_excel(file.path(resource_dir,"FullscreenV3.xlsx"))

D <- D %>%
  rename(
    Dispensed_well = `Dispensed\r\nwell`,
    Dispensed_row = `Dispensed\r\nrow`,
    Dispensed_col = `Dispensed\r\ncol`
  )

D <- D %>%
  mutate(condition = case_when(
    Dispensed_col == 1 & Dispensed_row %in% 1:8 ~ "Fluorescence",
    Dispensed_col == 2 & Dispensed_row %in% 1:8 ~ "Tween",
    Dispensed_col == 3 & Dispensed_row %in% 1:8 ~ "DMSO_1",
    TRUE ~ condition # keep the original value if none of the above match
  ))

D <- D %>%
  mutate(condition = case_when(
    'Dispensed\ col' == 2 & 'Dispensed\ row' %in% 1:8 ~ "Tween",
    TRUE ~ condition
  ))

D <- D %>%
  mutate(condition = case_when(
    
    TRUE ~ condition
  ))

D$condition <- ifelse(D$condition == "2 Fluids" & !is.na(D$`conc_SN-38`) & !is.na(D$conc_CHEK1), "SN38_CHEK1", D$condition)
D$condition <- ifelse(D$condition == "2 Fluids" & !is.na(D$conc_Lapatinib) & !is.na(D$conc_Binimetinib), "binimetinib_lapatinib", D$condition) 
D$condition <- ifelse(D$condition == "2 Fluids" & !is.na(D$conc_Alpelisib) & !is.na(D$conc_Binimetinib), "binimetinib_alpelisib", D$condition) 
D$condition <- ifelse(D$condition == "2 Fluids" & !is.na(D$conc_Lapatinib) & !is.na(D$conc_Alpelisib), "alpelisib_lapatinib", D$condition) 
D$condition <- ifelse(D$condition == "2 Fluids" & !is.na(D$conc_Vinorelbine) & !is.na(D$conc_Navitoclax), "navitoclax_vinorelbine", D$condition)
D$condition <- ifelse(D$condition == "3 Fluids", "vi_bi_la", D$condition) 

# select only the columns that start with "conc_" except "conc_condition"
conc_cols <- grep("^conc_", names(D), value = TRUE)[-11]

# create a new column "conc_condition" based on the max of conc_cols
D <- D %>%
  mutate(conc_condition = pmax(!!!syms(conc_cols), na.rm = TRUE))

# replace conc_condition with conc_Vinorelbine if it is a combination of conc_Lapatinib, conc_Binimetinib and conc_Vinorelbine
D <- D %>%
  mutate(conc_condition = ifelse(conc_Lapatinib > 0 & conc_Binimetinib > 0 & conc_Vinorelbine > 0,
                                 conc_Vinorelbine,
                                 conc_condition))

# # Rename some columns
# drug_screen <- drug_screen %>%
#   rename(condition = Fluidname,
#          conc_<drug> = Conc.<Um>.<drug>,
#          DMSO_pct = DMSO.%)