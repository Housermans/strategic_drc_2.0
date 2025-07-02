library(readxl)
library(dplyr)

setwd("D:/SURFdrive/Lsmabers/PhD/OPTIC_LS/Data/Export 16_04_2024")
OPTIC <- read.csv("Organoids_to_Predict_Treatment_r_export_20240416.csv", sep = ";", strip.white = TRUE)
OPTIC_1  <- OPTIC[OPTIC$f_evalu_evalu_please_sel.Organoid.was.established == 1, ]

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
home_dir <- dirname(dirname(script_dir))
resources_dir <- paste0(home_dir, "/Analyse 2.0/resources")
metadata <- read_xlsx(file.path(resources_dir, "STRATEGIC_PERSCO_Overview_Hartwig_WGS.xlsx"))
ID_STRATEGIC <- grep("^OPT379", ID_STRATEGIC, value = TRUE)
subset_STRATEGIC <- subset(OPTIC_1, Participant.Id %in% ID_STRATEGIC)

file_name <- "OPTIC_clinical_data.rds"
full_path <- file.path(resources_dir, file_name)
saveRDS(subset_STRATEGIC, file = full_path)
