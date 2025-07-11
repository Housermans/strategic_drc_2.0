library(readxl)
library(dplyr)

#version control v1: 21-06-24 updated export with OS data UMCU
#version control v1: 11-07-24 updated export with OS data Maastricht/Radboud

setwd("D:/SURFdrive/Lsmabers/PhD/OPTIC_LS/Data/Export_11_03_2025")
OPTIC <- read.csv("Organoids_to_Predict_Treatment_r_export_20250311.csv", sep = ";", strip.white = TRUE)

OPTIC$DA_Biopsy<- as.Date(OPTIC$DA_Biopsy, format = "%d-%m-%Y")
OPTIC$DA_Prim <- as.Date(OPTIC$DA_Prim, format = "%d-%m-%Y")
OPTIC$DA_Met <- as.Date(OPTIC$DA_Met, format = "%d-%m-%Y")
OPTIC$interval_prim_met <- OPTIC$DA_Met - OPTIC$DA_Prim
OPTIC$interval <- ifelse(OPTIC$interval>182, "Metachronous", "Synchronous")

OPTIC_1  <- OPTIC[OPTIC$f_evalu_evalu_please_sel.Organoid.was.established == 1, ]

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
home_dir <- dirname(dirname(script_dir))
resources_dir <- paste0(home_dir, "/Analyse 2.0/resources")
metadata <- read_xlsx(file.path(resources_dir, "STRATEGIC_PERSCO_Overview_Hartwig_WGS.xlsx"))
ID_STRATEGIC <- metadata$`Organoid line`
ID_STRATEGIC <- grep("^OPT379", ID_STRATEGIC, value = TRUE)
subset_STRATEGIC <- subset(OPTIC_1, Participant.Id %in% ID_STRATEGIC)

file_name <- "OPTIC_clinical_data.rds"
full_path <- file.path(resources_dir, file_name)
saveRDS(subset_STRATEGIC, file = full_path)

