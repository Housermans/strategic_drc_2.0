# Version Control

#----- Data loading----
library(dplyr)
library(readxl)
library(ggplot2)
library(openxlsx)
library(scales)   
library(drc)
library(nplr)
library(patchwork)
library(stringr)#str replace

rm(list=ls())

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# script.dir <- dirname(sys.frame(1)$ofile) 
home_dir <- dirname(script_dir)
raw_dir <- file.path(home_dir, "1_raw_files")
org_data_dir <- file.path(home_dir, "2_organoid_data")
QC_dir <- file.path(home_dir, "3_QC")
resource_dir <- file.path(home_dir, "resources")
plot_dir <- file.path(home_dir, "4_plot_data")
plot_output <- file.path(home_dir, "5_plot_output")
WGS_plot_dir <- file.path(home_dir, "7_WGS")
up_dir <- dirname(home_dir)
WGS_dir <- file.path(up_dir, "WGS")
Arne_driver_dir <- paste0(WGS_dir, "/analysis_files_arne/Drivers")
resources_dir <- paste0(up_dir, "/Analyse 2.0/resources")
overview <- read_excel(file.path(resource_dir, "Screening_overview.xlsx"))
analysis_dir <- dirname(script_dir)
metrics_dir <- file.path(analysis_dir, "6_metrics", "240607_QC_includes_RAS05")
R2_dir <- file.path(analysis_dir, "8_R2")

#PDO screen data (raw)
read_plot_data <- function(exp_id, organoid_name, AUC=FALSE) {
  
  exp_file_name <- list.files(file.path(plot_dir), pattern=paste0(exp_id, "_", organoid_name,"_plot_data.xlsx"))
    
  if (length(exp_file_name) == 0) {
    return(paste("ERROR: Combination of", exp_id, "and", organoid_name, "does not exist in this folder!"))
  }
  if (AUC) {
    print(paste("Reading", exp_id, organoid_name, "by individual reps"))
    organoid_data <- read_excel(file.path(plot_dir, exp_file_name), sheet="experimental_data_individual")
  } else {
    print(paste("Reading", exp_id, organoid_name))
    organoid_data <- read_excel(file.path(plot_dir, exp_file_name))
  }
  organoid_data
}

select_files_and_add_metadata <- function(select_on = "Analyse2", plot_condition = "all", selection_value = 1, metadata_fac=c("chemo_naive"), metadata_num=c("org_id"), AUC=FALSE) {
  Analyse <- overview %>% filter(get(select_on) == selection_value)
  d <- read_plot_data(Analyse[1, "STR_ID"], Analyse[1, "org_name"], AUC=AUC)
  for (n in 2:nrow(Analyse)) {
    r <- read_plot_data(Analyse[n, "STR_ID"], Analyse[n, "org_name"], AUC=AUC)
    d <- rbind(d, r)
  }
  d
}

PDO <- select_files_and_add_metadata(select_on = "Analyse2", selection_value = 1)

PDO$org_name <- sub("_Demi$", "", PDO$org_name)
PDO$org_name <- sub("_2$", "", PDO$org_name)
PDO$org_name <- sub("_3$", "", PDO$org_name)
PDO$org_name[grepl("^OPT", PDO$org_name)] <- 
  sub("OPT(\\d{2})(\\d+)", "OPT379-\\1\\2", PDO$org_name[grepl("^OPT", PDO$org_name)])
PDO$org_name[grepl("^OPT379", PDO$org_name)] <- 
  sub("OPT379-(\\d+)(\\d{2})$", "OPT379-\\1\\000\\2", PDO$org_name[grepl("^OPT379", PDO$org_name)])

#PDO screen data (AUC)
metrics <- read_xlsx(file.path(metrics_dir, "240607_QC_includes_RAS05_metrics_normalized.xlsx")) #PDO data drug screen metrics

metrics <- metrics %>% filter(metrics$org_name != "HUB-02-C2-89")
metrics <- metrics %>% filter(metrics$org_name != "RAS05_Demi")
metrics <- metrics %>% filter(metrics$org_name != "OPT0067")

metrics$org_name <- sub("_Demi$", "", metrics$org_name)
metrics$org_name <- sub("_2$", "", metrics$org_name)
metrics$org_name <- sub("_3$", "", metrics$org_name)
metrics$org_name[grepl("^OPT", metrics$org_name)] <- 
  sub("OPT(\\d{2})(\\d+)", "OPT379-\\1\\2", metrics$org_name[grepl("^OPT", metrics$org_name)])
metrics$org_name[grepl("^OPT379", metrics$org_name)] <- 
  sub("OPT379-(\\d+)(\\d{2})$", "OPT379-\\1\\000\\2", metrics$org_name[grepl("^OPT379", metrics$org_name)])
metrics <- metrics[-c(2,4,7,8)]

library(tidyr)
metrics_wide <- pivot_wider(metrics, names_from = condition, values_from = c(AUC, GRmax, GR50, xmid))

#WGS data
mat_HMF_treatment = readRDS(file.path(WGS_plot_dir,  "mat_HMF_treatment.rds"))
colnames(mat_HMF_treatment)

WGS <- mat_HMF_treatment [-c(62, 212:214, 216:226, 228:236, 238, 239, 241, 242, 246:299)]
WGS <- WGS %>% filter(WGS$org_name != "OPT0044")
WGS <- WGS %>% filter(WGS$org_name != "OPT0404")

WGS$org_name[grepl("^OPT", WGS$org_name)] <- 
  sub("OPT(\\d{2})(\\d+)", "OPT379-\\1\\2", WGS$org_name[grepl("^OPT", WGS$org_name)])
WGS$org_name[grepl("^OPT379", WGS$org_name)] <- 
  sub("OPT379-(\\d+)(\\d{2})$", "OPT379-\\1\\000\\2", WGS$org_name[grepl("^OPT379", WGS$org_name)])

columns_to_check <- c("FHIT","NAALADL2","CCSER1","PRKN","MACROD2","WWC3")
WGS$has_deep_deletion <- apply(WGS[, columns_to_check], 1, function(row) {
  if ('deep deletion' %in% row) {
    return('present')
  } else {
    return('absent')
  }
})

#merge
metrics_wide <- merge(metrics_wide, WGS[, c("org_name", "sampleId")], by = "org_name", all.x = TRUE)
PDO <- merge(PDO, WGS[, c("org_name", "sampleId")], by = "org_name", all.x = TRUE)
metrics_wide$org_name<- NULL
PDO$org_name<- NULL
PDO$STR_ID<- NULL

write.xlsx(WGS, file = file.path(R2_dir, "WGS_metadata_weekfrozen.xlsx"))
write.xlsx(PDO, file = file.path(R2_dir, "PDO_raw.xlsx"))
write.xlsx(metrics_wide, file = file.path(R2_dir, "PDO_metrics_wide.xlsx"))

R2_txt <-  read.table("r2_annotation_250104_jk.txt", header = TRUE, sep = "\t", row.names = NULL, fill = TRUE)

metrics_wide$organoid_name<- NULL
metrics_wide_WGS<- merge(metrics_wide, WGS, by = "sampleId")
metrics_wide_WGS$organoid_name<- NULL
metrics_wide_WGS$samplenames <- paste0("JR2406_", metrics_wide_WGS$sampleId)

R2_txt <- merge(R2_txt, metrics_wide_WGS, by="samplenames")
write.table(R2_txt, file = file.path(R2_dir, "r2_annotation_250106_ls.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
