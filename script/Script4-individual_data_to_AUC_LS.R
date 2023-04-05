library(dplyr)
library(readxl)
library(ggplot2)
library(openxlsx)
library(scales)   
library(drc)
library(nplr)
library(patchwork)

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
metrics_dir <- file.path(home_dir, "6_metrics")

overview <- read_excel(file.path(resource_dir, "Screening_overview.xlsx"))

read_plot_data <- function(exp_id, organoid_name, individual_reps=FALSE) {
  
  exp_file_name <- list.files(file.path(plot_dir), pattern=paste0(exp_id, "_", organoid_name,"_plot_data.xlsx"))
  
  if (length(exp_file_name) == 0) {
    return(paste("ERROR: Combination of", exp_id, "and", organoid_name, "does not exist in this folder!"))
  }
  if (individual_reps) {
    print(paste("Reading", exp_id, organoid_name, "by individual reps"))
    organoid_data <- read_excel(file.path(plot_dir, exp_file_name), sheet="experimental_data_individual")
  } else {
    print(paste("Reading", exp_id, organoid_name))
    organoid_data <- read_excel(file.path(plot_dir, exp_file_name))
  }
  organoid_data
}

select_files_and_add_metadata <- function(select_on = "Analyse", plot_condition = "all", selection_value = 1, metadata_fac=c("chemo_naive", "RASTRIC", "Passed_QC", "Tween_bad", "DMSO_bad", "Clin_benefit", "Analyse", "Analyse2"), metadata_num=c("RAS_pct_change"), individual_reps=FALSE) {
  if(select_on == "all") {
    Analyse <- overview %>% filter(Data_processed == 1)
  } else {
    Analyse <- overview %>% filter(get(select_on) == selection_value)
  }
  status <- Analyse %>% dplyr::select(STR_ID, org_name, all_of(metadata_fac), all_of(metadata_num))
  d <- read_plot_data(Analyse[1, "STR_ID"], Analyse[1, "org_name"], individual_reps=individual_reps)
  for (n in 2:nrow(Analyse)) {
    r <- read_plot_data(Analyse[n, "STR_ID"], Analyse[n, "org_name"], individual_reps=individual_reps)
    d <- rbind(d, r)
  }
  d <- left_join(d, status, by=c("STR_ID", "org_name"))
  d <- d %>% mutate(across(all_of(metadata_fac), .fns = ~factor(.,levels = c(0,1), labels = c("no", "yes"))))
  if (!plot_condition == "all") {
    d <- d %>% filter(condition == plot_condition)
  } 
  d
}

# TODO: get AUCs of organoid data
save_AUC <- function(exp_name, 
                     select_on="org_id", 
                     selection_value = "RAS25", 
                     selected_conditions="all" # of c("5-FU", "Oxaliplatin", "etc.")) 
) {
  
  if (!dir.exists(file.path(metrics_dir, exp_name))) {
    dir.create(file.path(metrics_dir, exp_name))
  }
  
  df_metrics <- tibble(condition = numeric(), 
              organoid = numeric(), 
              AUC_raw = numeric(), 
              AUC_fit_trapezoid = numeric(), 
              GRMax = numeric(), 
              GR50 = numeric()
  )
  
  # TO DO: make selection for individual metrics etc.
  d <- select_files_and_add_metadata(select_on = select_on, selection_value = selection_value)
  d_individual <- select_files_and_add_metadata(select_on = select_on, selection_value = selection_value, individual_reps=TRUE)
  
  # create an empty list to store the nplr models
  fit_list <- list()
  
  for (organoid in unique(d$org_name)) {
    organoid_data = d[d$org_name == organoid, ]
    for (condition in unique(d$condition)) {
      organoid_data$Max_Concentration_log <- log10(organoid_data$conc_condition)
      
      min_organoid = min(organoid_data$mean_GR)
      max_organoid = max(organoid_data$mean_GR)
      diff_organoid = max_organoid - min_organoid
      organoid_data$GR_prop <- convertToProp(organoid_data$mean_GR)
      
      # fit the nplr model and store it in the list
      fit <- nplr(organoid_data$conc_condition, organoid_data$GR_prop, 
                  useLog = TRUE,
                  LPweight = 0.25,
                  npars = "all",
                  method = "res",
                  silent = FALSE)
      
      fit_list[[organoid]] <- fit
    }
  }
  
  # create an empty data frame to store the fitted values
  dataframe_fit <- data.frame()
  
  # loop over the nplr models and extract the fitted values
  for (i in seq_along(fit_list)) {
    
    fit <- fit_list[[i]]
    
    dataframe_fit_i <- data.frame(getXcurve(fit), getYcurve(fit))
    colnames(dataframe_fit_i) <- c("Concentration_1", "y_fit")
    
    min_organoid = min(d[d$org_name == names(fit_list)[i], ]$mean_GR)
    max_organoid = max(d[d$org_name == names(fit_list)[i], ]$mean_GR)
    diff_organoid = max_organoid - min_organoid
    
    dataframe_fit_i$y_fit_original <- min_organoid + dataframe_fit_i$y_fit * diff_organoid
    
    dataframe_fit_i$condition <- condition
    dataframe_fit_i$org_name <- names(fit_list)[i]
    
    # append the data frame to the main one
    dataframe_fit <- rbind(dataframe_fit, dataframe_fit_i)
    
  }
  AUC_df <- tibble(condition = numeric(), 
              organoid = numeric(), 
              AUC_raw = numeric(), 
              AUC_fit_trapezoid = numeric(), 
              GRMax = numeric(), 
              GR50 = numeric()
  )
  
  
  d
}
d <- save_AUC("RAS25")

#   
#   

# }
#metrics

# organoid_data_AUC$GR_positive <- (((organoid_data_AUC$x)+1)/2) #GR range 0-1 for calculating AUC
# GRMax <- mean(subset(organoid_data, Max_Concentration == max(organoid_data$Max_Concentration))$GR)
# AUC_raw <- AUC(organoid_data_AUC$Group.1, organoid_data_AUC$GR_positive) #AUC raw
# AUC_fit <- getAUC(fit)
# AUC_fit_trapezoid <- AUC_fit$trapezoid #AUC fit
# AUC_fit_simpson <- AUC_fit$Simpson #AUC fit
# estim <- getEstimates(fit, .5)
# GR50_log <- format(estim$x, digits = 6, scientific = TRUE) #estimate at 0.5
# GR50_num <- as.numeric(GR50_log)
# GR50 <- log10(GR50_num)
# par <- getPar(fit)
# xmid_log <- (par$params$xmid)
# xmid<- as.numeric(xmid_log)

# DR_summary_local[nrow(DR_summary_local) + 1,]  <-  list(condition = condition, organoid = organoid, AUC_raw = AUC_raw, AUC_fit_trapezoid = AUC_fit_trapezoid, GRMax = GRMax, GR50 = GR50)

#checkGR<- getEstimates(fit, (1-((1-GRMax)*0.5)))
#check50_log <- format(esGR$x, digits = 6, scientific = TRUE)
#checkGR50_num <- as.numeric(GR50_log)
#checkGR50 <- log10(GR50_num)