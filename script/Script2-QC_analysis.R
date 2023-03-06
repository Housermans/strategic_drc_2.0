library(dplyr)
library(readxl)
library(ggplot2)
library(openxlsx)


rm(list=ls())

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# script.dir <- dirname(sys.frame(1)$ofile) 
home_dir <- dirname(script_dir)
raw_dir <- file.path(home_dir, "1_raw_files")
org_data_dir <- file.path(home_dir, "2_organoid_data")
QC_dir <- file.path(home_dir, "3_QC")
resource_dir <- file.path(home_dir, "resources")
plot_dir <- file.path(home_dir, "4_plot_data")

overview <- read_excel(file.path(resource_dir, "Screening_overview.xlsx"))

Get_QC_organoid <- function(exp_id, organoid_name) {
  # filter for a specific STR_ID and org_name
  selected_row <- filter(overview, STR_ID == exp_id, org_name == organoid_name)
  # get the filename from the selected row
  filename <- paste0(selected_row$Date_D0, "_", selected_row$STR_ID, "_",selected_row$org_name, ".xlsx")
  
  exp_file_name <- list.files(file.path(org_data_dir), pattern=paste0(selected_row$STR_ID, "_",selected_row$org_name))
  if (length(exp_file_name) == 0) {
    return(paste("ERROR: Combination of", exp_id, "and", organoid_name, "does not exist in this folder!"))
  } else {
    print(paste("Reading", exp_id, organoid_name))
  }
  # read the organoid data file
  organoid_data <- read_excel(file.path(org_data_dir, filename))
  
  # filter for control conditions
  control_df <- summarize(group_by(organoid_data, condition), mean_value = mean(Value), median_value = median(Value), sd_Value = sd(Value))
  control_df <- filter(control_df, condition %in% c("Tween", "DMSO_1", "D0_ctrl", "Fluorescence"))
  
  neg_ctrl1 <- filter(organoid_data, condition == "Navitoclax", conc_condition > 19)
  neg_ctrl2 <- filter(organoid_data, condition == "Oxaliplatin", conc_condition > 900)
  neg_ctrl <- rbind(neg_ctrl1, neg_ctrl2)
  neg_ctrl <- summarize(group_by(neg_ctrl, condition), mean_value = mean(Value), median_value = median(Value), sd_Value = sd(Value))
  
  control_df <- rbind(control_df, neg_ctrl)
  
  
  # calculate sem for each condition and conc_condition combination
  experiment_df <- summarize(group_by(organoid_data,
                                      condition,
                                      signif(conc_condition,4)),
                             mean_GR = mean(GR),
                             sem_GR = sd(GR) / sqrt(n()))
  
  # filter out control conditions
  experiment_df <- filter(experiment_df,
                          !condition %in% c("Tween",
                                             "DMSO_1",
                                             "D0_ctrl",
                                             "Fluorescence"))
  
  experiment_df <- rename(experiment_df, conc_condition = 'signif(conc_condition, 4)')
  
  
  # calculate the lower and upper bounds of the 95% confidence interval
  lower_bound_95 <- experiment_df$mean_GR - 1.96 * experiment_df$sem_GR
  upper_bound_95 <- experiment_df$mean_GR + 1.96 * experiment_df$sem_GR
  # add the new columns to the tibble
  experiment_df <- tibble::add_column(experiment_df, lower_bound_95 = lower_bound_95,
                              upper_bound_95 = upper_bound_95, org_name = organoid_name, STR_ID = exp_id)
  
  relocate(experiment_df, org_name:STR_ID, .before=condition)
  
  # Create a list of dataframes with names as worksheet names
  dfs <- list(experimental_data = experiment_df, controls_summary = control_df)
  
  QC <- read_excel(file.path(QC_dir, "QC_overview.xlsx"))
  
  control_vec <- as.vector(unlist(t(control_df[,2:4])), mode = "numeric")
  control_vec <- c(exp_id, organoid_name, control_vec)
  control_vec_df <- as.data.frame(t(control_vec))
  colnames(control_vec_df) <- colnames(QC)
  
  if (exp_id %in% QC$STR_ID & organoid_name %in% QC$org_name) {
    QC <- QC %>% filter (!(STR_ID == exp_id & org_name == organoid_name)) 
  } 
  QC <- rbind(QC, control_vec_df)
  
  # Write the list to an excel file
  write.xlsx(dfs, file = file.path(plot_dir, paste0(selected_row$STR_ID, "_" ,selected_row$org_name, "_plot_data.xlsx")))
  write.xlsx(QC, file = file.path(QC_dir, "QC_overview.xlsx"))
}

read_experiment <- function(exp_id) {
  experiment_orgs <- overview %>%
                              filter(STR_ID==exp_id)
  experiment_orgs <- unique(experiment_orgs$org_name)
  for (i in 1:length(experiment_orgs)) {
    print(paste("Opening:", experiment_orgs[i]))
    org <- Get_QC_organoid(exp_id, experiment_orgs[i])
    if (typeof(org) == "character") {
      print(org)
    }
  }
}

# read_experiment("STR23")
# all_exp <- unique(overview$STR_ID)
# lapply(all_exp, read_experiment)