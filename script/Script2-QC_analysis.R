library(dplyr)
library(readxl)
library(ggplot2)
library(openxlsx)
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
    print(paste("Summarizing conditions for plotting and QC", exp_id, organoid_name))
  }
  # read the organoid data file
  organoid_data <- read_excel(file.path(org_data_dir, filename))
  
  # filter for control conditions
  # control_df <- summarize(group_by(organoid_data, condition), mean_value = mean(Value), median_value = median(Value), sd_Value = sd(Value))
  # control_df <- filter(control_df, condition %in% c("Tween", "DMSO_1", "D0_ctrl", "Fluorescence"))
  # 
  # neg_ctrl1 <- filter(organoid_data, condition == "Navitoclax", conc_condition > 19)
  # if (length(neg_ctrl1) == 0) {neg_ctrl1 <- c(NA, NA, NA)}
  # neg_ctrl2 <- filter(organoid_data, condition == "Oxaliplatin", conc_condition > 900)
  # if (length(neg_ctrl2) == 0) {neg_ctrl1 <- c(NA, NA, NA)}
  # neg_ctrl <- rbind(neg_ctrl1, neg_ctrl2)
  # neg_ctrl <- summarize(group_by(neg_ctrl, condition), mean_value = mean(Value), median_value = median(Value), sd_Value = sd(Value))
  # 
  # control_df <- rbind(control_df, neg_ctrl)
  
  
  # calculate sem for each condition and conc_condition combination
  experiment_df <- filter(organoid_data, 
                          !condition %in% c("Tween",
                                            "DMSO_1",
                                            "D0_ctrl",
                                            "Fluorescence")) %>%
                    dplyr::select(Organoid, condition, conc_condition, GR) %>%
                    mutate(STR_ID = exp_id,
                           org_name = Organoid,
                           conc_condition = signif(conc_condition,4)) %>%
                    group_by(STR_ID, org_name, condition, conc_condition) %>%
                    summarise(mean_GR = mean(GR), 
                              sem_GR = sd(GR) / sqrt(n()))
                  
  no_avg_df <- filter(organoid_data, 
                      !condition %in% c("Tween", 
                                        "DMSO_1",
                                        "D0_ctrl",
                                        "Fluorescence")) %>%
              mutate(conc_condition = signif(conc_condition, 4), org_name = Organoid, STR_ID = exp_id) %>%
              dplyr::select(STR_ID, org_name, GR, condition, conc_condition) %>%
              arrange(condition)
  
  # calculate the lower and upper bounds of the 95% confidence interval
  lower_bound_95 <- experiment_df$mean_GR - 1.96 * experiment_df$sem_GR
  upper_bound_95 <- experiment_df$mean_GR + 1.96 * experiment_df$sem_GR
  # add the new columns to the tibble
  experiment_df <- tibble::add_column(experiment_df, 
                              lower_bound_95 = lower_bound_95,
                              upper_bound_95 = upper_bound_95)
  
  # Create a list of dataframes with names as worksheet names
  dfs <- list(experimental_data_average = experiment_df, experimental_data_individual = no_avg_df)
  
  QC <- read_excel(file.path(QC_dir, "QC_overview.xlsx"))
  # QC_wb <- loadWorkbook(file.path(QC_dir, "QC_overview.xlsx"))
  
  # control_vec <- as.vector(unlist(t(control_df[,2:4])), mode = "numeric")
  # control_vec <- c(exp_id, organoid_name, control_vec)
  # control_vec_df <- as.data.frame(t(control_vec))
  # colnames(control_vec_df) <- colnames(QC)
  
  # if (exp_id %in% QC$STR_ID & organoid_name %in% QC$org_name) {
  #   QC <- QC %>% filter (!(STR_ID == exp_id & org_name == organoid_name)) 
  # } 
  # QC <- rbind(QC, control_vec_df)
  
  # replace the sheet in the Excel file
  # removeWorksheet(QC_wb, sheet="Sheet 1")
  # new_sheet <- addWorksheet(QC_wb, sheet = "Sheet 1")
  # writeData(QC_wb, new_sheet, QC)
  
  # Write the list to an excel file
  write.xlsx(dfs, file = file.path(plot_dir, paste0(selected_row$STR_ID, "_" ,selected_row$org_name, "_plot_data.xlsx")))
  # saveWorkbook(QC_wb, file = file.path(QC_dir, "QC_overview.xlsx"), overwrite = TRUE)
  # write.xlsx(QC, file = file.path(QC_dir, "QC_overview.xlsx"), colWidths ="auto", numFmt="#.##0,00")
}

Plot_controls <- function(exp_id, organoid_name, set_binwidth = 300) {
  # filter for a specific STR_ID and org_name
  selected_row <- filter(overview, STR_ID == exp_id, org_name == organoid_name)
  # get the filename from the selected row
  filename <- paste0(selected_row$Date_D0, "_", selected_row$STR_ID, "_",selected_row$org_name, ".xlsx")
  
  # Looks for the experimental file based on the way it should be called and throws an error
  # if it can't find a file that matches the requirements
  exp_file_name <- list.files(file.path(org_data_dir), pattern=paste0(selected_row$STR_ID, "_",selected_row$org_name))
  if (length(exp_file_name) == 0) {
    return(paste("ERROR: Combination of", exp_id, "and", organoid_name, "does not exist in this folder!"))
  } else {
    print(paste("Plotting", exp_id, organoid_name, "control conditions"))
  }
  
  # read the organoid data file
  organoid_data <- read_excel(file.path(org_data_dir, filename))
  
  # Filter data by conditions
  D0_ctrl <- filter(organoid_data, condition == "D0_ctrl")
  Tween <- filter(organoid_data, condition == "Tween")
  DMSO_1 <- filter(organoid_data, condition == "DMSO_1")
  Fluorescence <- filter(organoid_data, condition == "Fluorescence")
  organoid_data_4 <- filter(organoid_data, condition %in% c("D0_ctrl", "Tween", "DMSO_1", "Fluorescence"))
  
  # Get the range of value_corr for D0_ctrl, Tween and DMSO_1
  x_range <- range(filter(organoid_data_4, condition != "Fluorescence")$value_corr)
  
  # Plot histograms for each condition with binwidth = 100 and adjusted x-axis
  p1 <- ggplot(D0_ctrl, aes(x = value_corr)) +
    geom_histogram(color = "black", fill = "blue", binwidth = set_binwidth) +
    labs(title = "Histogram of day 0 control values",
         x = "Fluorescence (A.U.)",
         y = "Count")
  
  p2 <- ggplot(Tween, aes(x = value_corr)) +
    geom_histogram(color = "black", fill = "red", binwidth = set_binwidth) +
    labs(title = "Histogram of day 5 0.3% Tween control values",
         x = "Fluorescence (A.U.)",
         y = "Count")
  
  p3 <- ggplot(DMSO_1, aes(x = value_corr)) +
    geom_histogram(color = "black", fill = "green", binwidth = set_binwidth) +
    labs(title = "Histogram of day 5 1% DMSO control fluorescence values",
         x = "Fluorescence (A.U.)",
         y = "Count")
  
  # Plot boxplots for all conditions
  p4 <- ggplot(organoid_data_4, aes(x = condition, y = value_corr)) +
    geom_boxplot(aes(color = condition), fill = "white") +
    geom_jitter(aes(color = condition), width = 0.1) +
    labs(title="Comparison of the different control conditions",
         x="control conditions",
         y="Fluorescence (A.U.)")
  
  # Arrange plots in a two by two figure
  p_all <- p1 + p2 + p3 + p4 + plot_layout(nrow=2)
  
  p_all <- p_all + plot_annotation(title = paste("Quality Controls", exp_id, organoid_name)) 
  p_all <- p_all +  theme(plot.title.position = "plot")

  ggsave(file.path(QC_dir, paste0(exp_id,"_", organoid_name, "_QC_plots.png")), p_all, width=4000, height=3000, units="px")
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
    } else {
      Plot_controls(exp_id, experiment_orgs[i])
    }
  }
}

### USE THIS TO READ OUT A SPECIFIC EXPERIMENT 
read_experiment("STR24A")
read_experiment("STR24B")

### USE THIS TO PROCESS ALL EXPERIMENTS
# all_exp <- unique(overview$STR_ID)
# lapply(all_exp, read_experiment)