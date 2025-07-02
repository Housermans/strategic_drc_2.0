## Versions
# versie 3: sucessvol aangepast dat een aantal concentraties van 5-FU en 
# oxaliplatin samen worden geplot in plaats van los van elkaar (bij: 0,0392 en 
# 0,0394 zijn nu allebei veranderd naar: 0,04)

#----library and data loading----

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

overview <- read_excel(file.path(resource_dir, "Screening_overview.xlsx"))

##---- plot functions ----
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

# Function to extract unique concentrations for specified conditions
extract_concentrations <- function(data, condition_name) {
  data %>%
    filter(condition == condition_name) %>%
    dplyr::select(conc_condition) %>%
    distinct() %>%
    arrange(conc_condition) %>%
    mutate(new_conc = NA)  # Placeholder for new concentration assignment
}

# Manually assign new concentrations for 5-FU and Oxaliplatin
assign_new_concentrations <- function(data) {
  # Extract unique concentrations for 5-FU and Oxaliplatin
  concentrations_5fu <- extract_concentrations(data, "5-FU")
  concentrations_oxaliplatin <- extract_concentrations(data, "Oxaliplatin")
  
  # Manually assign new concentrations for 5-FU
  concentrations_5fu$new_conc <- c(0.27, 0.27, 1.17, 1.17, 5.32, 5.32, 24.4, 24.4, 108, 108, 490, 490, 2213, 2213, 10000)
  
  # Manually assign new concentrations for Oxaliplatin
  concentrations_oxaliplatin$new_conc <- c(0.1, 0.1, 0.4, 0.4, 1.55, 1.55, 6.14, 6.14, 24.4, 24.4, 96.1, 96.1, 380, 1000, 1500)
  
  list(concentrations_5fu = concentrations_5fu, concentrations_oxaliplatin = concentrations_oxaliplatin)
}

# Function to update concentrations using the lookup table
update_concentrations <- function(data, lookup_table) {
  data %>%
    left_join(lookup_table, by = "conc_condition") %>%
    mutate(aggregate_conc = ifelse(!is.na(new_conc), new_conc, conc_condition)) %>%
    dplyr::select(-new_conc)
}

aggregate_by <- function(select_on = "Analyse", plot_condition = "all", selection_value = 1, aggregate_on = "chemo_naive") {
  # calculate sem for each condition and conc_condition combination
  d <- select_files_and_add_metadata(select_on="Passed_QC", plot_condition="all", individual_reps = T)
  
  # Assign new concentrations for 5-FU and Oxaliplatin
  new_conc_assignments <- assign_new_concentrations(d)
  concentrations_5fu <- new_conc_assignments$concentrations_5fu
  concentrations_oxaliplatin <- new_conc_assignments$concentrations_oxaliplatin
  
  # Update the original dataset for 5-FU and Oxaliplatin
  d_updated_5fu <- update_concentrations(d %>% filter(condition == "5-FU"), concentrations_5fu)
  d_updated_oxaliplatin <- update_concentrations(d %>% filter(condition == "Oxaliplatin"), concentrations_oxaliplatin)
  
  # Combine the updated datasets with the rest of the data
  d_rest <- d %>% filter(!condition %in% c("5-FU", "Oxaliplatin")) %>% mutate(aggregate_conc = conc_condition)
  d <- bind_rows(d_rest, d_updated_5fu, d_updated_oxaliplatin)
  
  experiment_df <- d %>%
    mutate(
      aggregate_var = aggregate_on,
      aggregate_var_status = get(aggregate_on)) %>%
    group_by(aggregate_var, aggregate_var_status, condition, aggregate_conc) %>%
    summarise(mean_GR = mean(GR),
              sem_GR = sd(GR) / sqrt(n()))
  
  # calculate the lower and upper bounds of the 95% confidence interval
  lower_bound_95 <- experiment_df$mean_GR - 1.96 * experiment_df$sem_GR
  upper_bound_95 <- experiment_df$mean_GR + 1.96 * experiment_df$sem_GR
  # add the new columns to the tibble
  experiment_df <- tibble::add_column(experiment_df,
                                      lower_bound_95 = lower_bound_95,
                                      upper_bound_95 = upper_bound_95)
}

plot_per_condition_aggregate <- function(df, condition, colorby="aggregate_var_status") {
  df <- df[df$condition == condition, ]
  title <- condition
  
  n1 <- length(unlist(unique(df[colorby])))               # Amount of default colors
  hex_codes1 <- hue_pal()(n1)                             # Identify hex codes
  
  # create an empty list to store the nplr models
  fit_list <- list()
  
  # loop over the unique organoids in the condition
  for (var in unique(df$aggregate_var_status)) {
    # print(paste("plotting summary of", df$aggregate_var, "on", condition))
    aggregate_data = df[df$aggregate_var_status == var, ]
    aggregate_data$Max_Concentration_log <- log10(aggregate_data$aggregate_conc)
    
    min_var = min(aggregate_data$mean_GR)
    max_var = max(aggregate_data$mean_GR)
    diff_var = max_var - min_var
    aggregate_data$GR_prop <- convertToProp(aggregate_data$mean_GR)
    
    # fit the nplr model and store it in the list
    fit <- nplr(aggregate_data$aggregate_conc, aggregate_data$GR_prop,
                useLog = TRUE,
                LPweight = 0.25,
                npars = "all",
                method = "res",
                silent = FALSE)
    
    fit_list[[var]] <- fit
  }
  
  # create an empty data frame to store the fitted values
  dataframe_fit <- data.frame()
  
  # loop over the nplr models and extract the fitted values
  for (i in seq_along(fit_list)) {
    fit <- fit_list[[i]]
    
    dataframe_fit_i <- data.frame(getXcurve(fit), getYcurve(fit))
    colnames(dataframe_fit_i) <- c("Concentration_1", "y_fit")
    df_filtered <- filter(df, aggregate_var_status == names(fit_list)[i])
    
    selection_vec <- df_filtered %>% pull(get(colorby))
    selection_var <- selection_vec[1]
    
    min_var = min(df[df$aggregate_var_status == names(fit_list)[i], ]$mean_GR)
    max_var = max(df[df$aggregate_var_status == names(fit_list)[i], ]$mean_GR)
    diff_var = max_var - min_var
    
    dataframe_fit_i$y_fit_original <- min_var + dataframe_fit_i$y_fit * diff_var
    
    dataframe_fit_i$condition <- condition
    dataframe_fit_i$aggregate_var_status <- names(fit_list)[i]
    dataframe_fit_i$colorby <- selection_var
    
    # append the data frame to the main one
    dataframe_fit <- rbind(dataframe_fit, dataframe_fit_i)
  }
  
  aggregate_var <- df$aggregate_var
  # plot all the vars with different colors based on their org_name
  p <- ggplot() +
    geom_point(data = df, aes(aggregate_conc, mean_GR, group=aggregate_var_status, color = get(colorby))) +
    geom_errorbar(data = df, aes(aggregate_conc, ymin = lower_bound_95, ymax = upper_bound_95, group=aggregate_var_status, color = get(colorby))) +
    geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_original, group=aggregate_var_status, color = colorby)) +
    labs (x= "Concentration (log10) μM", y="GR", main = condition, color = aggregate_var) +
    scale_x_log10(limits = c(min(aggregate_data$aggregate_conc),max(aggregate_data$aggregate_conc))) +
    theme_classic() +
    ylim(-1,1.5) +
    ggtitle(condition)
  
  return(p)
}

plot_multiple <- function(exp_name, select_on="Analyse", selection_value=1, aggregate_on="chemo_naive", color_numeric=FALSE) {
  if (!dir.exists(file.path(plot_output, exp_name))) {
    dir.create(file.path(plot_output, exp_name))
  }
  d <- aggregate_by(select_on = select_on, aggregate_on = aggregate_on)
  for (condition in unique(d$condition)) {
    p <- plot_per_condition_aggregate(d, condition)
    ggsave(file.path(plot_output, exp_name, paste0(condition,"_", exp_name, "_by_", aggregate_on ,"_plot.png")), p, width=2100, height=2100, units="px")
  }
}


# ---- function calls ----
plot_multiple("aggregatedv3", select_on="Passed_QC")



