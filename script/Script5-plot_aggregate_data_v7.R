## Versions
# versie 3: sucessvol aangepast dat een aantal concentraties van 5-FU en 
# oxaliplatin samen worden geplot in plaats van los van elkaar (bij: 0,0392 en 
# 0,0394 zijn nu allebei veranderd naar: 0,04)
#
# versie 4 (27-05-2024): aangepast dat de HMF data nu ook ingevoegd kan worden
# voor plots van de oxaliplatin resistance yes/no plots
#
# verise 5 (07-06-2024): aangepast dat hij pdfs maakt ipv pngs
#
# Versie 7 (19-09-2024): Kleurschema veranderd

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
up_dir <- dirname(home_dir)
WGS_dir <- file.path(up_dir, "WGS")
Arne_driver_dir <- paste0(WGS_dir, "/analysis_files_arne/Drivers")
resources_dir <- paste0(up_dir, "/Analyse 2.0/resources")

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

select_files_and_add_metadata <- function(select_on = "Analyse", plot_condition = "all", selection_value = 1, metadata_fac=c("chemo_naive", "RASTRIC", "Passed_QC", "Tween_bad", "DMSO_bad", "Clin_benefit", "Analyse", "Analyse2"), metadata_num=c("RAS_pct_change"), individual_reps=FALSE, WGS=FALSE) {
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
  if(WGS==TRUE){
    d$org_name_clean <- sub("_.*", "", d$org_name)
    d <- merge(d, mat_HMF_treatment, by = "org_name_clean") #merge WGS & drugscreen data)
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

aggregate_by <- function(select_on = "Analyse", plot_condition = "all", selection_value = 1, aggregate_on = "chemo_naive", WGS=FALSE) {
  # calculate sem for each condition and conc_condition combination
  d <- select_files_and_add_metadata(select_on="Passed_QC", plot_condition="all", individual_reps = T, WGS=WGS)
  
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

plot_multiple <- function(exp_name, select_on="Analyse", selection_value=1, aggregate_on="chemo_naive", color_numeric=FALSE, WGS=FALSE) {
  if (!dir.exists(file.path(plot_output, exp_name))) {
    dir.create(file.path(plot_output, exp_name))
  }
  d <- aggregate_by(select_on = select_on, aggregate_on = aggregate_on, WGS=WGS)
  for (condition in unique(d$condition)) {
    p <- plot_per_condition_aggregate(d, condition)
    ggsave(file.path(plot_output, exp_name, paste0(condition,"_", exp_name, "_by_", aggregate_on ,"_plot.pdf")), p, width=2100, height=2100, units="px")
  }
}

plot_single <- function(exp_name, condition = NULL, select_on="Passed_QC", selection_value=1, aggregate_on="chemo_naive", color_numeric=FALSE, WGS=FALSE) {
  if (!dir.exists(file.path(plot_output, exp_name))) {
    dir.create(file.path(plot_output, exp_name))
  }
  d <- aggregate_by(select_on = select_on, aggregate_on = aggregate_on, WGS=WGS)
  p <- plot_per_condition_aggregate(d, condition)
  ggsave(file.path(plot_output, exp_name, paste0(condition,"_", exp_name, "_by_", aggregate_on ,"_plot.pdf")), p, width=2100, height=2100, units="px")
}

# claude additions ----
# Updated calculate_spread function with NA handling
calculate_spread <- function(df, spread_type = "default") {
  df %>%
    group_by(aggregate_var, aggregate_var_status, condition, aggregate_conc) %>%
    summarise(
      mean_GR = mean(GR, na.rm = TRUE),
      lower = case_when(
        spread_type == "default" ~ mean_GR - (sd(GR, na.rm = TRUE) / sqrt(sum(!is.na(GR)))),
        spread_type == "sd" ~ mean_GR - sd(GR, na.rm = TRUE),
        spread_type == "quartile" ~ quantile(GR, 0.25, na.rm = TRUE),
        spread_type == "combined" ~ min(GR, na.rm = TRUE)
      ),
      upper = case_when(
        spread_type == "default" ~ mean_GR + (sd(GR, na.rm = TRUE) / sqrt(sum(!is.na(GR)))),
        spread_type == "sd" ~ mean_GR + sd(GR, na.rm = TRUE),
        spread_type == "quartile" ~ quantile(GR, 0.75, na.rm = TRUE),
        spread_type == "combined" ~ max(GR, na.rm = TRUE)
      ),
      lower_quartile = quantile(GR, 0.25, na.rm = TRUE),
      upper_quartile = quantile(GR, 0.75, na.rm = TRUE),
      sem = sd(GR, na.rm = TRUE) / sqrt(sum(!is.na(GR))),
      .groups = 'drop'
    )
}

# Updated plot_per_condition_aggregate function with error checking
# Updated plot_per_condition_aggregate function for numeric chemo_naive coding
# Updated plot_per_condition_aggregate function with error checking
plot_per_condition_aggregate <- function(df, condition, colorby = "aggregate_var_status", spread_type = "default") {
  df <- df[df$condition == condition, ]
  title <- condition
  
  # Check if dataframe is empty
  if (nrow(df) == 0) {
    warning("No data available for the specified condition.")
    return(NULL)
  }
  
  n1 <- length(unlist(unique(df[colorby])))
  hex_codes1 <- hue_pal()(n1)
  
  p <- ggplot(df, aes(x = aggregate_conc, y = mean_GR, group = aggregate_var_status, color = get(colorby))) +
    scale_x_log10(limits = c(min(df$aggregate_conc), max(df$aggregate_conc))) +
    labs(x = "Concentration (log10) uM", y = "GR", color = df$aggregate_var[1]) +
    scale_color_manual(values=c("#F8A29E", "#7FBFF5")) +
    theme_classic() +
    ylim(-1, 1.5) +
    ggtitle(condition)
  
  p <- ggplot(df, aes(x = aggregate_conc, y = mean_GR, group = aggregate_var_status, color = get(colorby))) +
    scale_x_log10(limits = c(min(df$aggregate_conc), max(df$aggregate_conc))) +
    labs(x = "Concentration (log10) uM", y = "GR", color = "Treatment Status") +
    theme_classic() +
    ylim(-1, 1.5) +
    ggtitle(condition) +
    scale_color_manual(values = c("#F8A29E", "#7FBFF5"), 
                       labels = c("Pretreated", "Chemonaive"),
                       name = "Treatment Status") +
    scale_fill_manual(values = c("#F8A29E", "#7FBFF5"), 
                      labels = c("Pretreated", "Chemonaive"),
                      name = "Treatment Status")
  
  if (spread_type == "default") {
    p <- p + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1)
    p <- p + geom_line()
    p <- p + geom_point()
  } else if (spread_type == "combined") {
    p <- p + geom_ribbon(aes(ymin = lower, ymax = upper, fill = get(colorby)), alpha = 0.1)
    p <- p + geom_ribbon(aes(ymin = lower_quartile, ymax = upper_quartile, fill = get(colorby)), alpha = 0.2)
    p <- p + geom_line()
    p <- p + geom_errorbar(aes(ymin = mean_GR - sem, ymax = mean_GR + sem), width = 0.1)
    p <- p + geom_point()
    p <- p + labs(fill = df$aggregate_var[1])
  } else {
    p <- p + geom_ribbon(aes(ymin = lower, ymax = upper, fill = get(colorby)), alpha = 0.3)
    p <- p + geom_line()
    p <- p + geom_point()
    p <- p + labs(fill = df$aggregate_var[1])
  }
  return(p)
}
# Updated aggregate_by function to ensure spread_type is passed correctly
aggregate_by <- function(select_on = "Analyse", plot_condition = "all", selection_value = 1, aggregate_on = "chemo_naive", WGS = FALSE, spread_type = "default") {
  d <- select_files_and_add_metadata(select_on = "Passed_QC", plot_condition = "all", individual_reps = T, WGS = WGS)
  
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
      aggregate_var_status = get(aggregate_on)
    )
  
  calculate_spread(experiment_df, spread_type)
}

# Modified plot_single function
plot_single <- function(exp_name, condition = NULL, select_on = "Passed_QC", selection_value = 1, aggregate_on = "chemo_naive", color_numeric = FALSE, WGS = FALSE, spread_type = "default") {
  if (!dir.exists(file.path(plot_output, exp_name))) {
    dir.create(file.path(plot_output, exp_name))
  }
  
  d <- aggregate_by(select_on = select_on, aggregate_on = aggregate_on, WGS = WGS, spread_type = spread_type)
  
  p <- plot_per_condition_aggregate(d, condition, spread_type = spread_type)
  
  file_name <- paste0(condition, "_", exp_name, "_by_", aggregate_on)
  if (spread_type != "default") {
    file_name <- paste0(file_name, "_", spread_type)
  }
  file_name <- paste0(file_name, "_plot.pdf")
  
  ggsave(file.path(plot_output, exp_name, file_name), p, width = 2100, height = 2100, units = "px")
}

# Example usage:
# Default (mean + whiskers):
# plot_single("aggregatedv5", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = '5-FU', WGS = FALSE)
# plot_single("aggregatedv5", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = '5-FU', WGS = FALSE, spread_type = "sd")
# plot_single("aggregatedv5", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = '5-FU', WGS = FALSE, spread_type = "quartile")
# plot_single("aggregatedv5", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = '5-FU', WGS = FALSE, spread_type = "combined")
# plot_single("aggregatedv5", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'Oxaliplatin', WGS = FALSE)
# plot_single("aggregatedv5", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'Oxaliplatin', WGS = FALSE, spread_type = "sd")
# plot_single("aggregatedv5", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'Oxaliplatin', WGS = FALSE, spread_type = "quartile")
# plot_single("aggregatedv5", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'Oxaliplatin', WGS = FALSE, spread_type = "combined")
# plot_single("aggregatedv5", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'SN-38', WGS = FALSE)
# plot_single("aggregatedv5", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'SN-38', WGS = FALSE, spread_type = "sd")
# plot_single("aggregatedv5", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'SN-38', WGS = FALSE, spread_type = "quartile")
# plot_single("aggregatedv5", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'SN-38', WGS = FALSE, spread_type = "combined")

plot_single("aggregatedv6", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = '5-FU', WGS = FALSE, spread_type = "quartile")
plot_single("aggregatedv6", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'Alpelisib', WGS = FALSE, spread_type = "quartile")
plot_single("aggregatedv6", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'Binimetinib', WGS = FALSE, spread_type = "quartile")
plot_single("aggregatedv6", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'CHEK1', WGS = FALSE, spread_type = "quartile")
plot_single("aggregatedv6", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'Lapatinib', WGS = FALSE, spread_type = "quartile")
plot_single("aggregatedv6", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'Navitoclax', WGS = FALSE, spread_type = "quartile")
plot_single("aggregatedv6", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'Oxaliplatin', WGS = FALSE, spread_type = "quartile")
plot_single("aggregatedv6", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'SN-38', WGS = FALSE, spread_type = "quartile")
plot_single("aggregatedv6", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'SN38_CHEK1', WGS = FALSE, spread_type = "quartile")
plot_single("aggregatedv6", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'Vinorelbine', WGS = FALSE, spread_type = "quartile")
plot_single("aggregatedv6", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'alpelisib_lapatinib', WGS = FALSE, spread_type = "quartile")
plot_single("aggregatedv6", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'binimetinib_alpelisib', WGS = FALSE, spread_type = "quartile")
plot_single("aggregatedv6", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'binimetinib_lapatinib', WGS = FALSE, spread_type = "quartile")
plot_single("aggregatedv6", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'navitoclax_vinorelbine', WGS = FALSE, spread_type = "quartile")
plot_single("aggregatedv6", select_on = "Passed_QC", aggregate_on = "chemo_naive", condition = 'vi_bi_la', WGS = FALSE, spread_type = "quartile")


# With SD spread:
# plot_single("aggregatedv5", select_on = "Passed_QC", aggregate_on = "response_5FU_chemo_naive", condition = '5-FU', WGS = TRUE, spread_type = "sd")

# With quartile spread:
# plot_single("aggregatedv5", select_on = "Passed_QC", aggregate_on = "response_5FU_chemo_naive", condition = '5-FU', WGS = TRUE, spread_type = "quartile")

####---- WGS plot data loading ----
#load and clean WGS data
ssig_PDO <- read.csv(file.path(paste0(WGS_dir, "/analysis_files_arne/MutSig"), "SBS_signature_contribution.txt"), sep = "\t")
dsig_PDO <- read.csv(file.path(paste0(WGS_dir, "/analysis_files_arne/MutSig"), "DBS_signature_contribution.txt"), sep = "\t")
names(dsig_PDO)[names(dsig_PDO) == "sample_ID"] <- "sampleId"
names(ssig_PDO)[names(ssig_PDO) == "SampleID"] <- "sampleId"
sig_PDO <- merge(ssig_PDO, dsig_PDO, by="sampleId")
sig_PDO <- sig_PDO[, c("sampleId", "SBS17a", "SBS17b", "SBS25", "SBS35", "DBS5")]

# Function to categorize each column
categorize_column <- function(col) {
  categories <- rep("middle", length(col))
  categories[col == 0] <- "low"
  
  # Get indices of top 3 values
  top_3_indices <- order(col, decreasing = TRUE)[1:3]
  categories[top_3_indices] <- "high"
  
  return(categories)
}

# Create a new data frame to store categorized columns
categorized_columns <- sig_PDO

# Apply the function to each column except the first one (sampleId)
for (col_name in names(sig_PDO)[-1]) {
  categorized_columns[[paste0(col_name, "_cat")]] <- categorize_column(sig_PDO[[col_name]])
}

categorized_columns <- categorized_columns %>%
  mutate(across(-sampleId, ~ifelse(grepl("Na", sampleId), NA, .)))

mat_PDO <- read.csv(file.path(Arne_driver_dir, "Driveroverview_rastric.txt"), sep = "\t")
metadata <- read_xlsx(file.path(resources_dir, "STRATEGIC_PERSCO_Overview_Hartwig_WGS.xlsx"))
metadata <- metadata %>% rename(sampleId = ProjectID)
metadata$sampleId <- gsub("MetNaive", "MetNa", metadata$sampleId)
metadata$sampleId <- gsub("MetPretreated", "MetPret", metadata$sampleId)
mat_sig_PDO <- merge(categorized_columns, mat_PDO, by = "sampleId") #merge meta data and WGS PDO data
mat_sig_PDO$sampleId <- gsub("T$", "", mat_sig_PDO$sampleId)
mat_HMF_treatment <- merge(mat_sig_PDO, metadata, by = "sampleId") #merge meta data and WGS PDO data

mat_HMF_treatment$`Organoid line` <- gsub("OPT379-(\\d{2})(\\d+)(\\d{2})", "OPT\\1\\3", mat_HMF_treatment$`Organoid line`, perl = TRUE)
mat_HMF_treatment <- mat_HMF_treatment[!grepl("tumor biopsy", mat_HMF_treatment$`Organoid line`), ]
names(mat_HMF_treatment)[names(mat_HMF_treatment) == "Organoid line"] <- "org_name_clean"

mat_HMF_treatment$response_5FU <- ifelse(mat_HMF_treatment$Palliative_5FU_response != "PD" & mat_HMF_treatment$`Sample type` == "pretreated", "pretreated, without PD", NA)

mat_HMF_treatment$response_SN38 <- ifelse(mat_HMF_treatment$Palliative_irinotecan_response != "PD" & mat_HMF_treatment$`Sample type` == "pretreated", "no PD",
                                          ifelse(mat_HMF_treatment$Palliative_irinotecan_response == "No" & mat_HMF_treatment$`Sample type` == "chemonaive", "resistant naive", NA))
mat_HMF_treatment$response_oxaliplatin<- ifelse(mat_HMF_treatment$Palliative_oxaliplatin_response != "PD" & mat_HMF_treatment$`Sample type` == "pretreated", "pretreated, without PD", NA)
mat_HMF_treatment$Adjuvant_received_chemonaive <- ifelse(mat_HMF_treatment$Adjuvant_received == "Yes" & mat_HMF_treatment$`Sample type` == "chemonaive", "Yes", mat_HMF_treatment$Adjuvant_received)
mat_HMF_treatment$Adjuvant_chemo_naive <- ifelse(mat_HMF_treatment$Adjuvant_received_chemonaive == "Yes","adjuvant", mat_HMF_treatment$`Sample type`)
mat_HMF_treatment$response_5FU_chemo_naive<- ifelse(!is.na(mat_HMF_treatment$response_5FU),mat_HMF_treatment$response_5FU, mat_HMF_treatment$`Sample type`)
mat_HMF_treatment$response_SN38_chemo_naive <- ifelse(!is.na(mat_HMF_treatment$response_SN38),mat_HMF_treatment$response_SN38, mat_HMF_treatment$`Sample type`)
mat_HMF_treatment$response_oxaliplatin_chemo_naive <- ifelse(!is.na(mat_HMF_treatment$response_oxaliplatin),mat_HMF_treatment$response_oxaliplatin, mat_HMF_treatment$`Sample type`)
mat_HMF_treatment$response_5FU_adjuvant_selection <- ifelse(mat_HMF_treatment$Adjuvant_chemo_naive == "chemonaive" | mat_HMF_treatment$Adjuvant_chemo_naive == "pretreated", 
                                                            mat_HMF_treatment$Adjuvant_chemo_naive, "NA")
mat_HMF_treatment$response_5FU_chemo_naive_selection <- ifelse(mat_HMF_treatment$response_5FU_chemo_naive == "chemonaive" | mat_HMF_treatment$response_5FU_chemo_naive == "pretreated", 
                                                               mat_HMF_treatment$response_5FU_chemo_naive, "NA")
mat_HMF_treatment$response_5FU_chemo_naive_adjuvant_selection <- ifelse((mat_HMF_treatment$response_5FU_adjuvant_selection == "chemonaive" & mat_HMF_treatment$response_5FU_chemo_naive_selection == "chemonaive") | 
                                                                          (mat_HMF_treatment$response_5FU_adjuvant_selection == "pretreated" & mat_HMF_treatment$response_5FU_chemo_naive_selection == "pretreated"),
                                                                        mat_HMF_treatment$response_5FU_chemo_naive, "NA")

# ---- function calls ----
plot_multiple("aggregatedv5", select_on="Passed_QC")
plot_single("aggregatedv5", select_on = "Passed_QC", aggregate_on = "response_5FU_chemo_naive", condition='5-FU', WGS=TRUE)
plot_single("aggregatedv5", select_on = "Passed_QC", aggregate_on = "response_oxaliplatin_chemo_naive", condition='Oxaliplatin', WGS=TRUE)

