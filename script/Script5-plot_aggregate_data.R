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

# aggregate_by <- function(select_on = "Analyse", plot_condition = "all", selection_value = 1, aggregate_on = "chemo_naive") {
#   # calculate sem for each condition and conc_condition combination
#   d <- select_files_and_add_metadata(select_on="Passed_QC", plot_condition="all", individual_reps = T)
#   
#   d <- d %>% mutate(conc_condition_rounded = signif(conc_condition, 4))
#   
#   unique_conditions <- unique(d$condition)
#   
#   for (cond in unique_conditions) {
#     condition_data <- d %>% filter(condition == cond)
#     unique_conc_conditions <- unique(condition_data$conc_condition_rounded)
#     
#     if (length(unique_conc_conditions) > 8) {
#       kmeans_result <- kmeans(condition_data$conc_condition_rounded, centers = 8)
#       representative_concentrations <- kmeans_result$centers
#       nearest_representative_index <- apply(as.matrix(dist(rbind(condition_data$conc_condition_rounded, representative_concentrations))), 1, which.min)
#       
#       # Update conc_condition_rounded values only for the current condition and non-NA indices
#       non_na_indices <- which(!is.na(d$conc_condition_rounded) & d$condition == cond)
#       d$conc_condition_rounded[non_na_indices] <- representative_concentrations[nearest_representative_index]
#     }
#   }
#   
#   experiment_df <- d %>%
#     mutate(
#       aggregate_var = aggregate_on, 
#       aggregate_var_status = get(aggregate_on)) %>%
#     group_by(aggregate_var, aggregate_var_status, condition, conc_condition_rounded) %>%
#     summarise(mean_GR = mean(GR), 
#               sem_GR = sd(GR) / sqrt(n()))
#   
#   # calculate the lower and upper bounds of the 95% confidence interval
#   lower_bound_95 <- experiment_df$mean_GR - 1.96 * experiment_df$sem_GR
#   upper_bound_95 <- experiment_df$mean_GR + 1.96 * experiment_df$sem_GR
#   # add the new columns to the tibble
#   experiment_df <- tibble::add_column(experiment_df, 
#                                       lower_bound_95 = lower_bound_95,
#                                       upper_bound_95 = upper_bound_95)
# }


aggregate_by <- function(select_on = "Analyse", plot_condition = "all", selection_value = 1, aggregate_on = "chemo_naive") {
    # calculate sem for each condition and conc_condition combination
    d <- select_files_and_add_metadata(select_on="Passed_QC", plot_condition="all", individual_reps = T)




    experiment_df <- d %>%
        mutate(
          aggregate_var = aggregate_on,
          aggregate_var_status = get(aggregate_on),
          conc_condition = signif(conc_condition, 4)) %>%
      group_by(aggregate_var, aggregate_var_status, condition, conc_condition) %>%
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
    aggregate_data$Max_Concentration_log <- log10(aggregate_data$conc_condition)

    min_var = min(aggregate_data$mean_GR)
    max_var = max(aggregate_data$mean_GR)
    diff_var = max_var - min_var
    aggregate_data$GR_prop <- convertToProp(aggregate_data$mean_GR)

    # fit the nplr model and store it in the list
    fit <- nplr(aggregate_data$conc_condition, aggregate_data$GR_prop,
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
    # selection_var <- df_filtered["chemo_naive"][[1]][1]


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
    geom_point(data = df, aes(conc_condition, mean_GR, group=aggregate_var_status, color = get(colorby))) +
    geom_errorbar(data = df, aes(conc_condition, ymin = lower_bound_95, ymax = upper_bound_95, group=aggregate_var_status, color = get(colorby))) +
    geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_original, group=aggregate_var_status, color = colorby)) +
    labs (x= "Concentration (log10) μM", y="GR", main = condition, color = aggregate_var) +
    scale_x_log10(limits = c(min(aggregate_data$conc_condition),max(aggregate_data$conc_condition))) +
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
    ggsave(file.path(plot_output, exp_name, paste0(condition,"_", exp_name, "_plot.png")), p, width=2100, height=2100, units="px")
  }
  
}

plot_multiple("aggregati0n", select_on="Passed_QC")





# plot_per_condition <- function(df, condition) {
#   df <- df[df$condition == condition, ]
#   title <- condition
#   
#   n1 <- length(unlist(unique(df$org_name)))               # Amount of default colors
#   hex_codes1 <- hue_pal()(n1)                             # Identify hex codes
#   
#   # create an empty list to store the nplr models
#   fit_list <- list()
#   
#   # loop over the unique organoids in the condition
#   for (organoid in unique(df$org_name)) {
#     print(paste("plotting", organoid, "on", condition, "by individual colors"))
#     organoid_data = df[df$org_name == organoid, ]
#     organoid_data$Max_Concentration_log <- log10(organoid_data$conc_condition)
#     
#     min_organoid = min(organoid_data$mean_GR)
#     max_organoid = max(organoid_data$mean_GR)
#     diff_organoid = max_organoid - min_organoid
#     organoid_data$GR_prop <- convertToProp(organoid_data$mean_GR)
#     
#     # fit the nplr model and store it in the list
#     fit <- nplr(organoid_data$conc_condition, organoid_data$GR_prop, 
#                 useLog = TRUE,
#                 LPweight = 0.25,
#                 npars = "all",
#                 method = "res",
#                 silent = FALSE)
#     
#     fit_list[[organoid]] <- fit
#     
#   }
#   
#   # create an empty data frame to store the fitted values
#   dataframe_fit <- data.frame()
#   
#   # loop over the nplr models and extract the fitted values
#   for (i in seq_along(fit_list)) {
#     
#     fit <- fit_list[[i]]
#     
#     dataframe_fit_i <- data.frame(getXcurve(fit), getYcurve(fit))
#     colnames(dataframe_fit_i) <- c("Concentration_1", "y_fit")
#     
#     min_organoid = min(df[df$org_name == names(fit_list)[i], ]$mean_GR)
#     max_organoid = max(df[df$org_name == names(fit_list)[i], ]$mean_GR)
#     diff_organoid = max_organoid - min_organoid
#     
#     dataframe_fit_i$y_fit_original <- min_organoid + dataframe_fit_i$y_fit * diff_organoid
#     
#     dataframe_fit_i$condition <- condition
#     dataframe_fit_i$org_name <- names(fit_list)[i]
#     
#     # append the data frame to the main one
#     dataframe_fit <- rbind(dataframe_fit, dataframe_fit_i)
#     
#   }
#   
#   # plot all the organoids with different colors based on their org_name
#   p <- ggplot() + 
#     geom_point(data = df, aes(conc_condition, mean_GR, color = org_name), size = 2) +
#     theme_classic() +
#     geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_original, color = org_name)) +
#     scale_color_manual(values = hex_codes1) + # use your custom colors here
#     labs (x= "Concentration (log10) μM", y="GR", main = condition, color = "Organoid") +  
#     scale_x_log10(limits = c(min(organoid_data$conc_condition),max(organoid_data$conc_condition))) +
#     ylim(-1,1.5) + ggtitle(condition)
#   return(p)
# }
# 
# plot_per_condition_numeric <- function(df, condition, colorby="RAS_pct_change") {
#   df <- df[df$condition == condition, ]
#   title <- condition
#   
#   # create an empty list to store the nplr models
#   fit_list <- list()
#   
#   # loop over the unique organoids in the condition
#   for (organoid in unique(df$org_name)) {
#     print(paste("plotting", organoid, "on", condition, "by", colorby))
#     organoid_data = df[df$org_name == organoid, ]
#     organoid_data$Max_Concentration_log <- log10(organoid_data$conc_condition)
#     
#     min_organoid = min(organoid_data$mean_GR)
#     max_organoid = max(organoid_data$mean_GR)
#     diff_organoid = max_organoid - min_organoid
#     organoid_data$GR_prop <- convertToProp(organoid_data$mean_GR)
#     
#     # fit the nplr model and store it in the list
#     fit <- nplr(organoid_data$conc_condition, organoid_data$GR_prop, 
#                 useLog = TRUE,
#                 LPweight = 0.25,
#                 npars = "all",
#                 method = "res",
#                 silent = FALSE)
#     
#     fit_list[[organoid]] <- fit
#     
#   }
#   
#   # create an empty data frame to store the fitted values
#   dataframe_fit <- data.frame()
#   
#   # loop over the nplr models and extract the fitted values
#   for (i in seq_along(fit_list)) {
#     
#     fit <- fit_list[[i]]
#     
#     dataframe_fit_i <- data.frame(getXcurve(fit), getYcurve(fit))
#     colnames(dataframe_fit_i) <- c("Concentration_1", "y_fit")
#     df_filtered <- filter(df, org_name == names(fit_list)[i])
#     
#     selection_vec <- df_filtered %>% pull(get(colorby))
#     selection_var <- selection_vec[1]
#     
#     min_organoid = min(df[df$org_name == names(fit_list)[i], ]$mean_GR)
#     max_organoid = max(df[df$org_name == names(fit_list)[i], ]$mean_GR)
#     diff_organoid = max_organoid - min_organoid
#     
#     dataframe_fit_i$y_fit_original <- min_organoid + dataframe_fit_i$y_fit * diff_organoid
#     
#     dataframe_fit_i$condition <- condition
#     dataframe_fit_i$org_name <- names(fit_list)[i]
#     dataframe_fit_i$colorby <- selection_var
#     
#     # append the data frame to the main one
#     dataframe_fit <- rbind(dataframe_fit, dataframe_fit_i)
#     
#   }
#   
#   # dataframe_fit$color_column <- dataframe_fit[[colorby]]
#   # Convert 'colorby' column to numeric
#   # df$colorby <- as.numeric(df$colorby)
#             
#   
#   
#   # plot all the organoids with different colors based on their org_name
#   p <- ggplot() + 
#     geom_point(data = df, aes(conc_condition, mean_GR, group=org_name, color = get(colorby))) +
#     geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_original, group=org_name, color = colorby)) +
#     labs (x= "Concentration (log10) μM", y="GR", main = condition, color = colorby) +  
#     scale_color_continuous(low="green", high= "red") +
#     scale_x_log10(limits = c(min(organoid_data$conc_condition),max(organoid_data$conc_condition))) +
#     theme_classic() +
#     ylim(-1,1.5) + 
#     ggtitle(condition)
#   
#   
#   
#   return(p)
# }
# 

# 
# plot_single <- function(exp_name, condition, select_on="Analyse", selection_value=1, colorby="chemo_naive", color_numeric=FALSE) {
#   if (!dir.exists(file.path(plot_output, exp_name))) {
#     dir.create(file.path(plot_output, exp_name))
#   }
#   d <- select_files_and_add_metadata(select_on = select_on, plot_condition = condition, selection_value = selection_value)
#   p <- plot_per_condition(d, condition)
#   if (color_numeric) {
#     q <- plot_per_condition_numeric(d, condition, colorby=colorby)
#   } else {
#     q <- plot_per_condition_factorial(d, condition, colorby=colorby)
#   }
#   p_and_q <- p + q
#   ggsave(file.path(plot_output, exp_name, paste0(condition,"_", exp_name, "_plots.png")), p_and_q, width=4200, height=2100, units="px")
# }
# 
# plot_selected_conditions <- function(exp_name, select_on="Passed_QC", selection_value = 1, colorby="chemo_naive", selected_conditions=c("5-FU", "Oxaliplatin"), color_numeric=FALSE) {
#   for (condition in selected_conditions) {
#     plot_single(exp_name = exp_name, condition = condition, select_on=select_on, selection_value=selection_value, colorby=colorby, color_numeric=color_numeric)
#   }
# }
# 
# dont_plot_selected_conditions <- function(exp_name, select_on="Passed_QC", selection_value = 1, colorby="chemo_naive",  unselected_conditions=c("5-FU", "Oxaliplatin"), color_numeric=FALSE) {
#   d <- select_files_and_add_metadata(select_on = select_on, selection_value = selection_value)
#   for (condition in unique(d$condition)) {
#     if (!condition %in% unselected_conditions) {
#       plot_single(exp_name = exp_name, condition = condition, select_on=select_on, selection_value=selection_value, colorby=colorby, color_numeric=color_numeric)
#     }
#   }
# }
