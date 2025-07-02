# Version Control
## V7
# Veranderd dat we ook pdfs krijgen van de plots
# Veranderd dat we ook van 1 responder meerdere behandelingen kunnen plotten.


#----- Data loading----
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
WGS_plot_dir <- file.path(home_dir, "7_WGS")
up_dir <- dirname(home_dir)
WGS_dir <- file.path(up_dir, "WGS")
Arne_driver_dir <- paste0(WGS_dir, "/analysis_files_arne/Drivers")
resources_dir <- paste0(up_dir, "/Analyse 2.0/resources")

overview <- read_excel(file.path(resource_dir, "Screening_overview.xlsx"))

#---- Plotting Functions----

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
  
  # Convert numeric columns after reading
  numeric_cols <- c("conc_condition", "mean_GR", "GR", "value_corr")
  for(col in numeric_cols) {
    if(col %in% names(organoid_data)) {
      organoid_data[[col]] <- as.numeric(organoid_data[[col]])
    }
  }
  
  organoid_data
}

select_files_and_add_metadata <- function(select_on = "Analyse", plot_condition = "all", selection_value = 1, metadata_fac=c("chemo_naive", "RASTRIC", "Passed_QC", "Tween_bad", "DMSO_bad", "Clin_benefit", "Analyse", "Analyse2", "OPTIC"), metadata_num=c("RAS_pct_change", "org_id", "TOT"), AUC=FALSE, WGS=FALSE, OPTIC=FALSE, OPTIC_only=FALSE) {
  if(select_on == "all") {
    Analyse <- overview %>% filter(Data_processed == 1)
  } else {
    Analyse <- overview %>% filter(get(select_on) == selection_value)
  }
  status <- Analyse %>% dplyr::select(STR_ID, org_name, org_id, all_of(metadata_fac), all_of(metadata_num), "Year", "Signature_platinum", "Signature_fluor")
  d <- read_plot_data(Analyse[1, "STR_ID"], Analyse[1, "org_name"], AUC=AUC)
  
  # Store column types before joining
  col_types <- sapply(d, class)
  numeric_cols <- names(d)[sapply(d, is.numeric)]
  
  # Read and combine the rest of the data
  for (n in 2:nrow(Analyse)) {
    r <- read_plot_data(Analyse[n, "STR_ID"], Analyse[n, "org_name"], AUC=AUC)
    d <- rbind(d, r)
  }
  d <- left_join(d, status, by=c("STR_ID", "org_name"))
  
  
  # Restore column types after join
  for(col in names(col_types)) {
    if(col %in% names(d)) {
      if(col_types[col] == "numeric") {
        d[[col]] <- suppressWarnings(as.numeric(d[[col]]))
      }
    }
  }
  
  
  # d <- d %>% filter((condition != "Vinorelbine" & condition != "vi_bi_la" & condition != "navitoclax_vinorelbine") | 
  #                     (condition == "Vinorelbine" & Year != 2024)| 
  #                     (condition == "vi_bi_la" & Year != 2024)| 
  #                     (condition == "navitoclax_vinorelbine" & Year != 2024))
  d <- d %>% mutate(across(all_of(metadata_fac), .fns = ~factor(.,levels = c(0,1), labels = c("no", "yes"))))
  if (!plot_condition == "all") {
    d <- d %>% filter(condition == plot_condition)
  } 
  if(WGS==TRUE){
    d$org_name <- sub("_.*", "", d$org_name)
    d <- merge(d, mat_HMF_treatment, by = "org_name") #merge WGS & drugscreen data)
  }else{
    if(OPTIC==TRUE){
    common_cols <- intersect(names(d), names(OPTIC_plotdata))
    d_subset <- d[, common_cols]
    d <- rbind(d_subset, OPTIC_plotdata)
    d <- transform(d, org_name_OPTIC = sapply(strsplit(as.character(d$org_name), "_"), `[`, 1))
    if(OPTIC_only==TRUE){
      d <- d %>% filter(STR_ID == "OPTIC")
    }else{
      d <- d
    }
    }else{
    d <- d
    }
  }
}

plot_per_condition <- function(df, condition, fitcurves = TRUE, viab = FALSE, custom_colors = NULL) {
  if(viab == TRUE) {
    df$mean_GR <- df$mean_viab
  }
  df <- df[df$condition == condition, ]
  title <- condition
  n1 <- length(unlist(unique(df$org_name)))
  hex_codes1 <- if (is.null(custom_colors)) hue_pal()(n1) else custom_colors
  
  # create an empty list to store the nplr models
  fit_list <- list()
  
  # loop over the unique organoids in the condition
  for (organoid in unique(df$org_name)) {
    print(paste("plotting", organoid, "on", condition, "by individual colors"))
    organoid_data = df[df$org_name == organoid, ]
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
  
  # create an empty data frame to store the fitted values
  dataframe_fit <- data.frame()
  
  # loop over the nplr models and extract the fitted values
  for (i in seq_along(fit_list)) {
    fit <- fit_list[[i]]
    
    dataframe_fit_i <- data.frame(getXcurve(fit), getYcurve(fit))
    colnames(dataframe_fit_i) <- c("Concentration_1", "y_fit")
    
    min_organoid = min(df[df$org_name == names(fit_list)[i], ]$mean_GR)
    max_organoid = max(df[df$org_name == names(fit_list)[i], ]$mean_GR)
    diff_organoid = max_organoid - min_organoid
    
    dataframe_fit_i$y_fit_original <- min_organoid + dataframe_fit_i$y_fit * diff_organoid
    
    dataframe_fit_i$condition <- condition
    dataframe_fit_i$org_name <- names(fit_list)[i]
    
    # append the data frame to the main one
    dataframe_fit <- rbind(dataframe_fit, dataframe_fit_i)
  }
  
  # plot all the organoids with different colors based on their org_name
  if(fitcurves == TRUE) {
    p <- ggplot() + 
      geom_point(data = df, aes(conc_condition, mean_GR, color = org_name), size = 2) +
      theme_classic() +
      geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_original, color = org_name)) +
      scale_color_manual(values = hex_codes1) +
      labs(x = expression(paste("Concentration (log"[10], " ", mu, "M)")),
           y = "GR", 
           main = condition, 
           color = "Organoid") +  
      scale_x_log10(limits = c(min(organoid_data$conc_condition), max(organoid_data$conc_condition))) +
      ylim(-1, 1.5) + 
      ggtitle(condition)
  } else {
    p <- ggplot() + 
      geom_point(data = df, aes(conc_condition, mean_GR, color = org_name), size = 2) +
      theme_classic() +
      geom_line(data = df, aes(conc_condition, mean_GR, color = org_name)) +
      scale_color_manual(values = hex_codes1) +
      labs(x = expression(paste("Concentration (log"[10], " ", mu, "M)")),
           y = "GR", 
           main = condition, 
           color = "Organoid") +  
      scale_x_log10(limits = c(min(organoid_data$conc_condition), max(organoid_data$conc_condition))) +
      ylim(-1, 1.5) + 
      ggtitle(condition)    
  }
  return(p)
}

plot_per_condition_factorial <- function(df, condition, colorby="chemo_naive", fitcurves = TRUE,viab = FALSE, custom_colors = NULL) {
  if(viab==TRUE){
    df$mean_GR <- df$mean_viab}
  df <- df[df$condition == condition, ]
  title <- condition
  
  n1 <- length(unlist(unique(df[colorby])))               # Amount of default colors
  hex_codes1 <- if (is.null(custom_colors)) hue_pal()(n1) else custom_colors
  
  # create an empty list to store the nplr models
  fit_list <- list()
  # loop over the unique organoids in the condition
  for (organoid in unique(df$org_name)) {
    print(paste("plotting", organoid, "on", condition, "by", colorby))
    organoid_data = df[df$org_name == organoid, ]
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
  
  # create an empty data frame to store the fitted values
  dataframe_fit <- data.frame()
  
  # loop over the nplr models and extract the fitted values
  for (i in seq_along(fit_list)) {
    
    fit <- fit_list[[i]]
    
    dataframe_fit_i <- data.frame(getXcurve(fit), getYcurve(fit))
    colnames(dataframe_fit_i) <- c("Concentration_1", "y_fit")
    df_filtered <- filter(df, org_name == names(fit_list)[i])
    # selection_var <- df_filtered["chemo_naive"][[1]][1]
    
    
    selection_vec <- df_filtered %>% pull(get(colorby))
    selection_var <- selection_vec[1]
    
    
    min_organoid = min(df[df$org_name == names(fit_list)[i], ]$mean_GR)
    max_organoid = max(df[df$org_name == names(fit_list)[i], ]$mean_GR)
    diff_organoid = max_organoid - min_organoid
    
    dataframe_fit_i$y_fit_original <- min_organoid + dataframe_fit_i$y_fit * diff_organoid
    
    dataframe_fit_i$condition <- condition
    dataframe_fit_i$org_name <- names(fit_list)[i]
    dataframe_fit_i$colorby <- selection_var
    
    # append the data frame to the main one
    dataframe_fit <- rbind(dataframe_fit, dataframe_fit_i)
    
  }
  
  # plot all the organoids with different colors based on their org_name
  if(fitcurves==TRUE){
    p <- ggplot() + 
      geom_point(data = df, aes(conc_condition, mean_GR, group=org_name, color = get(colorby))) +
      geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_original, group=org_name, color = colorby)) +
      scale_color_manual(values = hex_codes1) + # use your custom colors here
      labs (x= "Concentration (log10) μM", y="GR", main = condition, color = colorby) +  
      scale_x_log10(limits = c(min(organoid_data$conc_condition),max(organoid_data$conc_condition))) +
      theme_classic() +
      ylim(-1,1.5) + 
      ggtitle(condition)
  }else{
    p <- ggplot() + 
      geom_point(data = df, aes(conc_condition, mean_GR, group=org_name, color = get(colorby))) +
      geom_line(data = df, aes(conc_condition, mean_GR, group=org_name, color = get(colorby)))+
      labs (x= "Concentration (log10) μM", y="GR", main = condition, color = colorby) +  
      scale_x_log10(limits = c(min(organoid_data$conc_condition),max(organoid_data$conc_condition))) +
      theme_classic() +
      ylim(-1,1.5) + 
      # scale_color_manual(values = c(rainbow(length(unique(df[[colorby]])))), na.value = "lightgrey") +
      scale_color_manual(values = hex_codes1, na.value="lightgrey") + # use your custom colors here
      ggtitle(condition)    
  }
  
  return(p)
}

plot_per_condition_numeric <- function(df, condition, colorby="RAS_pct_change", fitcurves = TRUE, viab = FALSE) {
  if(viab==TRUE){
    df$mean_GR <- df$mean_viab}
  df <- df[df$condition == condition, ]
  title <- condition
  
  # create an empty list to store the nplr models
  fit_list <- list()
  
  # loop over the unique organoids in the condition
  for (organoid in unique(df$org_name)) {
    print(paste("plotting", organoid, "on", condition, "by", colorby))
    organoid_data = df[df$org_name == organoid, ]
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
  
  # create an empty data frame to store the fitted values
  dataframe_fit <- data.frame()
  
  # loop over the nplr models and extract the fitted values
  for (i in seq_along(fit_list)) {
    
    fit <- fit_list[[i]]
    
    dataframe_fit_i <- data.frame(getXcurve(fit), getYcurve(fit))
    colnames(dataframe_fit_i) <- c("Concentration_1", "y_fit")
    df_filtered <- filter(df, org_name == names(fit_list)[i])
    
    selection_vec <- df_filtered %>% pull(get(colorby))
    selection_var <- selection_vec[1]
    
    min_organoid = min(df[df$org_name == names(fit_list)[i], ]$mean_GR)
    max_organoid = max(df[df$org_name == names(fit_list)[i], ]$mean_GR)
    diff_organoid = max_organoid - min_organoid
    
    dataframe_fit_i$y_fit_original <- min_organoid + dataframe_fit_i$y_fit * diff_organoid
    
    dataframe_fit_i$condition <- condition
    dataframe_fit_i$org_name <- names(fit_list)[i]
    dataframe_fit_i$colorby <- selection_var
    
    # append the data frame to the main one
    dataframe_fit <- rbind(dataframe_fit, dataframe_fit_i)
    
  }
  
  # dataframe_fit$color_column <- dataframe_fit[[colorby]]
  # Convert 'colorby' column to numeric
  # df$colorby <- as.numeric(df$colorby)
            
  
  
  # plot all the organoids with different colors based on their org_name
  if(fitcurves==TRUE){
    p <- ggplot() + 
      geom_point(data = df, aes(conc_condition, mean_GR, group=org_name, color = get(colorby))) +
      geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_original, group=org_name, color = colorby)) +
      labs (x= "Concentration (log10) μM", y="GR", main = condition, color = colorby) +  
      scale_color_continuous(low="green", high= "red") +
      scale_x_log10(limits = c(min(organoid_data$conc_condition),max(organoid_data$conc_condition))) +
      theme_classic() +
      ylim(-1,1.5) + 
      ggtitle(condition)
  }else{
    p <- ggplot() + 
      geom_point(data = df, aes(conc_condition, mean_GR, group=org_name, color = get(colorby))) +
      geom_line(data = df, aes(conc_condition, mean_GR, group=org_name, color = get(colorby))) +
      labs (x= "Concentration (log10) μM", y="GR", main = condition, color = colorby) +  
      scale_color_continuous(low="green", high= "red") +
      scale_x_log10(limits = c(min(organoid_data$conc_condition),max(organoid_data$conc_condition))) +
      theme_classic() +
      ylim(-1,1.5) + 
      ggtitle(condition)
  }
  return(p)
}

plot_multiple <- function(exp_name, select_on="Passed_QC", selection_value=1, colorby="chemo_naive", color_numeric=FALSE, fitcurves=TRUE, WGS=FALSE, OPTIC=FALSE, OPTIC_only=TRUE, custom_colors=NULL) {
  if (!dir.exists(file.path(plot_output, exp_name))) {
    dir.create(file.path(plot_output, exp_name))
  }
  d <- select_files_and_add_metadata(select_on = select_on, selection_value = selection_value, WGS=WGS, OPTIC=OPTIC, OPTIC_only = OPTIC_only)
  for (condition in unique(d$condition)) {
    p <- plot_per_condition(d, condition, fitcurves=fitcurves)
    if (color_numeric) {
      q <- plot_per_condition_numeric(d, condition, colorby=colorby, fitcurves=fitcurves)
    } else {
      q <- plot_per_condition_factorial(d, condition, colorby=colorby, fitcurves=fitcurves, custom_colors=custom_colors)
    }
    p_and_q <- p + q
    ggsave(file.path(plot_output, exp_name, paste0(condition,"_", exp_name, "_plots.png")), p_and_q, width=3800, height=1250, units="px")
    ggsave(file.path(plot_output, exp_name, paste0(condition,"_", exp_name, "_plots.pdf")), p_and_q, width=3800, height=1250, units="px")
  }
}

plot_single <- function(exp_name, condition, select_on="Analyse", selection_value=1, colorby="chemo_naive", color_numeric=FALSE, fitcurves=TRUE, WGS=FALSE, OPTIC=FALSE, OPTIC_only=FALSE, custom_colors=NULL) {
  if (!dir.exists(file.path(plot_output, exp_name))) {
    dir.create(file.path(plot_output, exp_name))
  }
  d <- select_files_and_add_metadata(select_on = select_on, plot_condition = condition, selection_value = selection_value, WGS=WGS, OPTIC=OPTIC, OPTIC_only=OPTIC_only)
  p <- plot_per_condition(d, condition, fitcurves=fitcurves, custom_colors = custom_colors)
  if (color_numeric) {
    q <- plot_per_condition_numeric(d, condition, colorby=colorby, fitcurves=fitcurves)
  } else {
    q <- plot_per_condition_factorial(d, condition, colorby=colorby, fitcurves=fitcurves, custom_colors=custom_colors)
  }
  p_and_q <- p + q
  ggsave(file.path(plot_output, exp_name, paste0(condition,"_", exp_name, "_plots.png")), p_and_q, width=3800, height=1250, units="px")
}

plot_selected_conditions <- function(exp_name, select_on="Passed_QC", selection_value = 1, colorby="chemo_naive", selected_conditions=c("5-FU", "Oxaliplatin"), color_numeric=FALSE, fitcurves=TRUE, WGS=WGS, OPTIC=OPTIC, OPTIC_only=FALSE, custom_colors=NULL) {
  for (condition in selected_conditions) {
    plot_single(exp_name = exp_name, condition = condition, select_on=select_on, selection_value=selection_value, colorby=colorby, color_numeric=color_numeric, fitcurves=fitcurves, WGS=WGS, OPTIC=OPTIC, OPTIC_only = OPTIC_only, custom_colors=custom_colors)
  }
}

dont_plot_selected_conditions <- function(exp_name, select_on="Passed_QC", selection_value = 1, colorby="chemo_naive",  unselected_conditions=c("5-FU", "Oxaliplatin"), color_numeric=FALSE, fitcurves=fitcurves, custom_colors=NULL) {
  d <- select_files_and_add_metadata(select_on = select_on, selection_value = selection_value, WGS=WGS, OPTIC=OPTIC, OPTIC_only = OPTIC_only)
  for (condition in unique(d$condition)) {
    if (!condition %in% unselected_conditions) {
      plot_single(exp_name = exp_name, condition = condition, select_on=select_on, selection_value=selection_value, colorby=colorby, color_numeric=color_numeric, fitcurves=fitcurves, custom_colors=custom_colors)
    }
  }
}

# Function to plot with renamed lines
plot_renamed_multiple <- function(exp_name, select_on="Analyse", selection_value=1, rename_vector=NULL, fitcurves=TRUE, WGS=FALSE, OPTIC=FALSE, OPTIC_only=FALSE, custom_colors=NULL) {
  if (!dir.exists(file.path(plot_output, exp_name))) {
    dir.create(file.path(plot_output, exp_name))
  }
  d <- select_files_and_add_metadata(select_on = select_on, selection_value = selection_value, WGS=WGS, OPTIC=OPTIC, OPTIC_only=OPTIC_only)
  
  # Rename organoids if rename_vector is provided
  if (!is.null(rename_vector)) {
    d$org_name <- ifelse(d$org_name %in% names(rename_vector), 
                         rename_vector[d$org_name], 
                         d$org_name)
  }
  
  for (condition in unique(d$condition)) {
    p <- plot_per_condition(d, condition, fitcurves=fitcurves, custom_colors=custom_colors)
    ggsave(file.path(plot_output, exp_name, paste0(condition,"_", exp_name, "_renamed_plot.pdf")), p, width=8, height=6, units="in")
  }
}

# Function to plot with renamed lines
plot_renamed_lines <- function(exp_name, condition, select_on="Analyse", selection_value=1, rename_vector=NULL, fitcurves=TRUE, WGS=FALSE, OPTIC=FALSE, OPTIC_only=FALSE) {
  if (!dir.exists(file.path(plot_output, exp_name))) {
    dir.create(file.path(plot_output, exp_name))
  }
  d <- select_files_and_add_metadata(select_on = select_on, plot_condition = condition, selection_value = selection_value, WGS=WGS, OPTIC=OPTIC, OPTIC_only=OPTIC_only)
  
  # Rename organoids if rename_vector is provided
  if (!is.null(rename_vector)) {
    d$org_name <- ifelse(d$org_name %in% names(rename_vector), 
                         rename_vector[d$org_name], 
                         d$org_name)
  }
  
  p <- plot_per_condition(d, condition, fitcurves=fitcurves)
  
  ggsave(file.path(plot_output, exp_name, paste0(condition,"_", exp_name, "_renamed_plot.pdf")), p, width=10, height=6, units="in")
}

plot_multiple_conditions <- function(exp_name, conditions, organoid, select_on="Analyse2", rename_vector=NULL, fitcurves=TRUE) {
  if (!dir.exists(file.path(plot_output, exp_name))) {
    dir.create(file.path(plot_output, exp_name))
  }
  
  d <- select_files_and_add_metadata(select_on = select_on)
  
  # Filter for the specified organoid and conditions
  d_filtered <- d %>% 
    filter(org_name == organoid, condition %in% conditions)
  
  # Rename organoid if rename_vector is provided
  if (!is.null(rename_vector) && organoid %in% names(rename_vector)) {
    d_filtered$org_name <- rename_vector[organoid]
  }
  
  # Create the base plot
  p <- ggplot(d_filtered, aes(x = conc_condition, y = mean_GR, color = condition)) +
    geom_point(size = 2) +
    theme_classic() +
    labs(x = expression(paste("Concentration (log"[10], " ", mu, "M)")),
         y = "GR", 
         color = "Condition") +
    scale_x_log10() +
    ylim(-1, 1.5) +
    ggtitle(paste("Response of", unique(d_filtered$org_name)))
  
  if (fitcurves) {
    # Fit curves for each condition
    for (cond in conditions) {
      cond_data <- d_filtered %>% filter(condition == cond)
      if (nrow(cond_data) > 0) {
        fit <- nplr(cond_data$conc_condition, convertToProp(cond_data$mean_GR), 
                    useLog = TRUE, LPweight = 0.25, npars = "all", method = "res", silent = TRUE)
        
        x_curve <- getXcurve(fit)
        y_curve <- getYcurve(fit)
        
        # Convert y_curve back to original scale
        y_curve_original <- (max(cond_data$mean_GR) - min(cond_data$mean_GR)) * y_curve + min(cond_data$mean_GR)
        
        fitted_data <- data.frame(x = x_curve, y = y_curve_original, condition = cond)
        
        p <- p + geom_line(data = fitted_data, aes(x = x, y = y, color = condition))
      }
    }
  } else {
    p <- p + geom_line(aes(group = condition))
  }
  
  # Use a color-blind friendly palette
  p <- p + scale_color_brewer(palette = "Set2")
  
  # Save the plot
  ggsave(file.path(plot_output, exp_name, paste0(exp_name, "_plot.pdf")), p, width=8, height=6, units="in")
  
  return(p)
}
overview <- read_excel(file.path(resource_dir, "Screening_overview.xlsx"))


#------import WGS data PDO made in script 9-----
# mat_HMF_treatment <- readRDS(file.path(WGS_plot_dir, "mat_HMF_treatment.rds"))
# 
# View(mat_HMF_treatment)

####---- OPTIC plot data loading ----
# OPTIC_plotdata <- readRDS(file.path(resource_dir, "df_curve.rds"))
# OPTIC_plotdata <- OPTIC_plotdata %>% filter(treatment == "5-FU" | treatment == "Oxaliplatin"| treatment == "SN-38")
# OPTIC_plotdata <- OPTIC_plotdata %>% subset(!(treatment == "5-FU" & Medium == "WithoutNAC"))
# OPTIC_plotdata$org_name_clean <- gsub("OPT379-(\\d{2})(\\d+)(\\d{2})", "OPT\\1\\3", OPTIC_plotdata$Organoid_culture, perl = TRUE)
# STRATEGIC_overview <- overview
# STRATEGIC_overview <- transform(STRATEGIC_overview , org_name_clean = sapply(strsplit(as.character(STRATEGIC_overview $org_name), "_"), `[`, 1))
# OPTIC_plotdata <- OPTIC_plotdata[OPTIC_plotdata$org_name_clean %in% STRATEGIC_overview$org_name_clean, ]
# OPTIC_plotdata$org_name <- paste0(OPTIC_plotdata$org_name_clean, "_OPTIC")
# OPTIC_plotdata <- OPTIC_plotdata %>% rename(conc_condition = Concentration_1)
# OPTIC_plotdata <- OPTIC_plotdata %>% rename(mean_GR = GR_Intensity_mean_run_mean)
# OPTIC_plotdata <- OPTIC_plotdata %>% rename(condition = treatment)
# OPTIC_plotdata <- OPTIC_plotdata[, c("org_name", "condition", "conc_condition", "mean_GR")]
# OPTIC_plotdata <- distinct(OPTIC_plotdata)
# OPTIC_plotdata$STR_ID <- "OPTIC"
# OPTIC_plotdata$chemo_naive <- ifelse(OPTIC_plotdata$org_name=="OPT0067_OPTIC" | 
#                                       OPTIC_plotdata$org_name=="OPT0048_OPTIC" | 
#                                        OPTIC_plotdata$org_name=="OPT0424_OPTIC"| 
#                                        OPTIC_plotdata$org_name=="OPT0028_OPTIC", "no", "yes")
# 
# # plot_multiple("G12C_screens1", select_on="Analyse2", selection_value=2, fitcurves = FALSE)
# five_colors <- c("#00FFFF", "#00FF00", "#FF0000", "#8000FF", "#eeeeee")
# four_colors <- c("#00FFFF", "#00FF00", "#FF0000", "#8000FF")
# three_colors <- c("#00FFFF", "#00FF00", "#FF0000")
# three_colors2 <- c("#00FFFF", "#eeeeee", "#FF0000")
# three_colors3 <- c("#FF0000","#00FFFF", "#eeeeee")

####---- Regular Plotting Function Calls----
# rename_vector_g12c <- c(
#   "OPT0402_G12C" = "Wildtype KRAS",
#   "OPT0421_G12C" = "KRASG12C - untreated",
#   "RAS21_G12C" = "KRASG12C - refractory"
# )
# 
# plot_multiple("G12C_screen_3only", select_on="Analyse2", fitcurves = TRUE)
# plot_renamed_multiple("Renamed_G12C_screen", rename_vector = rename_vector_g12c, select_on="Analyse2", fitcurves = TRUE)
# 
# # Define the rename vector
# rename_vector_g12c <- c("OPT0421_G12C" = "KRASG12C - Untreated")
# 
# # 1. Plot single agent G12C inhibitors
# plot_multiple_conditions(
#   exp_name = "G12C_single_agents",
#   conditions = c("Sotorasib", "Adagrasib", "Divarasib"),
#   organoid = "OPT0421_G12C",
#   select_on = "Analyse2",
#   rename_vector = rename_vector_g12c,
#   fitcurves = FALSE
# )
# 
# # 2. Plot combinations with Lapatinib
# plot_multiple_conditions(
#   exp_name = "G12C_lapatinib_combos",
#   conditions = c("Sotorasib_Lapatinib", "Adagrasib_Lapatinib", "Divarasib_Lapatinib"),
#   organoid = "OPT0421_G12C",
#   select_on = "Analyse2",
#   rename_vector = rename_vector_g12c,
#   fitcurves = FALSE
# )
# 
# # 3. Plot combinations with Vinorelbine
# plot_multiple_conditions(
#   exp_name = "G12C_vinorelbine_combos",
#   conditions = c("Sotorasib_Lapatinib_Vinorelbine", "Adagrasib_Lapatinib_Vinorelbine", "Divarasib_Lapatinib_Vinorelbine"),
#   organoid = "OPT0421_G12C",
#   select_on = "Analyse2",
#   rename_vector = rename_vector_g12c,
#   fitcurves = FALSE
# )


overview <- read_excel(file.path(resource_dir, "Screening_overview.xlsx"))
# overview <- overview %>% 
#   mutate(Analyse = ifelse(STR_ID == "MHUI08", 1, 0))

plot_multiple("RASTRIC_by_response_pct", select_on="Analyse", fitcurves = FALSE, WGS=FALSE, colorby="RAS_pct_change", color_numeric = TRUE)
# plot_multiple("RASTRIC_by_TOT", select_on="Analyse", fitcurves = FALSE, WGS=FALSE, colorby="TOT", color_numeric = TRUE)

# plot_single("MHUI01", condition ="Vinorelbine", select_on="Analyse", fitcurves=FALSE)
# 
# plot_renamed_lines()
#plot_multiple("SMAD4", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "SMAD4")
#plot_multiple("adjuvant", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "Adjuvant_chemo_naive")
#plot_selected_conditions("response_SOC_chemo_naive", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "response_5FU_chemo_naive", selected_conditions = c("5-FU"))
#plot_multiple("Pairs", select_on="Pairs", fitcurves = FALSE, colorby = "org_id")
#plot_selected_conditions("OPTIC", select_on="OPTIC", fitcurves = FALSE, WGS=FALSE, OPTIC=TRUE, OPTIC_only=FALSE, colorby = "org_name_OPTIC", selected_conditions = c("5-FU"))
#plot_selected_conditions("OPTIC_chemonaive", select_on="OPTIC", fitcurves = FALSE, WGS=FALSE, OPTIC=TRUE, OPTIC_only=TRUE, colorby = "chemo_naive", selected_conditions = c("5-FU", "Oxaliplatin", "SN-38"))
#plot_selected_conditions("response_SOC_chemo_naive1", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "response_5FU_chemo_naive", selected_conditions = c("5-FU"), custom_colors=four_colors)
#plot_selected_conditions("response_SOC_chemo_naive3", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "response_5FU_chemo_naive_adjuvant_selection", selected_conditions = c("5-FU"), custom_colors=three_colors2)
#plot_selected_conditions("Signature", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "DBS5_cat", color_numeric=FALSE, selected_conditions = c("Oxaliplatin"), custom_colors = NULL)
# plot_multiple("Signature_SBS17_pretreated", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "SBS17_cat", color_numeric=FALSE, custom_colors = three_colors3)
# plot_multiple("Signature_SBSplatinum_pretreated", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "SBSplatinum_cat", color_numeric=FALSE, custom_colors = three_colors3)
# plot_multiple("Signature_SBS17a_pretreated", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "SBS17a_cat", color_numeric=FALSE, custom_colors = three_colors3)
# plot_multiple("Signature_SBS17threshold_pretreated", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "SBS17_threshold", color_numeric=FALSE, custom_colors = three_colors3)

#plot_multiple("viability", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE)
#plot_multiple("Time2", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "Year", color_numeric=FALSE, custom_colors=four_colors)
# 
# mat_HMF_treatment$org_name
# overv