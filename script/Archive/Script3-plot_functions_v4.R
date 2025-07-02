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
  organoid_data
}

select_files_and_add_metadata <- function(select_on = "Analyse", plot_condition = "all", selection_value = 1, metadata_fac=c("chemo_naive", "RASTRIC", "Passed_QC", "Tween_bad", "DMSO_bad", "Clin_benefit", "Analyse", "Analyse2", "OPTIC"), metadata_num=c("RAS_pct_change", "org_id"), AUC=FALSE, WGS=FALSE, OPTIC=FALSE, OPTIC_only=FALSE) {
  if(select_on == "all") {
    Analyse <- overview %>% filter(Data_processed == 1)
  } else {
    Analyse <- overview %>% filter(get(select_on) == selection_value)
  }
  status <- Analyse %>% dplyr::select(STR_ID, org_name, org_id, all_of(metadata_fac), all_of(metadata_num), "Year", "Signature_platinum", "Signature_fluor")
  d <- read_plot_data(Analyse[1, "STR_ID"], Analyse[1, "org_name"], AUC=AUC)
  for (n in 2:nrow(Analyse)) {
    r <- read_plot_data(Analyse[n, "STR_ID"], Analyse[n, "org_name"], AUC=AUC)
    d <- rbind(d, r)
    print(typeof(d$mean_GR))
  }
  d <- left_join(d, status, by=c("STR_ID", "org_name"))
  d <- d %>% mutate(across(all_of(metadata_fac), .fns = ~factor(.,levels = c(0,1), labels = c("no", "yes"))))
  if (!plot_condition == "all") {
    d <- d %>% filter(condition == plot_condition)
  } 
  if(WGS==TRUE){
    d$org_name_clean <- sub("_.*", "", d$org_name)
    d <- merge(d, mat_HMF_treatment, by = "org_name_clean") #merge WGS & drugscreen data)
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

plot_per_condition <- function(df, condition, fitcurves = TRUE, viab = FALSE, custom_colors=NULL) {
  if(viab==TRUE){
    df$mean_GR <- df$mean_viab}
  df <- df[df$condition == condition, ]
  title <- condition
  n1 <- length(unlist(unique(df$org_name)))               # Amount of default colors
  hex_codes1 <- hue_pal()(n1)  # Identify hex codes
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
  if(fitcurves==TRUE){
  p <- ggplot() + 
    geom_point(data = df, aes(conc_condition, mean_GR, color = org_name), size = 2) +
    theme_classic() +
    geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_original, color = org_name)) +
    scale_color_manual(values = hex_codes1) + # use your custom colors here
    labs (x= "Concentration (log10) μM", y="GR", main = condition, color = "Organoid") +  
    scale_x_log10(limits = c(min(organoid_data$conc_condition),max(organoid_data$conc_condition))) +
    ylim(-1,1.5) + ggtitle(condition)
  }else{
    p <- ggplot() + 
      geom_point(data = df, aes(conc_condition, mean_GR, color = org_name), size = 2) +
      theme_classic() +
      geom_line(data = df, aes(conc_condition, mean_GR, color = org_name)) +
      scale_color_manual(values = hex_codes1) + # use your custom colors here
      labs (x= "Concentration (log10) μM", y="GR", main = condition, color = "Organoid") +  
      scale_x_log10(limits = c(min(organoid_data$conc_condition),max(organoid_data$conc_condition))) +
      ylim(-1,1.5) + ggtitle(condition)    
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

overview <- read_excel(file.path(resource_dir, "Screening_overview.xlsx"))


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

mat_HMF_treatment$response_5FU <- ifelse(mat_HMF_treatment$Palliative_5FU_response != "PD" & mat_HMF_treatment$`Sample type` == "pretreated", "no PD",
                                         ifelse(mat_HMF_treatment$Palliative_5FU_response == "No" & mat_HMF_treatment$`Sample type` == "chemonaive", "resistant naive", NA))
mat_HMF_treatment$response_SN38 <- ifelse(mat_HMF_treatment$Palliative_irinotecan_response != "PD" & mat_HMF_treatment$`Sample type` == "pretreated", "no PD",
                                        ifelse(mat_HMF_treatment$Palliative_irinotecan_response == "No" & mat_HMF_treatment$`Sample type` == "chemonaive", "resistant naive", NA))
mat_HMF_treatment$response_oxaliplatin<- ifelse(mat_HMF_treatment$Palliative_oxaliplatin_response != "PD" & mat_HMF_treatment$`Sample type` == "pretreated", "no PD",
                                        ifelse(mat_HMF_treatment$Palliative_oxaliplatin_response == "No" & mat_HMF_treatment$`Sample type` == "chemonaive", "resistant naive", NA))
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

####---- OPTIC plot data loading ----
OPTIC_plotdata <- readRDS(file.path(resource_dir, "df_curve.rds"))
OPTIC_plotdata <- OPTIC_plotdata %>% filter(treatment == "5-FU" | treatment == "Oxaliplatin"| treatment == "SN-38")
OPTIC_plotdata <- OPTIC_plotdata %>% subset(!(treatment == "5-FU" & Medium == "WithoutNAC"))
OPTIC_plotdata$org_name_clean <- gsub("OPT379-(\\d{2})(\\d+)(\\d{2})", "OPT\\1\\3", OPTIC_plotdata$Organoid_culture, perl = TRUE)
STRATEGIC_overview <- overview
STRATEGIC_overview <- transform(STRATEGIC_overview , org_name_clean = sapply(strsplit(as.character(STRATEGIC_overview $org_name), "_"), `[`, 1))
OPTIC_plotdata <- OPTIC_plotdata[OPTIC_plotdata$org_name_clean %in% STRATEGIC_overview$org_name_clean, ]
OPTIC_plotdata$org_name <- paste0(OPTIC_plotdata$org_name_clean, "_OPTIC")
OPTIC_plotdata <- OPTIC_plotdata %>% rename(conc_condition = Concentration_1)
OPTIC_plotdata <- OPTIC_plotdata %>% rename(mean_GR = GR_Intensity_mean_run_mean)
OPTIC_plotdata <- OPTIC_plotdata %>% rename(condition = treatment)
OPTIC_plotdata <- OPTIC_plotdata[, c("org_name", "condition", "conc_condition", "mean_GR")]
OPTIC_plotdata <- distinct(OPTIC_plotdata)
OPTIC_plotdata$STR_ID <- "OPTIC"
OPTIC_plotdata$chemo_naive <- ifelse(OPTIC_plotdata$org_name=="OPT0067_OPTIC" | 
                                      OPTIC_plotdata$org_name=="OPT0048_OPTIC" | 
                                       OPTIC_plotdata$org_name=="OPT0424_OPTIC"| 
                                       OPTIC_plotdata$org_name=="OPT0028_OPTIC", "no", "yes")

####---- Regular Plotting Function Calls----
#plot_multiple("PassedQC7_nonfit", select_on="Passed_QC", fitcurves = FALSE)
#plot_multiple("SMAD4", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "SMAD4")
#plot_multiple("adjuvant", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "Adjuvant_chemo_naive")
# plot_selected_conditions("response_SOC_chemo_naive", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "response_5FU_chemo_naive", selected_conditions = c("5-FU"))
# #plot_selected_conditions("response_SOC_chemo_naive", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "response_oxaliplatin_chemo_naive", selected_conditions = c("Oxaliplatin"))
# #plot_selected_conditions("response_SOC_chemo_naive", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "response_SN-38_chemo_naive", selected_conditions = c("SN-38"))
# #plot_multiple("Pairs", select_on="Pairs", fitcurves = FALSE, colorby = "org_id")
# plot_selected_conditions("OPTIC", select_on="OPTIC", fitcurves = FALSE, WGS=FALSE, OPTIC=TRUE, OPTIC_only=FALSE, colorby = "org_name_OPTIC", selected_conditions = c("5-FU"))
# plot_selected_conditions("OPTIC_chemonaive", select_on="OPTIC", fitcurves = FALSE, WGS=FALSE, OPTIC=TRUE, OPTIC_only=TRUE, colorby = "chemo_naive", selected_conditions = c("5-FU", "Oxaliplatin", "SN-38"))

# plot_multiple("G12C_screens1", select_on="Analyse2", selection_value=2, fitcurves = FALSE)
five_colors <- c("#00FFFF", "#00FF00", "#FF0000", "#8000FF", "#eeeeee")
four_colors <- c("#00FFFF", "#00FF00", "#FF0000", "#8000FF")
three_colors <- c("#00FFFF", "#00FF00", "#FF0000")
three_colors2 <- c("#00FFFF", "#eeeeee", "#FF0000")

#plot_selected_conditions("response_SOC_chemo_naive1", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "response_5FU_chemo_naive", selected_conditions = c("5-FU"), custom_colors=four_colors)
#plot_selected_conditions("response_SOC_chemo_naive3", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "response_5FU_chemo_naive_adjuvant_selection", selected_conditions = c("5-FU"), custom_colors=three_colors2)
#plot_selected_conditions("Signature", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "DBS5_cat", color_numeric=FALSE, selected_conditions = c("Oxaliplatin"), custom_colors = NULL)
#plot_multiple("Signature_SBS17a_all", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE, colorby = "SBS17a_cat", color_numeric=FALSE, custom_colors = NULL)
#plot_multiple("viability", select_on="Passed_QC", fitcurves = FALSE, WGS=TRUE)
