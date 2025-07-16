### VERSION CONTROL
# Script11-RASTRIC_sensitivity_heatmap.R
# Adapted from Script10-individual_data_to_figure1_heatmap_v5.R
# 
# Purpose: Create heatmap showing organoid sensitivity to RASTRIC study combination (vi_bi_la) 
# and its components (Vinorelbine, Binimetinib, Lapatinib), arranged by sensitivity to vi_bi_la
# 
# Key modifications:
# - Focus on RASTRIC combination and components only
# - Arrange organoids by sensitivity to vi_bi_la
# - Show pretreatment status as color annotation (not split)
# - Show organoid names (RASxx, OPTxxxx)
# - No column splits between treatment types

library(tidyverse)
library(readxl)
library(ggplot2)
library(openxlsx)
library(scales)   
library(drc)
library(nplr)
library(DescTools) #AUC
library(cowplot)
library(ggpubr) #statcomparemeans
library(ComplexHeatmap)
library(circlize) # For color functions
library(lubridate) # For date handling

rm(list=ls())

# script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
script_dir <- "/Users/mhuismans/surfdrive/Promotie/12. STRATEGIC/Analyse 2.0/script/"
home_dir <- dirname(script_dir)
raw_dir <- file.path(home_dir, "1_raw_files")
org_data_dir <- file.path(home_dir, "2_organoid_data")
QC_dir <- file.path(home_dir, "3_QC")
resource_dir <- file.path(home_dir, "resources")
plot_dir <- file.path(home_dir, "4_plot_data")
plot_output <- file.path(home_dir, "5_plot_output")
metrics_dir <- file.path(home_dir, "6_metrics")
wgs_dir <- file.path(home_dir, "7_WGS")
metrics_plot_dir <- file.path(metrics_dir, "Plots metrics")
metrics_normalized_dir <- file.path(metrics_dir, "240607_QC_includes_RAS05")

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
    Analyse <- overview %>% dplyr::filter(Data_processed == 1) 
  } else {
    Analyse <- overview %>% dplyr::filter(get(select_on) == selection_value)
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

normalize_df <- function(df) {
  metrics_df_norm <- df
  metrics_df_norm = metrics_df_norm %>% dplyr::select(org_name,STR_ID,condition,chemo_naive,AUC,GRmax,AUC_low,GRmax_low,GR50,xmid)
  colnames(metrics_df_norm)
  
  metrics_df_norm <- metrics_df_norm %>%
    group_by(condition) %>%
    mutate_at(vars(5:last_col()), list(~(. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))))
  return(metrics_df_norm)
}

save_metrics <- function(exp_name, 
                      select_on="Passed_QC", 
                     selection_value = 1, 
                     selected_conditions="all",
                     normalize=TRUE) {

  if (!dir.exists(file.path(metrics_dir, exp_name))) {
    dir.create(file.path(metrics_dir, exp_name))
  }
    # create an empty data frame to store the fitted values
  dataframe_conditions <- data.frame()
  
  if (selected_conditions == "all") {
    selected_conditions <- c("5-FU", "vi_bi_la", "Oxaliplatin", "Lapatinib", "Binimetinib", "Alpelisib", "Navitoclax", "CHEK1", "SN-38" ,
          "alpelisib_lapatinib",  "binimetinib_lapatinib", "SN38_CHEK1", "alpelisib_lapatinib", "Vinorelbine", "binimetinib_alpelisib", "navitoclax_vinorelbine")
  }  
  for (condition in selected_conditions) {
    d <- select_files_and_add_metadata(select_on = select_on, selection_value = selection_value, plot_condition = condition)
    d_individual <- select_files_and_add_metadata(select_on = select_on, selection_value = selection_value, individual_reps=TRUE)

    # create an empty list to store the nplr models
    fit_list <- list()
    
    for (organoid in unique(d$org_name)) {
      
      organoid_data = d[d$org_name == organoid, ]
      
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
    
    dataframe_fit <- data.frame()
    dataframe_fit_metrics <- data.frame()
    
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
      
      #get metrics from fit and create dataframe
      estim <- getEstimates(fit, .5)
      GR50_log <- format(estim$x, digits = 6, scientific = TRUE) #estimate at 0.5
      GR50_num <- as.numeric(GR50_log)
      GR50 <- log10(GR50_num)
      par <- getPar(fit) 
      xmid_log <- (par$params$xmid)
      xmid<- as.numeric(xmid_log)
      
      dataframe_fit_m <- data.frame(GR50, xmid)
      dataframe_fit_m$org_name <- names(fit_list)[i]
      dataframe_fit_metrics <- rbind(dataframe_fit_metrics, dataframe_fit_m)
    }

    concentration_filters <- list(
      `5-FU` = c(489.5, 490.1, 489.8, 108.2, 108.1,108.4, 24.29, 24.36, 23.94, 5.316, 5.296, 5.297), 
      Oxaliplatin = c(0.1, 0.10100, 0.09937, 0.39540, 0.39450, 0.3947, 1.57000, 1.53300, 1.514, 6.13800, 6.132, 6.1930), 
      `SN-38` =  c(0.0069300, 0.020790, 0.064350, 0.1980000), 
      SN38_CHEK1 =  c(0.0069300, 0.020790, 0.064350, 0.1980000),
      CHEK1 = c(0.038610, 0.118800, 0.336600, 0.999900),
      Alpelisib = c(0.57420, 1.88100, 6.13800, 20.00000),
      Lapatinib = c(0.57420, 1.88100, 6.13800, 20.00000),
      Binimetinib = c(0.17330,0.57420, 1.88100, 6.13800), 
      lapatinib_binimetinib = c(0.17330,0.57420, 1.88100, 6.13800),
      binimetinib_alpelisib = c(0.17330,0.57420, 1.88100, 6.13800), 
      lapatinib_alpelisib = c(0.57420, 1.88100, 6.13800, 20.00000),
      Vinorelbine = c(0.001485, 0.003713, 0.010400, 0.028710),
      navitoclax_vinorelbine = c(0.01980,  0.07920,  0.32670, 1.28700),
      Navitoclax = c(0.07920, 0.32670, 1.28700, 5.24700),
      vi_bi_la = c(0.001485, 0.003713, 0.010400, 0.028710)
    )
    
    concentration_filters2 <- list(
      Binimetinib = c(0.00495, 0.01485, 0.05445, 0.1733, 0.5742, 1.881),
      vi_bi_la = c(0.000495, 0.001485, 0.003713, 0.010400, 0.028710, 0.0792)
    )
    
    filtered_data <- d %>%
      group_by(condition) %>%
      filter(conc_condition %in% concentration_filters[[unique(condition)]]) %>%
      ungroup()
    
    d <- d %>%
      group_by(condition) %>%
      anti_join(data.frame(conc_condition = unlist(concentration_filters2))) %>%
      ungroup()
    
    #fill df with metrics
    metric <- d %>% dplyr::select(conc_condition, mean_GR, org_name, STR_ID, condition, chemo_naive) %>%
      mutate(mean_GR_positive=(mean_GR +1)/2)  %>%
      group_by(org_name, STR_ID, condition) %>%
      mutate(AUC = AUC(log10(conc_condition), mean_GR_positive)) %>%
      group_by(org_name, condition) %>%
      mutate(max_index = which.max(conc_condition),
             GRmax = mean_GR[max_index]) %>%
      dplyr::select(-max_index) %>%
      group_by(org_name, STR_ID, condition, AUC, GRmax) %>% 
      summarise(AUC = mean(AUC))

      #dataframe for metrics only low concentrations
    metric_lowconc <- filtered_data %>% dplyr::select(conc_condition, mean_GR, org_name, STR_ID, condition, chemo_naive) %>%
      mutate(mean_GR_positive=(mean_GR +1)/2)  %>%
      group_by(org_name, condition, chemo_naive) %>%
      mutate(AUC_low = AUC(log10(conc_condition), mean_GR_positive)) %>% 
      dplyr::select(-condition) %>%
      mutate(max_index = which.max(conc_condition),
             GRmax_low = mean_GR[max_index]) %>%
      dplyr::select(-max_index) %>%
      group_by(org_name, AUC_low, GRmax_low, chemo_naive) %>%
      summarise(AUC_low = mean(AUC_low))

    #merge dataframe non-fitted metrics with fitted metrics
    metric_metriclowconc <- merge(metric, metric_lowconc, by="org_name") 
    
    #merge dataframe non-fitted metrics with fitted metrics
    dataframe_metrics <- merge(metric_metriclowconc, dataframe_fit_metrics, by="org_name") 
    
    #merge dataframe conditions
    dataframe_conditions <- rbind(dataframe_conditions, dataframe_metrics)
    
  }
  
  # save the file
  write.xlsx(dataframe_conditions, file.path(metrics_dir, exp_name, paste0(exp_name, "_metrics.xlsx")))
  
  # create and save the normalized file if desired
  if(normalize) {
    df_return <- normalize_df(dataframe_conditions)
    write.xlsx(df_return, file.path(metrics_dir, exp_name, paste0(exp_name, "_metrics_normalized.xlsx")))
  } else {
    df_return <- dataframe_conditions
  }
  # returns normalized df if Normalize is True (standard parameter in function)  
  return(df_return)
}

# Function to check annotation completeness
check_annotation_completeness <- function(data, annotation_cols = c("chemo_naive", "origin_location", "MAPK_status", "Sidedness", "OS")) {
  
  cat("\n=== ANNOTATION COMPLETENESS CHECK ===\n")
  
  # Get unique organoids
  unique_orgs <- data %>% distinct(org_id, org_name)
  
  # Check each annotation column
  completeness_report <- data.frame()
  
  for (col in annotation_cols) {
    if (col %in% colnames(data)) {
      # Get unique values for this annotation
      unique_vals <- data %>% 
        distinct(org_id, org_name, !!sym(col)) %>%
        pull(!!sym(col))
      
      # Count missing values
      missing_count <- sum(is.na(unique_vals) | unique_vals == "Unknown" | unique_vals == "")
      total_count <- length(unique_vals)
      
      # Find organoids with missing annotations
      missing_orgs <- data %>%
        distinct(org_id, org_name, !!sym(col)) %>%
        filter(is.na(!!sym(col)) | !!sym(col) == "Unknown" | !!sym(col) == "") %>%
        pull(org_name)
      
      completeness_report <- rbind(completeness_report, data.frame(
        Annotation = col,
        Total_Organoids = total_count,
        Missing_Count = missing_count,
        Completeness_Percent = round((total_count - missing_count) / total_count * 100, 1),
        Missing_Organoids = paste(missing_orgs, collapse = ", ")
      ))
      
      cat(sprintf("%s: %d/%d complete (%.1f%%)\n", col, total_count - missing_count, total_count, 
                  (total_count - missing_count) / total_count * 100))
      
      if (length(missing_orgs) > 0) {
        cat("  Missing for: ", paste(missing_orgs, collapse = ", "), "\n")
      }
    } else {
      cat(sprintf("%s: Column not found in data\n", col))
    }
  }
  
  # Overall completeness
  total_orgs <- nrow(unique_orgs)
  cat(sprintf("\nTotal unique organoids: %d\n", total_orgs))
  
  # Check if any organoids have completely missing annotations
  org_annotation_status <- data %>%
    distinct(org_id, org_name, .keep_all = TRUE) %>%
    mutate(
      missing_annotations = rowSums(is.na(dplyr::select(., all_of(annotation_cols))) | 
                                   dplyr::select(., all_of(annotation_cols)) == "Unknown" | 
                                   dplyr::select(., all_of(annotation_cols)) == "")
    )
  
  orgs_with_missing <- org_annotation_status %>%
    filter(missing_annotations > 0) %>%
    arrange(desc(missing_annotations))
  
  if (nrow(orgs_with_missing) > 0) {
    cat("\nOrganoids with missing annotations:\n")
    for (i in 1:nrow(orgs_with_missing)) {
      org <- orgs_with_missing[i, ]
      cat(sprintf("  %s (%s): %d missing annotations\n", 
                  org$org_name, org$org_id, org$missing_annotations))
    }
  } else {
    cat("\nAll organoids have complete annotations!\n")
  }
  
  return(completeness_report)
}

# Function to create RASTRIC sensitivity heatmap
make_RASTRIC_sensitivity_heatmap <- function(exp_name, data, excluded_organoids = NULL, output_filename = "RASTRIC_sensitivity_heatmap") {
  
  # Define RASTRIC combination and components
  rastric_conditions <- c("vi_bi_la", "Vinorelbine", "Binimetinib", "Lapatinib")
  
  # Filter out excluded organoids if specified
  if (!is.null(excluded_organoids) && length(excluded_organoids) > 0) {
    cat("Filtering out excluded organoids:", paste(excluded_organoids, collapse = ", "), "\n")
    data <- data %>% filter(!org_name %in% excluded_organoids)
    cat("Remaining organoids after exclusion:", nrow(data %>% distinct(org_id)), "\n")
  }
  
  # Prepare data for RASTRIC conditions only, and keep annotation variables
  data_select <- data %>% 
    filter(condition %in% rastric_conditions) %>%
    dplyr::select(org_id, org_name, condition, AUC, chemo_naive, origin_location, MAPK_status, Sidedness, OS) %>% 
    group_by(org_id, condition) %>%
    summarise(
      AUC = mean(AUC, na.rm = TRUE),
      org_name = first(org_name),
      chemo_naive = first(chemo_naive),
      origin_location = first(origin_location),
      MAPK_status = first(MAPK_status),
      Sidedness = first(Sidedness),
      OS = first(OS),
      .groups = 'drop'
    )

  # Add OS_category as in the official script
  data_select <- data_select %>%
    mutate(OS_category = cut(OS,
      breaks = c(-Inf, 12, 24, 36, 48, 60, Inf),
      labels = c("Less than 1 year",
                 "Between 1 and 2 years",
                 "Between 2 and 3 years",
                 "Between 3 and 4 years",
                 "Between 4 and 5 years",
                 "More than 5 years"),
      right = FALSE))

  # Reshape data for heatmap - separate chemo_naive from the matrix
  data_reshape <- tidyr::spread(data_select,
                                key = condition,
                                value = AUC)

  # Remove duplicated org_id rows if any
  if(any(duplicated(data_reshape$org_id))) {
    warning("Duplicated org_id rows found after spreading. Keeping only the first occurrence for each.")
    data_reshape <- data_reshape[!duplicated(data_reshape$org_id), ]
  }
  
  # Remove rows with NA org_id and set row names
  data_reshape <- data_reshape %>% filter(!is.na(org_id))
  rownames(data_reshape) <- data_reshape$org_id

  # Separate chemo_naive and annotation variables from the matrix data
  annotation_vars <- data_reshape %>% dplyr::select(chemo_naive, origin_location, MAPK_status, Sidedness, OS_category)
  # Rename vi_bi_la column to RASTRIC for display purposes
  colnames(data_reshape)[colnames(data_reshape) == "vi_bi_la"] <- "RASTRIC"
  # Reorder columns: RASTRIC, Vinorelbine, Lapatinib, Binimetinib
  condition_cols <- c("RASTRIC", "Vinorelbine", "Lapatinib", "Binimetinib")
  data_matrix_only <- data_reshape[, condition_cols]
  # Ensure all columns are numeric
  data_matrix_only <- data.matrix(data_matrix_only)
  data_matrix <- as.matrix(data_matrix_only)

  # Sort organoids by sensitivity to RASTRIC (ascending order = most sensitive first)
  vi_bi_la_order <- order(data_matrix[, "RASTRIC"], na.last = TRUE)
  data_matrix_sorted <- data_matrix[vi_bi_la_order, ]
  annotation_sorted <- annotation_vars[vi_bi_la_order, ]

  # Rename vi_bi_la column to RASTRIC for display purposes
  colnames(data_matrix_sorted)[colnames(data_matrix_sorted) == "vi_bi_la"] <- "RASTRIC"
  
  # Reorder columns: RASTRIC, Vinorelbine, Lapatinib, Binimetinib
  data_matrix_sorted <- data_matrix_sorted[, c("RASTRIC", "Vinorelbine", "Lapatinib", "Binimetinib")]

  # Set up color schemes as in the official script with fallback colors
  row_anno_colors <- list(
    `Organoid before SOC treatment` = c("yes" = "#2ECC71", "no" = "#E74C3C", "Unknown" = "lightgrey"),
    `Origin of biopsy for organoid` = c("Liver" = "#D35400", "Skin" = "#FAD7A0", "Lymph node" = "#2980B9", "Lung" = "#A3E4D7", 
                                       "Soft tissue" = "#8B4513", "Abdominal wall" = "#CD853F", "Primary" = "#FF6347", 
                                       "Pleural" = "#20B2AA", "Duodenum" = "#9370DB", "Unknown" = "lightgrey"),
    `MAPK mutational status` = c("RAS/BRAF-wildtype" = "#DCF763", "RAS-mutant" = "#8E44AD", "KRAS wildtype amplification" = "#E67E22", "Unknown" = "lightgrey"),
    Sidedness = c("Right-sided (coecum-transverse colon)" = "#E74C3C", 
                  "Left-sided (splenic flexure-sigmoid)" = "#2ECC71", 
                  "Rectum (rectosigmoid/rectal)" = "#F39C12", 
                  "Unknown" = "lightgrey"),
    `Overall Survival` = c(
      "Not reached" = "white",
      "Less than 1 year" = "#7B241C",         
      "Between 1 and 2 years" = "#C0392B",    
      "Between 2 and 3 years" = "#E74C3C",    
      "Between 3 and 4 years" = "#EC7063",    
      "Between 4 and 5 years" = "#F5B7B1",    
      "More than 5 years" = "#FADBD8",
      "Unknown" = "lightgrey")
  )

  # Debug: Print unique values for each annotation column
  cat("\nDebug annotation values:\n")
  cat("chemo_naive:", paste(unique(annotation_sorted$chemo_naive), collapse=", "), "\n")
  cat("origin_location:", paste(unique(annotation_sorted$origin_location), collapse=", "), "\n")
  cat("MAPK_status:", paste(unique(annotation_sorted$MAPK_status), collapse=", "), "\n")
  cat("Sidedness:", paste(unique(annotation_sorted$Sidedness), collapse=", "), "\n")
  cat("OS_category:", paste(unique(annotation_sorted$OS_category), collapse=", "), "\n")
  
  # Check row names
  cat("Row names check:\n")
  cat("Number of rows:", nrow(data_matrix_sorted), "\n")
  cat("Row names set:", !is.null(rownames(data_matrix_sorted)), "\n")
  cat("Sample row names:", head(rownames(data_matrix_sorted), 5), "\n")

  # Prepare annotation data frame for ComplexHeatmap
  row_annotation <- annotation_sorted %>%
    dplyr::rename(
      `Organoid before SOC treatment` = chemo_naive,
      `Origin of biopsy for organoid` = origin_location,
      `MAPK mutational status` = MAPK_status,
      `Overall Survival` = OS_category
    )

  # Create row annotation
  ra <- rowAnnotation(
    df = row_annotation,
    col = row_anno_colors,
    annotation_name_rot = 25,
    annotation_legend_param = list(border = TRUE),
    simple_anno_size = unit(1, "cm")
  )

  # Create directory if it doesn't exist
  if (!dir.exists(file.path(metrics_dir, exp_name))) {
    dir.create(file.path(metrics_dir, exp_name), recursive = TRUE)
  }

  # Ensure row names are properly set
  if(is.null(rownames(data_matrix_sorted))) {
    warning("Row names are NULL, setting them to organoid IDs")
    rownames(data_matrix_sorted) <- rownames(data_reshape)[vi_bi_la_order]
  }
  
  # Create heatmap with dynamic color scaling based on actual data range
  h <- Heatmap(data_matrix_sorted,
               name = "AUC",
               left_annotation = ra,
               col = colorRamp2(c(min(data_matrix_sorted, na.rm = TRUE), 
                                  mean(c(min(data_matrix_sorted, na.rm = TRUE), 
                                         max(data_matrix_sorted, na.rm = TRUE))), 
                                max(data_matrix_sorted, na.rm = TRUE)), 
                                c("blue", "white", "red")),
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_row_names = TRUE,
               row_names_side = "left",
               row_names_gp = gpar(fontsize = 8),  # Smaller font size
               show_column_names = TRUE,
               column_names_rot = 25,  # Same angle as annotations
               column_names_side = "bottom",
               border = FALSE,
               heatmap_legend_param = list(title = "AUC", 
                                           at = c(min(data_matrix_sorted, na.rm = TRUE), 
                                                 max(data_matrix_sorted, na.rm = TRUE)), 
                                           labels = c("sensitive", "resistant")),
               width = unit(4, "cm"),  # Same width as annotations
               height = unit(12, "cm"))

  # Save the heatmap
  pdf(file.path(metrics_dir, exp_name, paste0(output_filename, ".pdf")), 
      width = 12, height = 10)
  draw(h)
  dev.off()

  return(h)
}

# Main execution
# Set experiment name
exp_name_current <- "250714_RASTRIC_sensitivity"

# Load data
table1 <- read_excel(file.path(home_dir, "STRATEGIC_table1", "Combined_table1.xlsx"))

# Load the RNA database for legacy sample annotations (as in previous scripts)
RNA <- readRDS(file.path(wgs_dir, "PDO_RNA_screen.rds"))

# Load the RASTRIC metrics data (contains all organoids)
rastric_metrics <- read_excel(file.path(metrics_dir, "250714_STRATEGIC_AND_RASTRIC_QC", "250714_STRATEGIC_AND_RASTRIC_QC_metrics.xlsx"))

# Load screening overview to identify RASTRIC organoids
overview <- read_excel(file.path(resource_dir, "Screening_overview.xlsx"))

# Get all organoids that have Passed_QC_RASTRIC = 1
rastric_organoids <- overview %>%
  filter(Passed_QC_RASTRIC == 1) %>%
  dplyr::select(STR_ID, org_id, org_name) %>%
  distinct(org_id, .keep_all = TRUE)

# Debug: Check if target organoids are in rastric_organoids
cat("\n=== CHECKING RASTRIC ORGANOIDS FILTER ===\n")
target_orgs <- c("RAS40_mhui2", "RAS42_mhui4", "RAS44_mhui5", "RAS45_mhui4", "RAS53_mhui")
for(org in target_orgs) {
  found <- rastric_organoids %>% filter(org_name == org)
  if(nrow(found) > 0) {
    cat(sprintf("  %s: Found in rastric_organoids (org_id=%s)\n", org, found$org_id[1]))
  } else {
    cat(sprintf("  %s: NOT FOUND in rastric_organoids\n", org))
  }
}

# Load RASTRIC clinical data for annotation
rastric_clinical <- read_delim(file.path(resource_dir, "Gegevens_RASTRIC_uptodate.csv"), delim=",")
rastric_mutations <- read_delim(file.path(resource_dir, "RASTRIC_RAS_mutaties_volledig.csv"), delim=";")

# Process RASTRIC clinical data to match table1 format
rastric_annotations <- rastric_clinical %>%
  mutate(
    organoid_name = str_replace(SUBJID, "RAS0", "RAS"),
    # Use localization directly as sidedness (already in correct format)
    Sidedness = localization,
    # Use Location directly as origin_location (already in correct format)
    origin_location = Location,
    # Calculate age from birth year and registration date (regisdat is already a date)
    birthday = as.Date(dob), 
    Age = as.integer(floor(interval(birthday, regisdat) / years(1))), 
    Sex = ifelse(sex == 1, "Male", "Female"),
    # Calculate OS from first palliative therapy to death
    first_pall_date = as.POSIXct(first_Pall_therapy_date),
    death_date = dmy(overlijden),
    OS = ifelse(!is.na(first_pall_date) & !is.na(death_date), 
                as.integer(floor(interval(first_pall_date, death_date) / months(1))), 
                NA_real_),
    # All RASTRIC patients are pretreated (not chemo-naive)
    chemo_naive = "no"
  ) %>% 
  # Join with mutation data to get MAPK status
  left_join(rastric_mutations, by="SUBJID") %>%
  mutate(
    # Determine MAPK status based on RAS mutation type
    MAPK_status = case_when(
      ras_mutation_type == "KRAS" ~ "RAS-mutant",
      ras_mutation_type == "NRAS" ~ "RAS-mutant", 
      is.na(ras_mutation_type) ~ "RAS-mutant", # Default for RASTRIC patients
      TRUE ~ "RAS-mutant"
    ),
    # Map origin_location to expected annotation values
    origin_location = case_when(
      str_detect(origin_location, "Liver") ~ "Liver",
      str_detect(origin_location, "Skin") & !str_detect(origin_location, "only") ~ "Skin",
      str_detect(origin_location, "Lymph node") ~ "Lymph node",
      str_detect(origin_location, "Lung") ~ "Lung",
      str_detect(origin_location, "Soft tissue") ~ "Soft tissue",
      str_detect(origin_location, "Abdominal wall") ~ "Abdominal wall",
      str_detect(origin_location, "Primary tumour") ~ "Primary",
      str_detect(origin_location, "Pleural") ~ "Pleural",
      str_detect(origin_location, "Duodenum") ~ "Duodenum",
      str_detect(origin_location, "skin only") ~ NA_character_,
      TRUE ~ origin_location
    )
  ) %>%
  {cat("Unique mapped origin_location values in rastric_annotations:\n"); print(unique(.$origin_location)); .} %>%
  dplyr::select(organoid_name, Age, Sex, Sidedness, origin_location, MAPK_status, chemo_naive, OS)

# Filter RASTRIC metrics data for RASTRIC conditions only
rastric_conditions <- c("vi_bi_la", "Vinorelbine", "Binimetinib", "Lapatinib")

# Get mapping from org_name to org_id from overview
org_mapping <- overview %>%
  dplyr::select(org_id, org_name) %>%
  distinct(org_id, .keep_all = TRUE)

# --- Use org_id for annotation joins ---
rastric_data <- rastric_metrics %>%
  filter(condition %in% rastric_conditions) %>%
  dplyr::select(org_name, condition, AUC, STR_ID) %>%
  left_join(org_mapping, by="org_name") %>%
  left_join(rastric_annotations, by=c("org_id" = "organoid_name"))

# Combine data using the RNA database approach (as in previous scripts)
fig1_df_legacy <- RNA %>% 
  left_join(table1, by="organoid_name") %>%
  mutate(
    Palliative_duration.y = case_when(
      sampleId == "MetPret25" ~ 25,
      .default = Palliative_duration.y)
  )

# Add org_id mapping to legacy data
fig1_df_legacy <- fig1_df_legacy %>%
  left_join(org_mapping, by="org_name")

# Filter legacy data for RASTRIC conditions only and select matching columns
fig1_df_legacy <- fig1_df_legacy %>%
  filter(condition %in% rastric_conditions) %>%
  dplyr::select(org_id, org_name, condition, AUC, chemo_naive, origin_location, MAPK_status, Sidedness, OS)

# Ensure RASTRIC data has the same column structure
rastric_data_aligned <- rastric_data %>%
  dplyr::select(org_id, org_name, condition, AUC, chemo_naive, origin_location, MAPK_status, Sidedness, OS)

# Combine legacy and RASTRIC data
fig1_df <- rbind(fig1_df_legacy, rastric_data_aligned)

# Ensure we have unique organoid-condition combinations using org_id
fig1_df <- fig1_df %>%
  distinct(org_id, condition, AUC, .keep_all = TRUE)

# EXCLUDE SPECIFIC ORGANOIDS
excluded_organoids <- c("OPT0067", "HUB-02-C2-89", "OPT0426", "OPT0426_Demi")
cat("\n=== EXCLUDING ORGANOIDS ===\n")
cat("Excluding organoids:", paste(excluded_organoids, collapse = ", "), "\n")

# Check if these organoids exist in the data
existing_excluded <- fig1_df %>%
  filter(org_name %in% excluded_organoids) %>%
  distinct(org_name, org_id)

if (nrow(existing_excluded) > 0) {
  cat("Found organoids to exclude:\n")
  print(existing_excluded)
  
  # Remove excluded organoids
  fig1_df <- fig1_df %>%
    filter(!org_name %in% excluded_organoids)
  
  cat("Removed", nrow(existing_excluded), "organoids from analysis\n")
} else {
  cat("No excluded organoids found in the data\n")
}

# CHECK ANNOTATION COMPLETENESS BEFORE CREATING HEATMAP
cat("\n=== PRE-HEATMAP ANNOTATION CHECK ===\n")
annotation_check <- check_annotation_completeness(fig1_df)

# Also check annotation completeness for organoids that will be in the final heatmap
cat("\n=== ANNOTATION CHECK FOR FINAL HEATMAP ORGANOIDS ===\n")
rastric_conditions <- c("vi_bi_la", "Vinorelbine", "Binimetinib", "Lapatinib")
final_heatmap_data <- fig1_df %>%
  filter(condition %in% rastric_conditions) %>%
  filter(!org_name %in% excluded_organoids) %>%
  distinct(org_id, org_name, .keep_all = TRUE)

final_annotation_check <- check_annotation_completeness(final_heatmap_data)

# Debug information
cat("Columns in table1:", paste(colnames(table1), collapse=", "), "\n")
cat("Columns in rastric_data:", paste(colnames(rastric_data), collapse=", "), "\n")
cat("Columns in fig1_df:", paste(colnames(fig1_df), collapse=", "), "\n")

# Debug: Check data ranges and organoid counts
cat("\nData validation before heatmap creation:\n")
cat("Total organoids in fig1_df:", nrow(fig1_df %>% distinct(org_id)), "\n")
cat("Organoids with vi_bi_la data:", nrow(fig1_df %>% filter(condition == "vi_bi_la") %>% distinct(org_id)), "\n")
cat("AUC range in fig1_df:", range(fig1_df$AUC, na.rm = TRUE), "\n")
cat("Sample organoid IDs:", head(fig1_df$org_id, 10), "\n")
cat("Sample organoid names:", head(fig1_df$org_name, 10), "\n")

# Create RASTRIC sensitivity heatmap
rastric_heatmap <- make_RASTRIC_sensitivity_heatmap(exp_name_current, fig1_df, excluded_organoids = excluded_organoids, "RASTRIC_sensitivity_heatmap")

# Print summary statistics
cat("RASTRIC sensitivity heatmap created successfully!\n")
cat("Number of organoids included:", nrow(fig1_df %>% filter(condition == "vi_bi_la") %>% distinct(org_name)), "\n")
cat("Number of RASTRIC organoids:", nrow(rastric_organoids), "\n")
cat("RASTRIC organoids:", paste(rastric_organoids$org_name, collapse=", "), "\n")
cat("Output saved to:", file.path(metrics_dir, exp_name_current), "\n")

# Debug information
cat("\nDebug information:\n")
cat("Columns in rastric_metrics:", paste(colnames(rastric_metrics), collapse=", "), "\n")
cat("Columns in rastric_annotations:", paste(colnames(rastric_annotations), collapse=", "), "\n")
cat("Number of RASTRIC conditions found:", nrow(rastric_metrics %>% filter(condition %in% rastric_conditions)), "\n")

# Additional RASTRIC data validation
cat("\nRASTRIC data validation:\n")
cat("Number of RASTRIC patients with OS data:", sum(!is.na(rastric_annotations$OS)), "\n")
cat("RASTRIC sidedness distribution:\n")
print(table(rastric_annotations$Sidedness))
cat("RASTRIC origin locations:\n")
print(table(rastric_annotations$origin_location))
cat("RASTRIC MAPK status distribution:\n")
print(table(rastric_annotations$MAPK_status))

# DEBUG: Trace specific RASTRIC organoids
cat("\n=== DEBUGGING SPECIFIC RASTRIC ORGANOIDS ===\n")
target_rastric <- c("RAS40", "RAS42", "RAS44", "RAS45", "RAS53")

cat("\n1. Checking rastric_annotations (clinical data):\n")
for (org in target_rastric) {
  org_data <- rastric_annotations %>% filter(organoid_name == org)
  if (nrow(org_data) > 0) {
    cat(sprintf("  %s: Found with %d rows\n", org, nrow(org_data)))
  } else {
    cat(sprintf("  %s: NOT FOUND\n", org))
  }
}

cat("\n2. Checking rastric_metrics (drug screen data):\n")
for (org in target_rastric) {
  org_data <- rastric_metrics %>% filter(org_name == org)
  if (nrow(org_data) > 0) {
    cat(sprintf("  %s: Found with %d rows, conditions: %s\n", 
                org, nrow(org_data), 
                paste(unique(org_data$condition), collapse=", ")))
  } else {
    cat(sprintf("  %s: NOT FOUND\n", org))
  }
}

cat("\n3. Checking overview (screening data):\n")
for (org in target_rastric) {
  org_data <- overview %>% filter(org_name == org)
  if (nrow(org_data) > 0) {
    cat(sprintf("  %s: Found with org_id=%s, Passed_QC_RASTRIC=%s\n", 
                org, org_data$org_id[1], org_data$Passed_QC_RASTRIC[1]))
  } else {
    cat(sprintf("  %s: NOT FOUND\n", org))
  }
}

cat("\n4. Checking rastric_data (after joining):\n")
for (org in target_rastric) {
  org_data <- rastric_data %>% filter(org_name == org)
  if (nrow(org_data) > 0) {
    cat(sprintf("  %s: Found with %d rows\n", org, nrow(org_data)))
    # Check if annotations are present
    missing_annos <- sum(is.na(org_data$chemo_naive) | org_data$chemo_naive == "Unknown")
    cat(sprintf("    Missing annotations: %d/%d\n", missing_annos, nrow(org_data)))
  } else {
    cat(sprintf("  %s: NOT FOUND\n", org))
  }
}

cat("\n5. Checking final fig1_df:\n")
for (org in target_rastric) {
  org_data <- fig1_df %>% filter(org_name == org)
  if (nrow(org_data) > 0) {
    cat(sprintf("  %s: Found with %d rows\n", org, nrow(org_data)))
  } else {
    cat(sprintf("  %s: NOT FOUND\n", org))
  }
}

cat("\n6. Checking final heatmap data:\n")
final_data <- fig1_df %>%
  filter(condition %in% rastric_conditions) %>%
  filter(!org_name %in% excluded_organoids) %>%
  distinct(org_id, org_name, .keep_all = TRUE)

for (org in target_rastric) {
  org_data <- final_data %>% filter(org_name == org)
  if (nrow(org_data) > 0) {
    cat(sprintf("  %s: Found with %d rows\n", org, nrow(org_data)))
  } else {
    cat(sprintf("  %s: NOT FOUND\n", org))
  }
}

cat("\n7. Checking RNA dataset (legacy data):\n")
for (org in target_rastric) {
  org_data <- RNA %>% filter(organoid_name == org)
  if (nrow(org_data) > 0) {
    cat(sprintf("  %s: Found with %d rows, conditions: %s\n", 
                org, nrow(org_data), 
                paste(unique(org_data$condition), collapse=", ")))
  } else {
    cat(sprintf("  %s: NOT FOUND\n", org))
  }
}

cat("\n8. Checking which data source each organoid comes from:\n")
for (org in target_rastric) {
  # Check if in legacy data
  legacy_data <- fig1_df_legacy %>% filter(org_name == org)
  rastric_data_check <- rastric_data_aligned %>% filter(org_name == org)
  
  if (nrow(legacy_data) > 0) {
    cat(sprintf("  %s: From LEGACY data (%d rows)\n", org, nrow(legacy_data)))
  } else if (nrow(rastric_data_check) > 0) {
    cat(sprintf("  %s: From RASTRIC data (%d rows)\n", org, nrow(rastric_data_check)))
  } else {
    cat(sprintf("  %s: NOT FOUND in either source\n", org))
  }
} 