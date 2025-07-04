# Script voor vergelijking van organoids OPT0032_2, OPT0067 en RAS34_3
# Versie 1.0 (2024-12-19)
# Doel: Vergelijken van drie specifieke organoids met curves en metrics

#----- Libraries -----
library(dplyr)
library(readxl)
library(ggplot2)
library(openxlsx)
library(scales)   
library(drc)
library(nplr)
library(patchwork)
library(DescTools) #AUC
library(cowplot)
library(ggpubr) #statcomparemeans
library(pheatmap)
library(tidyr) #spread

rm(list=ls())

#----- Directory setup -----
# script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
script_dir <- "/Users/mhuismans/surfdrive/Promotie/12. STRATEGIC/Analyse 2.0/script"

home_dir <- dirname(script_dir)
raw_dir <- file.path(home_dir, "1_raw_files")
org_data_dir <- file.path(home_dir, "2_organoid_data")
QC_dir <- file.path(home_dir, "3_QC")
resource_dir <- file.path(home_dir, "resources")
plot_dir <- file.path(home_dir, "4_plot_data")
plot_output <- file.path(home_dir, "5_plot_output")
metrics_dir <- file.path(home_dir, "6_metrics")
metrics_normalized_dir <- file.path(metrics_dir, "240607_QC_includes_RAS05")

#----- Load overview data -----
overview <- read_excel(file.path(resource_dir, "Screening_overview.xlsx"))

#----- Functions from Script3 -----
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

#----- Define organoids to compare -----
organoids_to_compare <- c("OPT0032_2", "OPT0067", "RAS34_3")

#----- Define organoid labels for plotting -----
organoid_labels <- c(
  "OPT0032_2" = "Metastasis organoid before SOC treatment",
  "OPT0067" = "Metastasis organoid during SOC treatment", 
  "RAS34_3" = "Metastasis organoid after SOC treatment"
)

#----- Get data for selected organoids -----
get_organoid_data <- function(organoid_names) {
  # Filter overview for selected organoids
  selected_orgs <- overview %>% 
    filter(org_name %in% organoid_names)
  
  if(nrow(selected_orgs) == 0) {
    stop("No data found for the selected organoids")
  }
  
  # Read data for each organoid
  all_data <- data.frame()
  
  for(i in 1:nrow(selected_orgs)) {
    org_data <- read_plot_data(selected_orgs$STR_ID[i], selected_orgs$org_name[i])
    if(is.character(org_data)) {
      warning(paste("Could not read data for:", selected_orgs$STR_ID[i], selected_orgs$org_name[i]))
      next
    }
    
    # Add metadata
    org_data$org_id <- selected_orgs$org_id[i]
    org_data$chemo_naive <- selected_orgs$chemo_naive[i]
    org_data$RASTRIC <- selected_orgs$RASTRIC[i]
    
    all_data <- rbind(all_data, org_data)
  }
  
  return(all_data)
}

#----- Plot functions -----
plot_organoid_comparison <- function(df, condition, fitcurves = TRUE, custom_colors = NULL) {
  
  # Filter for specific condition
  df_condition <- df[df$condition == condition, ]
  
  if(nrow(df_condition) == 0) {
    warning(paste("No data found for condition:", condition))
    return(NULL)
  }
  
  # Apply organoid labels
  df_condition$org_name_label <- organoid_labels[df_condition$org_name]
  
  n_orgs <- length(unique(df_condition$org_name))
  hex_codes <- if (is.null(custom_colors)) hue_pal()(n_orgs) else custom_colors
  
  # Fit curves if requested
  if(fitcurves) {
    fit_list <- list()
    dataframe_fit <- data.frame()
    
    for (organoid in unique(df_condition$org_name)) {
      print(paste("Fitting curve for", organoid, "on", condition))
      organoid_data = df_condition[df_condition$org_name == organoid, ]
      
      if(nrow(organoid_data) < 3) {
        warning(paste("Insufficient data points for", organoid, "on", condition))
        next
      }
      
      organoid_data$Max_Concentration_log <- log10(organoid_data$conc_condition)
      
      min_organoid = min(organoid_data$mean_GR, na.rm = TRUE)
      max_organoid = max(organoid_data$mean_GR, na.rm = TRUE)
      diff_organoid = max_organoid - min_organoid
      organoid_data$GR_prop <- convertToProp(organoid_data$mean_GR)
      
      # Fit the nplr model
      tryCatch({
        fit <- nplr(organoid_data$conc_condition, organoid_data$GR_prop, 
                    useLog = TRUE,
                    LPweight = 0.25,
                    npars = "all",
                    method = "res",
                    silent = TRUE)
        
        fit_list[[organoid]] <- fit
        
        # Extract fitted values
        dataframe_fit_i <- data.frame(getXcurve(fit), getYcurve(fit))
        colnames(dataframe_fit_i) <- c("Concentration_1", "y_fit")
        
        dataframe_fit_i$y_fit_original <- min_organoid + dataframe_fit_i$y_fit * diff_organoid
        dataframe_fit_i$condition <- condition
        dataframe_fit_i$org_name <- organoid
        dataframe_fit_i$org_name_label <- organoid_labels[organoid]
        
        dataframe_fit <- rbind(dataframe_fit, dataframe_fit_i)
        
      }, error = function(e) {
        warning(paste("Could not fit curve for", organoid, "on", condition, ":", e$message))
      })
    }
  }
  
  # Create plot
  if(fitcurves && nrow(dataframe_fit) > 0) {
    # Zorg dat de factor levels van org_name_label in de juiste volgorde staan
    df_condition$org_name_label <- factor(df_condition$org_name_label, levels = organoid_labels[organoids_to_compare])
    dataframe_fit$org_name_label <- factor(dataframe_fit$org_name_label, levels = organoid_labels[organoids_to_compare])
    p <- ggplot() + 
      geom_point(data = df_condition, aes(conc_condition, mean_GR, color = org_name_label), size = 3) +
      theme_classic() +
      geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_original, color = org_name_label), size = 1.2) +
      scale_color_manual(values = hex_codes, labels = organoid_labels[organoids_to_compare]) +
      labs(x = expression(paste("Concentration (log"[10], " ", mu, "M)")),
           y = "GR", 
           title = paste("Drug response curves:", condition), 
           color = "Organoid") +  
      scale_x_log10() +
      ylim(-1, 1.5) + 
      theme(plot.title = element_text(hjust = 0.5, size = 14),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10))
  } else {
    df_condition$org_name_label <- factor(df_condition$org_name_label, levels = organoid_labels[organoids_to_compare])
    p <- ggplot() + 
      geom_point(data = df_condition, aes(conc_condition, mean_GR, color = org_name_label), size = 3) +
      theme_classic() +
      geom_line(data = df_condition, aes(conc_condition, mean_GR, color = org_name_label), size = 1.2) +
      scale_color_manual(values = hex_codes, labels = organoid_labels[organoids_to_compare]) +
      labs(x = expression(paste("Concentration (log"[10], " ", mu, "M)")),
           y = "GR", 
           title = paste("Drug response curves:", condition), 
           color = "Organoid") +  
      scale_x_log10() +
      ylim(-1, 1.5) + 
      theme(plot.title = element_text(hjust = 0.5, size = 14),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10))
  }
  
  return(p)
}

#----- Plot SOC curves los, maar samen in één PDF -----
plot_soc_curves_grid <- function(df, soc_conditions = c("5-FU", "Oxaliplatin", "SN-38"), custom_colors = NULL) {
  plots <- list()
  for (cond in soc_conditions) {
    p <- plot_organoid_comparison(df, cond, fitcurves = TRUE, custom_colors = custom_colors)
    if (!is.null(p)) {
      plots[[cond]] <- p
    }
  }
  if (length(plots) > 0) {
    combined <- wrap_plots(plots, ncol = 1)
    return(combined)
  } else {
    return(NULL)
  }
}

#----- Metrics: boxplot met alle organoids grijs, drie interesse gekleurd+gelabeld -----
plot_metrics_highlight <- function(metrics_data, selected_conditions = c("5-FU", "Oxaliplatin", "SN-38"), custom_colors = NULL) {
  all_metrics <- read_excel(file.path(metrics_normalized_dir, "240607_QC_includes_RAS05_metrics_normalized.xlsx"))
  metrics_to_plot <- c("AUC", "GR50")
  plot_list <- list()
  
  # Labels met nieuwe regel
  organoid_labels_newline <- c(
    "OPT0032_2" = "Metastasis organoid before\nSOC treatment",
    "OPT0067" = "Metastasis organoid during\nSOC treatment",
    "RAS34_3" = "Metastasis organoid after\nSOC treatment"
  )
  
  # Define colors if not provided - gebruik de originele organoid namen voor de kleuren
  if(is.null(custom_colors)) {
    custom_colors <- setNames(hue_pal()(3), organoids_to_compare)
  }
  
  for (metric in metrics_to_plot) {
    if (metric %in% colnames(all_metrics)) {
      for (cond in selected_conditions) {
        # Data voor alle organoids
        all_cond <- all_metrics %>% filter(condition == cond)
        # Data voor de drie van interesse
        interest_cond <- all_cond %>% filter(org_name %in% organoids_to_compare)
        # Voeg labels toe
        interest_cond$org_name_label <- organoid_labels_newline[interest_cond$org_name]
        # Rest grijs
        all_cond$highlight <- ifelse(all_cond$org_name %in% organoids_to_compare, "highlight", "background")
        all_cond$org_name_label <- ifelse(all_cond$org_name %in% organoids_to_compare, organoid_labels_newline[all_cond$org_name], NA)
        # Plot
        p <- ggplot(all_cond, aes(x = "", y = !!sym(metric))) +
          geom_boxplot(outlier.shape = NA, fill = "grey90", color = "grey60") +
          geom_jitter(data = all_cond[all_cond$highlight == "background",], aes(x = "", y = !!sym(metric)), color = "grey", width = 0.15, size = 2, alpha = 0.5) +
          geom_point(data = interest_cond, aes(x = "", y = !!sym(metric), color = org_name_label), size = 4) +
          scale_color_manual(values = setNames(custom_colors, organoid_labels_newline[organoids_to_compare]), na.value = "grey") +
          geom_text(data = interest_cond, aes(x = 1.15, y = !!sym(metric), label = org_name_label, color = org_name_label), hjust = 0, size = 3.5, show.legend = FALSE) +
          labs(title = paste(metric, "-", cond), x = NULL, y = metric, color = "Organoid") +
          theme_classic() +
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_text(size = 10), axis.text = element_text(size = 9), legend.position = "none")
        plot_list[[paste(metric, cond, sep = "_")]] <- p
      }
    }
  }
  return(plot_list)
}

#----- Load metrics data -----
load_metrics_data <- function() {
  metrics_file <- file.path(metrics_normalized_dir, "240607_QC_includes_RAS05_metrics_normalized.xlsx")
  
  if(!file.exists(metrics_file)) {
    stop("Metrics file not found:", metrics_file)
  }
  
  # Read metrics data
  metrics_data <- read_excel(metrics_file)
  
  # Filter for our organoids
  metrics_filtered <- metrics_data %>% 
    filter(org_name %in% organoids_to_compare)
  
  return(metrics_filtered)
}

#----- Main execution -----
main <- function() {
  print("Starting organoid comparison analysis...")
  organoid_data <- get_organoid_data(organoids_to_compare)
  if(nrow(organoid_data) == 0) stop("No data loaded for the selected organoids")
  print(paste("Loaded data for", length(unique(organoid_data$org_name)), "organoids"))
  print(paste("Available conditions:", paste(unique(organoid_data$condition), collapse = ", ")))
  output_dir <- file.path(plot_output, "organoid_comparison_OPT0032_2_OPT0067_RAS34_3")
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Kleuren consistent houden
  organism_colors <- setNames(hue_pal()(3), organoid_labels[organoids_to_compare])

  # SOC curves los, samen in 1 PDF
  soc_grid <- plot_soc_curves_grid(organoid_data, custom_colors = organism_colors)
  if(!is.null(soc_grid)) {
    ggsave(file.path(output_dir, "SOC_curves_grid.pdf"), soc_grid, width = 8, height = 12)
  }

  # Metrics highlight
  metrics_data <- load_metrics_data()
  metrics_plots <- plot_metrics_highlight(metrics_data, custom_colors = organism_colors)
  if(length(metrics_plots) > 0) {
    for(plot_name in names(metrics_plots)) {
      ggsave(file.path(output_dir, paste0("metrics_highlight_", plot_name, ".pdf")), metrics_plots[[plot_name]], width = 6, height = 6)
    }
    combined_metrics <- wrap_plots(metrics_plots, ncol = 3)
    ggsave(file.path(output_dir, "combined_metrics_highlight.pdf"), combined_metrics, width = 18, height = 12)
  }

  # Overige output (optioneel)
  summary_data <- organoid_data %>%
    group_by(org_name, condition) %>%
    summarise(n_points = n(), mean_GR = mean(mean_GR, na.rm = TRUE), sd_GR = sd(mean_GR, na.rm = TRUE), min_GR = min(mean_GR, na.rm = TRUE), max_GR = max(mean_GR, na.rm = TRUE), .groups = 'drop')
  write.csv(summary_data, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)
  print(paste("Analysis complete! Results saved to:", output_dir))
  return(list(organoid_data = organoid_data))
}

#----- Run analysis -----
results <- main() 
