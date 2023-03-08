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

Analyse <- overview %>% filter(Analyse == 1)
read_plot_data <- function(exp_id, organoid_name) {
  
  exp_file_name <- list.files(file.path(plot_dir), pattern=paste(exp_id, organoid_name, sep="_"))
    
  if (length(exp_file_name) == 0) {
    return(paste("ERROR: Combination of", exp_id, "and", organoid_name, "does not exist in this folder!"))
  }
  print(paste("Reading", exp_id, organoid_name))
  organoid_data <- read_excel(file.path(plot_dir, exp_file_name))
  
  organoid_data
}

read_experiment <- function() {
  Analyse <- overview %>% filter(Analyse == 1)
  status <- Analyse %>% dplyr::select(STR_ID, org_name, chemo_naive, RASTRIC)
  d <- read_plot_data(Analyse[1, "STR_ID"], Analyse[1, "org_name"])
  for (n in 2:nrow(Analyse)) {
    r <- read_plot_data(Analyse[n, "STR_ID"], Analyse[n, "org_name"])
    d <- rbind(d, r)
  }
  d <- left_join(d, status, by=c("STR_ID", "org_name"))
  d$chemo_naive <- as.factor(d$chemo_naive)
  levels(d$chemo_naive) <- c("no", "yes")
  d$RASTRIC <- as.factor(d$RASTRIC)
  levels(d$RASTRIC) <- c("no", "yes")
  d
}
d <- read_experiment()

# organoids = unique(d$org_name)
# conditions = unique(d$condition)

plot_per_condition <- function(df, condition) {
  df <- df[df$condition == condition, ]
  title <- condition
  
  n1 <- length(unlist(unique(d$org_name)))               # Amount of default colors
  hex_codes1 <- hue_pal()(n1)                             # Identify hex codes
  
  # create an empty list to store the nplr models
  fit_list <- list()
  
  # loop over the unique organoids in the condition
  for (organoid in unique(df$org_name)) {
    print(paste("plotting", organoid))
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
  p <- ggplot() + 
    geom_point(data = df, aes(conc_condition, mean_GR, color = org_name), size = 2) +
    theme_classic() +
    geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_original, color = org_name)) +
    scale_color_manual(values = hex_codes1) + # use your custom colors here
    labs (x= "Concentration (log10) uM", y="GR", main = condition, color = "Organoid") +  
    scale_x_log10(limits = c(min(organoid_data$conc_condition),max(organoid_data$conc_condition))) +
    ylim(-1,1.5) + ggtitle(condition)
  return(p)
}

plot_per_condition_factorial <- function(df, condition, colorby="chemo_naive") {
  df <- df[df$condition == condition, ]
  title <- condition
  
  n1 <- length(unlist(unique(d[colorby])))               # Amount of default colors
  hex_codes1 <- hue_pal()(n1)                             # Identify hex codes
  
  # create an empty list to store the nplr models
  fit_list <- list()
  
  # loop over the unique organoids in the condition
  for (organoid in unique(df$org_name)) {
    print(paste("plotting", organoid))
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
    selection_var <- df_filtered["chemo_naive"][[1]][1]
    
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
  p <- ggplot() + 
    geom_point(data = df, aes(conc_condition, mean_GR, group=org_name, color = !!sym(colorby)), size = 2) +
    geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_original, group=org_name, color = colorby)) +
    scale_color_manual(values = hex_codes1) + # use your custom colors here
    labs (x= "Concentration (log10) uM", y="GR", main = condition, color = colorby) +  
    scale_x_log10(limits = c(min(organoid_data$conc_condition),max(organoid_data$conc_condition))) +
    theme_classic() +
    ylim(-1,1.5) + 
    ggtitle(condition)
  return(p)
}

plot_multiple <- function(exp_name) {
  for (condition in unique(d$condition)) {
    p <- plot_per_condition(d, condition)
    q <- plot_per_condition_factorial(d, condition, colorby="chemo_naive")
    p_and_q <- p + q
    ggsave(file.path(plot_output, paste0(condition,"_", exp_name, "_plots.png")), p_and_q, width=4200, height=2100, units="px")
  }
}




# plot_organoid_with_data <- function(df, condition, n=1, save=F) {
#   df <- df[df$condition == condition, ]
#   title <- condition
#   
#   organoid = organoids[n]
#   print(paste("plotting", organoid))
#   organoid_data = df[df$org_name == organoid, ]
#   organoid_data$Max_Concentration_log <- log10(organoid_data$conc_condition)
#   organoid_data$organoid_color <- hex_codes1[n]
#   # Add drug concentration 0 
#   # conc_zero <- organoid_data[nrow(organoid_data) + 1,]
#   # conc_zero$org_name <- organoid
#   # conc_zero$mean_GR <- 1
#   # conc_zero$condition <- title
#   # conc_zero$Concentration <- min(organoid_data$Max_Concentration, na.rm = TRUE)
#   # conc_zero$Max_Concentration <- (min(organoid_data$Max_Concentration, na.rm = TRUE)/10)
#   # conc_zero$Max_Concentration_log <- (min(organoid_data$Max_Concentration_log) -1)
#   # organoid_data[nrow(organoid_data) + 1,] = conc_zero
#   
#   min_organoid = min(organoid_data$mean_GR)
#   max_organoid = max(organoid_data$mean_GR)
#   diff_organoid = max_organoid - min_organoid
#   organoid_data$GR_prop <- convertToProp(organoid_data$mean_GR)
#   
#   fit <- nplr(organoid_data$conc_condition, organoid_data$GR_prop, 
#               useLog = TRUE, #should conc values be log10 transformed
#               LPweight = 0.25, 
#               npars = "all", #number to specify the number of parameters to use in the model; "all" --> the logistic model will be tested with 2 to 5 parameters, best option will be returned
#               method = "res", #model optimized using sum squared errors (using weighting method, different options)
#               silent = FALSE) #should warning messages be silenced
#   
#  
#   dataframe_fit <- data.frame(getXcurve(fit), getYcurve(fit))
#   colnames(dataframe_fit) <- c("Concentration_1", "y_fit")
#   dataframe_fit$y_fit_original <- min_organoid + dataframe_fit$y_fit * diff_organoid
#   # dataframe_fit$y_fit<- ((dataframe_fit$y_fit)*2)-1
#   
#   organoid_data$organoid_color <- organoids[n]
#   dataframe_fit$organoid_color <- organoids[n]
#   
#   #metrics
#   
#   # organoid_data_AUC$GR_positive <- (((organoid_data_AUC$x)+1)/2) #GR range 0-1 for calculating AUC
#   # GRMax <- mean(subset(organoid_data, Max_Concentration == max(organoid_data$Max_Concentration))$GR)
#   # AUC_raw <- AUC(organoid_data_AUC$Group.1, organoid_data_AUC$GR_positive) #AUC raw
#   # AUC_fit <- getAUC(fit)
#   # AUC_fit_trapezoid <- AUC_fit$trapezoid #AUC fit
#   # AUC_fit_simpson <- AUC_fit$Simpson #AUC fit
#   # estim <- getEstimates(fit, .5)
#   # GR50_log <- format(estim$x, digits = 6, scientific = TRUE) #estimate at 0.5
#   # GR50_num <- as.numeric(GR50_log)
#   # GR50 <- log10(GR50_num)
#   # par <- getPar(fit) 
#   # xmid_log <- (par$params$xmid)
#   # xmid<- as.numeric(xmid_log)
#   
#   # DR_summary_local[nrow(DR_summary_local) + 1,]  <-  list(condition = condition, organoid = organoid, AUC_raw = AUC_raw, AUC_fit_trapezoid = AUC_fit_trapezoid, GRMax = GRMax, GR50 = GR50)
#   
#   #checkGR<- getEstimates(fit, (1-((1-GRMax)*0.5)))
#   #check50_log <- format(esGR$x, digits = 6, scientific = TRUE) 
#   #checkGR50_num <- as.numeric(GR50_log)
#   #checkGR50 <- log10(GR50_num)
#   if (n==1) {
#     p <- ggplot() + 
#       geom_point(data = organoid_data, aes(conc_condition, mean_GR, color = organoid_color), size = 2) +
#       theme_classic() +
#       geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_original, color = organoid_color)) +
#       #geom_vline(xintercept = checkGR50, color = hex_codes1[1]) +
#       labs (x= expression(paste("Concentration (log10) ", mu, "M")), y="GR", main = condition, colour = "Organoid") +  
#       scale_x_log10(limits = c(min(organoid_data$conc_condition),max(organoid_data$conc_condition))) +
#       ylim(-1,1.5) + ggtitle(condition)
#   } else {
#     p <- p +
#       geom_point(data = organoid_data, aes(conc_condition, mean_GR, color = organoid_color), size = 2) +
#       geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_original, color = organoid_color)) 
#   }
#   p
# }