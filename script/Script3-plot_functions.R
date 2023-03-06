library(dplyr)
library(readxl)
library(ggplot2)
library(openxlsx)
library(scales)   
library(drc)
library(nplr)

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

passed_QC <- overview %>% filter(Passed_QC == 1)
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
  passed_QC <- overview %>% filter(Passed_QC == 1)
  status <- passed_QC %>% dplyr::select(STR_ID, org_name, chemo_naive, RASTRIC)
  d <- read_plot_data(passed_QC[1, "STR_ID"], passed_QC[1, "org_name"])
  for (n in 2:nrow(passed_QC)) {
    r <- read_plot_data(passed_QC[n, "STR_ID"], passed_QC[n, "org_name"])
    d <- rbind(d, r)
  }
  d <- left_join(d, status, by=c("STR_ID", "org_name"))
  d
}
d <- read_experiment()

organoids = unique(d$org_name)


n1 <- length(organoids)                                 # Amount of default colors
hex_codes1 <- hue_pal()(n1)                             # Identify hex codes
hex_codes1                                              # Print hex codes to console
# "#F8766D" "#7CAE00" "#00BFC4" "#C77CFF"

plot_organoid_condition <- function(df, condition, n=1, save=F) {
  df <- df[df$condition == condition, ]
  title <- condition
  
  organoid = organoids[n]
  print(paste("plotting", organoid))
  organoid_data = df[df$org_name == organoid, ]
  organoid_data$Max_Concentration_log <- log10(organoid_data$conc_condition)
  organoid_data$organoid_color <- hex_codes1[n]
  # Add drug concentration 0 
  # conc_zero <- organoid_data[nrow(organoid_data) + 1,]
  # conc_zero$org_name <- organoid
  # conc_zero$mean_GR <- 1
  # conc_zero$condition <- title
  # conc_zero$Concentration <- min(organoid_data$Max_Concentration, na.rm = TRUE)
  # conc_zero$Max_Concentration <- (min(organoid_data$Max_Concentration, na.rm = TRUE)/10)
  # conc_zero$Max_Concentration_log <- (min(organoid_data$Max_Concentration_log) -1)
  # organoid_data[nrow(organoid_data) + 1,] = conc_zero
  
  min_organoid = min(organoid_data$mean_GR)
  max_organoid = max(organoid_data$mean_GR)
  diff_organoid = max_organoid - min_organoid
  organoid_data$GR_prop <- convertToProp(organoid_data$mean_GR)
  
  fit <- nplr(organoid_data$conc_condition, organoid_data$GR_prop, 
              useLog = TRUE, #should conc values be log10 transformed
              LPweight = 0.25, 
              npars = "all", #number to specify the number of parameters to use in the model; "all" --> the logistic model will be tested with 2 to 5 parameters, best option will be returned
              method = "res", #model optimized using sum squared errors (using weighting method, different options)
              silent = FALSE) #should warning messages be silenced
  
 
  dataframe_fit <- data.frame(getXcurve(fit), getYcurve(fit))
  colnames(dataframe_fit) <- c("Concentration_1", "y_fit")
  dataframe_fit$y_fit_original <- min_organoid + dataframe_fit$y_fit * diff_organoid
  # dataframe_fit$y_fit<- ((dataframe_fit$y_fit)*2)-1
  
  organoid_data$organoid_color <- organoids[n]
  dataframe_fit$organoid_color <- organoids[n]
  
  #metrics
  
  # organoid_data_AUC$GR_positive <- (((organoid_data_AUC$x)+1)/2) #GR range 0-1 for calculating AUC
  # GRMax <- mean(subset(organoid_data, Max_Concentration == max(organoid_data$Max_Concentration))$GR)
  # AUC_raw <- AUC(organoid_data_AUC$Group.1, organoid_data_AUC$GR_positive) #AUC raw
  # AUC_fit <- getAUC(fit)
  # AUC_fit_trapezoid <- AUC_fit$trapezoid #AUC fit
  # AUC_fit_simpson <- AUC_fit$Simpson #AUC fit
  # estim <- getEstimates(fit, .5)
  # GR50_log <- format(estim$x, digits = 6, scientific = TRUE) #estimate at 0.5
  # GR50_num <- as.numeric(GR50_log)
  # GR50 <- log10(GR50_num)
  # par <- getPar(fit) 
  # xmid_log <- (par$params$xmid)
  # xmid<- as.numeric(xmid_log)
  
  # DR_summary_local[nrow(DR_summary_local) + 1,]  <-  list(condition = condition, organoid = organoid, AUC_raw = AUC_raw, AUC_fit_trapezoid = AUC_fit_trapezoid, GRMax = GRMax, GR50 = GR50)
  
  #checkGR<- getEstimates(fit, (1-((1-GRMax)*0.5)))
  #check50_log <- format(esGR$x, digits = 6, scientific = TRUE) 
  #checkGR50_num <- as.numeric(GR50_log)
  #checkGR50 <- log10(GR50_num)
  if (n==1) {
    p <- ggplot() + 
      geom_point(data = organoid_data, aes(conc_condition, mean_GR, color = organoid_color), size = 2) +
      theme_classic() +
      geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_original, color = organoid_color)) +
      #geom_vline(xintercept = checkGR50, color = hex_codes1[1]) +
      labs (x= expression(paste("Concentration (log10) ", mu, "M")), y="GR", main = condition, colour = "Organoid") +  
      scale_x_log10(limits = c(min(organoid_data$conc_condition),max(organoid_data$conc_condition))) +
      ylim(-1,1.5) + ggtitle(condition)
  } else {
    p <- p +
      geom_point(data = organoid_data, aes(conc_condition, mean_GR, color = organoid_color), size = 2) +
      geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_original, color = organoid_color)) 
  }
  p
}
p <- plot_organoid_condition(df, "5-FU", 1)
p <- plot_organoid_condition(df, "5-FU", 2)
p <- plot_organoid_condition(df, "5-FU", 3)
p <- plot_organoid_condition(df, "5-FU", 4)
p <- plot_organoid_condition(df, "5-FU", 5)
p <- plot_organoid_condition(df, "5-FU", 6)
p <- plot_organoid_condition(df, "5-FU", 7)
p <- plot_organoid_condition(df, "5-FU", 8)
p <- plot_organoid_condition(df, "5-FU", 9)
p <- plot_organoid_condition(df, "5-FU", 10)
p <- plot_organoid_condition(df, "5-FU", 11)
p <- plot_organoid_condition(df, "5-FU", 12)
p <- plot_organoid_condition(df, "5-FU", 13)
p <- plot_organoid_condition(df, "5-FU", 14)

show(p)

plot_multiple <- function(df, condition, save=F) {
  for (i in 1:length(organoids)) {
    p <- plot_organoid_condition(df, condition, n=i, save=save)
    p
  }
  show(p)
}

plot_multiple(d, "5-FU")

p <- plot_organoid_condition(d,"5-FU")

# plot all organoid responses to one condition
plot_condition_drcs <- function(df, condition) {
  df = d
  condition = "5-FU"
  # plaatselijke dataframe met de resultaten van de conditie, zodat je deze later bijelkaar kan zetten
  # DR_summary_local <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("condition", "organoid", "AUC_raw", "AUC_fit_trapezoid", "GRMax", "GR50"))))
  
  df <- df[df$condition == condition, ]
  title <- condition
 
  organoid = organoids[1]
  organoid_data = df[df$org_name == organoid, ]
  organoid_data$Max_Concentration_log <- log10(organoid_data$conc_condition)
  organoid_data$organoid_color <- hex_codes1[1]
  # Add drug concentration 0 
  # conc_zero <- organoid_data[nrow(organoid_data) + 1,]
  # conc_zero$org_name <- organoid
  # conc_zero$mean_GR <- 1
  # conc_zero$condition <- title
  # conc_zero$Concentration <- min(organoid_data$Max_Concentration, na.rm = TRUE)
  # conc_zero$Max_Concentration <- (min(organoid_data$Max_Concentration, na.rm = TRUE)/10)
  # conc_zero$Max_Concentration_log <- (min(organoid_data$Max_Concentration_log) -1)
  # organoid_data[nrow(organoid_data) + 1,] = conc_zero
  
  min_organoid = min(organoid_data$mean_GR)
  max_organoid = max(organoid_data$mean_GR)
  diff_organoid = max_organoid - min_organoid
  organoid_data$GR_prop <- convertToProp(organoid_data$mean_GR)
  
  fit <- nplr(organoid_data$conc_condition, organoid_data$GR_prop, 
              useLog = TRUE, #should conc values be log10 transformed
              LPweight = 0.25, 
              npars = "all", #number to specify the number of parameters to use in the model; "all" --> the logistic model will be tested with 2 to 5 parameters, best option will be returned
              method = "res", #model optimized using sum squared errors (using weighting method, different options)
              silent = FALSE) #should warning messages be silenced
  
  min_organoid_lower = min(organoid_data$lower_bound_95)
  max_organoid_lower = max(organoid_data$lower_bound_95)
  diff_organoid_lower = max_organoid_lower - min_organoid_lower
  organoid_data$GR_prop_lower <- convertToProp(organoid_data$lower_bound_95)
  
  fit_lower <- nplr(organoid_data$conc_condition, organoid_data$GR_prop_lower, 
              useLog = TRUE, #should conc values be log10 transformed
              LPweight = 0.25, 
              npars = "all", #number to specify the number of parameters to use in the model; "all" --> the logistic model will be tested with 2 to 5 parameters, best option will be returned
              method = "res", #model optimized using sum squared errors (using weighting method, different options)
              silent = FALSE) #should warning messages be silenced
  
  min_organoid_upper = min(organoid_data$upper_bound_95)
  max_organoid_upper = max(organoid_data$upper_bound_95)
  diff_organoid_upper = max_organoid_upper - min_organoid_upper
  organoid_data$GR_prop_upper <- convertToProp(organoid_data$upper_bound_95)
  
  fit_upper <- nplr(organoid_data$conc_condition, organoid_data$GR_prop_upper, 
                    useLog = TRUE, #should conc values be log10 transformed
                    LPweight = 0.25, 
                    npars = "all", #number to specify the number of parameters to use in the model; "all" --> the logistic model will be tested with 2 to 5 parameters, best option will be returned
                    method = "res", #model optimized using sum squared errors (using weighting method, different options)
                    silent = FALSE) #should warning messages be silenced
  
  
  dataframe_fit <- data.frame(getXcurve(fit), getYcurve(fit), getYcurve(fit_lower), getYcurve(fit_upper))
  colnames(dataframe_fit) <- c("Concentration_1", "y_fit", "y_fit_lower", "y_fit_upper")
  dataframe_fit$y_fit_original <- min_organoid + dataframe_fit$y_fit * diff_organoid
  dataframe_fit$y_fit_lower_original <- min_organoid_lower + dataframe_fit$y_fit_lower * diff_organoid_lower
  dataframe_fit$y_fit_upper_original <- min_organoid_upper + dataframe_fit$y_fit_upper * diff_organoid_upper
    #dataframe_fit$y_fit<- ((dataframe_fit$y_fit)*2)-1
  
  organoid_data$organoid_color <- organoids[1]
  dataframe_fit$organoid_color <- organoids[1]
  
  #metrics
   
  # organoid_data_AUC$GR_positive <- (((organoid_data_AUC$x)+1)/2) #GR range 0-1 for calculating AUC
  # GRMax <- mean(subset(organoid_data, Max_Concentration == max(organoid_data$Max_Concentration))$GR)
  # AUC_raw <- AUC(organoid_data_AUC$Group.1, organoid_data_AUC$GR_positive) #AUC raw
  # AUC_fit <- getAUC(fit)
  # AUC_fit_trapezoid <- AUC_fit$trapezoid #AUC fit
  # AUC_fit_simpson <- AUC_fit$Simpson #AUC fit
  # estim <- getEstimates(fit, .5)
  # GR50_log <- format(estim$x, digits = 6, scientific = TRUE) #estimate at 0.5
  # GR50_num <- as.numeric(GR50_log)
  # GR50 <- log10(GR50_num)
  # par <- getPar(fit) 
  # xmid_log <- (par$params$xmid)
  # xmid<- as.numeric(xmid_log)
  
  # DR_summary_local[nrow(DR_summary_local) + 1,]  <-  list(condition = condition, organoid = organoid, AUC_raw = AUC_raw, AUC_fit_trapezoid = AUC_fit_trapezoid, GRMax = GRMax, GR50 = GR50)
  
  #checkGR<- getEstimates(fit, (1-((1-GRMax)*0.5)))
  #check50_log <- format(esGR$x, digits = 6, scientific = TRUE) 
  #checkGR50_num <- as.numeric(GR50_log)
  #checkGR50 <- log10(GR50_num)
  
  plot <- ggplot() + geom_point(data = organoid_data, aes(conc_condition, mean_GR, color = organoid_color), size = 2) +
    theme_classic() +
    geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_original, color = organoid_color)) + 
    geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_lower_original, color = organoid_color)) + 
    geom_line(data = dataframe_fit, aes(x = 10^Concentration_1, y = y_fit_upper_original, color = organoid_color)) + 
    #geom_vline(xintercept = checkGR50, color = hex_codes1[1]) +
    labs (x= expression(paste("Concentration (log10) ", mu, "M")), y="GR", main = condition,
          colour = "Organoid") +  
    scale_x_log10(limits = c(min(organoid_data$conc_condition),max(organoid_data$conc_condition))) +
    ylim(-1,1.5) + ggtitle(condition)

  
  for (i in 2:length(organoids)) {
    organoid = organoids[i]
    organoid_data = df[df$org_name == organoid, ]
    organoid_data$Max_Concentration_log <- log10(organoid_data$Max_Concentration)
    
    # Add drug concentration 0 
    conc_zero <- organoid_data[nrow(organoid_data) + 1,]
    conc_zero$org_name <- organoid
    conc_zero$GR <- 1
    conc_zero$condition <- title
    conc_zero$Concentration <- min(organoid_data$Max_Concentration, na.rm = TRUE)
    conc_zero$Max_Concentration <- (min(organoid_data$Max_Concentration, na.rm = TRUE)/10)
    conc_zero$Max_Concentration_log <- (min(organoid_data$Max_Concentration_log) -1)
    organoid_data[nrow(organoid_data) + 1,] = conc_zero
    
    min_organoid = min(organoid_data$GR)
    max_organoid = max(organoid_data$GR)
    diff_organoid = max_organoid - min_organoid
    organoid_data$GR_prop <- convertToProp(organoid_data$GR)
    
    fit <- nplr(organoid_data$Max_Concentration, organoid_data$GR_prop, 
                useLog = TRUE, #should conc values be log10 transformed
                LPweight = 0.25, 
                npars = "all", #number to specify the number of parameters to use in the model; "all" --> the logistic model will be tested with 2 to 5 parameters, best option will be returned
                method = "res", #model optimized using sum squared errors (using weighting method, different options)
                silent = FALSE) #should warning messages be silenced
    
    dataframe_fit <- data.frame(getXcurve(fit), getYcurve(fit))
    colnames(dataframe_fit) <- c("Concentration_1", "y_fit")
    dataframe_fit$y_fit_original <- min_organoid + dataframe_fit$y_fit * diff_organoid
    #dataframe_fit$y_fit<- ((dataframe_fit$y_fit)*2)-1
    
    dataframe_fit$organoid_color <- organoids[i]
    organoid_data$organoid_color <- organoids[i]
    
    #metrics
    organoid_data_AUC <- aggregate(organoid_data$GR, list(organoid_data$Max_Concentration), FUN=mean) #calculate mean AUC 
    organoid_data_AUC$GR_positive <- (((organoid_data_AUC$x)+1)/2) #GR range 0-1 for calculating AUC
    GRMax <- mean(subset(organoid_data, Max_Concentration == max(organoid_data$Max_Concentration))$GR)
    AUC_raw <- AUC(organoid_data_AUC$Group.1, organoid_data_AUC$GR_positive) #AUC raw
    AUC_fit <- getAUC(fit)
    AUC_fit_trapezoid <- AUC_fit$trapezoid #AUC fit
    AUC_fit_simpson <- AUC_fit$Simpson #AUC fit
    estim <- getEstimates(fit, .5)
    GR50_log <- format(estim$x, digits = 6, scientific = TRUE) #estimate at 0.5
    GR50_num <- as.numeric(GR50_log)
    GR50 <- log10(GR50_num)
    par <- getPar(fit) 
    xmid_log <- (par$params$xmid)
    xmid<- as.numeric(xmid_log)
    
    plot <- plot + 
      geom_point(data = organoid_data, aes(Max_Concentration_log, GR, color = organoid_color), size = 2) +
      geom_line(data = dataframe_fit, aes(x = Concentration_1, y = y_fit_original, color = organoid_color))  +
      geom_line(data = dataframe_fit, aes(x = Concentration_1, y = y_fit_original, color = organoid_color)) + 
      geom_hline(yintercept = GRMax, color = hex_codes1[i]) +
      geom_vline(xintercept = GR50, color = hex_codes1[i]) 
    
    DR_summary_local[nrow(DR_summary_local) + 1,]  <-  list(condition = condition, organoid = organoid, AUC_raw = AUC_raw, AUC_fit_trapezoid = AUC_fit_trapezoid, GRMax = GRMax, GR50 = GR50)
  }
  ggsave(filename=paste0(experiment_name, "_",  condition, ".png"), plot)  
  DR_summary_local
}


# # Plotten van de curves met ggplot2
# ggplot(experiment_df,
#        aes(x = conc_condition,
#            y = mean_GR)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = lower_bound_95,
#                     ymax = upper_bound_95)) +
#   facet_wrap(~condition) +
#   labs(x = "Concentratie",
#        y = "Gemiddelde GR") +
#   stat_smooth(method=drm, fct=LL.4(), se=FALSE)
# 
# # Maak een lege lijst om de fit objecten op te slaan
# fit_list <- list()
# 
# # Loop over de unieke waarden van condition
# for (cond in unique(experiment_df$condition)) {
#   # Maak een subset van de data voor elke condition
#   sub_data <- experiment_df[experiment_df$condition == cond, ]
#   # Pas de drm functie toe op de subset data
#   fit <- drc::drm(mean_GR ~ conc_condition, data = sub_data, fct = LL.4())
#   # Sla het fit object op in de lijst met de naam van de condition
#   fit_list[[cond]] <- fit
# }
# 
# # Bekijk de lijst met fit objecten
# fit_list
# 
# # Write a function that makes a plot for a given condition
# plot_condition <- function(condition) {
#   
#   # Subset the data frame by condition using filter()
#   df <- filter(experiment_df,
#                condition == condition)
#   
#   # Make a plot with log scale x axis and adjusted range
#   p <- ggplot(df,
#               aes(x = conc_condition,
#                   y = mean_GR)) +
#     geom_point() +
#     geom_errorbar(aes(ymin = lower_bound_95,
#                       ymax = upper_bound_95)) +
#     labs(x = "Concentration",
#          y = "Mean GR",
#          title = paste("Condition", condition)) +
#     stat_smooth(method=drm,fct=LL.4(),se=FALSE) +
#     
#     # Add log scale transformation and limits
#     scale_x_log10(limits = c(0.001,100))
#   
#   # Return the plot
#   return(p)
# }
# 
# # Loop over the unique values of condition
# for (c in unique(experiment_df$condition)) {
#   
#   # Clear the plot device
#   dev.off()
#   
#   # Call the function and assign the result to p_c
#   p_c <- plot_condition(c)
#   
#   # Print out or save the plot as desired
#   print(p_c)
#   ggsave(p_c,file.path(plot_dir, paste("p_",c,"png",sep="")))
# }