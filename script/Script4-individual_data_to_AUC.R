library(dplyr)
library(readxl)
library(ggplot2)
library(openxlsx)
library(scales)   
library(drc)
library(nplr)
# library(patchwork)
library(DescTools) #AUC
library(cowplot)
library(Cairo)

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
metrics_dir <- file.path(home_dir, "6_metrics")
metrics_plot_dir <- file.path(metrics_dir, "Plots metrics")

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

normalize_df <- function(df) {
  metrics_df_norm <- df
  metrics_df_norm = metrics_df_norm %>% dplyr::select(org_name,condition,chemo_naive,AUC,GRmax,AUC_low,GRmax_low,GR50,xmid)
  
  #normalize per condition IN PROGRESS
  for(i in 4:ncol(metrics_df_norm)) {       
    metrics_df_norm[i] <- (metrics_df_norm[i]- min(metrics_df_norm[i], na.rm=TRUE))/(max(metrics_df_norm[i], na.rm=TRUE)-min(metrics_df_norm[i], na.rm=TRUE))
  }
  
  metrics_df_norm <- metrics_df_norm %>%
    group_by(condition) %>%
    mutate_at(vars(4:last_col()), list(~(. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))))
  return(metrics_df_norm)
}

save_metrics <- function(exp_name, 
                     select_on="Passed_QC", 
                     selection_value = 1, 
                     selected_conditions="all", # of c("5-FU", "Oxaliplatin", "etc.")) 
                     normalize=TRUE) 
{
  
    # create an empty data frame to store the fitted values
  dataframe_conditions <- data.frame()
  
  #DIT WERKT WEL
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
    
    #fill df with metrics
    metric <- d %>% dplyr::select(conc_condition, mean_GR, org_name, condition, chemo_naive) %>%
      mutate(mean_GR_positive=(mean_GR +1)/2)  %>% #GR range 0-1 for calculating AUC
      group_by(org_name, condition) %>%
      mutate(AUC = AUC(log10(conc_condition), mean_GR_positive)) %>%
      group_by(org_name, condition) %>%
      mutate(max_index = which.max(conc_condition),
             GRmax = mean_GR[max_index]) %>%
      dplyr::select(-max_index) %>%
      group_by(org_name, condition, AUC, GRmax) %>% 
      summarise(AUC = mean(AUC))
    
      #dataframe for metrics only low concentrations
    metric_lowconc <- d %>% dplyr::select(conc_condition, mean_GR, org_name, condition, chemo_naive) %>%
      mutate(mean_GR_positive=(mean_GR +1)/2)  %>% #GR range 0-1 for calculating AUC
      group_by(org_name, condition, chemo_naive) %>%
      slice_min(conc_condition, n = 4) %>%
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
    write.xlsx(dataframe_conditions, file.path(metrics_dir, exp_name, paste0(exp_name, "_metrics_normalized.xlsx")))
  } else {
    df_return <- dataframe_conditions
  }
  # returns normalized df if Normalize is True (standard parameter in function)  
  return(df_return)
}

# Use Mann-Whitney U-test (non-parametric alternative to 2-sample t-test)
test_difference <- function(df, parameter, filter_condition, factor="chemo_naive") {
  test_metric <- filter(df, condition == filter_condition)
  print(head(test_metric))
  results <- wilcox.test(test_metric[[parameter]] ~ test_metric[[factor]], data=test_metric) 
  # results_chemo_naive$p.value <- round(results_chemo_naive$p.value, digits = 3)
  return(results)
}
# use as follows: 
# test_difference(df, parameter = "AUC_low", filter_condition = "vi_bi_la")

make_boxplots <- function(data){
  boxplot_AUC_low <- ggplot(data, aes(y=AUC_low, x=chemo_naive, color = chemo_naive)) + 
    facet_grid(~condition) +
    geom_boxplot(size =0.6, outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.1))+
    #stat_summary(fun = mean, geom = "text", col = "red", vjust = +2, aes(label = paste("Mean:", round(..y.., digits = 2))))+
    #annotate("text", x = 1, y = 1.7, label = substitute(paste(italic('p'), " = ...")), size =3)+
    labs(title = "", y=expression(GR[AUC]*" low conc."), x = "Chemo-exposure") + 
    theme_classic() 
  boxplot_AUC_low
  
  boxplot_AUC <- ggplot(data, aes(y=AUC, x=chemo_naive, color = chemo_naive)) + 
    facet_grid(~condition) +
    geom_boxplot(size =0.6, outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.1))+
    #stat_summary(fun = mean, geom = "text", col = "red", vjust = +2, aes(label = paste("Mean:", round(..y.., digits = 2))))+
    #annotate("text", x = 1, y = 1.7, label = substitute(paste(italic('p'), " = ...")), size =3)+
    labs(title = "", y=expression(GR[AUC]), x = "Chemo-exposure") + 
    theme_classic() 
  boxplot_AUC
  
  boxplot_GRmax_low <- ggplot(data, aes(y=GRmax_low, x=chemo_naive, color = chemo_naive)) + 
    facet_grid(~condition) +
    geom_boxplot(size =0.6, outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.1))+
    #stat_summary(fun = mean, geom = "text", col = "red", vjust = +2, aes(label = paste("Mean:", round(..y.., digits = 2))))+
    #annotate("text", x = 1, y = 1.7, label = substitute(paste(italic('p'), " = ...")), size =3)+
    labs(title = "", y=expression(GR[max]* "low conc."), x = "Chemo-exposure") + 
    theme_classic() 
  boxplot_GRmax_low
  
  boxplot_GRmax <- ggplot(data, aes(y=GRmax, x=chemo_naive, color = chemo_naive)) + 
    facet_grid(~condition) +
    geom_boxplot(size =0.6, outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.1))+
    #stat_summary(fun = mean, geom = "text", col = "red", vjust = +2, aes(label = paste("Mean:", round(..y.., digits = 2))))+
    #annotate("text", x = 1, y = 1.7, label = substitute(paste(italic('p'), " = ...")), size =3)+
    labs(title = "", y=expression(GR[max]), x = "Chemo-exposure") + 
    theme_classic() 
  boxplot_GRmax
  
  boxplot_GR0 <- ggplot(data, aes(y=GR50, x=chemo_naive, color = chemo_naive)) + 
    facet_grid(~condition) +
    geom_boxplot(size =0.6, outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.1))+
    #stat_summary(fun = mean, geom = "text", col = "red", vjust = +2, aes(label = paste("Mean:", round(..y.., digits = 2))))+
    #annotate("text", x = 1, y = 1.7, label = substitute(paste(italic('p'), " = ...")), size =3)+
    labs(title = "", y=expression(GR["0-intercept"]), x = "Chemo-exposure") + 
    theme_classic() 
  boxplot_GR0
  
  grid_boxplots <- plot_grid(
    boxplot_AUC,
    boxplot_AUC_low,
    boxplot_GRmax,
    boxplot_GRmax_low,
    boxplot_GR0,
    ncol = 1,
    nrow = 5,
    labels = c(''))
  
  tiff(file = file.path(metrics_plot_dir, "boxplots_metrics_QCpassed_norm.tiff"),
       units="in",
       width = (27),
       height = (15),
       res=300) 
  print(grid_boxplots)
  dev.off()

}

df <- save_metrics("QC_normalized")
make_boxplots(df)

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