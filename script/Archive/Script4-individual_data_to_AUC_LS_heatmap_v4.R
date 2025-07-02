library(dplyr)
library(readxl)
library(ggplot2)
library(openxlsx)
library(scales)   
library(drc)
library(nplr)
library(DescTools) #AUC
library(cowplot)
library(ggpubr) #statcomparemeans

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
metrics_normalized_dir <- file.path(metrics_dir, "240506_QC")

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
    Analyse <- overview %>% filter(get(select_on) == selection_value) #default waarop hij overview filtert is Analyse, als je all invult selecteert hij op Data_processed ipv Analyse
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
  #normalize per condition IN PROGRESS
  for(i in 5:ncol(metrics_df_norm)) {       
    metrics_df_norm[i] <- (metrics_df_norm[i]- min(metrics_df_norm[i], na.rm=TRUE))/(max(metrics_df_norm[i], na.rm=TRUE)-min(metrics_df_norm[i], na.rm=TRUE))
  }
  
  metrics_df_norm <- metrics_df_norm %>%
    group_by(condition) %>%
    mutate_at(vars(5:last_col()), list(~(. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))))
  return(metrics_df_norm)
}

save_metrics <- function(exp_name, 
                     select_on="Passed_QC", 
                     selection_value = 1, 
                     selected_conditions="all", # of c("5-FU", "Oxaliplatin", "etc.")) 
                     normalize=TRUE) 
{

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

    #fill df with metrics
    metric <- d %>% dplyr::select(conc_condition, mean_GR, org_name, STR_ID, condition, chemo_naive) %>%
      mutate(mean_GR_positive=(mean_GR +1)/2)  %>% #GR range 0-1 for calculating AUC
      group_by(org_name, STR_ID, condition) %>%
      mutate(AUC = AUC(log10(conc_condition), mean_GR_positive)) %>%
      group_by(org_name, condition) %>%
      mutate(max_index = which.max(conc_condition),
             GRmax = mean_GR[max_index]) %>%
      dplyr::select(-max_index) %>%
      group_by(org_name, STR_ID, condition, AUC, GRmax) %>% 
      summarise(AUC = mean(AUC))
   
    concentration_filters <- list(
      `5-FU` = c(489.5, 490.1, 489.8, 108.2, 108.1,108.4, 24.29, 24.36, 23.94, 5.316, 5.296, 5.297), 
      Oxaliplatin = c(0.1, 0.10100, 0.09937, 0.39540, 0.39450, 0.3947, 1.57000, 1.53300, 1.514, 6.13800, 6.132, 6.1930), 
      `SN-38` =  c(0.0069300, 0.020790, 0.064350, 0.1980000), 
      SN38_CHEK1 =  c(0.0069300, 0.020790, 0.064350, 0.1980000),
      CHEK1 = c(0.038610, 0.118800, 0.336600, 0.999900),
      Alpelisib = c(0.57420, 1.88100, 6.13800, 20.00000),#1 conc hoger, geen effect
      Lapatinib = c(0.57420, 1.88100, 6.13800, 20.00000), #1 conc hoger, geen effect
      Binimetinib = c(0.17330,0.57420, 1.88100, 6.13800), 
      lapatinib_binimetinib = c(0.17330,0.57420, 1.88100, 6.13800),
      binimetinib_alpelisib = c(0.17330,0.57420, 1.88100, 6.13800), 
      lapatinib_alpelisib = c(0.57420, 1.88100, 6.13800, 20.00000),#1 conc hoger, geen effect
      Vinorelbine = c(0.001485, 0.003713, 0.010400, 0.028710), # 3 conc lager (overleg!)
      navitoclax_vinorelbine = c(0.01980,  0.07920,  0.32670, 1.28700), #3 conc lager (overleg!)
      Navitoclax = c(0.07920, 0.32670, 1.28700, 5.24700), # 2 conc hoger, geen effect
      vi_bi_la = c(0.001485, 0.003713, 0.010400, 0.028710) # 3 conc lager (overleg!)
    )
    
    filtered_data <- d %>%
      group_by(condition) %>%
      filter(conc_condition %in% concentration_filters[[unique(condition)]]) %>%
      ungroup()

      #dataframe for metrics only low concentrations
    metric_lowconc <- filtered_data %>% dplyr::select(conc_condition, mean_GR, org_name, STR_ID, condition, chemo_naive) %>%
      mutate(mean_GR_positive=(mean_GR +1)/2)  %>% #GR range 0-1 for calculating AUC
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

make_boxplots <- function(exp_name, data){
  boxplot_AUC_low <- ggplot(data, aes(y=AUC_low, x=chemo_naive, color = chemo_naive)) + 
    facet_grid(~condition) +
    geom_boxplot(size =0.6, outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.1))+
    #stat_summary(fun = mean, geom = "text", col = "red", vjust = +2, aes(label = paste("Mean:", round(..y.., digits = 2))))+
    #annotate("text", x = 1, y = 1.7, label = substitute(paste(italic('p'), " = ...")), size =3)+
    labs(title = "", y=expression(GR[AUC]*" low conc."), x = "Chemo-exposure") + 
    theme_classic() +
    stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE),
                       label = "p.signif", hide.ns = FALSE, label.y = 0.97)
  boxplot_AUC_low
  
  boxplot_AUC <- ggplot(data, aes(y=AUC, x=chemo_naive, color = chemo_naive)) + 
    facet_grid(~condition) +
    geom_boxplot(size =0.6, outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.1))+
    #stat_summary(fun = mean, geom = "text", col = "red", vjust = +2, aes(label = paste("Mean:", round(..y.., digits = 2))))+
    #annotate("text", x = 1, y = 1.7, label = substitute(paste(italic('p'), " = ...")), size =3)+
    labs(title = "", y=expression(GR[AUC]), x = "Chemo-exposure") + 
    theme_classic() +
    stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE),
                       label = "p.signif", hide.ns = FALSE, label.y = 0.97)
  boxplot_AUC
  
  boxplot_GRmax_low <- ggplot(data, aes(y=GRmax_low, x=chemo_naive, color = chemo_naive)) + 
    facet_grid(~condition) +
    geom_boxplot(size =0.6, outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.1))+
    #stat_summary(fun = mean, geom = "text", col = "red", vjust = +2, aes(label = paste("Mean:", round(..y.., digits = 2))))+
    #annotate("text", x = 1, y = 1.7, label = substitute(paste(italic('p'), " = ...")), size =3)+
    labs(title = "", y=expression(GR[max]* "low conc."), x = "Chemo-exposure") + 
    theme_classic() +
    stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE),
                       label = "p.signif", hide.ns = FALSE, label.y = 0.97)
  boxplot_GRmax_low
  
  boxplot_GRmax <- ggplot(data, aes(y=GRmax, x=chemo_naive, color = chemo_naive)) + 
    facet_grid(~condition) +
    geom_boxplot(size =0.6, outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.1))+
    #stat_summary(fun = mean, geom = "text", col = "red", vjust = +2, aes(label = paste("Mean:", round(..y.., digits = 2))))+
    #annotate("text", x = 1, y = 1.7, label = substitute(paste(italic('p'), " = ...")), size =3)+
    labs(title = "", y=expression(GR[max]), x = "Chemo-exposure") + 
    theme_classic() +
    stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE),
                       label = "p.signif", hide.ns = FALSE, label.y = 0.97)
  boxplot_GRmax
  
  boxplot_GR0 <- ggplot(data, aes(y=GR50, x=chemo_naive, color = chemo_naive)) + 
    facet_grid(~condition) +
    geom_boxplot(size =0.6, outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.1))+
    #stat_summary(fun = mean, geom = "text", col = "red", vjust = +2, aes(label = paste("Mean:", round(..y.., digits = 2))))+
    #annotate("text", x = 1, y = 1.7, label = substitute(paste(italic('p'), " = ...")), size =3)+
    labs(title = "", y=expression(GR["0-intercept"]), x = "Chemo-exposure") + 
    theme_classic() +
    stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE),
                       label = "p.signif", hide.ns = FALSE, label.y = 0.97)
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
  
  tiff(file = file.path(metrics_dir, exp_name, "boxplots_metrics_QCpassed_norm.tiff"),
       units="in",
       width = (27),
       height = (15),
       res=300) 
  print(grid_boxplots)

}

df <- save_metrics("240506_QC")
make_boxplots("240506_QC", df)
#make_boxplots("240506_QC", norm_df2)
norm_df2 <- read_excel(file.path(metrics_normalized_dir, "240506_QC_metrics_normalized.xlsx")) #if you do not want to run the full script

#make heatmap AUC----------------------
library(pheatmap)
library(tidyr) #spread

all <- c("5-FU", "vi_bi_la", "Oxaliplatin", "Lapatinib", "Binimetinib", "Alpelisib", "Navitoclax", "CHEK1", "SN-38" ,
         "alpelisib_lapatinib",  "binimetinib_lapatinib", "SN38_CHEK1", "alpelisib_lapatinib", "Vinorelbine", "binimetinib_alpelisib", "navitoclax_vinorelbine")
targeted <-  c("vi_bi_la" ,"Lapatinib", "Binimetinib","Alpelisib","Navitoclax",
               "Vinorelbine","binimetinib_alpelisib", "navitoclax_vinorelbine")
chemo <- c("5-FU", "Oxaliplatin", "SN-38")
chek <- c("5-FU", "Oxaliplatin", "SN-38", "vi_bi_la" ,"Lapatinib", "Binimetinib","Alpelisib","Navitoclax",
          "Vinorelbine","binimetinib_alpelisib", "navitoclax_vinorelbine")
singles <- c("5-FU", "Oxaliplatin", "Lapatinib", "Binimetinib", "Alpelisib", "Navitoclax", "CHEK1", "SN-38", "Vinorelbine")

make_heatmap <- function(exp_name, data, selection=NULL){
  selection <- selection

  data <- data %>% filter(condition %in% selection)
  
  data_select = data %>% dplyr::select("org_name","condition", "chemo_naive", "AUC")
  data_select <- data_select %>% distinct(org_name, condition, chemo_naive, AUC, .keep_all = TRUE)
  data_reshape2 <- spread(data_select,                                  
                          key = condition,
                          value = AUC)  
  data_reshape2 <- as.data.frame(data_reshape2)
  rownames(data_reshape2) <- data_reshape2$org_name #reassigning row names
  data_reshape2$org_name <- NULL
  heatmap_annotation <- data_reshape2[ c("chemo_naive")]
  heatmap_annotation <- heatmap_annotation %>% 
    rename(
      "Chemo naive" = chemo_naive)
  # Create a named vector for annotation colors
  ann_colors <- c("no" = "red", "yes" = "green")  # Assuming "0" is not naive and "1" is naive
  
  data_reshape2$chemo_naive <- NULL
  data_reshape2$STR_ID <- NULL
  pheatmap_AUC <- pheatmap(data_reshape2,
                     cluster_cols = TRUE, clustering_distance_cols  = "euclidean",
                     cluster_rows = TRUE, clustering_distance_rows  = "euclidean",
                     annotation_row = heatmap_annotation,
                     legend_breaks = c(0, 1),
                     legend_labels = c("sensitive", "resistant"),
                     annotation_colors = list("Chemo naive" = ann_colors),
                     angle_col = "315",
                     row_names_rot = 90,
                     color = colorRampPalette(c("blue", "white", "red"))(50),
                     main = "AUC")
  
  ggsave(file = file.path(metrics_dir, exp_name, "AUC_heatmap.pdf"), plot = pheatmap_AUC, width = 7, height = 8)
  
  data_select = data %>% dplyr::select("org_name","condition", "chemo_naive", "AUC_low")
  data_select <- data_select %>% distinct(org_name, condition, chemo_naive, AUC_low, .keep_all = TRUE)
  data_reshape2 <- spread(data_select,                                  
                          key = condition,
                          value = AUC_low)  
  data_reshape2 <- as.data.frame(data_reshape2)
  rownames(data_reshape2) <- data_reshape2$org_name #reassigning row names
  data_reshape2$org_name <- NULL
  heatmap_annotation <- data_reshape2[ c("chemo_naive")]
  heatmap_annotation <- heatmap_annotation %>% 
    rename(
      "Chemo naive" = chemo_naive)
  data_reshape2$chemo_naive <- NULL
  data_reshape2$STR_ID <- NULL
  pheatmap_AUClow <- pheatmap(data_reshape2,
                       cluster_cols = TRUE, clustering_distance_cols  = "euclidean",
                       cluster_rows = TRUE, clustering_distance_rows  = "euclidean",
                       annotation_row = heatmap_annotation,
                       legend_breaks = c(0, 1),
                       legend_labels = c("sensitive", "resistant"),
                       annotation_colors = list("Chemo naive" = ann_colors),
                       angle_col = "315",
                       row_names_rot = 90,
                       color = colorRampPalette(c("blue", "white", "red"))(50),
                       main = "AUC_low")
  
  ggsave(file = file.path(metrics_dir, exp_name, "AUC_low_heatmap.pdf"), plot = pheatmap_AUClow, width = 7, height = 8)
}

make_heatmap("240513_QC_single", norm_df2, selection=singles)

