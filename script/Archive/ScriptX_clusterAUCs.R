library(dplyr)
library(tidyr)
library(readxl)
library(openxlsx)
library(ggplot2)
library(cowplot)

rm(list=ls())

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# script.dir <- dirname(sys.frame(1)$ofile) 
home_dir <- dirname(script_dir)
metrics_dir <- file.path(home_dir, "6_metrics")
metrics_normalized_dir <- file.path(metrics_dir, "QC_normalized")
resource_dir <- file.path(home_dir, "resources")

data <- read_excel(file.path(metrics_normalized_dir, "QC_normalized_metrics_normalized.xlsx"))
cluster_file <- read_excel(file.path(resource_dir, "Treatment_Clusters.xlsx"))

cluster_long <- cluster_file %>%
  gather(key = "DrugNum", value = "condition", -ClusterName)

result <- data %>%
  left_join(cluster_long, by = "condition") %>%
  select(-DrugNum) %>%
  mutate(ClusterName = replace_na(ClusterName, "Unknown"))

drug_count <- cluster_long %>%
  filter(!is.na(condition)) %>%
  group_by(ClusterName) %>%
  summarise(drug_count = n())


result <- result %>%
  group_by(org_name, chemo_naive, ClusterName) %>%
  summarise(
    AUC_sum = sum(AUC, na.rm = TRUE),
    GRmax_sum = sum(GRmax, na.rm = TRUE),
    AUC_low_sum = sum(AUC_low, na.rm = TRUE),
    GRmax_low_sum = sum(GRmax_low, na.rm = TRUE),
    GR50_sum = sum(GR50, na.rm = TRUE),
    xmid_sum = sum(xmid, na.rm = TRUE)
  )

result_summary_with_count <- result %>%
  left_join(drug_count, by = "ClusterName")

result_normalized <- result_summary_with_count %>%
  mutate(
    AUC_normalized = AUC_sum / drug_count,
    GRmax_normalized = GRmax_sum / drug_count,
    AUC_low_normalized = AUC_low_sum / drug_count,
    GRmax_low_normalized = GRmax_low_sum / drug_count,
    GR50_normalized = GR50_sum / drug_count,
    xmid_normalized = xmid_sum / drug_count
  ) %>%
  select(-c(AUC_sum, GRmax_sum, AUC_low_sum, GRmax_low_sum, GR50_sum, xmid_sum, drug_count))


write.xlsx(result_normalized, file.path(metrics_normalized_dir,paste0("QC_normalized_metrics_normalized_clustered_normalized.xlsx")))

make_boxplots <- function(data){
  boxplot_AUC_low <- ggplot(data, aes(y=AUC_low_normalized, x=chemo_naive, color=chemo_naive)) + 
    facet_grid( ~ClusterName) +
    geom_boxplot(size =0.6, outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.1))+
    #stat_summary(fun = mean, geom = "text", col = "red", vjust = +2, aes(label = paste("Mean:", round(..y.., digits = 2))))+
    #annotate("text", x = 1, y = 1.7, label = substitute(paste(italic('p'), " = ...")), size =3)+
    labs(title = "", y=expression(GR[AUC]*" low conc."), x = "Chemo-exposure") + 
    theme_classic() 
  boxplot_AUC_low
  
  boxplot_AUC <- ggplot(data, aes(y=AUC_normalized, x=chemo_naive, color=chemo_naive)) + 
    facet_grid( ~ClusterName) +
    geom_boxplot(size =0.6, outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.1))+
    #stat_summary(fun = mean, geom = "text", col = "red", vjust = +2, aes(label = paste("Mean:", round(..y.., digits = 2))))+
    #annotate("text", x = 1, y = 1.7, label = substitute(paste(italic('p'), " = ...")), size =3)+
    labs(title = "", y=expression(GR[AUC]), x = "Chemo-exposure") + 
    theme_classic() 
  boxplot_AUC
  
  boxplot_GRmax_low <- ggplot(data, aes(y=GRmax_low_normalized, x=chemo_naive, color=chemo_naive)) + 
    facet_grid( ~ClusterName) +
    geom_boxplot(size =0.6, outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.1))+
    #stat_summary(fun = mean, geom = "text", col = "red", vjust = +2, aes(label = paste("Mean:", round(..y.., digits = 2))))+
    #annotate("text", x = 1, y = 1.7, label = substitute(paste(italic('p'), " = ...")), size =3)+
    labs(title = "", y=expression(GR[max]* "low conc."), x = "Chemo-exposure") + 
    theme_classic() 
  boxplot_GRmax_low
  
  boxplot_GRmax <- ggplot(data, aes(y=GRmax_normalized, x=chemo_naive, color=chemo_naive)) + 
    facet_grid( ~ClusterName) +
    geom_boxplot(size =0.6, outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.1))+
    #stat_summary(fun = mean, geom = "text", col = "red", vjust = +2, aes(label = paste("Mean:", round(..y.., digits = 2))))+
    #annotate("text", x = 1, y = 1.7, label = substitute(paste(italic('p'), " = ...")), size =3)+
    labs(title = "", y=expression(GR[max]), x = "Chemo-exposure") + 
    theme_classic() 
  boxplot_GRmax
  
  boxplot_GR0 <- ggplot(data, aes(y=GR50_normalized, x=chemo_naive, color=chemo_naive)) + 
    facet_grid( ~ClusterName) +
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
  
  tiff(file = file.path(metrics_normalized_dir, "boxplots_metrics_QCpassed_norm_clustered.tiff"),
       units="in",
       width = (15),
       height = (15),
       res=300) 
  print(grid_boxplots)
  dev.off()
  
}
make_boxplots(result_normalized)
