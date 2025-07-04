# Version control
# v14 added color scale

#library(colorspace) #diverging colors
#q6 <- diverging_hcl(10, palette = "Berlin") 
#"#E0E1E0" "#D8E7F6" "#7FBFF5" "#005682" "#241211" "#4B201D" "#7F3C37" "#F8A29E" "#F15D3F"
# Voor boxplots gebruik ik "#7FBFF5" en "#F8A29E"

#---- load libraries and data ----
# Load dplyr and readxl packages for data manipulation and reading Excel files
library(dplyr)
library(data.table)
library(ggplot2)
library(readxl)
library(RColorBrewer)
library(pheatmap)
library(ggpubr) #stat_compare_means
library(gridExtra) #gridarrang

rm(list=ls()) #Clear existing objects in the workspace

# Define necessary directories for accessing and saving files
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
home_dir <- dirname(dirname(script_dir))
WGS_dir <- file.path(home_dir, "WGS")
WGS_data_dir <- file.path(WGS_dir, "Data_HMF")
resource_dir <- file.path(home_dir, "Analyse 2.0/resources")
WGS_plot_dir <- file.path(home_dir, "Analyse 2.0/7_WGS")

#------import DEPMAP data (validation set)-----
plot_depmap_facet <- function(file_paths, titles, facet_labels) {
  df_list <- lapply(seq_along(file_paths), function(i) {
    df <- read_xlsx(file_paths[i])
    df[df == "NA"] <- NA
    df <- na.omit(df)
    colnames(df) <- c("model", "cell_line_name", "AUC", "mutation")
    df$mutation <- ifelse(df$mutation == 0, "absent", "present")
    df$AUC <- as.numeric(df$AUC)
    df$facet_label <- facet_labels[i]  # Toevoegen van een facet variabele
    return(df)
  })
  
  df_combined <- do.call(rbind, df_list)  # Combineer de datasets
  
  plot <- ggplot(df_combined, aes(y = AUC, x = mutation, color = mutation)) + 
    geom_jitter(width = 0.2, size = 1) +
    stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5) +
    #stat_compare_means(method = "t.test", aes(label = ..p.format..), label.y = max(df_combined$AUC) * 0.95) + 
    labs(title = "DEPMAP Analysis", y = "AUC irinotecan", x = "Deep Deletion CSF") + 
    theme_classic() +
    scale_color_manual(values = c("#F8A29E", "#7FBFF5")) +
    facet_grid(. ~ facet_label, scales = "free_y") +  # Facet toevoegen
    theme(legend.position = "none")
  
  return(plot)
}

# Gebruik de functie met meerdere bestanden
file_paths <- c(file.path(resource_dir, "AUC_irinotecan_deepdeletions_CRC.xlsx"),
                file.path(resource_dir, "AUC_irinotecan_deepdeletions.xlsx"))
facet_labels <- c("CRC cell lines", "cancer cell lines")
boxplot_depmap_facet <- plot_depmap_facet(file_paths, titles = facet_labels, facet_labels = facet_labels)

# Opslaan als PDF
ggsave(file = sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_DEPMAP_faceted"), plot = boxplot_depmap_facet, width = 3, height = 2.5)

corplot_depmap <- function(file_path, title) {
  df <- read_xlsx(file_path)
  df[df == "NA"] <- NA
  df <- na.omit(df)
  colnames(df) <- c("model", "cell_line_name", "AUC", "exp_cor", "exp")
  df$exp_cor <- as.numeric(df$exp_cor)
  df$AUC <- as.numeric(df$AUC)
  
  plot <- ggplot(df, aes(y = exp_cor, x = AUC)) + 
    geom_point(size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "#E0E1E0") +
    labs(title = title, y = "Expression", x = "AUC 5-FU") + 
    stat_cor(method = "spearman", size = 3, color = "black")+
    theme_classic() +
    theme(legend.position = "none")
  return(plot)
}

boxplot_depmap_5FU_L1CAM <- corplot_depmap(file.path(resource_dir, "AUC_5FU_L1CAM.xlsx"), "DEPMAP CRC cell lines")
boxplot_depmap_5FU_L1CAM 

#------import WGS data Hartwig (validation set)-----

# Import the data
HMF_treatment <- read.csv(file.path(WGS_data_dir,"Hartwig_treatments_DR279_clean.tsv"), sep = "\t") #Hartwig data (validation set) treatment status
drivers <- read.csv(file.path(WGS_data_dir,"linx_driver_gene_data_prefilter.tsv"), sep = "\t")   #Hartwig data (validation set) drivers
purple <- read.csv(file.path(WGS_data_dir,"purple_purity_data_prefilter.tsv"), sep = "\t")   #Hartwig data (validation set) ploidy, purity, TMB/TML

# Clean Hartwig data validation set
HMF_treatment_filter <- HMF_treatment %>% filter(hadTreatment != "unknown")
HMF_treatment_filter <- HMF_treatment_filter %>% filter(!(standardOfCare == "no" & hadTreatment == "yes"))
HMF_treatment_filter$chemotherapy <- ifelse(is.na(HMF_treatment_filter$standardOfCare), "chemonaive", "pretreated")
HMF_treatment_filter$hadTreatment <- NULL
HMF_treatment_filter$standardOfCare <- NULL

# Merge Hartwig data (validation set) into 'combined' dataframe
HMF_treatment_filter <- HMF_treatment_filter %>% dplyr::rename(sampleId = patientId) #change patientId to sampleId for merge
purple_treatment <- merge(HMF_treatment_filter, purple, by = "sampleId") #merge dataframes
purple_treatment <- purple_treatment %>% filter(msStatus == "MSS") 
combined <- merge(purple_treatment, drivers, by = "sampleId") #merge dataframes

result <- anti_join(HMF_treatment_filter, purple, by = "sampleId") #samples not present in both datasets


#------WGS data Hartwig plots-----
# Make boxplots for TMB/TML/WG disruptions in Hartwig data (validation set)
boxplot_tml <- ggplot(purple_treatment, aes(y = tml, x = chemotherapy, color=chemotherapy)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "tml", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "t.test", label = "p.format")+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)+
  scale_color_manual(values = c("#7FBFF5", "#F8A29E"))+
  theme(legend.position = "none")+
  stat_summary(fun = median, geom = "text", aes(label = format(..y..,)), vjust = -2, color = "black")
boxplot_tml

purple_treatment %>%
  group_by(chemotherapy) %>%
  summarise(count = n())

boxplot_tmb <- ggplot(purple_treatment, aes(y = tmbPerMb, x = chemotherapy, color=chemotherapy)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "tmb per Mb", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "t.test", label = "p.format")+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)+
  scale_color_manual(values = c("#7FBFF5", "#F8A29E"))+
  theme(legend.position = "none")+
  stat_summary(fun = median, geom = "text", aes(label = format(..y..,)), vjust = -2, color = "black")
boxplot_tmb

boxplot_svtmb <- ggplot(purple_treatment, aes(y = svTumorMutationalBurden, x = chemotherapy, color=chemotherapy)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "sv tmb", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "t.test", label = "p.format")+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)+
  scale_color_manual(values = c("#7FBFF5", "#F8A29E"))+
  theme(legend.position = "none")+
  stat_summary(fun = median, geom = "text", aes(label = format(..y..,)), vjust = -2, color = "black")
boxplot_svtmb

#contingency_table for whole genome duplications in Hartwig data (validation set)
contingency_table <- table(purple_treatment$wholeGenomeDuplication, purple_treatment$chemotherapy == "pretreated")
dimnames(contingency_table) <- list(wholeGenomeDuplication = c("FALSE", "TRUE"), chemotherapy_pretreated = c("FALSE", "TRUE"))
fisher_result <- fisher.test(contingency_table) #not significant
contingency_table

boxplot_ploidy <- ggplot(purple_treatment, aes(y = ploidy, x = chemotherapy, color=chemotherapy)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "ploidy", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "t.test", label = "p.format")+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)+
  scale_color_manual(values = c("#7FBFF5", "#F8A29E"))+
  theme(legend.position = "none")+
  stat_summary(fun = median, geom = "text", aes(label = format(..y..,)), vjust = -2, color = "black")
boxplot_ploidy

grid_plot <- grid.arrange(grobs = c(list(boxplot_tml, boxplot_tmb, boxplot_svtmb, boxplot_ploidy)), ncol=4)
ggsave(file = sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_HMF_tmb_etc"), plot = grid_plot, width = 12, height = 3)

#------import WGS data PDO made in script 9-----
mat_HMF_treatment <- readRDS(file.path(WGS_plot_dir, "mat_HMF_treatment.rds")) #OPT0044 en OPT0404 niet hier in
sig_PDO_categorized <- readRDS(file.path(WGS_plot_dir, "sig_PDO_categorized.rds"))

# Define a function to filter dataframes based on the 'treatment' column
filter_treatment <- function(df) {
  df %>% filter(condition %in% c("5-FU", "Oxaliplatin", "SN-38"))}

# List of dataframes
dataframes <- list(
  PDO_WGS_screen = readRDS(file.path(WGS_plot_dir, "PDO_WGS_screen.rds")),
  PDO_WGS_screen_top5 = readRDS(file.path(WGS_plot_dir, "PDO_WGS_screen_top5.rds")),
  PDO_WGS_screen_pretreatedfiltered = readRDS(file.path(WGS_plot_dir, "PDO_WGS_screen_pretreatedfiltered.rds")),
  PDO_WGS_screen_pretreatedonly = readRDS(file.path(WGS_plot_dir, "PDO_WGS_screen_pretreatedonly.rds")),
  PDO_metadata_screen = readRDS(file.path(WGS_plot_dir, "PDO_metadata_screen.rds")),
  PDO_WGS_screen_top5_chemonaive = readRDS(file.path(WGS_plot_dir, "PDO_WGS_screen_top5_chemonaive.rds")),
  PDO_WGS_screen_top5_pretreated = readRDS(file.path(WGS_plot_dir, "PDO_WGS_screen_top5_pretreated.rds"))
)

filtered_dataframes <- lapply(dataframes, filter_treatment)# Apply the filtering function to each dataframe
list2env(filtered_dataframes, .GlobalEnv)# Assign each filtered dataframe back to the global environment with the same names
ls()# Verify the dataframes are now available in the global environment

PDO_WGS_screen_top5_all <- rbind(PDO_WGS_screen_top5_chemonaive, PDO_WGS_screen_top5_pretreated)

PDO_WGS_screen_chemonaivefiltered <- PDO_WGS_screen %>% 
  filter(`Sample type`=="chemonaive")

OPTIC_WGS_SBS <- read.table("D:/SURFdrive/Lsmabers/PhD/OPTIC_LS/Analyse/Plots/WGS_OPTIC/SBS_DBS_signature_contribution.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
OPTIC_WGS_SBS$DBS5[OPTIC_WGS_SBS$SampleID == ""]


grid_plot <- grid.arrange(grobs = c(list(boxplot_tml, boxplot_tmb, boxplot_svtmb, boxplot_ploidy)), ncol=4)
ggsave(file = sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_HMF_tmb_etc"), plot = grid_plot, width = 12, height = 3)

#------rebuttal-----
PDO_WGS_screen_pair <- PDO_WGS_screen %>%
  filter(org_name == "OPT0032" | org_name == "RAS34")
new_row <- PDO_WGS_screen_pair[1, ]
new_row[,] <- NA
new_row$tml <- 87
new_row$tmbPerMb <-7.5
new_row$org_name <- "OPT0067"
PDO_WGS_screen_pair <- rbind(PDO_WGS_screen_pair, new_row)
#The tumor mutational load represents the total number of somatic missense variants acrossthe whole genome of the tumor.
#The tumor mutational burden score represents the number of all somatic variants across the whole genome of the tumor per Mb. 

pair_DBS <- ggplot(subset(PDO_WGS_screen_pair, org_name != "OPT0067"), aes(x = org_name, y = DBS5)) +
  geom_bar(stat = "identity", fill = "#7FBFF5") +
  labs(title = "DBS5\nduring treatment", x = "PDO", y = "Abs. contribution") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

long_data <- pivot_longer(PDO_WGS_screen_pair, cols = c(tml, tmbPerMb),
                          names_to = "metric", values_to = "value")
long_data$metric <- factor(long_data$metric, levels = c("tml", "tmbPerMb"))
long_data$metric <- recode(long_data$metric,
                           "tmbPerMb" = "tmb per Mb")
# Facet plot
pair_tml <- ggplot(long_data, aes(x = org_name, y = value, group = 1)) +
  geom_line(color = "black") +
  geom_point(color = "black") +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +  # 2 kolommen, naast elkaar
  labs(title = "Tml & tmb during treatment", x = "PDO", y = "value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pair <- grid.arrange(pair_DBS, pair_tml, ncol = 2, widths = c(1, 2))  # naast elkaar
ggsave(file = sprintf("%s/%s.pdf", WGS_plot_dir, "pair"), plot = pair, width = 6, height = 3)

#------WGS data PDO plots-----

#first test normality
PDO_WGS_screen %>%
  group_by(condition, chemo_naive) %>%
  summarise(
    W = shapiro.test(AUC)$statistic,        # Teststatistiek
    p_value = shapiro.test(AUC)$p.value     # P-waarde
  )
#also test equal variances
bartlett.test(AUC ~ condition, data = PDO_WGS_screen)

# Function to create a boxplot
make_boxplot <- function(data, x_var, y_var, comparisons, x_labels = NULL, x_lab = NULL, y_lab = NULL, 
                         label_y = 0.9, facet="chemo_condition") {
  plot <- ggplot(data, aes_string(x = x_var, y = y_var, color = x_var)) +
    geom_jitter(width = 0.2, size = 1) +
    labs(x = x_lab, y = y_lab) + 
    theme_classic() +
    stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5) +
    stat_compare_means(comparisons = comparisons, method = "t.test", aes(sprintf("p = %.3f", ..p..)), p.adj = "none", label.y = label_y) +
    scale_x_discrete(labels = x_labels)+
    scale_color_manual(values = c("#F8A29E", "#7FBFF5"))+
    stat_summary(fun = median, geom = "text", aes(label = format(..y..,)), vjust = -2, color = "black")
  if (facet == "condition") {
    plot <- plot + facet_grid(. ~ condition, scales = "free_y")
  } else if (facet == "chemo_naive") {
    plot <- plot + facet_grid(. ~ chemo_naive, scales = "free_y")
  } else if (facet == "chemo_condition") {
    plot <- plot + facet_grid(. ~ chemo_condition, scales = "free_y")
  } else if(facet == "no") {
    plot <- plot }
  return(plot)
}

# Function to create a correlation plot
make_corplot <- function(data, x_var, y_var, x_lab = NULL, y_lab = NULL, annotate_stats = FALSE, 
                         facet="chemo_condition", log_y = FALSE) {
  if (log_y) {
    data <- data %>%
      mutate_at(vars(y_var), log)
    y_lab <- paste0(y_lab, "(log)")
  }
  
  plot <- ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point(size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "#E0E1E0") +
    labs(x = x_lab, y = y_lab) + 
    theme_classic()
  if (facet == "condition") {
    plot <- plot + facet_grid(. ~ condition, scales = "free")
  } else if (facet == "chemo_naive") {
    plot <- plot + facet_grid(. ~ chemo_naive, scales = "free")
  } else if (facet == "chemo_condition") {
    plot <- plot + facet_grid(. ~ chemo_condition, scales = "free")
  } else if(facet == "no") {
    plot <- plot }
  if (annotate_stats) {
    plot <- plot + stat_cor(method = "spearman", size = 3, color = "black")
  }
  return(plot)
}

# new combined variable for facetting on both condition and chemonaive status
PDO_WGS_screen$chemo_condition <- interaction(PDO_WGS_screen$chemo_naive, PDO_WGS_screen$condition)
PDO_WGS_screen_top5$chemo_condition <- interaction(PDO_WGS_screen_top5$chemo_naive, PDO_WGS_screen_top5$condition)
PDO_WGS_screen_top5_all$chemo_condition <- interaction(PDO_WGS_screen_top5_all$chemo_naive, PDO_WGS_screen_top5_all$condition)

# Generate boxplots for genes of interest
genes <- c("SMAD4", "MACROD2", "PRKN", "KRAS", "PIK3CA", "TCF7L2", "PTEN", "FBXW7", "CCSER1", "FLCN", "FHIT")
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_gene_AUC"), width = 5*2, height = 2 * length(genes))
plots <- lapply(genes, function(gene) {
  make_boxplot(PDO_WGS_screen, x_var = sprintf("!is.na(%s)", gene), y_var = "AUC",
               comparisons = list(c("TRUE", "FALSE")), x_labels = c("Wt", "Mut"), x_lab = gene, y_lab = "AUC")
})
grid.arrange(grobs = plots, nrow = length(genes))
dev.off()

# Variables of interest
WGS_variable <- c("tml", "tmbPerMb", "svTumorMutationalBurden")
WGS_variable_status <- c("tmlStatus", "tmbStatus")

# Generate boxplots and corplots for TML, TMB, etc.
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_tmbetc"), width = 5, height = 2 * (2 * length(WGS_variable) + length(WGS_variable_status))) 
plots <- lapply(WGS_variable, function(var) {
  make_boxplot(PDO_WGS_screen_pretreatedfiltered, x_var = "AUC_category", y_var = var,
               comparisons = list(c("resistant", "sensitive")), x_lab = "PDO", y_lab = var, facet = "condition")
})
plots_status <- lapply(WGS_variable_status, function(var) {
  make_boxplot(PDO_WGS_screen_pretreatedfiltered, x_var = var, y_var = "AUC", 
               comparisons = list(c("HIGH", "LOW")), x_lab = var, y_lab = "AUC", label_y = 0.9, facet = "condition")
})
corplots <- lapply(WGS_variable, function(var) {
  make_corplot(PDO_WGS_screen_pretreatedfiltered, x_var = "AUC", y_var = var, 
               x_lab = "AUC", y_lab = var, annotate_stats = TRUE, log_y = TRUE, facet = "condition")
})
grid.arrange(grobs = c(plots, plots_status, corplots), nrow = 2 * length(WGS_variable) + length(WGS_variable_status))
dev.off()

# Generate boxplots and corplots for TML, TMB, etc.
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_tmbetc_all"), width = 5, height = 2 * (2 * length(WGS_variable) + length(WGS_variable_status))) 
plots <- lapply(WGS_variable, function(var) {
  make_boxplot(PDO_WGS_screen, x_var = "AUC_category", y_var = var,
               comparisons = list(c("resistant", "sensitive")), x_lab = "PDO", y_lab = var, facet = "condition")
})
plots_status <- lapply(WGS_variable_status, function(var) {
  make_boxplot(PDO_WGS_screen, x_var = var, y_var = "AUC", 
               comparisons = list(c("HIGH", "LOW")), x_lab = var, y_lab = "AUC", label_y = 0.9, facet = "condition")
})
corplots <- lapply(WGS_variable, function(var) {
  make_corplot(PDO_WGS_screen, x_var = "AUC", y_var = var, 
               x_lab = "AUC", y_lab = var, annotate_stats = TRUE, log_y = TRUE, facet = "condition")
})
grid.arrange(grobs = c(plots, plots_status, corplots), nrow = 2 * length(WGS_variable) + length(WGS_variable_status))
dev.off()

#Variables of interest
WGS_variable_karyotype <- colnames(PDO_WGS_screen)[startsWith(colnames(PDO_WGS_screen), "chr")]
WGS_variable_karyotype_top5 <- colnames(PDO_WGS_screen_top5)[startsWith(colnames(PDO_WGS_screen_top5), "chr")]

PDO_WGS_screen_pheatmap <- subset(PDO_WGS_screen, select = c("sampleId", "chemo_naive", "condition", WGS_variable_karyotype))
PDO_WGS_screen_pheatmap_5FU <- PDO_WGS_screen_pheatmap %>%
  filter(condition == "5-FU")
PDO_WGS_screen_pheatmap_5FU_ordered <- PDO_WGS_screen_pheatmap_5FU[order(PDO_WGS_screen_pheatmap_5FU$chemo_naive), ]

rownames(PDO_WGS_screen_pheatmap_5FU_ordered) <- PDO_WGS_screen_pheatmap_5FU_ordered$sampleId
PDO_WGS_screen_pheatmap_5FU_ordered$sampleId <- NULL
heatmap_annotation <- PDO_WGS_screen_pheatmap_5FU_ordered[ c("chemo_naive")]
colnames(heatmap_annotation)[1]<- "Treatment status" 
PDO_WGS_screen_pheatmap_5FU_ordered<- PDO_WGS_screen_pheatmap_5FU_ordered %>% dplyr::select(-chemo_naive, -condition)

heatmap_annotation$`Treatment status`[heatmap_annotation$`Treatment status` == "yes"] <- "chemonaive"
heatmap_annotation$`Treatment status`[heatmap_annotation$`Treatment status` == "no"] <- "pretreated"

pheatmap_karyotype_5FU <- pheatmap(PDO_WGS_screen_pheatmap_5FU_ordered,
                     cluster_cols = FALSE,
                     cluster_rows = FALSE,
                     #legend_breaks = c(0, 1),
                     legend_labels = c("sensitive", "resistant"),
                     color = colorRampPalette(c("#8B8BFF","white", "#FFD0D0", "#FFB9B9", "#FFA2A2", "#FF8B8B", "#FF0000", "#FF0000", "#FF0000"))(50),
                     main = "         Organoid sensitivity (AUC)\n         per treatment",
                     angle_col = "315",
                     row_names_rot = 90,
                     annotation_row = heatmap_annotation,
                     annotation_colors = list("Treatment status" = c("chemonaive" = "#7FBFF5", "pretreated" = "#F8A29E")))

my_color_fun <- colorRampPalette(c("blue","#FFFFFF", "#FF0000"))

pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "pheatmap_karyotype"), width = 10, height = 7) 
pheatmap_karyotype_5FU
dev.off()

# Generate boxplots for Palliative oxaliplatin PD
PDO_metadata_screen$Palliative_oxaliplatin_PD <- ifelse(PDO_metadata_screen$Palliative_oxaliplatin_response_all == "PD" & PDO_metadata_screen$chemo_naive == "no", "PD",
                                                   ifelse(PDO_metadata_screen$chemo_naive == "no", "no PD", NA))
PDO_metadata_screen$Palliative_oxaliplatin_PD <- ifelse(is.na(PDO_metadata_screen$Palliative_oxaliplatin_response_all) & PDO_metadata_screen$chemo_naive == "no", "no PD",
                                                   PDO_metadata_screen$Palliative_oxaliplatin_PD)
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_PD_oxaliplatin_all"), width = 5, height = 2)
plots_PD_oxaliplatin <- lapply(c("Palliative_oxaliplatin_PD"), function(var) {
  make_boxplot(PDO_metadata_screen, x_var = var, y_var = "AUC", 
               comparisons = list(c("PD", "no PD")), x_lab = var, y_lab = "AUC", label_y = 0.9, facet = "condition")
})
grid.arrange(grobs = plots_PD_oxaliplatin, nrow = 1)
dev.off()

# Generate boxplots for Palliative oxaliplatin PD
PDO_WGS_screen_pretreatedfiltered$Palliative_oxaliplatin_PD <- ifelse(PDO_WGS_screen_pretreatedfiltered$Palliative_oxaliplatin_response_all == "PD" & PDO_WGS_screen_pretreatedfiltered$chemo_naive == "no", "PD",
                                                        ifelse(PDO_WGS_screen_pretreatedfiltered$chemo_naive == "no", "no PD", NA))
PDO_WGS_screen_pretreatedfiltered$Palliative_oxaliplatin_PD <- ifelse(is.na(PDO_WGS_screen_pretreatedfiltered$Palliative_oxaliplatin_response_all) & PDO_WGS_screen_pretreatedfiltered$chemo_naive == "no", "no PD",
                                                                      PDO_WGS_screen_pretreatedfiltered$Palliative_oxaliplatin_PD)
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_PD_oxaliplatin"), width = 5, height = 2)
plots_PD_oxaliplatin <- lapply(c("Palliative_oxaliplatin_PD"), function(var) {
  make_boxplot(PDO_WGS_screen_pretreatedfiltered, x_var = var, y_var = "AUC", 
               comparisons = list(c("PD", "no PD")), x_lab = "Clinical oxaliplatin response", y_lab = "AUC", label_y = 0.9, facet = "condition")
})
grid.arrange(grobs = plots_PD_oxaliplatin, nrow = 1)
dev.off()

# Generate boxplots for signature variables
signature_variables <- c("SBS17a_cat", "SBS17b_cat", "has_deep_deletion", "SBS2_cat", "SBS13_cat")
signature_variables_cont <- c("SBS17a", "SBS17b", "SBS35", "SBS17", "SBS2", "SBS13", "DBS5", "SBS25")
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_signatures"), width = 5, 
    height = 2 * length(signature_variables)+2*2*length(signature_variables_cont))
plots_signature <- lapply(signature_variables, function(var) {
  make_boxplot(PDO_WGS_screen_pretreatedfiltered, x_var = var, y_var = "AUC", 
               comparisons = list(c("present", "absent")), x_lab = var, y_lab = "AUC", label_y = 0.9, facet ="condition")
})
plots_signature_cont <- lapply(signature_variables_cont, function(var) {
  make_boxplot(PDO_WGS_screen_pretreatedfiltered, x_var = "AUC_category", y_var = var, 
               comparisons = list(c("resistant", "sensitive")), x_lab = "PDO", y_lab = var, facet ="condition")
})
corplots_signature_cont <- lapply(signature_variables_cont, function(var) {
  make_corplot(PDO_WGS_screen_pretreatedfiltered, x_var = "AUC", y_var = var, 
               x_lab = "PDO", y_lab = var, annotate_stats = TRUE, facet ="condition")
})
grid.arrange(grobs = c(plots_signature, plots_signature_cont, corplots_signature_cont), 
             nrow = length(signature_variables)+length(signature_variables_cont)+length(corplots_signature_cont))
dev.off()

# Generate boxplots for signature variables chemonaive
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_signatures_chemonaive"), width = 5, 
    height = 2 * length(signature_variables)+2*2*length(signature_variables_cont))
plots_signature <- lapply(signature_variables, function(var) {
  make_boxplot(PDO_WGS_screen_chemonaivefiltered, x_var = var, y_var = "AUC", 
               comparisons = list(c("present", "absent")), x_lab = var, y_lab = "AUC", label_y = 0.9, facet ="condition")
})
plots_signature_cont <- lapply(signature_variables_cont, function(var) {
  make_boxplot(PDO_WGS_screen_chemonaivefiltered, x_var = "AUC_category", y_var = var, 
               comparisons = list(c("resistant", "sensitive")), x_lab = "PDO", y_lab = var, facet ="condition")
})
corplots_signature_cont <- lapply(signature_variables_cont, function(var) {
  make_corplot(PDO_WGS_screen_chemonaivefiltered, x_var = "AUC", y_var = var, 
               x_lab = "PDO", y_lab = var, annotate_stats = TRUE, facet ="condition")
})
grid.arrange(grobs = c(plots_signature, plots_signature_cont, corplots_signature_cont), 
             nrow = length(signature_variables)+length(signature_variables_cont)+length(corplots_signature_cont))
dev.off()

# Generate boxplots for signature variables chemonaive vs pretreated
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_signatures_all"), width = 3, 
    height = 2 * length(signature_variables_cont))
plots_signature <- lapply(signature_variables_cont, function(var) {
  make_boxplot(mat_HMF_treatment, x_var = "`Sample type`", y_var = var, 
               comparisons = list(c("chemonaive", "pretreated")), x_lab = "Treatment status", y_lab = var, label_y = 0.9, facet = "no")
})
grid.arrange(grobs = c(plots_signature), 
             nrow = length(signature_variables_cont))
dev.off()

# Generate boxplots for signature variables top 5
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_signatures_top5"), width = 5, 
    height = 2 * length(signature_variables)+2*length(signature_variables_cont))
plots_signature <- lapply(signature_variables, function(var) {
  make_boxplot(PDO_WGS_screen_top5_pretreated, x_var = var, y_var = "AUC", facet ="condition", 
               comparisons = list(c("present", "absent")), x_lab = var, y_lab = "AUC top 5", label_y = 0.9)
})
plots_signature_cont <- lapply(signature_variables_cont, function(var) {
  make_boxplot(PDO_WGS_screen_top5_pretreated, x_var = "AUC_category", y_var = var, facet ="condition",
               comparisons = list(c("resistant", "sensitive")), x_lab = "PDO top 5", y_lab = var)
})
grid.arrange(grobs = c(plots_signature, plots_signature_cont), 
             nrow = length(signature_variables)+length(signature_variables_cont))
dev.off()

# Generate boxplots for signature variables
signature_variables <- c("SBS2_cat", "SBS13_cat")
signature_variables_cont <- c("SBS2", "SBS13")
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_signatures_all"), width = 5, 
    height = 2 * length(signature_variables)+2*2*length(signature_variables_cont))
plots_signature <- lapply(signature_variables, function(var) {
  make_boxplot(PDO_WGS_screen_pretreatedfiltered, x_var = var, y_var = "AUC", 
               comparisons = list(c("present", "absent")), x_lab = var, y_lab = "AUC", label_y = 0.9, facet ="condition")
})
plots_signature_cont <- lapply(signature_variables_cont, function(var) {
  make_boxplot(PDO_WGS_screen_pretreatedfiltered, x_var = "AUC_category", y_var = var, 
               comparisons = list(c("resistant", "sensitive")), x_lab = "PDO", y_lab = var, facet ="condition")
})
corplots_signature_cont <- lapply(signature_variables_cont, function(var) {
  make_corplot(PDO_WGS_screen_pretreatedfiltered, x_var = "AUC", y_var = var, 
               x_lab = "PDO", y_lab = var, annotate_stats = TRUE, facet ="condition")
})
grid.arrange(grobs = c(plots_signature, plots_signature_cont, corplots_signature_cont), 
             nrow = length(signature_variables)+length(signature_variables_cont)+length(corplots_signature_cont))
dev.off()

# chemonaive and pretreated together 
plots_signature <- lapply(signature_variables, function(var) {
make_boxplot(PDO_WGS_screen, x_var = var, y_var = "AUC", comparisons = list(c("present", "absent")), 
             x_lab = var, y_lab = "AUC", label_y = 0.9, facet="condition")})

# Check deep deletions in chemonaive group
PDO_WGS_screen %>%
  filter(has_deep_deletion == "present" & chemo_naive == "yes") %>%
  dplyr::select(sampleId) %>%
  distinct() #MetNa03 MetNa11 MetNa16 MetNa17 #MetNa23


PDO_WGS_screen %>%
  filter(has_deep_deletion == "present" & chemo_naive == "no") %>%
  dplyr::select(sampleId) %>%
  distinct()
#44%

PDO_WGS_screen_cont <- PDO_WGS_screen %>%
  dplyr::select(sampleId, chemo_naive, has_deep_deletion) %>%
  distinct()

# Create a contingency table
contingency_table <- table(PDO_WGS_screen_cont$chemo_naive, PDO_WGS_screen_cont$has_deep_deletion)
print(contingency_table)
chisq.test(contingency_table)

#Variables of interest
WGS_variable <- c("tml", "tmbPerMb", "svTumorMutationalBurden")

# Generate corplots for signatures vs treatment cycles
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "plots_tmb_sig"), width = 16, height = 2*(length(WGS_variable))) 
boxplots_SBS17a <- lapply(WGS_variable, function(var) {
  make_boxplot(mat_HMF_treatment, x_var = "SBS17a_cat", y_var = var,
               comparisons = list(c("absent", "present")), x_lab = "signature", y_lab = var, facet="no")})
boxplots_SBS17b <- lapply(WGS_variable, function(var) {
  make_boxplot(mat_HMF_treatment, x_var = "SBS17b_cat", y_var = var,
               comparisons = list(c("absent", "present")), x_lab = "signature", y_lab = var, facet="no")})
boxplots_DBS5 <- lapply(WGS_variable, function(var) {
  make_boxplot(mat_HMF_treatment, x_var = "DBS5_cat", y_var = var,
               comparisons = list(c("absent", "present")), x_lab = "signature", y_lab = var, facet="no")})
boxplots_SBS35 <- lapply(WGS_variable, function(var) {
  make_boxplot(mat_HMF_treatment, x_var = "SBS35_cat", y_var = var,
               comparisons = list(c("absent", "present")), x_lab = "signature", y_lab = var, facet="no")})

corplots_SBS17a <- lapply(WGS_variable, function(var) {
  make_corplot(mat_HMF_treatment, x_var = "SBS17a", y_var = var, 
               x_lab = "SBS17a", y_lab = var, annotate_stats=TRUE, facet="no")})
corplots_SBS17b <- lapply(WGS_variable, function(var) {
  make_corplot(mat_HMF_treatment, x_var = "SBS17b", y_var = var, 
               x_lab = "SBS17b", y_lab = var, annotate_stats=TRUE, facet="no")})
corplots_DBS5 <- lapply(WGS_variable, function(var) {
  make_corplot(mat_HMF_treatment, x_var = "DBS5", y_var = var, 
               x_lab = "DBS5", y_lab = var, annotate_stats=TRUE, facet="no")})
corplots_SBS35 <- lapply(WGS_variable, function(var) {
  make_corplot(mat_HMF_treatment, x_var = "SBS35", y_var = var, 
               x_lab = "SBS35", y_lab = var, annotate_stats=TRUE, facet="no")})
grid.arrange(grobs = c(boxplots_SBS17a, corplots_SBS17a, boxplots_SBS17b,corplots_SBS17b,
                       boxplots_DBS5, corplots_DBS5, boxplots_SBS35, corplots_SBS35), 
             nrow = 4, ncol = 8)
dev.off()


# Generate corplots for signatures vs treatment cycles
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "plots_tmb_sig_pretreated"), width = 16, height = 2*(length(WGS_variable))) 
boxplots_SBS17a <- lapply(WGS_variable, function(var) {
  make_boxplot(PDO_WGS_screen_pretreatedfiltered, x_var = "SBS17a_cat", y_var = var,
               comparisons = list(c("absent", "present")), x_lab = "signature", y_lab = var, facet="no")})
boxplots_SBS17b <- lapply(WGS_variable, function(var) {
  make_boxplot(PDO_WGS_screen_pretreatedfiltered, x_var = "SBS17b_cat", y_var = var,
               comparisons = list(c("absent", "present")), x_lab = "signature", y_lab = var, facet="no")})
boxplots_DBS5 <- lapply(WGS_variable, function(var) {
  make_boxplot(PDO_WGS_screen_pretreatedfiltered, x_var = "DBS5_cat", y_var = var,
               comparisons = list(c("absent", "present")), x_lab = "signature", y_lab = var, facet="no")})
boxplots_SBS35 <- lapply(WGS_variable, function(var) {
  make_boxplot(PDO_WGS_screen_pretreatedfiltered, x_var = "SBS35_cat", y_var = var,
               comparisons = list(c("absent", "present")), x_lab = "signature", y_lab = var, facet="no")})

corplots_SBS17a <- lapply(WGS_variable, function(var) {
  make_corplot(PDO_WGS_screen_pretreatedfiltered, x_var = "SBS17a", y_var = var, 
               x_lab = "SBS17a", y_lab = var, annotate_stats=TRUE, facet="no")})
corplots_SBS17b <- lapply(WGS_variable, function(var) {
  make_corplot(PDO_WGS_screen_pretreatedfiltered, x_var = "SBS17b", y_var = var, 
               x_lab = "SBS17b", y_lab = var, annotate_stats=TRUE, facet="no")})
corplots_DBS5 <- lapply(WGS_variable, function(var) {
  make_corplot(PDO_WGS_screen_pretreatedfiltered, x_var = "DBS5", y_var = var, 
               x_lab = "DBS5", y_lab = var, annotate_stats=TRUE, facet="no")})
corplots_SBS35 <- lapply(WGS_variable, function(var) {
  make_corplot(PDO_WGS_screen_pretreatedfiltered, x_var = "SBS35", y_var = var, 
               x_lab = "SBS35", y_lab = var, annotate_stats=TRUE, facet="no")})
grid.arrange(grobs = c(boxplots_SBS17a, corplots_SBS17a, boxplots_SBS17b,corplots_SBS17b,
                       boxplots_DBS5, corplots_DBS5, boxplots_SBS35, corplots_SBS35), 
             nrow = 4, ncol = 8)
dev.off()



#clean cycles data, make it integer
mat_HMF_treatment$Palliative_oxaliplatin_cycles<- ifelse(mat_HMF_treatment$Palliative_oxaliplatin_cycles == "No", 0, mat_HMF_treatment$Palliative_oxaliplatin_cycles)
columns_to_replace <- c("Palliative_5FU_cycles", "Palliative_oxaliplatin_cycles", "Palliative_irinotecan_cycles")
for (col in columns_to_replace) {  mat_HMF_treatment[[col]][mat_HMF_treatment[[col]] == "NA"] <- NA}
mat_HMF_treatment[columns_to_replace] <- lapply(mat_HMF_treatment[columns_to_replace], as.integer)

mat_HMF_treatment$Palliative_cycles <- mat_HMF_treatment$Palliative_irinotecan_cycles+
  mat_HMF_treatment$Palliative_oxaliplatin_cycles+
  mat_HMF_treatment$Palliative_5FU_cycles #make 1 variable with all cycles

PDO_WGS_screen_pretreatedfiltered[columns_to_replace] <- lapply(PDO_WGS_screen_pretreatedfiltered[columns_to_replace], as.integer)

PDO_WGS_screen_pretreatedfiltered$Palliative_cycles <- PDO_WGS_screen_pretreatedfiltered$Palliative_irinotecan_cycles+
  PDO_WGS_screen_pretreatedfiltered$Palliative_oxaliplatin_cycles+
  PDO_WGS_screen_pretreatedfiltered$Palliative_5FU_cycles #make 1 variable with all cycles

#Variables of interest
WGS_variable <- colnames(mat_HMF_treatment)[endsWith(colnames(mat_HMF_treatment), "cycles")]

# Generate corplots for signatures vs treatment cycles
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "plots_cycles_sig"), width = 16, height = 2.5*length(WGS_variable)) 
boxplots_SBS17a <- lapply(WGS_variable, function(var) {
  make_boxplot(mat_HMF_treatment, x_var = "SBS17a_cat", y_var = var,
               comparisons = list(c("absent", "present")), x_lab = "signature", y_lab = var, facet="no")})
boxplots_SBS17b <- lapply(WGS_variable, function(var) {
  make_boxplot(mat_HMF_treatment, x_var = "SBS17b_cat", y_var = var,
               comparisons = list(c("absent", "present")), x_lab = "signature", y_lab = var, facet="no")})
boxplots_DBS5 <- lapply(WGS_variable, function(var) {
  make_boxplot(mat_HMF_treatment, x_var = "DBS5_cat", y_var = var,
               comparisons = list(c("absent", "present")), x_lab = "signature", y_lab = var, facet="no")})
boxplots_SBS35 <- lapply(WGS_variable, function(var) {
  make_boxplot(mat_HMF_treatment, x_var = "SBS35_cat", y_var = var,
               comparisons = list(c("absent", "present")), x_lab = "signature", y_lab = var, facet="no")})
boxplots_SBS25 <- lapply(WGS_variable, function(var) {
  make_boxplot(mat_HMF_treatment, x_var = "SBS25_cat", y_var = var,
               comparisons = list(c("absent", "present")), x_lab = "signature", y_lab = var, facet="no")})


corplots_SBS17a <- lapply(WGS_variable, function(var) {
  make_corplot(mat_HMF_treatment, x_var = "SBS17a", y_var = var, 
               x_lab = "SBS17a", y_lab = var, annotate_stats=TRUE, facet="no")})
corplots_SBS17b <- lapply(WGS_variable, function(var) {
  make_corplot(mat_HMF_treatment, x_var = "SBS17b", y_var = var, 
               x_lab = "SBS17b", y_lab = var, annotate_stats=TRUE, facet="no")})
corplots_DBS5 <- lapply(WGS_variable, function(var) {
  make_corplot(mat_HMF_treatment, x_var = "DBS5", y_var = var, 
               x_lab = "DBS5", y_lab = var, annotate_stats=TRUE, facet="no")})
corplots_SBS35 <- lapply(WGS_variable, function(var) {
  make_corplot(mat_HMF_treatment, x_var = "SBS35", y_var = var, 
               x_lab = "SBS35", y_lab = var, annotate_stats=TRUE, facet="no")})
corplots_SBS25 <- lapply(WGS_variable, function(var) {
  make_corplot(mat_HMF_treatment, x_var = "SBS25", y_var = var, 
               x_lab = "SBS25", y_lab = var, annotate_stats=TRUE, facet="no")})
grid.arrange(grobs = c(boxplots_SBS17a, corplots_SBS17a, boxplots_SBS17b,corplots_SBS17b,
                       boxplots_DBS5, corplots_DBS5, boxplots_SBS35, corplots_SBS35, boxplots_SBS25, corplots_SBS25), 
             nrow = 5, ncol = 8)
dev.off()

# Generate corplots treatment cycles vs sensitivity & TMB
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "plots_cycles_sens_tmb"), width = 8, height = 1.5*length(WGS_variable)) 
corplots_cycles_tml <- lapply(WGS_variable, function(var) {
  make_corplot(mat_HMF_treatment, x_var = "tml", y_var = var,annotate_stats=TRUE,
               x_lab = "tml", y_lab = var, facet="no")})
corplots_cycles_tmbPerMb <- lapply(WGS_variable, function(var) {
  make_corplot(mat_HMF_treatment, x_var = "tmbPerMb", y_var = var,annotate_stats=TRUE,
               x_lab = "tmbPerMb", y_lab = var, facet="no")})
corplots_cycles_svTumorMutationalBurden <- lapply(WGS_variable, function(var) {
  make_corplot(mat_HMF_treatment, x_var = "svTumorMutationalBurden", y_var = var,annotate_stats=TRUE,
               x_lab = "svTumorMutationalBurden", y_lab = var, facet="no")})

grid.arrange(grobs = c(corplots_cycles_tml, corplots_cycles_tmbPerMb, corplots_cycles_svTumorMutationalBurden), 
             ncol = 4)
dev.off()

corplots_cycles_sens <- lapply(WGS_variable, function(var) {
  make_corplot(PDO_WGS_screen_pretreatedfiltered, x_var = "AUC", y_var = var,annotate_stats=TRUE,
               x_lab = "AUC", y_lab = var, facet="condition")}) #geen correlaties
corplots_cycles_sens

# Make boxplots for TMB/TML/WG disruptions in PDO data (validation set)
boxplot_tml <- ggplot(mat_HMF_treatment, aes(y = tml, x = `Sample type`, color=`Sample type`)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "tml", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "t.test", label = "p.format", label.y = 1)+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)+
  scale_color_manual(values = c("#7FBFF5", "#F8A29E"))+
  theme(legend.position = "none")+
  stat_summary(fun = median, geom = "text", aes(label = format(..y..,)), vjust = -2, color = "black")
boxplot_tml

#ggsave(file=sprintf("%s/%s.pdf",WGS_plot_dir, "boxplot_tml"), width = (3), height = (3)) #save boxplot TLM
#dev.off()

boxplot_tmb <- ggplot(mat_HMF_treatment, aes(y = tmbPerMb, x = `Sample type`, color=`Sample type`)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "tmb per Mb", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "t.test", label = "p.format", label.y = 0.1)+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)+
  scale_color_manual(values = c("#7FBFF5", "#F8A29E"))+
  theme(legend.position = "none")+
  stat_summary(fun = median, geom = "text", aes(label = format(..y..,)), vjust = -2, color = "black")
boxplot_tmb

#ggsave(file=sprintf("%s/%s.pdf",WGS_plot_dir, "boxplot_tmb"), width = (3), height = (3)) #save boxplot TMB
#dev.off()

boxplot_svtmb <- ggplot(mat_HMF_treatment, aes(y = svTumorMutationalBurden, x = `Sample type`, color=`Sample type`)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "sv tmb", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "t.test", label = "p.format", label.y = 1)+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)+
  scale_color_manual(values = c("#7FBFF5", "#F8A29E"))+
  theme(legend.position = "none")+
  stat_summary(fun = median, geom = "text", aes(label = format(..y..,)), vjust = -2, color = "black")
boxplot_svtmb

#ggsave(file=sprintf("%s/%s.pdf",WGS_plot_dir, "boxplot_svtmb"), width = (3), height = (3)) #save boxplot svTMB
#dev.off()

#contingency_table for whole genome duplications in Hartwig data (validation set)
contingency_table <- table(mat_HMF_treatment$wholeGenomeDuplication, mat_HMF_treatment$`Sample type` == "pretreated")
dimnames(contingency_table) <- list(wholeGenomeDuplication = c("FALSE", "TRUE"), chemotherapy_pretreated = c("FALSE", "TRUE"))
fisher_result <- fisher.test(contingency_table) #not significant
contingency_table

boxplot_ploidy <- ggplot(mat_HMF_treatment, aes(y = ploidy, x = `Sample type`, color=`Sample type`)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "ploidy", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "t.test", label = "p.format", label.y = 1)+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)+
  scale_color_manual(values = c("#7FBFF5", "#F8A29E"))+
  theme(legend.position = "none")+
  stat_summary(fun = median, geom = "text", aes(label = format(..y..,)), vjust = -2, color = "black")
boxplot_ploidy

#ggsave(file=sprintf("%s/%s.pdf",WGS_plot_dir, "boxplot_ploidy"), width = (3), height = (3)) #save boxplot ploidy
#dev.off()

mat_HMF_treatment$treatment <- mat_HMF_treatment$`Sample type`

generate_boxplot <- function(data, y_var, x_var, y_label, x_label) {
  plot <- ggplot(data, aes_string(y = y_var, x = x_var, color=x_var)) + 
    geom_jitter(width = 0.2, size = 0.5) +
    labs(title = "", y = y_label, x = x_label) + 
    theme_classic() +
    scale_y_continuous(trans = "log10") +
    stat_compare_means(method = "t.test", label = "p.format") +
    stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)+
    scale_color_manual(values = c("#7FBFF5", "#F8A29E"))+
    theme(legend.position = "none")+
    stat_summary(fun = median, geom = "text", aes(label = format(..y..,)), vjust = -2, color = "black")
  return(plot)
}

# Generate the individual plots
boxplot_tml <- generate_boxplot(mat_HMF_treatment, "tml", "treatment", "tml", "chemotherapy")
boxplot_tmb <- generate_boxplot(mat_HMF_treatment, "tmbPerMb", "treatment", "tmb per Mb", "chemotherapy")
boxplot_svtmb <- generate_boxplot(mat_HMF_treatment, "svTumorMutationalBurden", "treatment", "sv tmb", "chemotherapy")
boxplot_ploidy <- generate_boxplot(mat_HMF_treatment, "ploidy", "treatment", "ploidy", "chemotherapy")

grid_plot <- grid.arrange(grobs = c(list(boxplot_tml, boxplot_tmb, boxplot_svtmb, boxplot_ploidy)), ncol=4)
ggsave(file = sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_PDO_tmb_etc"), plot = grid_plot, width = 12, height = 3)

total_chemonaive <- sum(PDO_WGS_screen$`Sample type` == "chemonaive")
count_condition_met <- sum(PDO_WGS_screen$has_deep_deletion == "present" & PDO_WGS_screen$`Sample type` == "chemonaive")
percentage <- (count_condition_met / total_chemonaive) * 100
print(percentage)

total_pretreated <- sum(PDO_WGS_screen$`Sample type` == "pretreated")
count_condition_met <- sum(PDO_WGS_screen$has_deep_deletion == "present" & PDO_WGS_screen$`Sample type` == "pretreated")
percentage <- (count_condition_met / total_pretreated) * 100
print(percentage)

data <- matrix(c(
  5, 12,  # Chemonaive met deletie, Chemonaive zonder deletie
  8, 10   # Pretreated met deletie, Pretreated zonder deletie
), nrow = 2, byrow = TRUE)

# Zet de matrix om in een tabel met rijnamen en kolomnamen
colnames(data) <- c("Deletion", "No Deletion")
rownames(data) <- c("Chemonaive", "Pretreated")

# Bekijk de tabel
print(data)
# Voer de chi-kwadraat test uit
test_result <- chisq.test(data)

# Bekijk de resultaten
print(test_result)
