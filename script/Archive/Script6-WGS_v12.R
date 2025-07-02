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
WGS_plot_dir <- file.path(home_dir, "Analyse 2.0/7_WGS")

#------WGS data Hartwig (validation set)-----

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
HMF_treatment_filter <- HMF_treatment_filter %>% rename(sampleId = patientId) #change patientId to sampleId for merge
purple_treatment <- merge(HMF_treatment_filter, purple, by = "sampleId") #merge dataframes
purple_treatment <- purple_treatment %>% filter(msStatus == "MSS") 
combined <- merge(purple_treatment, drivers, by = "sampleId") #merge dataframes

result <- anti_join(HMF_treatment_filter, purple, by = "sampleId") #samples not present in both datasets


#------WGS data Hartwig plots-----
# Make boxplots for TMB/TML/WG disruptions in Hartwig data (validation set)
boxplot_tml <- ggplot(purple_treatment, aes(y = tml, x = chemotherapy)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "tml", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "wilcox.test", label = "p.format")+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)
boxplot_tml

boxplot_tmb <- ggplot(purple_treatment, aes(y = tmbPerMb, x = chemotherapy)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "tmb per Mb", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "wilcox.test", label = "p.format")+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)
boxplot_tmb

boxplot_svtmb <- ggplot(purple_treatment, aes(y = svTumorMutationalBurden, x = chemotherapy)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "sv tmb", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "wilcox.test", label = "p.format")+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)
boxplot_svtmb

#contingency_table for whole genome duplications in Hartwig data (validation set)
contingency_table <- table(purple_treatment$wholeGenomeDuplication, purple_treatment$chemotherapy == "pretreated")
dimnames(contingency_table) <- list(wholeGenomeDuplication = c("FALSE", "TRUE"), chemotherapy_pretreated = c("FALSE", "TRUE"))
fisher_result <- fisher.test(contingency_table) #not significant
contingency_table

boxplot_ploidy <- ggplot(purple_treatment, aes(y = ploidy, x = chemotherapy)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "ploidy", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "wilcox.test", label = "p.format")+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)
boxplot_ploidy

grid_plot <- grid.arrange(grobs = c(list(boxplot_tml, boxplot_tmb, boxplot_svtmb, boxplot_ploidy)), ncol=4)
ggsave(file = sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_HMF_tmb_etc"), plot = grid_plot, width = 12, height = 3)

#------import WGS data PDO made in script 9-----
mat_HMF_treatment <- readRDS(file.path(WGS_plot_dir, "mat_HMF_treatment.rds"))
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

#------WGS data PDO plots-----
# Function to create a boxplot
make_boxplot <- function(data, x_var, y_var, comparisons, x_labels = NULL, x_lab = NULL, y_lab = NULL, label_y = 1.1, facet="chemo_condition") {
  plot <- ggplot(data, aes_string(x = x_var, y = y_var, color = x_var)) +
    geom_jitter(width = 0.2, size = 1) +
    labs(x = x_lab, y = y_lab) + 
    theme_classic() +
    stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5) +
    stat_compare_means(comparisons = comparisons, method = "t.test", aes(sprintf("p = %.3f", ..p..)), p.adj = "none", label.y = label_y) +
    scale_x_discrete(labels = x_labels)
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
make_corplot <- function(data, x_var, y_var, x_lab = NULL, y_lab = NULL, annotate_stats = FALSE, facet="chemo_condition", log_y = FALSE) {
  if (log_y) {
    data <- data %>%
      mutate_at(vars(y_var), log)
    y_lab <- paste0(y_lab, "(log)")
  }
  
  plot <- ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point(size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
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
    plot <- plot + stat_cor(method = "spearman", size = 3, color = "red")
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

#Variables of interest
WGS_variable_karyotype <- colnames(PDO_WGS_screen)[startsWith(colnames(PDO_WGS_screen), "chr")]
WGS_variable_karyotype_top5 <- colnames(PDO_WGS_screen_top5)[startsWith(colnames(PDO_WGS_screen_top5), "chr")]

PDO_WGS_screen_pheatmap <- subset(PDO_WGS_screen, select = c("sampleId", "chemo_naive", "AUC", "condition", WGS_variable_karyotype))
PDO_WGS_screen_pheatmap_5FU <- PDO_WGS_screen_pheatmap %>%
  filter(condition == "5-FU")
PDO_WGS_screen_pheatmap_5FU_ordered <- PDO_WGS_screen_pheatmap_5FU[order(PDO_WGS_screen_pheatmap_5FU$chemo_naive, PDO_WGS_screen_pheatmap_5FU$AUC), ]

rownames(PDO_WGS_screen_pheatmap_5FU_ordered) <- PDO_WGS_screen_pheatmap_5FU_ordered$sampleId
PDO_WGS_screen_pheatmap_5FU_ordered$sampleId <- NULL
heatmap_annotation <- PDO_WGS_screen_pheatmap_5FU_ordered[ c("chemo_naive", "AUC")]
PDO_WGS_screen_pheatmap_5FU_ordered<- PDO_WGS_screen_pheatmap_5FU_ordered %>% select(-chemo_naive, -AUC, -condition)

pheatmap_karyotype_5FU <- pheatmap(PDO_WGS_screen_pheatmap_5FU_ordered,
                     cluster_cols = FALSE,
                     cluster_rows = FALSE,
                     #legend_breaks = c(0, 1),
                     legend_labels = c("sensitive", "resistant"),
                     #color = colorRampPalette(c("#7FBFF5","#005F8F","#012033","#321614","#8C423D","#F8A29E"))(50),
                     main = "         Organoid sensitivity (AUC)\n         per treatment",
                     angle_col = "315",
                     row_names_rot = 90,
                     annotation_row = heatmap_annotation,)

pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "pheatmap_karyotype_5FU"), width = 10, height = 7) 
pheatmap_karyotype_5FU
dev.off()

# Generate boxplots for Palliative oxaliplatin PD
PDO_metadata_screen$Palliative_oxaliplatin_PD <- ifelse(PDO_metadata_screen$Palliative_oxaliplatin_response_all == "PD" & PDO_metadata_screen$chemo_naive == "no", "PD",
                                                   ifelse(PDO_metadata_screen$chemo_naive == "no", "no PD", NA))
PDO_metadata_screen$Palliative_oxaliplatin_PD <- ifelse(is.na(PDO_metadata_screen$Palliative_oxaliplatin_response_all) & PDO_metadata_screen$chemo_naive == "no", "no PD",
                                                   PDO_metadata_screen$Palliative_oxaliplatin_PD)
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_PD_oxaliplatin"), width = 19, height = 2)
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
               comparisons = list(c("PD", "no PD")), x_lab = var, y_lab = "AUC", label_y = 0.9, facet = "condition")
})
grid.arrange(grobs = plots_PD_oxaliplatin, nrow = 1)
dev.off()

# Generate boxplots for signature variables
signature_variables <- c("SBS17a_cat", "SBS17b_cat", "has_deep_deletion", "SBS2_cat", "SBS13_cat")
signature_variables_cont <- c("SBS17a", "SBS17b", "SBS35", "SBS17", "SBS2", "SBS13")
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
  select(sampleId) %>%
  distinct() #MetNa03 MetNa11 MetNa16 MetNa17

#clean cycles data, make it integer
mat_HMF_treatment$Palliative_oxaliplatin_cycles<- ifelse(mat_HMF_treatment$Palliative_oxaliplatin_cycles == "No", 0, mat_HMF_treatment$Palliative_oxaliplatin_cycles)
columns_to_replace <- c("Palliative_5FU_cycles", "Palliative_oxaliplatin_cycles", "Palliative_irinotecan_cycles")
for (col in columns_to_replace) {  mat_HMF_treatment[[col]][mat_HMF_treatment[[col]] == "NA"] <- NA}
mat_HMF_treatment[columns_to_replace] <- lapply(mat_HMF_treatment[columns_to_replace], as.integer)

mat_HMF_treatment$Palliative_cycles <- mat_HMF_treatment$Palliative_irinotecan_cycles+
  mat_HMF_treatment$Palliative_oxaliplatin_cycles+
  mat_HMF_treatment$Palliative_5FU_cycles #make 1 variable with all cycles

#Variables of interest
WGS_variable <- c("tml", "tmbPerMb", "svTumorMutationalBurden")

# Generate corplots for signatures vs treatment cycles
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "plots_tmb_sig"), width = 16, height = 2*length(WGS_variable)) 
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
             nrow = 4, ncol = 6)
dev.off()

# Generate corplots for signatures vs treatment cycles
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "plots_cycles_sig"), width = 16, height = 2*length(WGS_variable)) 
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
             nrow = length(WGS_variable), ncol = 8)
dev.off()



#------Final figures-----

# Make boxplots for TMB/TML/WG disruptions in Hartwig data (validation set)
boxplot_tml <- ggplot(mat_HMF_treatment, aes(y = tml, x = `Sample type`)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "tml", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 1)+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)
boxplot_tml

ggsave(file=sprintf("%s/%s.pdf",WGS_plot_dir, "boxplot_tml"), width = (3), height = (3)) #save boxplot TLM
dev.off()

boxplot_tmb <- ggplot(mat_HMF_treatment, aes(y = tmbPerMb, x = `Sample type`)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "tmb per Mb", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 0.1)+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)
boxplot_tmb

ggsave(file=sprintf("%s/%s.pdf",WGS_plot_dir, "boxplot_tmb"), width = (3), height = (3)) #save boxplot TMB
dev.off()

boxplot_svtmb <- ggplot(mat_HMF_treatment, aes(y = svTumorMutationalBurden, x = `Sample type`)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "sv tmb", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 1)+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)
boxplot_svtmb

ggsave(file=sprintf("%s/%s.pdf",WGS_plot_dir, "boxplot_svtmb"), width = (3), height = (3)) #save boxplot svTMB
dev.off()

#contingency_table for whole genome duplications in Hartwig data (validation set)
contingency_table <- table(mat_HMF_treatment$wholeGenomeDuplication, mat_HMF_treatment$`Sample type` == "pretreated")
dimnames(contingency_table) <- list(wholeGenomeDuplication = c("FALSE", "TRUE"), chemotherapy_pretreated = c("FALSE", "TRUE"))
fisher_result <- fisher.test(contingency_table) #not significant
contingency_table

boxplot_ploidy <- ggplot(mat_HMF_treatment, aes(y = ploidy, x = `Sample type`)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "ploidy", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 1)+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)
boxplot_ploidy

ggsave(file=sprintf("%s/%s.pdf",WGS_plot_dir, "boxplot_ploidy"), width = (3), height = (3)) #save boxplot ploidy
dev.off()

mat_HMF_treatment$treatment <- mat_HMF_treatment$`Sample type`

generate_boxplot <- function(data, y_var, x_var, y_label, x_label) {
  plot <- ggplot(data, aes_string(y = y_var, x = x_var)) + 
    geom_jitter(width = 0.2, size = 0.5) +
    labs(title = "", y = y_label, x = x_label) + 
    theme_classic() +
    scale_y_continuous(trans = "log10") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)
  
  return(plot)
}

# Generate the individual plots
boxplot_tml <- generate_boxplot(mat_HMF_treatment, "tml", "treatment", "tml", "chemotherapy")
boxplot_tmb <- generate_boxplot(mat_HMF_treatment, "tmbPerMb", "treatment", "tmb per Mb", "chemotherapy")
boxplot_svtmb <- generate_boxplot(mat_HMF_treatment, "svTumorMutationalBurden", "treatment", "sv tmb", "chemotherapy")
boxplot_ploidy <- generate_boxplot(mat_HMF_treatment, "ploidy", "treatment", "ploidy", "chemotherapy")

grid_plot <- grid.arrange(grobs = c(list(boxplot_tml, boxplot_tmb, boxplot_svtmb, boxplot_ploidy)), ncol=4)
ggsave(file = sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_PDO_tmb_etc"), plot = grid_plot, width = 12, height = 3)


