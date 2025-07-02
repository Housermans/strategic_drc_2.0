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
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 1)+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)
boxplot_tml

ggsave(file=sprintf("%s/%s.pdf",WGS_plot_dir, "boxplot_tml"), width = (3), height = (3)) #save boxplot TLM
dev.off()

boxplot_tmb <- ggplot(purple_treatment, aes(y = tmbPerMb, x = chemotherapy)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "tmb per Mb", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 0.1)+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)
boxplot_tmb

ggsave(file=sprintf("%s/%s.pdf",WGS_plot_dir, "boxplot_tmb"), width = (3), height = (3)) #save boxplot TMB
dev.off()

boxplot_svtmb <- ggplot(purple_treatment, aes(y = svTumorMutationalBurden, x = chemotherapy)) + 
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
contingency_table <- table(purple_treatment$wholeGenomeDuplication, purple_treatment$chemotherapy == "pretreated")
dimnames(contingency_table) <- list(wholeGenomeDuplication = c("FALSE", "TRUE"), chemotherapy_pretreated = c("FALSE", "TRUE"))
fisher_result <- fisher.test(contingency_table) #not significant
contingency_table

boxplot_ploidy <- ggplot(purple_treatment, aes(y = ploidy, x = chemotherapy)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "ploidy", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)
boxplot_ploidy

ggsave(file=sprintf("%s/%s.pdf",WGS_plot_dir, "boxplot_ploidy"), width = (3), height = (3)) #save boxplot ploidy
dev.off()

# Prep merged (combined) Hartwig data (validation set)  for table with driver genes
overview <- setNames(data.frame(matrix(ncol = length(unique(combined$gene)), nrow = length(unique(combined$sampleId)))),unique(combined$gene))
row.names(overview) <- unique(combined$sampleId)
overview <- overview %>% tibble::rownames_to_column(var = "sampleId") 

#make table with driver genes


#------Heavy function (do not use!!!)----
# for (geneID in unique(combined$gene)){
#   for(sampleID in unique(combined$sampleId)){
#     
#     tt <- combined  %>% dplyr::filter(sampleId==sampleID) %>% dplyr::filter(gene==geneID) 
#     if(nrow(tt)>1){
#       tt <- tt %>% dplyr::filter(isCanonical=="true")}
#     
#     if(nrow(tt)==0){
#       next
#     }else{
#       
#       simplemuts<- tt[c("missense","nonsense","splice","inframe","frameshift")]
#       print(tt)
#       if(tt$likelihoodMethod=="DEL"){
#         overview[overview$sampleId==sampleID,geneID] <- "deep deletion"}
#       if(tt$likelihoodMethod=="AMP" & rowSums(simplemuts)==0){
#         overview[overview$sampleId==sampleID,geneID] <- "amplification"}
#       if(tt$likelihoodMethod=="AMP" & simplemuts$missense>1){
#         overview[overview$sampleId==sampleID,geneID] <- "amplification/missense"}
#       if(tt$likelihoodMethod=="AMP" & simplemuts$nonsense>1){
#         overview[overview$sampleId==sampleID,geneID] <- "amplification/nonsense"}
#       if(tt$likelihoodMethod=="AMP" & simplemuts$splice>1){
#         overview[overview$sampleId==sampleID,geneID] <- "amplification/splice"}
#       if(tt$likelihoodMethod=="AMP" & simplemuts$inframe>1){
#         overview[overview$sampleId==sampleID,geneID] <- "amplification/inframe"}
#       if(tt$likelihoodMethod=="AMP" & simplemuts$frameshift>1){
#         overview[overview$sampleId==sampleID,geneID] <- "amplification/frameshift"}
#       if(tt$driver=="MUTATION"){
#         simplemuts <- simplemuts %>% select_if(colSums(.) != 0)
#         if(rowSums(simplemuts)==1){
#           overview[overview$sampleId==sampleID,geneID] <- sprintf(names(simplemuts))
#         }else{
#           if(ncol(simplemuts)>1){
#             if(ncol(simplemuts)==2){
#               overview[overview$sampleId==sampleID,geneID] <- paste(names(simplemuts[1]),names(simplemuts[2]),sep = "/")
#             }else{overview[overview$sampleId==sampleID,geneID] <- paste(names(simplemuts[2]),names(simplemuts[3]),sep = "/")}
#           }
#           if(ncol(simplemuts)==1){
#             overview[overview$sampleId==sampleID,geneID] <- paste(names(simplemuts[1]),names(simplemuts[1]),sep = "/")
#           }
#         }
#       }
#     }
#   }
# }


#------import WGS data PDO made in script 9-----
mat_HMF_treatment <- readRDS(file.path(WGS_plot_dir, "mat_HMF_treatment.rds"))
sig_PDO_categorized <- readRDS(file.path(WGS_plot_dir, "sig_PDO_categorized.rds"))
PDO_WGS_screen <- readRDS(file.path(WGS_plot_dir, "PDO_WGS_screen.rds"))
PDO_WGS_screen_top5 <- readRDS(file.path(WGS_plot_dir, "PDO_WGS_screen_top5.rds"))
PDO_WGS_screen_pretreatedfiltered <- readRDS(file.path(WGS_plot_dir, "PDO_WGS_screen_pretreatedfiltered.rds"))
PDO_WGS_screen_pretreatedonly <- readRDS(file.path(WGS_plot_dir, "PDO_WGS_screen_pretreatedonly.rds"))
PDO_metadata_screen <- readRDS(file.path(WGS_plot_dir, "PDO_metadata_screen.rds"))
PDO_WGS_screen_top5_chemonaive <- readRDS(file.path(WGS_plot_dir, "PDO_WGS_screen_top5_chemonaive.rds"))
PDO_WGS_screen_top5_pretreated <- readRDS(file.path(WGS_plot_dir, "PDO_WGS_screen_top5.rds"))
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
    plot <- plot + facet_grid(. ~ condition, scales = "free")
  } else if (facet == "chemo_naive") {
    plot <- plot + facet_grid(. ~ chemo_naive, scales = "free")
  } else if (facet == "chemo_condition") {
    plot <- plot + facet_grid(. ~ chemo_condition, scales = "free")
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
    plot <- plot + stat_cor(method = "pearson", size = 3, color = "red")
  }
  return(plot)
}

# new combined variable for facetting on both condition and chemonaive status
PDO_WGS_screen$chemo_condition <- interaction(PDO_WGS_screen$chemo_naive, PDO_WGS_screen$condition)
PDO_WGS_screen_top5$chemo_condition <- interaction(PDO_WGS_screen_top5$chemo_naive, PDO_WGS_screen_top5$condition)

# Generate boxplots for genes of interest
genes <- c("SMAD4", "MACROD2", "PRKN", "KRAS", "PIK3CA", "TCF7L2", "PTEN", "FBXW7", "CCSER1", "FLCN", "FHIT")
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_gene_AUC"), width = 19*2, height = 2 * length(genes))
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
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_tmbetc"), width = 19*2, height = 2 * (2 * length(WGS_variable) + length(WGS_variable_status))) 
plots <- lapply(WGS_variable, function(var) {
  make_boxplot(PDO_WGS_screen, x_var = "AUC_category", y_var = var,
               comparisons = list(c("resistant", "sensitive")), x_lab = "PDO", y_lab = var)
})
plots_status <- lapply(WGS_variable_status, function(var) {
  make_boxplot(PDO_WGS_screen, x_var = var, y_var = "AUC", 
               comparisons = list(c("HIGH", "LOW")), x_lab = var, y_lab = "AUC", label_y = 0.9)
})
corplots <- lapply(WGS_variable, function(var) {
  make_corplot(PDO_WGS_screen, x_var = "AUC", y_var = var, 
               x_lab = "AUC", y_lab = var, annotate_stats = TRUE, log_y = TRUE)
})
grid.arrange(grobs = c(plots, plots_status, corplots), nrow = 2 * length(WGS_variable) + length(WGS_variable_status))
dev.off()

#Variables of interest
WGS_variable_karyotype <- colnames(PDO_WGS_screen)[startsWith(colnames(PDO_WGS_screen), "chr")]
WGS_variable_karyotype_top5 <- colnames(PDO_WGS_screen_top5)[startsWith(colnames(PDO_WGS_screen_top5), "chr")]

# Generate boxplots and corplots for karyotype
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "corplots_karyotype"), width = 2*19, height = 2*length(WGS_variable_karyotype)) 
corplots_karyotype <- lapply(WGS_variable_karyotype, function(var) {
  make_corplot(PDO_WGS_screen, x_var = "AUC", y_var = var, 
               x_lab = "AUC", y_lab = paste0(var, "(log)"), annotate_stats=TRUE)
})
grid.arrange(grobs = c(corplots_karyotype), nrow = length(WGS_variable_karyotype))
dev.off()

pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "corplots_karyotype_top5"), width = 2*19, height = 2*length(WGS_variable_karyotype_top5)) 
corplots_karyotype <- lapply(WGS_variable_karyotype_top5, function(var) {
  make_corplot(PDO_WGS_screen_top5, x_var = "AUC", y_var = var, 
               x_lab = "AUC", y_lab = paste0(var, "(log)"), annotate_stats=TRUE)
})
grid.arrange(grobs = c(corplots_karyotype), nrow = length(WGS_variable_karyotype))
dev.off()


# Generate boxplots for Palliative oxaliplatin PD
PDO_metadata_screen$Palliative_oxaliplatin_PD <- ifelse(PDO_metadata_screen$Palliative_oxaliplatin_response_all == "PD" & PDO_metadata_screen$chemo_naive == "no", "PD",
                                                   ifelse(PDO_metadata_screen$chemo_naive == "no", "no PD", NA))
PDO_metadata_screen$Palliative_oxaliplatin_PD <- ifelse(is.na(PDO_metadata_screen$Palliative_oxaliplatin_response_all) & PDO_metadata_screen$chemo_naive == "no", "no PD",
                                                   PDO_metadata_screen$Palliative_oxaliplatin_PD)
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_PD_oxaliplatin"), width = 19, height = 2)
plots_PD_oxaliplatin <- lapply(c("Palliative_oxaliplatin_PD"), function(var) {
  make_boxplot(PDO_metadata_screen, x_var = var, y_var = "AUC", 
               comparisons = list(c("PD", "no PD")), x_lab = var, y_lab = "AUC", label_y = 0.9)
})
grid.arrange(grobs = plots_PD_oxaliplatin, nrow = 1)
dev.off()

# Generate boxplots for signature variables
signature_variables <- c("SBS17a_cat", "SBS17b_cat", "SBS35_cat", "SBS17_threshold", "has_deep_deletion", "SBS2_cat", "SBS13_cat")
signature_variables_cont <- c("SBS17a", "SBS17b", "SBS35", "SBS17", "SBS2", "SBS13")
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "boxplots_signatures"), width = 2*19, height = 2 * length(signature_variables)+length(signature_variables_cont))
plots_signature <- lapply(signature_variables, function(var) {
  make_boxplot(PDO_WGS_screen, x_var = var, y_var = "AUC", 
               comparisons = list(c("present", "absent")), x_lab = var, y_lab = "AUC", label_y = 0.9)
})
plots_signature_cont <- lapply(signature_variables_cont, function(var) {
  make_boxplot(PDO_WGS_screen, x_var = "AUC_category", y_var = var, 
               comparisons = list(c("resistant", "sensitive")), x_lab = "PDO", y_lab = var)
})
grid.arrange(grobs = c(plots_signature, plots_signature_cont), nrow = length(signature_variables)+length(signature_variables_cont))
dev.off()

PDO_WGS_screen_pretreatedonly$org_name[PDO_WGS_screen_pretreatedonly$DBS5 == 0] #RAS05 heeft geen DBS5, niet goed gecalld?

# chemonaive and pretreated together 
plots_signature <- lapply(signature_variables, function(var) {
make_boxplot(PDO_WGS_screen, x_var = var, y_var = "AUC", comparisons = list(c("present", "absent")), 
             x_lab = var, y_lab = "AUC", label_y = 0.9, facet="condition")})
 
# Check deep deletions in chemonaive group
PDO_WGS_screen %>%
  filter(has_deep_deletion == "present" & chemo_naive == "yes") %>%
  select(sampleId) %>%
  distinct() #MetNa03 MetNa11 MetNa16 MetNa17

PDO_WGS_screen_chemonaivefiltered<- PDO_WGS_screen %>%
  filter(chemo_naive == "yes")

# Generate boxplots for adjuvant treatment
make_boxplot(PDO_WGS_screen, x_var = "Adjuvant_received", y_var = "AUC", 
             comparisons = list(c("Yes", "No")), x_lab = "Adjuvant", y_lab = "AUC", label_y = 0.9)

#clean cycles data, make it integer
mat_HMF_treatment$Palliative_oxaliplatin_cycles<- ifelse(mat_HMF_treatment$Palliative_oxaliplatin_cycles == "No", 0, mat_HMF_treatment$Palliative_oxaliplatin_cycles)
columns_to_replace <- c("Palliative_5FU_cycles", "Palliative_oxaliplatin_cycles", "Palliative_irinotecan_cycles")
for (col in columns_to_replace) {  mat_HMF_treatment[[col]][mat_HMF_treatment[[col]] == "NA"] <- NA}
mat_HMF_treatment[columns_to_replace] <- lapply(mat_HMF_treatment[columns_to_replace], as.integer)

mat_HMF_treatment$Palliative_cycles <- mat_HMF_treatment$Palliative_irinotecan_cycles+
  mat_HMF_treatment$Palliative_oxaliplatin_cycles+
  mat_HMF_treatment$Palliative_5FU_cycles #make 1 variable with all cycles

#Variables of interest
WGS_variable <- colnames(mat_HMF_treatment)[endsWith(colnames(mat_HMF_treatment), "cycles")]

# Generate corplots for signatures vs treatment cycles
pdf(sprintf("%s/%s.pdf", WGS_plot_dir, "plots_cycles_sig"), width = 20, height = 2*length(WGS_variable)) 
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

# Check signatures in chemonaive group
PDO_WGS_screen %>%
  filter(SBS17a_cat == "present" & chemo_naive == "yes") %>%
  select(sampleId) %>%
  distinct() #MetNa02

PDO_WGS_screen %>%
  filter(SBS17b_cat == "present" & chemo_naive == "yes") %>%
  select(sampleId) %>%
  distinct() #MetNa03 MetNa11 MetNa16 MetNa17

boxplots_SBS17b <- lapply(WGS_variable, function(var) {
  make_boxplot(mat_HMF_treatment, x_var = "SBS17b_cat", y_var = var,
               comparisons = list(c("absent", "present")), x_lab = "signature", y_lab = var, facet="no")})


# vanaf hier geen toelichting, deze plots gaan we wss niet gebruiken
# Making df for PDO drivers----------------------
plotting_data <- mat_PDO
row.names(plotting_data) <- plotting_data$sampleId
plotting_data$sampleId <- NULL

modify_gene_names <- function(gene_names = names(plotting_data),df=plotting_data){
  
  percent <- function(x, digits = 0, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  }
  
  data <- df
  n_samples<-nrow(data)
  data[!is.na(data)] <- 1
  data[is.na(data)] <- 0
  data <- mutate_all(data, function(x) as.numeric(as.character(x)))
  percentage <- percent(colSums(as.matrix(data), na.rm = FALSE)/n_samples)
  
  df2 <- data.frame(Gene=names(df),
                    percentage_ID=percentage)
  df2$percentage_ID_2 <- with(df2, paste("(",percentage_ID,")", sep=""))
  cols<-c("Gene","percentage_ID_2")
  df2$unite <- apply( df2[ , cols ] , 1 , paste , collapse = " " )
  return(df2$unite)
}

gene_names<-modify_gene_names(plotting_data)

# Make a new variable that indicates if patiens are pretreated based on sample Ids
plotting_data <- tibble::rownames_to_column(plotting_data,var = "sampleId")
row.names(plotting_data) <- plotting_data$sampleId
plotting_data$chemotherapy <- ifelse(grepl("Pret", plotting_data$sampleId), "pretreated", 
                                     ifelse(grepl("Na", plotting_data$sampleId), "chemonaive", NA))
plotting_data$sampleId <- NULL #remove Sample Id variable, so data can be reshaped
plotting_data$chemotherapy <- NULL #remove chemotherapy variable, so data can be reshaped

# Reshape data to make a dataframe suitable for a geom_tile plot with driver categories per allele
make_dt <- function(plotting_data){
  m <- plotting_data
  names(m)<-seq(1, ncol(m), by=1)
  m <- as.matrix(m)
  
  dt <- data.table(reshape2::melt(m))
  dt <- dt[, strsplit(as.character(value), "/"), by=list(Var1, Var2)]  # this expands "X/Y/Z" into three rows
  dt[, shift:=(1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(Var1, Var2)]
  dt[, height:=1/.N, by=list(Var1, Var2)]
  dt$test=dt$Var2+dt$shift
  dt$test2=as.numeric(format(dt$test,digits =1))
  dt$height2=as.numeric(1)
  return(dt)
}

dt <- make_dt(plotting_data)

colorpalette_drivers = c("amplification" =  '#b3de69',
                         "deep deletion" =  'coral1',
                         "INDEL" =  'skyblue3',
                         "Synonymous" =  '#fb8072',
                         "NA" =  'white',
                         "missense" =  '#fdb462',
                         "nonsense" =  '#ffffb3',
                         "splice" =  '#8dd3c7',
                         "inframe" =  '#d9d9d9',
                         "frameshift" = 'skyblue',
                         "Nonsense" = 'plum3' ,
                         "grey" = 'gray80')

plot <- ggplot() + 
  geom_tile(data = dt,aes(Var1,y=test, fill=V1, height=height),
            color="grey", size=0.1) + xlab('sample ID') + ylab('Gene')+
  scale_y_continuous(breaks = seq(1, ncol(plotting_data), by=1), labels = gene_names)+
  scale_fill_manual(values=colorpalette_drivers,na.value = "white") +
  guides(fill=guide_legend(title="Effect"))+
  theme(axis.title.y=element_text(size=5,vjust=3),
        axis.text.y=element_text(size=5),
        axis.title.x=element_text(size=5),
        axis.text.x = element_text(colour="grey20",size=5,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.ticks.y = element_line(size=0.5),
        strip.text.x=element_text(size=5),
        strip.text.y=element_text(size=5,vjust=2),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = 'grey20',linewidth = 0.5),
        strip.background=element_blank(),
        legend.text=element_text(size=5),
        legend.title=element_text(size=5))

plot

geom_tile(data = dt,aes(Var1,y=test2, height=height2),
          color="gray46", size=0.1,alpha=0)

# Convert Var2 to factor with gene names
longformat <- function(dt, fiveplus){
  dt$Var2 <- factor(dt$Var2, levels = unique(dt$Var2), labels = gene_names)
  dt$percentage <- as.numeric(gsub(".*\\((\\d+)%\\).*", "\\1", dt$Var2))
  dt$Gene_Name <- sub("\\s.*", "", dt$Var2)
  if (fiveplus == TRUE) {
    dt <- dt[dt$Gene_Name %in% PDO_fiveplus, ]
  } else {
    dt <- dt[dt$percentage >= 5, ]
  }
  
  dt_allele <- dt %>%
    group_by(Var1, Var2) %>%
    mutate(
      allele = if(n() == 2) c(1, 2) else 1
    ) %>%
    ungroup()
  
  # Filter rows with only one V1 per Var1 and Var2
  single_v1 <- dt_allele %>%
    group_by(Var1, Var2) %>%
    filter(n() == 1) %>%
    ungroup()
  
  # Duplicate these rows and assign allele 2
  additional_rows <- single_v1 %>%
    mutate(allele = 2) 
  
  # Combine the original data with the additional rows
  dt_allele2 <- bind_rows(dt_allele, additional_rows)
  
  # Arrange the data to calculate row numbers for allele 1
  ordered_rows_allele1 <- dt_allele2 %>%
    filter(allele == 1) %>%
    arrange(Var2, V1)
  
  # Calculate row numbers for allele 1
  ordered_rows_allele1 <- ordered_rows_allele1 %>%
    group_by(Var2, allele) %>%
    mutate(row_number_allele1 = row_number() - min(row_number()) + 1) %>%
    ungroup()
  
  dt_rownumbers <- dt_allele2 %>%
    left_join(ordered_rows_allele1 %>% select(Var1, Var2, row_number_allele1), by = c("Var1", "Var2")) %>%
    mutate(row_number = row_number_allele1)
  
  dt_rownumbers$facet_variable <- paste(dt_rownumbers$Var2, dt_rownumbers$allele, sep = " / ")
  return(dt_rownumbers)
}

dt_rownumbers <- longformat(dt, fiveplus = FALSE) # Make a long file suitable for bargraph
PDO_fiveplus <- unique(dt_rownumbers$Gene_Name) # Select driver genes with 5+% frequency

dt_rownumbers <- dt_rownumbers[dt_rownumbers$Gene_Name != "HLA.C" ,] # Remove gene not present in WGS data
dt_rownumbers$facet_variable <- as.factor(dt_rownumbers$facet_variable) # Make factor for bargraph

# Relevel in order of frequency
levels_new <- c("TP53 (93%) / 1", "TP53 (93%) / 2",
                "APC (85%) / 1", "APC (85%) / 2",
                "KRAS (63%) / 1", "KRAS (63%) / 2",
                "PIK3CA (29%) / 1", "PIK3CA (29%) / 2",
                "SMAD4 (24%) / 1", "SMAD4 (24%) / 2",
                "PRKN (17%) / 1", "PRKN (17%) / 2",
                "TCF7L2 (17%) / 1", "TCF7L2 (17%) / 2",
                "MACROD2 (12%) / 1", "MACROD2 (12%) / 2",
                "PTEN (12%) / 1", "PTEN (12%) / 2",
                "BRAF (12%) / 1", "BRAF (12%) / 2",
                "FLCN (10%) / 1", "FLCN (10%) / 2",
                "FHIT (10%) / 1", "FHIT (10%) / 2",
                "CCSER1 (10%) / 1", "CCSER1 (10%) / 2",
                "CDKN2A (7%) / 1", "CDKN2A (7%) / 2",
                "NAALADL2 (7%) / 1", "NAALADL2 (7%) / 2",
                "FBXW7 (7%) / 1", "FBXW7 (7%) / 2",
                "ZFP36L2 (7%) / 1", "ZFP36L2 (7%) / 2",
                "CCND2 (5%) / 1", "CCND2 (5%) / 2",
                "CDK12 (5%) / 1", "CDK12 (5%) / 2",
                "SMAD3 (5%) / 1", "SMAD3 (5%) / 2",
                "SEMG2 (5%) / 1", "SEMG2 (5%) / 2",
                "ELF3 (5%) / 1", "ELF3 (5%) / 2",
                "GATA3 (5%) / 1", "GATA3 (5%) / 2",
                "STS (5%) / 1", "STS (5%) / 2",
                "MTOR (5%) / 1", "MTOR (5%) / 2",
                "RNF43 (5%) / 1", "RNF43 (5%) / 2",
                "KMT2C (5%) / 1", "KMT2C (5%) / 2",
                "STAG2 (5%) / 1", "STAG2 (5%) / 2",
                "CDX2 (5%) / 1", "CDX2 (5%) / 2", 
                "FLT3 (5%) / 1", "FLT3 (5%) / 2",
                "GNAS (5%) / 1", "GNAS (5%) / 2",
                "MAP3K21 (5%) / 1", "MAP3K21 (5%) / 2")

dt_rownumbers$facet_variable <- factor(dt_rownumbers$facet_variable, levels = levels_new)

barplot_genes <- ggplot(dt_rownumbers, aes(x = factor(row_number), fill = factor(V1))) +
  geom_bar(stat = "count", position = "stack", color = "transparent") + 
  facet_grid(facet_variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Frequence", y = "Count", title = "Driver mutations PDO") +
  theme_minimal() +
  scale_fill_manual(values=colorpalette_drivers,na.value = "white") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),  # Remove text on the y-axis
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        strip.text.y = element_text(size = 5, angle = 0))  # Adjust strip text properties

barplot_genes <- barplot_genes + 
  theme(panel.spacing.y = unit(-0.2, "lines"))

print(barplot_genes)

ggsave(file=sprintf("%s/%s.pdf", WGS_plot_dir, "barplot_genes_PDO"), width = (5), height = (5)) 
dev.off()


#making df for HMF drivers----------------------
#write.table(mat, file.path(WGS_plot_dir, "Driveroverview.txt"), sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE) #safe table
mat <- read.csv(file.path(WGS_plot_dir, "Driveroverview.txt"), sep = "\t")
mat_HMF_treatment <- merge(mat, HMF_treatment_filter, by = "sampleId")

#create table with data
mat_HMF_treatment$MACROD2_no <- ifelse(is.na(mat_HMF_treatment$MACROD2), "No deletion", "deep deletion")
mat_HMF_treatment$PRKN_no <- ifelse(is.na(mat_HMF_treatment$PRKN), "No deletion", "deep deletion")
counts <- table(mat_HMF_treatment$MACROD2_no, mat_HMF_treatment$chemotherapy)
counts_chemo <- table(mat_HMF_treatment$chemotherapy)
count_perc <- apply(counts, 2, function(x) {round(x * 100 / sum(x, na.rm = TRUE), 1)})

fisher.test(counts)
fisher_result

pdf(sprintf("%s/%s.pdf",WGS_plot_dir, "MACROD2_bargraph") , useDingbats = F, width = 3.5, height = 3)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(count_perc, col=brewer.pal(3, "Pastel2") , border="white")
text(x = c(0.7, 1.9), y=c(100, 100), paste(c(count_perc[1], count_perc[3]), "%"), cex=1, pos = 3)
par(mar=c(4.1, 2.1, 2.1, 2.1), xpd=FALSE)

legend("topright", legend = c("Deep\ndeletion", "No\ndeletion") , 
       col = c("#B3E2CD", "#FDCDAC", "#CBD5E8") , "#B3E2CD", "#FDCDAC", "#CBD5E8" ,
       bty = "n", pch=20 , pt.cex = 3, cex = 1, horiz = FALSE, inset = c(-1, 0.0))

dev.off() # Close PDF device

row.names(mat) <- mat$sampleId #sample ID as row name
mat$sampleId <- NULL 

select_driver_genes <- names(mat)[nrow(mat)- colSums(is.na(mat))>0] #filter driver genes with low frequency in data
mat<-dplyr::select(mat,select_driver_genes)
driver_names_order <-names(sort(colSums(!is.na(mat)),decreasing = T))
mat<-dplyr::select(mat,rev(driver_names_order))

plotting_data <- tibble::rownames_to_column(mat)
names(plotting_data)[1]<-("sampleId")

row.names(plotting_data) <- plotting_data$sampleId
plotting_data$sampleId <- NULL

gene_names<-modify_gene_names(plotting_data)

plotting_data <- tibble::rownames_to_column(plotting_data,var = "sampleId")
plotting_data <- merge(plotting_data, HMF_treatment_filter, by = "sampleId")
row.names(plotting_data) <- plotting_data$sampleId
plotting_data$sampleId <- NULL

plotting_data <- plotting_data[order(plotting_data$chemotherapy), ]
chemotherapy <- plotting_data$chemotherapy
chemotherapy <- as.data.frame(chemotherapy)
row_names <- rownames(plotting_data)
rownames(chemotherapy) <- row_names
plotting_data$chemotherapy <- NULL
plotting_data$EGFR_inhibition <- NULL

levels <- c("missense", "frameshift", "missense/nonsense", "splice", 
            "deep deletion", "missense/inframe", "nonsense", "inframe", "amplification", 
            "missense/missense", "nonsense/inframe", "missense/frameshift", "nonsense/frameshift", "inframe/frameshift",
            "frameshift/frameshift", "nonsense/splice", "nonsense/nonsense", "missense/splice", "splice/frameshift",
            "splice/splice")

factor <- lapply(plotting_data, function(x) {
  factor(x, levels = levels)
})

numeric <- as.data.frame(lapply(factor, as.numeric))
rownames(numeric) <- row_names

colors <- colorRampPalette(brewer.pal(8, "Spectral"))(12)

pheatmap <- pheatmap(numeric,
                     legend_labels = levels,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     main = "Frequent alterations in HMF CRC data",
                     angle_col = 315,
                     row_names_rot = 90,
                     annotation_row = chemotherapy,
                     fontsize_row = 5,
                     fontsize_col = 8,
                     legend_breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20),
                     color = colorRampPalette(brewer.pal(8, "Spectral"))(20))

ggsave(file=sprintf("%s/%s.pdf",WGS_plot_dir, "heatmap_drivers"), pheatmap, width = (8), height = (13)) 
dev.off()

dt <- make_dt(plotting_data)

plot <- ggplot() + 
  geom_tile(data = dt,aes(Var1,y=test, fill=V1, height=height),
            color="grey", size=0.1) + xlab('sample ID') + ylab('Gene')+
  scale_y_continuous(breaks = seq(1, ncol(plotting_data), by=1), labels = gene_names)+
  scale_fill_manual(values=colorpalette_drivers,na.value = "white") +
  guides(fill=guide_legend(title="Effect"))+
  theme(axis.title.y=element_text(size=5,vjust=3),
        axis.text.y=element_text(size=5),
        axis.title.x=element_text(size=5),
        axis.text.x = element_text(colour="grey20",size=5,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.ticks.y = element_line(size=0.5),
        strip.text.x=element_text(size=5),
        strip.text.y=element_text(size=5,vjust=2),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = 'grey20',linewidth = 0.5),
        strip.background=element_blank(),
        legend.text=element_text(size=5),
        legend.title=element_text(size=5))

geom_tile(data = dt,aes(Var1,y=test2, height=height2),
          color="gray46", size=0.1,alpha=0)

dt_rownumbers <- longformat(dt, fiveplus = TRUE)
dt_rownumbers$facet_variable <- as.factor(dt_rownumbers$facet_variable)

levels_new_HMF <- c("TP53 (84%) / 1", "TP53 (84%) / 2",
                    "APC (86%) / 1", "APC (86%) / 2",
                    "KRAS (48%) / 1", "KRAS (48%) / 2",
                    "PIK3CA (17%) / 1", "PIK3CA (17%) / 2",
                    "SMAD4 (22%) / 1", "SMAD4 (22%) / 2",
                    "PRKN (8%) / 1", "PRKN (8%) / 2",
                    "TCF7L2 (16%) / 1", "TCF7L2 (16%) / 2",
                    "MACROD2 (6%) / 1", "MACROD2 (6%) / 2",
                    "PTEN (4%) / 1", "PTEN (4%) / 2",
                    "BRAF (6%) / 1", "BRAF (6%) / 2",
                    "FLCN (1%) / 1", "FLCN (1%) / 2",
                    "FHIT (6%) / 1", "FHIT (6%) / 2",
                    "CCSER1 (7%) / 1", "CCSER1 (7%) / 2",
                    "CDKN2A (0%) / 1", "CDKN2A (0%) / 2",
                    "NAALADL2 (1%) / 1", "NAALADL2 (1%) / 2",
                    "FBXW7 (8%) / 1", "FBXW7 (8%) / 2",
                    "ZFP36L2 (12%) / 1", "ZFP36L2 (12%) / 2",
                    "CCND2 (3%) / 1", "CCND2 (3%) / 2",
                    "CDK12 (1%) / 1", "CDK12 (1%) / 2",
                    "SMAD3 (5%) / 1", "SMAD3 (5%) / 2",
                    "SEMG2 (3%) / 1", "SEMG2 (3%) / 2",
                    "ELF3 (1%) / 1", "ELF3 (1%) / 2",
                    "GATA3 (1%) / 1", "GATA3 (1%) / 2",
                    "STS (1%) / 1", "STS (1%) / 2",
                    "MTOR (2%) / 1", "MTOR (2%) / 2",
                    "RNF43 (2%) / 1", "RNF43 (2%) / 2",
                    "KMT2C (6%) / 1", "KMT2C (6%) / 2",
                    "STAG2 (2%) / 1", "STAG2 (2%) / 2",
                    "CDX2 (4%) / 1", "CDX2 (4%) / 2", 
                    "FLT3 (5%) / 1", "FLT3 (5%) / 2",
                    "GNAS (5%) / 1", "GNAS (5%) / 2",
                    "MAP3K21 (5%) / 1", "MAP3K21 (5%) / 2")

dt_rownumbers$facet_variable <- factor(dt_rownumbers$facet_variable, levels = levels_new_HMF)

barplot_genes <- ggplot(dt_rownumbers, aes(x = factor(row_number), fill = factor(V1))) +
  geom_bar(stat = "count", position = "stack", color = "transparent") + 
  facet_grid(facet_variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Frequence", y = "Count", title = "Driver mutations validation dataset") +
  theme_minimal() +
  scale_fill_manual(values=colorpalette_drivers,na.value = "white") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),  # Remove text on the y-axis
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        strip.text.y = element_text(size = 5, angle = 0))  # Adjust strip text properties

barplot_genes <- barplot_genes + 
  theme(panel.spacing.y = unit(-0.2, "lines"))

print(barplot_genes)

ggsave(file=sprintf("%s/%s.pdf", WGS_plot_dir, "barplot_genes"), width = (5), height = (5)) 
dev.off()






