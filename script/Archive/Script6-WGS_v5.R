# Load dplyr and readxl packages for data manipulation and reading Excel files
library(dplyr)
library(readxl)
library(data.table)
library(ggplot2)
library(readxl)
library(RColorBrewer)
library(pheatmap)
library(ggpubr) #stat_compare_means

# Clear existing objects in the workspace
rm(list=ls())

# Define necessary directories for accessing and saving files
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
home_dir <- dirname(dirname(script_dir))
WGS_dir <- file.path(home_dir, "WGS")
WGS_data_dir <- file.path(WGS_dir, "Data_HMF")
WGS_plot_dir <- file.path(home_dir, "Analyse 2.0/7_WGS")

# Read the tsv
HMF_treatment <- read.csv(file.path(WGS_data_dir,"Hartwig_treatments_DR279_clean.tsv"), sep = "\t")   
drivers <- read.csv(file.path(WGS_data_dir,"linx_driver_gene_data_prefilter.tsv"), sep = "\t")   
purple <- read.csv(file.path(WGS_data_dir,"purple_purity_data_prefilter.tsv"), sep = "\t")   

HMF_treatment_filter <- HMF_treatment %>% filter(hadTreatment != "unknown")
HMF_treatment_filter <- HMF_treatment_filter %>% filter(!(standardOfCare == "no" & hadTreatment == "yes"))
HMF_treatment_filter$chemotherapy <- ifelse(is.na(HMF_treatment_filter$standardOfCare), "chemonaive", "pretreated")
HMF_treatment_filter$hadTreatment <- NULL
HMF_treatment_filter$standardOfCare <- NULL
  
HMF_treatment_filter <- HMF_treatment_filter %>% rename(sampleId = patientId) #change patientId to sampleId for merge
purple_treatment <- merge(HMF_treatment_filter, purple, by = "sampleId") #merge dataframes
purple_treatment <- purple_treatment %>% filter(msStatus == "MSS") 
combined <- merge(purple_treatment, drivers, by = "sampleId") #merge dataframes

result <- anti_join(HMF_treatment_filter, purple, by = "sampleId")

boxplot_tml <- ggplot(purple_treatment, aes(y = tml, x = chemotherapy)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "tml", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 1)+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)
boxplot_tml

ggsave(file=sprintf("%s/%s.pdf",WGS_plot_dir, "boxplot_tml"), width = (3), height = (3))
dev.off()

boxplot_tmb <- ggplot(purple_treatment, aes(y = tmbPerMb, x = chemotherapy)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "tmb per Mb", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 0.1)+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)
boxplot_tmb

ggsave(file=sprintf("%s/%s.pdf",WGS_plot_dir, "boxplot_tmb"), width = (3), height = (3)) 
dev.off()

boxplot_svtmb <- ggplot(purple_treatment, aes(y = svTumorMutationalBurden, x = chemotherapy)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "sv tmb", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 1)+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)
boxplot_svtmb

ggsave(file=sprintf("%s/%s.pdf",WGS_plot_dir, "boxplot_svtmb"), width = (3), height = (3)) 
dev.off()

contingency_table <- table(purple_treatment$wholeGenomeDuplication, purple_treatment$chemotherapy == "pretreated")
dimnames(contingency_table) <- list(wholeGenomeDuplication = c("FALSE", "TRUE"), chemotherapy_pretreated = c("FALSE", "TRUE"))
fisher_result <- fisher.test(contingency_table)
contingency_table

boxplot_ploidy <- ggplot(purple_treatment, aes(y = ploidy, x = chemotherapy)) + 
  geom_jitter(width = 0.2, size=0.5) +
  labs(title = "", y = "ploidy", x = "chemotherapy") + 
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  stat_summary(fun.data = mean_cl_normal, geom = "crossbar", width = 0.5)
boxplot_ploidy

ggsave(file=sprintf("%s/%s.pdf",WGS_plot_dir, "boxplot_ploidy"), width = (3), height = (3)) 
dev.off()

overview <- setNames(data.frame(matrix(ncol = length(unique(combined$gene)), nrow = length(unique(combined$sampleId)))),unique(combined$gene))
row.names(overview) <- unique(combined$sampleId)
overview <- overview %>% tibble::rownames_to_column(var = "sampleId") 

#make table with driver genes
#for (geneID in unique(combined$gene)){
  for(sampleID in unique(combined$sampleId)){
    
    tt <- combined  %>% dplyr::filter(sampleId==sampleID) %>% dplyr::filter(gene==geneID) 
    if(nrow(tt)>1){
      tt <- tt %>% dplyr::filter(isCanonical=="true")}
    
    if(nrow(tt)==0){
      next
    }else{
      
      simplemuts<- tt[c("missense","nonsense","splice","inframe","frameshift")]
      print(tt)
      if(tt$likelihoodMethod=="DEL"){
        overview[overview$sampleId==sampleID,geneID] <- "deep deletion"}
      if(tt$likelihoodMethod=="AMP" & rowSums(simplemuts)==0){
        overview[overview$sampleId==sampleID,geneID] <- "amplification"}
      if(tt$likelihoodMethod=="AMP" & simplemuts$missense>1){
        overview[overview$sampleId==sampleID,geneID] <- "amplification/missense"}
      if(tt$likelihoodMethod=="AMP" & simplemuts$nonsense>1){
        overview[overview$sampleId==sampleID,geneID] <- "amplification/nonsense"}
      if(tt$likelihoodMethod=="AMP" & simplemuts$splice>1){
        overview[overview$sampleId==sampleID,geneID] <- "amplification/splice"}
      if(tt$likelihoodMethod=="AMP" & simplemuts$inframe>1){
        overview[overview$sampleId==sampleID,geneID] <- "amplification/inframe"}
      if(tt$likelihoodMethod=="AMP" & simplemuts$frameshift>1){
        overview[overview$sampleId==sampleID,geneID] <- "amplification/frameshift"}
      if(tt$driver=="MUTATION"){
        simplemuts <- simplemuts %>% select_if(colSums(.) != 0)
        if(rowSums(simplemuts)==1){
          overview[overview$sampleId==sampleID,geneID] <- sprintf(names(simplemuts))
        }else{
          if(ncol(simplemuts)>1){
            if(ncol(simplemuts)==2){
              overview[overview$sampleId==sampleID,geneID] <- paste(names(simplemuts[1]),names(simplemuts[2]),sep = "/")
            }else{overview[overview$sampleId==sampleID,geneID] <- paste(names(simplemuts[2]),names(simplemuts[3]),sep = "/")}
          }
          if(ncol(simplemuts)==1){
            overview[overview$sampleId==sampleID,geneID] <- paste(names(simplemuts[1]),names(simplemuts[1]),sep = "/")
          }
        }
      }
    }
  }
}

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

# Close PDF device
dev.off()

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

geom_tile(data = dt,aes(Var1,y=test2, height=height2),
          color="gray46", size=0.1,alpha=0)

# Convert Var2 to factor with gene names
dt$Var2 <- factor(dt$Var2, levels = unique(dt$Var2), labels = gene_names)
dt$percentage <- as.numeric(gsub(".*\\((\\d+)%\\).*", "\\1", dt$Var2))
dt$Gene_Name <- sub("\\s.*", "", dt$Var2)
dt <- dt[dt$Gene_Name %in% PDO_genes_fiveplus, ]

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
unique(dt_rownumbers$facet_variable )
unique(PDO_genes_fiveplus)

levels(dt_rownumbers$facet_variable) <- c("APC (86%) / 1", "APC (86%) / 2", "TP53 (84%) / 1", "TP53 (84%) / 2", "KRAS (48%) / 1", "KRAS (48%) / 2", "SMAD4 (22%) / 1", "SMAD4 (22%) / 2", "FAT4 (17%) / 1", "FAT4 (17%) / 2", "PIK3CA (17%) / 1", "PIK3CA (17%) / 2", "TCF7L2 (16%) / 1", "TCF7L2 (16%) / 2", "ZFP36L2 (12%) / 1", "ZFP36L2 (12%) / 2", "FBXW7 (8%) / 1", "FBXW7 (8%) / 2", "SOX9 (8%) / 1", "SOX9 (8%) / 2", "PRKN (8%) / 1", "PRKN (8%) / 2", "AMER1 (7%) / 1", "AMER1 (7%) / 2", "CCSER1 (7%) / 1", "CCSER1 (7%) / 2", "ZFHX3 (7%) / 1", "ZFHX3 (7%) / 2", "ATM (9%) / 1", "ATM (9%) / 2", "BRAF (6%) / 1", "BRAF (6%) / 2", "ARID1A (6%) / 1", "ARID1A (6%) / 2", "ERBB4 (6%) / 1", "ERBB4 (6%) / 2", "FHIT (6%) / 1", "FHIT (6%) / 2", "GNAS (5%) / 1", "GNAS (5%) / 2", "GRIN2A (6%) / 1", "GRIN2A (6%) / 2", "KMT2C (6%) / 1", "KMT2C (6%) / 2", "MACROD2 (6%) / 1", "MACROD2 (6%) / 2", "MAP3K21 (5%) / 1", "MAP3K21 (5%) / 2", "NRAS (6%) / 1", "NRAS (6%) / 2", "PTPN13 (6%) / 1", "PTPN13 (6%) / 2", "SMAD3 (5%) / 1", "SMAD3 (5%) / 2", "BRCA2 (5%) / 1", "BRCA2 (5%) / 2", "TG (6%) / 1", "TG (6%) / 2")
levels(dt_rownumbers$facet_variable) <-  c("APC (86%) / 1", "APC (86%) / 2", "TP53 (84%) / 1", "TP53 (84%) / 2", "KRAS (48%) / 1", "KRAS (48%) / 2", "SMAD4 (22%) / 1", "SMAD4 (22%) / 2", "TCF7L2 (16%) / 1", "TCF7L2 (16%) / 2", "PIK3CA (17%) / 1", "PIK3CA (17%) / 2", "ZFP36L2 (12%) / 1", "ZFP36L2 (12%) / 2", "FBXW7 (8%) / 1", "FBXW7 (8%) / 2", "PRKN (8%) / 1", "PRKN (8%) / 2", "CCSER1 (7%) / 1", "CCSER1 (7%) / 2", "FHIT (6%) / 1", "FHIT (6%) / 2", "MACROD2 (6%) / 1", "MACROD2 (6%) / 2", "KMT2C (6%) / 1", "KMT2C (6%) / 2", "BRAF (6%) / 1", "BRAF (6%) / 2", "GNAS (5%) / 1", "GNAS (5%) / 2", "SMAD3 (5%) / 1", "SMAD3 (5%) / 2", "FLT3 (5%) / 1", "FLT3 (5%) / 2", "MAP3K21 (5%) / 1", "MAP3K21 (5%) / 2")

dt_rownumbers$facet_variable <- factor(dt_rownumbers$facet_variable, levels = levels(dt_rownumbers$facet_variable))

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

#import WGS data PDOs and metadata
Arne_driver_dir <- paste0(WGS_dir, "/analysis_files_arne/Drivers")
resources_dir <- paste0(home_dir, "/Analyse 2.0/resources")
metrics_dir <- paste0(home_dir, "/Analyse 2.0/6_metrics/240104_QC")

mat_PDO <- read.csv(file.path(Arne_driver_dir, "Driveroverview_rastric.txt"), sep = "\t")
metrics <- read_xlsx(file.path(metrics_dir, "240104_QC_metrics_normalized.xlsx"))
metadata <- read_xlsx(file.path(resources_dir, "STRATEGIC_PERSCO_Overview_Hartwig_WGS.xlsx"))
metadata <- metadata %>% rename(sampleId = ProjectID)

#rename sample IDs for merging
metadata$sampleId <- gsub("MetNaive", "MetNa", metadata$sampleId)
metadata$sampleId <- gsub("MetPretreated", "MetPret", metadata$sampleId)
mat_PDO$sampleId <- gsub("T$", "", mat_PDO$sampleId)

mat_HMF_treatment <- merge(mat_PDO, metadata, by = "sampleId") #merge meta data and WGS PDO data
mat_HMF_treatment$sampleId
View(metrics)

unique(metrics$org_name)

#select 1 from double screens, adjust later!
condition <- !(grepl("^OPT0024_WNT_high$|^OPT0024_WNT_low$|^OPT0042_2$|^RAS36|^RAS37|^RAS38|^RAS39", metrics$org_name))
filtered_metrics <- metrics[condition, ]
filtered_metrics$org_name <- gsub("_Demi|_2", "", filtered_metrics$org_name)

# Check the filtered data
unique(filtered_metrics$org_name)






plotting_data <- mat_PDO
row.names(plotting_data) <- plotting_data$sampleId
plotting_data$sampleId <- NULL

gene_names<-modify_gene_names(plotting_data)

plotting_data <- tibble::rownames_to_column(plotting_data,var = "sampleId")
row.names(plotting_data) <- plotting_data$sampleId
plotting_data$chemotherapy <- ifelse(grepl("Pret", plotting_data$sampleId), "pretreated", 
                                     ifelse(grepl("Na", plotting_data$sampleId), "chemonaive", NA))
plotting_data$sampleId <- NULL
plotting_data$chemotherapy <- NULL

factor <- lapply(plotting_data, function(x) {
  factor(x, levels = levels)
})

numeric <- as.data.frame(lapply(factor, as.numeric))
rownames(numeric) <- row_names

colors <- colorRampPalette(brewer.pal(8, "Spectral"))(12)

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

# Convert Var2 to factor with gene names
dt$Var2 <- factor(dt$Var2, levels = unique(dt$Var2), labels = gene_names)
dt$percentage <- as.numeric(gsub(".*\\((\\d+)%\\).*", "\\1", dt$Var2))
dt <- dt[dt$percentage >= 5, ]
dt$gene_name  <- sub("\\s.*", "", dt$Var2)
PDO_genes_fiveplus<- unique(dt$gene_name)

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
levels(dt_rownumbers$facet_variable) <-  c("TP53 (93%) / 1", "TP53 (93%) / 2", "APC (85%) / 1", "APC (85%) / 2", "KRAS (63%) / 1", "KRAS (63%) / 2", "PIK3CA (29%) / 1", "PIK3CA (29%) / 2", "PRKN (17%) / 1", "PRKN (17%) / 2", "FHIT (10%) / 1", "FHIT (10%) / 2", "MACROD2 (12%) / 1", "MACROD2 (12%) / 2", "BRAF (12%) / 1", "BRAF (12%) / 2", "PTEN (12%) / 1", "PTEN (12%) / 2", "CCSER1 (10%) / 1", "CCSER1 (10%) / 2", "CDKN2A (7%) / 1", "CDKN2A (7%) / 2", "NAALADL2 (7%) / 1", "NAALADL2 (7%) / 2", "FLCN (10%) / 1", "FLCN (10%) / 2", "ZFP36L2 (7%) / 1", "ZFP36L2 (7%) / 2", "FBXW7 (7%) / 1", "FBXW7 (7%) / 2", "STS (5%) / 1", "STS (5%) / 2", "GATA3 (5%) / 1", "GATA3 (5%) / 2", "GNAS (5%) / 1", "GNAS (5%) / 2", "SEMG2 (5%) / 1", "SEMG2 (5%) / 2", "RNF43 (5%) / 1", "RNF43 (5%) / 2", "KMT2C (5%) / 1", "KMT2C (5%) / 2", "CDK12 (5%) / 1", "CDK12 (5%) / 2", "MTOR (5%) / 1", "MTOR (5%) / 2", "MAP3K21 (5%) / 1", "MAP3K21 (5%) / 2", "SMAD3 (5%) / 1", "SMAD3 (5%) / 2", "CDX2 (5%) / 1", "CDX2 (5%) / 2", "CCND2 (5%) / 1", "CCND2 (5%) / 2", "ELF3 (5%) / 1", "ELF3 (5%) / 2", "FLT3 (5%) / 1", "FLT3 (5%) / 2", "HLA.C (5%) / 1", "HLA.C (5%) / 2")
dt_rownumbers$facet_variable <- factor(dt_rownumbers$facet_variable, levels = levels(dt_rownumbers$facet_variable))

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
