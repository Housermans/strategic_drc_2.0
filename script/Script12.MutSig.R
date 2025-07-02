#!/usr/bin/env Rscript
library(GenomicRanges)
library(VariantAnnotation)
library(ggplot2)
library(scatterpie)
ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)
library(devtools)
library(readxl)
library(reshape2)
library(nlme)
library(stringr)
library(ggforce)
library(ggrepel)
library(ggalluvial)
library("MutationalPatterns")
library(ggpubr)
library(dplyr)
library(gridExtra)

# Define necessary directories for accessing and saving files
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
home_dir <- dirname(dirname(script_dir))
WGS_dir <- file.path(home_dir, "WGS")
WGS_plot_dir <- file.path(home_dir, "Analyse 2.0/7_WGS")
WGS_analysis_files_arne2_dir <- file.path(WGS_dir, "analysis_files_arne2/GCRuns_organoidWGSdata")

# Import all .purple.somatic.vcf.gz files in the WGS folder
dirs<-Sys.glob(file.path(WGS_analysis_files_arne2_dir))
sampleslist=list.files(dirs, pattern = ".purple.somatic.vcf.gz$",full.names=TRUE, recursive = TRUE)

# Remove contaminated samples, MSI and 2nd line etc.
sampleslist <- sampleslist[!grepl("MetNa22T", sampleslist)] #MSI
sampleslist <- sampleslist[!grepl("MetNa07T", sampleslist)] #2nd line
sampleslist <- sampleslist[!grepl("MetNa06T", sampleslist)] #contaminated 
sampleslist <- sampleslist[!grepl("MetPret08T", sampleslist)] #contaminated 
sampleslist <- sampleslist[!grepl("MetPret16T", sampleslist)] #contaminated 
sampleslist <- sampleslist[!grepl("MetPret17T", sampleslist)] #contaminated 
sampleslist <- sampleslist[!grepl("MetPret21T", sampleslist)] #niet in STRATEGIC
sampleslist <- sampleslist[!grepl("MetPret22T", sampleslist)] #niet in STRATEGIC
sampleslist <- sampleslist[!grepl("MetPret23T", sampleslist)] #niet in STRATEGIC
sampleslist <- sampleslist[!grepl("MetPret24T", sampleslist)] #niet in STRATEGIC
sampleslist <- sampleslist[!grepl("MetPret26T", sampleslist)] #dubbel met 0032 
sampleslist <- sampleslist[!grepl("RAS32T", sampleslist)] #niet in STRATEGIC
sampleslist <- sampleslist[!grepl("RAS35T", sampleslist)] #niet in STRATEGIC
sampleslist <- sampleslist[!grepl("RAS36T", sampleslist)] #niet in STRATEGIC
sampleslist <- sampleslist[!grepl("RAS37T", sampleslist)] #niet in STRATEGIC
sampleslist <- sampleslist[!grepl("RAS38T", sampleslist)] #niet in STRATEGIC
sampleslist <- sampleslist[!grepl("RAS39T", sampleslist)] #niet in STRATEGIC
sampleslist <- sampleslist[!grepl("RAS40T", sampleslist)] #niet in STRATEGIC

# Remove samples not screened.
sampleslist <- sampleslist[!grepl("MetNa12T", sampleslist)]
sampleslist <- sampleslist[!grepl("MetNa15T", sampleslist)]

vcf_files_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
vcf_files_names <- gsub("\\-","\\.",vcf_files_names)
length(vcf_files_names)

VCF_GRList = GenomicRanges::GRangesList()
VAFlist = list()
for(i in 1:length(sampleslist)){
  # Read vcf files with VariantAnnotation
  vcf_object = readVcf(sampleslist[i], "hg38")
  print(vcf_files_names[i])
  #seqlevels(vcf_object) <- paste("chr", seqlevels(vcf_object), sep="")
  vcf_object = vcf_object[ rowRanges(vcf_object)$FILTER == "PASS",]
  vcf_object = vcf_object[ elementNROWS(rowRanges(vcf_object)$ALT) == 1,]
  length(vcf_object)
  
  uniq_df = data.frame(ID =names(rowRanges(vcf_object)),
                       VAF =info(vcf_object)$PURPLE_AF,
                       SampleID = vcf_files_names[i])
  uniq_df$VAF[uniq_df$VAF>1] <- 1
  
  vafplot <- gghistogram(
    uniq_df, x = "VAF", y = "..density..",binwidth = 0.01,
    add = "mean", rug = TRUE,
    fill = "SampleID",
    alette = c("#00AFBB"),
    add_density = TRUE)
  
  
  info(vcf_object)$PURPLE_AF[is.na(info(vcf_object)$PURPLE_AF)] <- 0
  #vcf_object = vcf_object[ info(vcf_object)$PURPLE_AF > 0.3,]
  #seqlevelsStyle(vcf_object) = "UCSC"
  VCF_GRList[[i]] = granges(vcf_object)
  VAFlist[[i]] = vafplot
}

names(VCF_GRList) = vcf_files_names
names(VAFlist) = vcf_files_names
ggarrange(plotlist=VAFlist, widths = c(2,2,2,2))

current_names <- names(VCF_GRList)
new_names <- gsub("Met|a|ret|T", "", current_names)
names(VCF_GRList) <- new_names

#save VAF plot
dirpath=file.path(WGS_plot_dir, "/")
plot_name="VAFplot"
pdf(sprintf("%s%s.pdf",dirpath, plot_name) , useDingbats = F, width = 15, height = 15) 
ggarrange(plotlist=VAFlist, widths = c(2,2,2,2))
dev.off()

VCF_GRList

# Extract mutation matrices
snv_grl <- get_mut_type(VCF_GRList, type = "snv",predefined_dbs_mbs = TRUE)
indel_grl <- get_mut_type(VCF_GRList, type = "indel",predefined_dbs_mbs = TRUE)
dbs_grl <- get_mut_type(VCF_GRList, type = "dbs",predefined_dbs_mbs = TRUE)
mbs_grl <- get_mut_type(VCF_GRList, type = "mbs",predefined_dbs_mbs = TRUE)

mut_mat <- mut_matrix(vcf_list = snv_grl, ref_genome = ref_genome)
indel_grl <- get_indel_context(indel_grl, ref_genome= ref_genome)
indel_counts <- count_indel_contexts(indel_grl)
dbs_grl <- get_dbs_context(dbs_grl)
dbs_counts <- count_dbs_contexts(dbs_grl)

# make combined df with SBS & DBS mutational load
mut_load <- data.frame(SBS = as.data.frame(rowSums(as.data.frame(t(mut_mat)))),
                       indel = as.data.frame(rowSums(as.data.frame(t(indel_counts)))),
                       DBS = as.data.frame(rowSums(as.data.frame(t(dbs_counts)))))
colnames(mut_load) <- c("SBSmuts","indelmuts","DBSmuts")
mut_load <- tibble::rownames_to_column(mut_load, "sample_ID")

# save table SBS & DBS mutational load
write.table(mut_load , file = file.path(WGS_plot_dir, "mut_load.txt"), sep="\t", col.names = TRUE, quote = F,row.names = FALSE)

# plot 96 trinucleotide profile matrix
plot_96_profile(mut_mat)
plot_indel_contexts(indel_counts[,6:8], condensed = TRUE)

dirpath=file.path(WGS_plot_dir, "/")
plot_name="SBS_Mutation_profiles_PASS_variants"
pdf(sprintf("%s%s.pdf",dirpath,plot_name) , useDingbats = F, width = 8, height = 30) 
plot_96_profile(mut_mat)
dev.off()

plot_name="INDEL_Mutation_profiles_PASS_variants"
pdf(sprintf("%s%s.pdf",dirpath,plot_name) , useDingbats = F, width = 15, height = 30) 
plot_indel_contexts(indel_counts)
dev.off()

#mut_load <- left_join(mut_load,SV_load,by=("sample_name_R"))
#mut_load$mut_load <- rowSums(mut_load[c("SBSmuts","indelmuts","DBSmuts")])

colorpalette_sbs = c( '#8dd3c7',
                               '#ffffb3',
                               '#bebada',
                               '#fb8072',
                               '#996ffb',
                               '#fdb462',
                               '#b3de69',
                               '#fccde5',
                               '#d9d9d9',
                               '#ff1417' ,
                               '#ff6611' ,
                               '#c4ff00' ,
                               '#ff8844' ,
                               '#ffee55' ,
                               '#ffff99' ,
                               '#78FA37' ,
                               '#aacc22' ,
                               '#bbdd77' ,
                               '#c8cf82' ,
                               '#92a77e' ,
                               '#5599ee' ,
                               '#0088cc' ,
                               '#226688' ,
                               '#175279' ,
                               '#557777' ,
                               '#ddbb33' ,
                               '#d3a76d' ,
                               '#a9834b' ,
                               '#aa6688',
                               '#767676',
                               '#458B00' ,
                               '#D2691E' ,
                               '#6495ED' ,
                               '#A2CD5A' ,
                               '#CD3333' ,
                               '#7AC5CD' ,
                               '#009ACD' ,
                               '#CD2626' ,
                               '#FFB90F' ,
                               '#76EEC6' ,
                               '#EEB422' ,
                               '#97FFFF' ,
                               '#E9967A' ,
                               '#ff994b',
                               '#463ec0',
                               '#88c928',
                               '#80b1d3',
                               '#68b1c0',
                               '#e34bd9',
                               '#106b00',
                               '#d10073',
                               '#98d76a',
                               '#6b3a9d',
                               '#d5c94e',
                               '#0072e2',
                               '#ff862c',
                               '#31528d',
                               '#d7003a',
                               '#ff4791',
                               '#01837a',
                               '#ff748a',
                               '#777700',
                               '#ff86be',
                               '#4a5822',
                               '#ffabe4',
                               '#6a4e03',
                               '#c6c0fb',
                               '#ffb571',
                               '#873659',
                               '#dea185',
                               '#a0729d')

# make SBS dataframe
signatures_SBS = get_known_signatures(muttype = "snv")
strict_refit <- fit_to_signatures_strict(mut_mat, as.matrix(signatures_SBS), max_delta = 0.004)
contribution_SBS <- as.data.frame(t(strict_refit$fit_res$contribution)) %>% tibble::rownames_to_column("sample_ID")
plot_contribution(strict_refit$fit_res$contribution,as.matrix(signatures_SBS),mode = "relative",palette = colorpalette_sbs,coord_flip = T)

# make SBS contr plot relative & absolute
plot_name="SBS_singatures_relative"
pdf(sprintf("%s%s.pdf",dirpath,plot_name) , useDingbats = F, width = 8, height = 5) 
plot_contribution(strict_refit$fit_res$contribution,as.matrix(signatures_SBS),mode = "relative",palette = colorpalette_sbs,coord_flip = T)
dev.off()
plot_name="SBS_singatures_absolute"
pdf(sprintf("%s%s.pdf",dirpath,plot_name) , useDingbats = F, width = 8, height = 5) 
plot_contribution(strict_refit$fit_res$contribution,as.matrix(signatures_SBS),mode = "absolute",palette = colorpalette_sbs,coord_flip = T)
dev.off()

# make dataframe SBS_singatures
tt <- as.data.frame(t(strict_refit$fit_res$contribution)) %>% tibble::rownames_to_column(., "SampleID") %>% mutate(pretreated = ifelse(substr(SampleID,1,5) == "MetNa", "No","Yes"))
df_cluster_SBS_sample = as.data.frame(t(strict_refit$fit_res$contribution))%>% tibble::rownames_to_column(., "SampleID")
df_cluster_SBS = as.data.frame(t(strict_refit$fit_res$contribution))

#save table
write.table(tt , file = file.path(WGS_plot_dir, "SBS_signature_contribution.txt"), sep="\t", col.names = TRUE, quote = F,row.names = FALSE)

tt$pretreated <- factor(tt$pretreated, levels = c("Yes","No")) # New variable for treatment status added to tt dataframe

dirpathfigures <- file.path(WGS_plot_dir, "/")
signatures <- colnames(tt)[substr(colnames(tt),1,3)=="SBS"]

# create a plot with boxplots of all signatures
output_pdf <- file.path(WGS_plot_dir, "boxplots_all_SBS_signatures_grid.pdf")
plot_list <- list()
signatures <- colnames(tt)[substr(colnames(tt),1,3) == "SBS"]
for (i in signatures) {
  print(i)
  ts <- ggboxplot(tt, x = "pretreated", y = sprintf("%s", i),
                  fill = "white", add = "jitter", color = "pretreated") +
    stat_compare_means(method = "wilcox.test") +
    xlab("") +
    ylab(sprintf("%s Signature mutation contribution", i)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none")
  plot_list[[i]] <- ts
}
pdf(output_pdf, useDingbats = FALSE, width = 25, height = 20)  # Adjust width and height as needed
do.call(grid.arrange, c(plot_list, ncol = 10))  # Adjust ncol for the number of columns in the grid
dev.off()

# cluster by signatures
df_cluster = as.data.frame(t(strict_refit$fit_res$contribution))

# cluster samples based on eucledian distance between relative contribution
dista <- hclust(dist(df_cluster), method = "complete")

# order samples according to clustering
sample_order = rownames(df_cluster)[dista$order]

df_cluster$sample = rownames(df_cluster)
#df_cluster$sample = c("FR11123336")
df_cluster_melt = melt(df_cluster,id=c("sample"))
colnames(df_cluster_melt) = c("Sample","Signature","Contribution")
df_cluster_melt <- within(df_cluster_melt, Sample <- factor(Sample,levels=sample_order))

# make SBS contr plot with clustering on SBS 
clustbysignatures <- ggplot(data = df_cluster_melt, aes(x = factor(Sample), y = Contribution, fill = factor(Signature),order = Sample))+ 
  geom_bar(stat = "identity") + #position="fill",
  # white background
  theme_bw() +
  coord_flip()+
  # no gridlines
  scale_fill_manual(values = colorpalette_sbs,name= "Signature") +
  labs(x = "Sample", y = "Absolute contribution \n (no. mutations)") +
  guides(fill=guide_legend(reverse=FALSE), 
         colour=guide_legend(reverse=TRUE)) +
  theme(legend.position = "right",
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title=element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(colour="grey20",size=18,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=11,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=18,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=18,angle=90,hjust=.5,vjust=.5,face="plain"),
        legend.text=element_text(size=12))

plot_name="SBS_signature_absolute_clustbysignatures"
pdf(sprintf("%s%s.pdf",dirpath,plot_name) , useDingbats = F, width = 10, height = 7) 
plot(clustbysignatures)
dev.off()

# cluster by TMB
TMBall = as.data.frame(t(strict_refit$fit_res$contribution))
sample_order = names(sort(rowSums(TMBall)))

#TMBall$mutload = rowSums(TMBall)

TMBall <- TMBall %>% tibble::rownames_to_column("sample_ID")
TMBall_melt = melt(TMBall)
colnames(TMBall_melt) = c("Sample","Signature", "Contribution")

#df_cluster$sample = rownames(df_cluster)
#df_cluster$sample = c("FR11123336")
#df_cluster_melt = melt(df_cluster,id=c("sample"))
#colnames(df_cluster_melt) = c("Sample","Signature","Contribution")
TMBall_melt <- within(TMBall_melt, Sample <- factor(Sample,levels=rev(sample_order)))

# make SBS contr plot with clustering on TMB
ggplot(data = TMBall_melt, aes(x = factor(Sample), y = Contribution, fill = factor(Signature),order = Sample))+ 
  geom_bar(stat = "identity") + #position="fill",
  # white background
  theme_bw() +
  coord_flip()+
  # no gridlines
  scale_fill_manual(values = colorpalette_sbs,name= "Signature") +
  labs(x = "Sample", y = "Absolute contribution \n (no. mutations)") +
  guides(fill=guide_legend(reverse=FALSE), 
         colour=guide_legend(reverse=TRUE)) +
  theme(legend.position = "right",
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title=element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(colour="grey20",size=18,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=11,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=18,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=18,angle=90,hjust=.5,vjust=.5,face="plain"),
        legend.text=element_text(size=12))

#make new SBS df and select chemo signatures
df_cluster = as.data.frame(t(strict_refit$fit_res$contribution))
df_cluster = df_cluster[c("SBS17a","SBS17b","SBS25","SBS35")]

# cluster samples based on eucledian distance between relative contribution
dista <- hclust(dist(df_cluster), method = "complete")

# order samples according to clustering
sample_order = rownames(df_cluster)[dista$order]

df_cluster$sample = rownames(df_cluster)
#df_cluster$sample = c("FR11123336")
df_cluster_melt = melt(df_cluster,id=c("sample"))
colnames(df_cluster_melt) = c("Sample","Signature","Contribution")
df_cluster_melt <- within(df_cluster_melt, Sample <- factor(Sample,levels=df_cluster$sample))

# make chemo SBS contr plot with clustering on sample ID 
signatures_contribution_chemo <- ggplot(data = df_cluster_melt, aes(x = factor(Sample), y = Contribution, fill = factor(Signature),order = Sample))+ 
  geom_bar(stat = "identity") + #position="fill",
  # white background
  theme_bw() +
  coord_flip()+
  # no gridlines
  scale_fill_manual(values = colorpalette_sbs,name= "Signature") +
  labs(x = "Sample", y = "Absolute contribution \n (no. mutations)") +
  guides(fill=guide_legend(reverse=FALSE), 
         colour=guide_legend(reverse=TRUE)) +
  theme(legend.position = "right",
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title=element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(colour="grey20",size=18,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=11,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=18,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=18,angle=90,hjust=.5,vjust=.5,face="plain"),
        legend.text=element_text(size=12))

plot_name="SBS_signature_chemo"
pdf(sprintf("%s%s.pdf",dirpathfigures,plot_name) , useDingbats = F, width = 6, height = 8) 
plot(signatures_contribution_chemo)
dev.off()

# save table SBS as column names, sample ID as row names
plot_name="mutSign_contributions_SBS_PASSvariants"
write.table(strict_refit$fit_res$contribution, file = sprintf("%s%s.tsv",dirpath,plot_name), sep="\t", col.names = TRUE, row.names = TRUE)

# explain treatment signatures the increase in TMB
TMBall = as.data.frame(t(strict_refit$fit_res$contribution))
sample_order = names(sort(rowSums(TMBall)))
TMBall$mutload = rowSums(TMBall) #sum off al signatures

TMBall$mutload_wo_treatmentmutations <- (TMBall$mutload - TMBall$SBS17a - TMBall$SBS17b - TMBall$SBS25 - TMBall$SBS35) #sum minus chemo sign
TMBall <- TMBall %>% tibble::rownames_to_column(., "SampleID") %>% mutate(pretreated = ifelse(substr(SampleID,1,5) == "MetNa", "No","Yes"))
TMBall$pretreated <- factor(TMBall$pretreated, levels = c("Yes","No"))

# boxplot mut load without chemo treatment mutations
ggboxplot(TMBall, x = "pretreated", y = "mutload_wo_treatmentmutations",
                fill = "white",add = "jitter",color="pretreated")+  #fill = "variable", #palette = c("#00AFBB", "#E7B800", "#FC4E07")
  #stat_compare_means(aes(label = ..p.signif..))+
  stat_compare_means(method = "wilcox.test")+
  xlab(c(""))+
  ylab("mut load w/o treatment mutations")+
  #scale_y_continuous(limits = c(0, 0.5))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")

# make INDEL dataframe
signatures_indel = get_known_signatures(muttype = "indel")
strict_refit <- fit_to_signatures_strict(indel_counts, signatures_indel, max_delta = 0.004)

# plot INDEL contribution absulote
acontribution_indel <- as.data.frame(t(strict_refit$fit_res$contribution)) %>% tibble::rownames_to_column("sample_ID")
plot_name="INDEL_signature_contributions_PASS_variants"
pdf(sprintf("%s%s.pdf",dirpathfigures,plot_name) , useDingbats = F, width = 8, height = 5) 
plot_contribution(strict_refit$fit_res$contribution,as.matrix(signatures_SBS),coord_flip = TRUE,mode = "absolute")
dev.off()

plot_name="mutSign_contributions_indel_PASSvariants"
write.table(strict_refit$fit_res$contribution, file = sprintf("%s%s.tsv", WGS_plot_dir ,plot_name), sep="\t", col.names = TRUE, row.names = TRUE,quote = F)

tt <- as.data.frame(t(strict_refit$fit_res$contribution)) %>% tibble::rownames_to_column(., "SampleID") %>% mutate(pretreated = ifelse(substr(SampleID,1,5) == "MetNa", "No","Yes"))
tt$pretreated <- factor(tt$pretreated, levels = c("Yes","No"))

# create a plot with boxplots of all INDELs
output_pdf <- file.path(WGS_plot_dir, "boxplots_all_INDEL_signatures_grid.pdf")
plot_list <- list()
signatures <- colnames(tt)[substr(colnames(tt),1,2)=="ID"]
for (i in signatures) {
  print(i)
  ts <- ggboxplot(tt, x = "pretreated", y = sprintf("%s", i),
                  fill = "white", add = "jitter", color = "pretreated") +
    stat_compare_means(method = "wilcox.test") +
    xlab("") +
    ylab(sprintf("%s Signature mutation contribution", i)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none")
  plot_list[[i]] <- ts
}
pdf(output_pdf, useDingbats = FALSE, width = 20, height = 7)  # Adjust width and height as needed
do.call(grid.arrange, c(plot_list, ncol = 9))  # Adjust ncol for the number of columns in the grid
dev.off()

# make DBS dataframe
signatures_DBS = get_known_signatures(muttype = "dbs")
strict_refit <- fit_to_signatures_strict(dbs_counts, signatures_DBS, max_delta = 0.004)
contribution_DBS <- as.data.frame(t(strict_refit$fit_res$contribution)) %>% tibble::rownames_to_column("sample_ID")

# plot DBS contribution absolute
plot_name="mutSign_contributions_DBS_PASSvariants"
pdf(sprintf("%s%s.pdf", dirpathfigures,plot_name) , useDingbats = F, width = 8, height = 5) 
plot_contribution(strict_refit$fit_res$contribution,signatures_DBS,coord_flip = TRUE,mode = "absolute")
dev.off()

tt <- as.data.frame(t(strict_refit$fit_res$contribution)) %>% tibble::rownames_to_column(., "SampleID") %>% mutate(pretreated = ifelse(substr(SampleID,1,5) == "MetNa", "No","Yes"))
tt$pretreated <- factor(tt$pretreated, levels = c("Yes","No"))

# create a plot with boxplots of all DBS
output_pdf <- file.path(WGS_plot_dir, "boxplots_all_DBS_signatures_grid.pdf")
plot_list <- list()
signatures <- colnames(tt)[substr(colnames(tt),1,3)=="DBS"]
for (i in signatures) {
  print(i)
  ts <- ggboxplot(tt, x = "pretreated", y = sprintf("%s", i),
                  fill = "white", add = "jitter", color = "pretreated") +
    stat_compare_means(method = "wilcox.test") +
    xlab("") +
    ylab(sprintf("%s Signature mutation contribution", i)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none")
  plot_list[[i]] <- ts
}
pdf(output_pdf, useDingbats = FALSE, width = 20, height = 7)  # Adjust width and height as needed
do.call(grid.arrange, c(plot_list, ncol = 6))  # Adjust ncol for the number of columns in the grid
dev.off()

df_cluster_DBS_sample = as.data.frame(t(strict_refit$fit_res$contribution))%>% tibble::rownames_to_column(., "SampleID")
df_cluster_DBS = as.data.frame(t(strict_refit$fit_res$contribution))

ttt <- cbind(df_cluster_SBS,df_cluster_DBS) # Combine SBS & DBS dataframe
ttt_sample <- cbind(df_cluster_SBS_sample,df_cluster_DBS_sample) # Combine SBS & DBS dataframe

#save table
write.table(ttt_sample , file = file.path(WGS_plot_dir, "SBS_DBS_signature_contribution.txt"), sep="\t", col.names = TRUE, quote = F,row.names = FALSE)

df_cluster = ttt[c("SBS17a","SBS17b","SBS25","SBS35","DBS5")] #select chemo treatment signatures

# cluster samples based on eucledian distance between relative contribution
dista <- hclust(dist(df_cluster), method = "complete")

# order samples according to clustering
sample_order = rownames(df_cluster)[dista$order]

# plot clustering
dhc <- as.dendrogram(dista) 
plot(dhc)

ttt <- ttt[c("SBS17a","SBS17b","SBS25","SBS35","DBS5")]
ttt <- ttt %>% tibble::rownames_to_column("sample")
df_cluster_melt = melt(ttt,id=c("sample"))
colnames(df_cluster_melt) = c("Sample","Signature","Contribution")
df_cluster_melt <- df_cluster_melt %>% dplyr::mutate(stype = ifelse(substr(Signature,1,3)=="DBS","DBS","SBS"))

colorpalette_SBS_DBS_contribution_chemo = c("#F15D3F",
                                            '#EEB3B0',
                                            '#A0A9D2', 
                                            "#FCF6B7",
                                             '#7FBFF5')
              
SBS_DBS_contribution_chemo<- ggplot(data = df_cluster_melt, aes(x = factor(Sample), y = Contribution, fill = factor(Signature),order = Sample))+ 
  geom_bar(stat = "identity") + #position="fill",
  # white background
  theme_bw() +
  facet_wrap(~stype, scales = "free")+
  coord_flip()+
  # no gridlines
  scale_fill_manual(values = colorpalette_SBS_DBS_contribution_chemo,name= "Signature") +
  labs(x = "Sample", y = "Absolute contribution \n (no. mutations)") +
  guides(fill=guide_legend(reverse=FALSE), 
         colour=guide_legend(reverse=TRUE)) +
  theme(legend.position = "right",
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title=element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(colour="grey20",size=18,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=11,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=18,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=18,angle=90,hjust=.5,vjust=.5,face="plain"),
        legend.text=element_text(size=12))

plot_name="SBS_DBS_contribution_chemo"
pdf(sprintf("%s%s.pdf", dirpathfigures ,plot_name) , useDingbats = F, width = 10, height = 8) 
plot(SBS_DBS_contribution_chemo)
dev.off()
