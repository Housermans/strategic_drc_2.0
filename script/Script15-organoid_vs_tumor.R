#!/usr/bin/env Rscript
library(GenomicRanges)
library(VariantAnnotation)
library(dplyr)
library(ggplot2)
library(scatterpie)
library(GenomicRanges)
library(VariantAnnotation)
library(dplyr)
library(ggplot2)
library(scatterpie)
ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(devtools)
library(readxl)
library(reshape2)
library(nlme)
library(stringr)
library(ggforce)
library(ggrepel)
library(ggalluvial)
library(MutationalPatterns)
library(ggpubr)
library(data.table)

# Define necessary directories for accessing and saving files
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
home_dir <- dirname(dirname(script_dir))
WGS_dir <- file.path(home_dir, "WGS")
WGS_plot_dir <- file.path(home_dir, "Analyse 2.0/7_WGS")
WGS_analysis_files_arne2_dir <- file.path(WGS_dir, "analysis_files_arne2/GCRuns_organoidWGSdata")

#metadata
metadata <- as.data.frame(read_xlsx(file.path(home_dir,"WGS/analysis_files_arne2/analysis_update/organoid_vs_tumor/metadata_organoid_vs_tumor.xlsx"), col_names = T,sheet = "Sheet1"))
metadata <- metadata %>% dplyr::filter(Match=="Yes")

#organoid data
dirs<-Sys.glob(file.path(WGS_analysis_files_arne2_dir, "/"))
sampleslist=list.files(dirs, pattern = ".purple.somatic.vcf.gz$",full.names=TRUE, recursive = TRUE)
vcf_files_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
vcf_files_names <- gsub("\\-","\\.",vcf_files_names)
length(vcf_files_names)
sampleslist <- sampleslist[vcf_files_names %in% metadata$GCrunID]

#tissue data
dirs<-Sys.glob("L:/lsmabers/OPTIC WGS data/somatic~/tmp/DR-OPTIC-250109")
sampleslist=list.files(dirs, pattern = ".purple.somatic.vcf.gz$",full.names=TRUE, recursive = TRUE)
vcf_files_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
vcf_files_names <- gsub("\\-","\\.",vcf_files_names)
length(vcf_files_names)
sampleslist <- sampleslist[vcf_files_names %in% metadata$Hartwig_ID]

#Load datasets
vcf_files_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
vcf_files_names <- gsub("\\-","\\.",vcf_files_names)
length(vcf_files_names)
VCF_GRList = GenomicRanges::GRangesList()
VAFlist = list()
subclonal_dt = setNames(data.frame(matrix(ncol = 3, nrow = 0)),c("sampleID","clonalmuts","subclonalmuts"))

for(i in 1:length(sampleslist)){
  # Read vcf files with VariantAnnotation
  vcf_object = readVcf(sampleslist[i], "hg38")
  print(vcf_files_names[i])
  samplename <- vcf_files_names[i]
  
  #samplename <- metadata %>% dplyr::filter(GCrunID==samplename)%>% dplyr::pull(Organoid_ID)
  samplename <- metadata %>% dplyr::filter(Hartwig_ID==samplename)%>% dplyr::pull(Organoid_ID)
  
  
  #seqlevels(vcf_object) <- paste("chr", seqlevels(vcf_object), sep="")
  vcf_object = vcf_object[ rowRanges(vcf_object)$FILTER == "PASS",]
  vcf_object = vcf_object[ elementNROWS(rowRanges(vcf_object)$ALT) == 1,]
  length(vcf_object)
  
  uniq_df = data.frame(ID =names(rowRanges(vcf_object)),
                       VAF =info(vcf_object)$PURPLE_AF,
                       SUBCL=info(vcf_object)$SUBCL,
                       SampleID = samplename)
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
  clonalDT=setNames(data.frame(matrix(ncol = 3, nrow = 1)),c("sampleID","clonalmuts","subclonalmuts"))
  clonalDT$sampleID <- samplename
  clonalDT$clonalmuts <- uniq_df %>% dplyr::filter(SUBCL<=0.8) %>% nrow()
  clonalDT$subclonalmuts <- uniq_df %>% dplyr::filter(SUBCL>0.8) %>% nrow()
  subclonal_dt <- rbind(subclonal_dt,clonalDT)
  
}
names(VCF_GRList) = vcf_files_names
names(VAFlist) = vcf_files_names
ggarrange(plotlist=VAFlist, widths = c(2,2,2,2))
organoid <- subclonal_dt
hartwig <- subclonal_dt

organoid$tissue <- "organoid"
hartwig$tissue <- "tumor_tissue"
organoid$clonality_organoid <- organoid$clonalmuts/(organoid$subclonalmuts+organoid$clonalmuts)
hartwig$clonality_tumor <- hartwig$clonalmuts/(hartwig$subclonalmuts+hartwig$clonalmuts)

clonality_matrix <- left_join(organoid[c("sampleID","clonality_organoid")],hartwig[c("sampleID","clonality_tumor")],by="sampleID")
clonality_plot<- ggplot() +
  geom_point(data=clonality_matrix,aes(x=clonality_organoid, y=clonality_tumor),colour="black",size=5) +
  xlab("Clonality PDO") + ylab("Clonality tumor") +
  theme_bw(base_size=22) +scale_y_continuous(limits = c(0, 1))+ scale_x_continuous(limits = c(0, 1))
setwd(WGS_plot_dir)
ggsave(filename=paste('tumor_PDO_clonality.pdf'), clonality_plot , width =4, height = 4)

#organoid data
dirs<-Sys.glob(file.path(WGS_analysis_files_arne2_dir, "/"))
sampleslist=list.files(dirs, pattern = "purple.purity.tsv$",full.names=TRUE, recursive = TRUE)
vcf_files_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
vcf_files_names <- gsub("\\-","\\.",vcf_files_names)
length(sampleslist)
sampleslist <- sampleslist[vcf_files_names %in% metadata$GCrunID]
sampleslist_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
datalist = lapply(sampleslist, function(x){read.table(file=x,sep="\t", header=T)})
names(datalist) <- sampleslist_names
datalist_organoid <- do.call(rbind, datalist)
colnames(datalist_organoid) <- paste0(colnames(datalist_organoid),"_organoid")

#tissue data
dirs<-Sys.glob("L:/lsmabers/OPTIC WGS data/somatic~/tmp/DR-OPTIC-250109")
sampleslist=list.files(dirs, pattern = "purple.purity.tsv$",full.names=TRUE, recursive = TRUE)
vcf_files_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
vcf_files_names <- gsub("\\-","\\.",vcf_files_names)
length(sampleslist)

sampleslist <- sampleslist[vcf_files_names %in% metadata$Hartwig_ID]
sampleslist_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
datalist = lapply(sampleslist, function(x){read.table(file=x,sep="\t", header=T)})
names(datalist) <- sampleslist_names
datalist_hartwig <- do.call(rbind, datalist)
colnames(datalist_hartwig) <- paste0(colnames(datalist_hartwig),"_tissue")

datalist_organoid <- datalist_organoid %>% tibble::rownames_to_column("GCrunID")
datalist_organoid <- left_join(datalist_organoid,metadata,by="GCrunID")
datalist_hartwig <- datalist_hartwig%>% tibble::rownames_to_column("Hartwig_ID")

# Functie om scatterplots te maken
make_plot <- function(df1, df2, var1, var2, xlab_text, ylab_text, lim, scale_factor = 1) {
  data <- left_join(df1[c("Hartwig_ID", var1)], df2[c("Hartwig_ID", var2)], by="Hartwig_ID") %>%
    mutate(across(c(var1, var2), ~ . * scale_factor))
  
  ggplot(data) +
    geom_point(aes_string(x=var1, y=var2), colour="black", size=5) +
    xlab(xlab_text) + ylab(ylab_text) + theme_bw(base_size=22) +
    scale_x_continuous(limits=c(0, lim)) + scale_y_continuous(limits=c(0, lim))
}

# Plots genereren
purity_plot <- make_plot(datalist_organoid, datalist_hartwig, "purity_organoid", "purity_tissue", "Purity PDO (%)", "Purity tumor (%)", 100, 100)
ploidy_plot <- make_plot(datalist_organoid, datalist_hartwig, "diploidProportion_organoid", "diploidProportion_tissue", "Diploid proportion PDO", "Diploid proportion tumor", 0.2)
tml_plot <- make_plot(datalist_organoid, datalist_hartwig, "tml_organoid", "tml_tissue", "tml PDO", "tml tumor", 200)
tmb_plot <- make_plot(datalist_organoid, datalist_hartwig, "tmbPerMb_organoid", "tmbPerMb_tissue", "tmb per Mb PDO", "tmb per Mb tumor", 15)
svtmb_plot <- make_plot(datalist_organoid, datalist_hartwig, "svTumorMutationalBurden_organoid", "svTumorMutationalBurden_tissue", "sv tmb PDO", "sv tmb tumor", 500)

library(cowplot)
grid <- plot_grid(
  purity_plot,
  tml_plot,
  tmb_plot,
  svtmb_plot,
  nrow = 2)
setwd(WGS_plot_dir)
ggsave(filename=paste('tumor_PDO.pdf'), grid , width =8, height = 8)

#drivers tissue
dirs<-Sys.glob("L:/lsmabers/OPTIC WGS data/somatic~/tmp/DR-OPTIC-250109")
sampleslist=list.files(dirs, pattern = "catalog.somatic.tsv$",full.names=TRUE, recursive = TRUE)
vcf_files_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
vcf_files_names <- gsub("\\-","\\.",vcf_files_names)
length(sampleslist)
sampleslist <- sampleslist[vcf_files_names %in% metadata$Hartwig_ID]
length(sampleslist)
sampleslist_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
datalist = lapply(sampleslist, function(x){read.table(file=x,sep="\t", header=T)})
names(datalist) <- sampleslist_names
datalist1 <- do.call(rbind, datalist)
#rownames(datalist1) <- gsub('-','',rownames(datalist1))
datalist2 <- datalist1 %>% tibble::rownames_to_column(var = "Hartwig_ID") %>% tidyr::separate(Hartwig_ID, c("Hartwig_ID", "mut_nr") )
datalist2 <- datalist2  %>% dplyr::filter(driverLikelihood>0.5)
datalist2 <- datalist2 %>% 
  rename(
    sampleId=Hartwig_ID
  )
#start processing
overview <- setNames(data.frame(matrix(ncol = length(unique(datalist2$gene)), nrow = length(unique(datalist2$sampleId)))),unique(datalist2$gene))
row.names(overview) <- unique(datalist2$sampleId)
overview <- overview %>% tibble::rownames_to_column(var = "sampleId") 

for (geneID in unique(datalist2$gene)){
  for(sampleID in unique(datalist2$sampleId)){
    
    tt <- datalist2  %>% dplyr::filter(sampleId==sampleID) %>% dplyr::filter(gene==geneID) 
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

mat <- overview
row.names(mat) <- mat$sampleId
mat$sampleId <- NULL
#select_driver_genes <- names(mat)[nrow(mat)- colSums(is.na(mat))>1]
#mat<-dplyr::select(mat,select_driver_genes)

driver_names_order <-names(sort(colSums(!is.na(mat)),decreasing = T))
mat<-dplyr::select(mat,rev(driver_names_order))
#names(mat)

mat <- tibble::rownames_to_column(mat)
names(mat)[1]<-("sampleId")
mat %>% 
  dplyr::arrange(factor(APC, levels = c("splice","frameshift","nonsense","nonsense/nonsense","nonsense/frameshift")),
                 factor(TP53, levels = c("inframe","missense","splice","frameshift","deep deletion","missense/missense","missense/nonsense","missense/splice","nonsense")),
                 factor(KRAS,levels = c("missense","missense/missense","amplification")),) -> plotting_data
mat -> plotting_data
levels(as.factor(mat$APC))

# Subset the mat dataframe to include only frequent drivers
tumor_drivers <- colnames(mat)
gene_order <- colnames(mat)
mat <- mat[, tumor_drivers]

row.names(plotting_data) <- plotting_data$sampleId
plotting_data$sampleId <- NULL
#na_count <- sapply(plotting_data, function(y) sum(length(which(is.na(y)))))
#plotting_data <- plotting_data[na_count!=2]
#gene_names_2<-names(plotting_data)
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
  
  
  #df2 <- left_join(df2, dplyr::select(cancer_genes,Gene,Gene_MoA), by = "Gene")
  #df2 <- df2[!duplicated(df2), ] 
  #df2 %>%
  #  dplyr::mutate(MOA = ifelse(Gene_MoA == "LoF", " ", 
  #                             ifelse(Gene_MoA == "Act", "*","NA"))) -> df2
  df2$percentage_ID_2 <- with(df2, paste("(",percentage_ID,")", sep=""))
  #df2$Gene_2 <- with(df2, paste(MOA,Gene, sep=""))
  cols<-c("Gene","percentage_ID_2")
  df2$unite <- apply( df2[ , cols ] , 1 , paste , collapse = " " )
  return(df2$unite)
}
#gene_names<-names(plotting_data)
#gene_names<-names(plotting_data)
gene_names<-modify_gene_names(plotting_data)

plotting_data <- tibble::rownames_to_column(plotting_data,var = "sampleId")
row.names(plotting_data) <- plotting_data$sampleId
plotting_data$sampleId <- NULL

plotting_data[plotting_data=="NA/NA"]<-NA
plotting_data[plotting_data=="deep deletion/deep deletion"]<-"deep deletion"

plotting_data <- plotting_data[ c("OPTC01030042T","OPTC01080002T","OPTC01080013T","OPTC01080019T","OPTC01030048T","OPTC01080023T"), ]


names(plotting_data)<-seq(1, ncol(plotting_data), by=1)
m <- as.matrix(plotting_data)
dt <- data.table(melt(m))
dt <- dt[, strsplit(as.character(value), "/"), by=list(Var1, Var2)]  # this expands "X/Y/Z" into three rows
dt[, shift:=(1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(Var1, Var2)]
dt[, height:=1/.N, by=list(Var1, Var2)]
dt$test=dt$Var2+dt$shift
dt$test2=as.numeric(format(dt$test,digits =0))
dt$height2=as.numeric(1)

colorpalette_drivers = c("amplification" =  '#b3de69',
                         "deep deletion" =  '#F15D3F',
                         "INDEL" =  'skyblue3',
                         "Synonymous" =  '#fb8072',
                         "NA" =  'white',
                         "missense" =  '#EEB2B0',
                         "nonsense" =  '#7FBFF5',
                         "splice" =  '#A0A9D2',
                         "inframe" =  '#E0E1E0',
                         "frameshift" = '#FCF6B5',
                         "Nonsense" = 'plum3' ,
                         "grey" = 'gray80')


plot <- ggplot() + 
  geom_tile(data = dt,aes(Var1,y=test, fill=V1, height=height),
            color="grey", size=0.1) + xlab('sample ID') + ylab('Gene')+
  scale_y_continuous(breaks = seq(1, ncol(plotting_data), by=1), labels = gene_names)+
  scale_fill_manual(values=colorpalette_drivers,na.value = "white") +
  guides(fill=guide_legend(title="Effect"))+
  theme(axis.title.y=element_text(size=20,vjust=3),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=20),
        axis.text.x = element_text(colour="grey20",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.ticks.y = element_line(size=0.5),
        strip.text.x=element_text(size=15),
        strip.text.y=element_text(size=15,vjust=2),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = 'grey20',size = 0.5),
        strip.background=element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))

geom_tile(data = dt,aes(Var1,y=test2, height=height2),
          color="gray46", size=0.1,alpha=0)
plot(plot)
ggsave(file=sprintf("%s/%s.pdf", WGS_plot_dir, "drivers_and_suppressors_tumor"),plot, width = 6, height = 10)
#rm(list = ls())

dirs<-Sys.glob(file.path(WGS_analysis_files_arne2_dir, "/"))
##dirs<-Sys.glob("/Users/avanhoeck/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/201106_HMFregWIDE_FR22825413_FR22825420_WIDE01011135/purple/")
sampleslist=list.files(dirs, pattern = "catalog.somatic.tsv$",full.names=TRUE, recursive = TRUE)
vcf_files_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
vcf_files_names <- gsub("\\-","\\.",vcf_files_names)
sampleslist <- sampleslist[vcf_files_names %in% metadata$GCrunID]
length(sampleslist)
sampleslist_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
datalist = lapply(sampleslist, function(x){read.table(file=x,sep="\t", header=T)})
names(datalist) <- sampleslist_names
datalist1 <- do.call(rbind, datalist)
#rownames(datalist1) <- gsub('-','',rownames(datalist1))
datalist2 <- datalist1 %>% tibble::rownames_to_column(var = "GCrunID") %>% tidyr::separate(GCrunID, c("GCrunID", "mut_nr") )
datalist2 <- left_join(datalist2,metadata[c("Hartwig_ID","GCrunID")],by="GCrunID")
datalist2 <- datalist2  %>% dplyr::filter(driverLikelihood>0.5)
datalist2$GCrunID <- NULL
datalist2 <- datalist2 %>% 
  rename(
    sampleId=Hartwig_ID
  )

#start processing
overview <- setNames(data.frame(matrix(ncol = length(unique(datalist2$gene)), nrow = length(unique(datalist2$sampleId)))),unique(datalist2$gene))
row.names(overview) <- unique(datalist2$sampleId)
overview <- overview %>% tibble::rownames_to_column(var = "sampleId") 

for (geneID in unique(datalist2$gene)){
  for(sampleID in unique(datalist2$sampleId)){
    
    tt <- datalist2  %>% dplyr::filter(sampleId==sampleID) %>% dplyr::filter(gene==geneID) 
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

mat <- overview
row.names(mat) <- mat$sampleId
mat$sampleId <- NULL
#select_driver_genes <- names(mat)[nrow(mat)- colSums(is.na(mat))>1]
#mat<-dplyr::select(mat,select_driver_genes)

driver_names_order <-names(sort(colSums(!is.na(mat)),decreasing = T))
mat<-dplyr::select(mat,rev(driver_names_order))
#names(mat)

mat <- tibble::rownames_to_column(mat)
names(mat)[1]<-("sampleId")
mat %>% 
  dplyr::arrange(factor(APC, levels = c("splice","frameshift","nonsense","nonsense/nonsense","nonsense/frameshift")),
                 factor(TP53, levels = c("inframe","missense","splice","frameshift","deep deletion","missense/missense","missense/nonsense","missense/splice","nonsense")),
                 factor(KRAS,levels = c("missense","missense/missense","amplification")),) -> plotting_data
mat -> plotting_data
levels(as.factor(mat$APC))

# Subset the mat dataframe to include only frequent drivers
plotting_data<-dplyr::select(mat,(gene_order))

row.names(plotting_data) <- plotting_data$sampleId
plotting_data$sampleId <- NULL
#na_count <- sapply(plotting_data, function(y) sum(length(which(is.na(y)))))
#plotting_data <- plotting_data[na_count!=2]
#gene_names_2<-names(plotting_data)
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
  
  
  #df2 <- left_join(df2, dplyr::select(cancer_genes,Gene,Gene_MoA), by = "Gene")
  #df2 <- df2[!duplicated(df2), ] 
  #df2 %>%
  #  dplyr::mutate(MOA = ifelse(Gene_MoA == "LoF", " ", 
  #                             ifelse(Gene_MoA == "Act", "*","NA"))) -> df2
  df2$percentage_ID_2 <- with(df2, paste("(",percentage_ID,")", sep=""))
  #df2$Gene_2 <- with(df2, paste(MOA,Gene, sep=""))
  cols<-c("Gene","percentage_ID_2")
  df2$unite <- apply( df2[ , cols ] , 1 , paste , collapse = " " )
  return(df2$unite)
}
#gene_names<-names(plotting_data)
#gene_names<-names(plotting_data)
gene_names<-modify_gene_names(plotting_data)

plotting_data <- tibble::rownames_to_column(plotting_data,var = "sampleId")
row.names(plotting_data) <- plotting_data$sampleId
plotting_data$sampleId <- NULL

plotting_data[plotting_data=="NA/NA"]<-NA
plotting_data[plotting_data=="deep deletion/deep deletion"]<-"deep deletion"

plotting_data <- plotting_data[ c("OPTC01030042T","OPTC01080002T","OPTC01080013T","OPTC01080019T","OPTC01030048T","OPTC01080023T"), ]


names(plotting_data)<-seq(1, ncol(plotting_data), by=1)
m <- as.matrix(plotting_data)
dt <- data.table(melt(m))
dt <- dt[, strsplit(as.character(value), "/"), by=list(Var1, Var2)]  # this expands "X/Y/Z" into three rows
dt[, shift:=(1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(Var1, Var2)]
dt[, height:=1/.N, by=list(Var1, Var2)]
dt$test=dt$Var2+dt$shift
dt$test2=as.numeric(format(dt$test,digits =0))
dt$height2=as.numeric(1)

colorpalette_drivers = c("amplification" =  '#b3de69',
                         "deep deletion" =  '#F15D3F',
                         "INDEL" =  'skyblue3',
                         "Synonymous" =  '#fb8072',
                         "NA" =  'white',
                         "missense" =  '#EEB2B0',
                         "nonsense" =  '#7FBFF5',
                         "splice" =  '#A0A9D2',
                         "inframe" =  '#E0E1E0',
                         "frameshift" = '#FCF6B5',
                         "Nonsense" = 'plum3' ,
                         "grey" = 'gray80')


plot <- ggplot() + 
  geom_tile(data = dt,aes(Var1,y=test, fill=V1, height=height),
            color="grey", size=0.1) + xlab('sample ID') + ylab('Gene')+
  scale_y_continuous(breaks = seq(1, ncol(plotting_data), by=1), labels = gene_names)+
  scale_fill_manual(values=colorpalette_drivers,na.value = "white") +
  guides(fill=guide_legend(title="Effect"))+
  theme(axis.title.y=element_text(size=20,vjust=3),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=20),
        axis.text.x = element_text(colour="grey20",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.ticks.y = element_line(size=0.5),
        strip.text.x=element_text(size=15),
        strip.text.y=element_text(size=15,vjust=2),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = 'grey20',size = 0.5),
        strip.background=element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))

geom_tile(data = dt,aes(Var1,y=test2, height=height2),
          color="gray46", size=0.1,alpha=0)
plot(plot)
ggsave(file=sprintf("%s/%s.pdf", WGS_plot_dir, "drivers_and_suppressors_PDO"),plot, width = 6, height = 10)
#rm(list = ls())












#signatures (lukt niet ivm ref genome)

#tissue data
dirs<-Sys.glob("L:/lsmabers/OPTIC WGS data/somatic~/tmp/DR-OPTIC-250109")
sampleslist=list.files(dirs, pattern = ".purple.somatic.vcf.gz$",full.names=TRUE, recursive = TRUE)
vcf_files_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
vcf_files_names <- gsub("\\-","\\.",vcf_files_names)
length(vcf_files_names)
sampleslist <- sampleslist[vcf_files_names %in% metadata$Hartwig_ID]

#Load datasets
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
plot_name="VAFplot_tumor"
pdf(sprintf("%s%s.pdf",dirpath, plot_name) , useDingbats = F, width = 15, height = 15) 
ggarrange(plotlist=VAFlist, widths = c(2,2,2,2))
dev.off()

# Extract mutation matrices
snv_grl <- get_mut_type(VCF_GRList, type = "snv",predefined_dbs_mbs = TRUE)
indel_grl <- get_mut_type(VCF_GRList, type = "indel",predefined_dbs_mbs = TRUE)
dbs_grl <- get_mut_type(VCF_GRList, type = "dbs",predefined_dbs_mbs = TRUE)
mbs_grl <- get_mut_type(VCF_GRList, type = "mbs",predefined_dbs_mbs = TRUE)

mut_mat <- mut_matrix(vcf_list = snv_grl, ref_genome = ref_genome)
