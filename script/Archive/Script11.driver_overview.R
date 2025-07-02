#load packages
library(GenomicRanges)
library(VariantAnnotation)
library(ggplot2)
library(reshape2)
library(BSgenome)

#available.genomes()[1:5]
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg19")
ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
library(BSgenome)
library(grid)
library(gridExtra)
library(VariantAnnotation)
library("readxl")
library(stringr)
library(openxlsx)
library(data.table)
library(forcats)
library(dplyr)

# Define necessary directories for accessing and saving files
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
home_dir <- dirname(dirname(script_dir))
WGS_dir <- file.path(home_dir, "WGS")
WGS_plot_dir <- file.path(home_dir, "Analyse 2.0/7_WGS")
WGS_analysis_files_arne2_dir <- file.path(WGS_dir, "analysis_files_arne2/GCRuns_organoidWGSdata")

# Import all catalog.somatic.tsv files in the WGS folder
dirs<-Sys.glob(file.path(WGS_analysis_files_arne2_dir, "/"))
sampleslist=list.files(dirs, pattern = "catalog.somatic.tsv$",full.names=TRUE, recursive = TRUE)
sampleslist_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
datalist = lapply(sampleslist, function(x){read.table(file=x,sep="\t", header=T)})
names(datalist) <- sampleslist_names
datalist1 <- do.call(rbind, datalist)

#assign mutation number to samples with multiple driver gene mutations
datalist2 <- datalist1 %>% tibble::rownames_to_column(var = "sampleId") %>% tidyr::separate(sampleId, c("sampleId", "mut_nr") ) #MetPret09T 1 mut (NA)
datalist2$mut_nr[is.na(datalist2$mut_nr)] <- 1 #assign nr 1 to MetPret09T with only 1 driver gene mutation
datalist2 <- datalist2  %>% dplyr::filter(driverLikelihood>0.5) #filter for driverLikelihood>0.5

#filter samples without contamination, MSI or 2nd line
datalist2 <- datalist2  %>% dplyr::filter(sampleId!="MetNa22T")  #MSI
datalist2 <- datalist2  %>% dplyr::filter(sampleId!="MetNa07T")  #2nd line
datalist2 <- datalist2  %>% dplyr::filter(sampleId!="MetNa06T") %>% 
  dplyr::filter(sampleId!="MetPret08T")%>% 
  dplyr::filter(sampleId!="MetPret16T")%>% 
  dplyr::filter(sampleId!="MetPret17T")             #contamination
datalist2 <- datalist2  %>% dplyr::filter(sampleId!="MetPret21T") %>% 
  dplyr::filter(sampleId!="MetPret22T")%>% 
  dplyr::filter(sampleId!="MetPret23T")%>% 
  dplyr::filter(sampleId!="MetPret24T")%>% 
  dplyr::filter(sampleId!="MetPret26T")%>% 
  dplyr::filter(sampleId!="RAS32T")%>% 
  dplyr::filter(sampleId!="RAS35T")%>% 
  dplyr::filter(sampleId!="RAS36T")%>% 
  dplyr::filter(sampleId!="RAS37T")%>% 
  dplyr::filter(sampleId!="RAS38T")%>% 
  dplyr::filter(sampleId!="RAS39T")%>% 
  dplyr::filter(sampleId!="RAS40T")            #biopsies and RASTRIC organoids

#save the raw driver gene file
write.table(datalist2 , file = file.path(WGS_plot_dir, "Driveroverview_raw.txt"), sep="\t", col.names = TRUE,quote = F, row.names = FALSE)

#change format of the  raw driver gene file datalist2 > overview
overview <- setNames(data.frame(matrix(ncol = length(unique(datalist2$gene)), nrow = length(unique(datalist2$sampleId)))),unique(datalist2$gene))
row.names(overview) <- unique(datalist2$sampleId)
overview <- overview %>% tibble::rownames_to_column(var = "sampleId") 

#make table with driver genes
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

colorpalette_drivers = c("amplification" =  '#B0D3DC',
                         "deep deletion" =  'coral1',
                         "INDEL" =  'skyblue3',
                         "Synonymous" =  '#fb8072',
                         "NA" =  'white',
                         "missense" =  '#EEB2B0',
                         "nonsense" =  '#A0A9D2',
                         "splice" =  '#b3de69',
                         "inframe" =  '#d9d9d9',
                         "frameshift" = '#FCF6B5',
                         "Nonsense" = 'plum3' ,
                         "grey" = 'gray80')
                         

#save the adapted driver gene file
mat <- overview 
write.table(mat , file = file.path(WGS_plot_dir, "Driveroverview.txt"), sep="\t", col.names = TRUE,quote = F, row.names = FALSE)

row.names(mat) <- mat$sampleId #sample Id as rowname
mat$sampleId <- NULL #remove sample Id column
select_driver_genes <- names(mat)[nrow(mat)- colSums(is.na(mat))>1]  #select driver genes with low frequency in data
mat<-dplyr::select(mat,select_driver_genes) #filter driver genes with low frequency in data

# order driver genes on frequency
driver_names_order <-names(sort(colSums(!is.na(mat)),decreasing = T)) 
mat<-dplyr::select(mat,rev(driver_names_order))
#names(mat)

mat <- tibble::rownames_to_column(mat) #sample Id as column
names(mat)[1]<-("sampleId") #sample Id as column

# re-order levels of frequently occuring drivers
mat %>% 
  dplyr::arrange(factor(APC, levels = c("splice","frameshift","nonsense","nonsense/nonsense","nonsense/frameshift")),
                 factor(TP53, levels = c("inframe","missense","splice","frameshift","deep deletion","missense/missense","missense/nonsense","missense/splice","nonsense")),
                 factor(KRAS,levels = c("missense","missense/missense","amplification")),) -> plotting_data
mat -> plotting_data
levels(as.factor(mat$APC))

row.names(plotting_data) <- plotting_data$sampleId #sample Id as rowname
plotting_data$sampleId <- NULL #remove sample Id column
gene_names_2<-names(plotting_data) #get driver gene names

#modify gene names: genes with percentages (frequency)
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

#rename double character values
plotting_data[plotting_data=="NA/NA"]<-NA 
plotting_data[plotting_data=="deep deletion/deep deletion"]<-"deep deletion" 

#plotting_data <- plotting_data[ c("5CarT","5InvT","p8AdeT","p8Car2T","p8Car3T","p8Inv4T","p9AdeT","p9Mid2T","p9RandT","10AdeT","10CarT","p10AdeNutT","p10AdeWntT","11CarT","11RandT","11InvT","14AdeLT","14CarRT","14InvRT","16CarT","16Inv1T","17AdeLNoWntT","17CarLPN6T","17CarRPN6T","17InvLPN6T"), ]

#format data to split colors (shift) in the ggplot
names(plotting_data)<-seq(1, ncol(plotting_data), by=1)
m <- as.matrix(plotting_data)
dt <- data.table(melt(m))
dt <- dt[, strsplit(as.character(value), "/"), by=list(Var1, Var2)]  # this expands "X/Y/Z" into three rows
dt[, shift:=(1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(Var1, Var2)]
dt[, height:=1/.N, by=list(Var1, Var2)]
dt$test=dt$Var2+dt$shift
dt[, test2 := as.numeric(round(test))]
dt$height2=as.numeric(1)

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

plot
  
#save_plot
ggsave(file=sprintf("%s/%s.pdf",WGS_plot_dir, "drivers_and_supressors"), plot, width = 13, height = 8.5)
dev.off()
      