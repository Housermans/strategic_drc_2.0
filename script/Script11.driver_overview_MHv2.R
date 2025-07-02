# Version control
# Version v2: added seperate chemonaive & pretreated plot 

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

transform_sample_id <- function(sample_id) {
  sapply(sample_id, function(id) {
    if (startsWith(id, "MetNa")) {
      return(paste0("N", sprintf("%02d", as.numeric(substr(id, 6, 7)))))
    } else if (startsWith(id, "MetPret")) {
      return(paste0("P", sprintf("%02d", as.numeric(substr(id, 8, 9)))))
    } else {
      return(id)  # Return unchanged if it doesn't match the patterns
    }
  })
}

# Import all catalog.somatic.tsv files in the WGS folder
dirs<-Sys.glob(file.path(WGS_analysis_files_arne2_dir, "/"))
sampleslist=list.files(dirs, pattern = "catalog.somatic.tsv$",full.names=TRUE, recursive = TRUE)
sampleslist_names = basename(sampleslist) %>% gsub(pattern = "\\..*$",replacement =  "")
datalist = lapply(sampleslist, function(x){read.table(file=x,sep="\t", header=T)})
names(datalist) <- sampleslist_names

# After reading the data
datalist1 <- do.call(rbind, datalist)
datalist1$sampleId <- transform_sample_id(rownames(datalist1))

# When creating datalist2
datalist2 <- datalist1 %>% 
  tibble::rownames_to_column(var = "original_sampleId") %>% 
  tidyr::separate(original_sampleId, c("sampleId", "mut_nr"), fill = "right") %>%
  mutate(sampleId = transform_sample_id(sampleId))

# Update the filtering steps
datalist2 <- datalist2 %>% 
  dplyr::filter(sampleId != "N22") %>%  # MSI
  dplyr::filter(sampleId != "N07") %>%  # 2nd line
  dplyr::filter(!sampleId %in% c("N06", "P08", "P16", "P17")) %>%  # contamination
  dplyr::filter(!sampleId %in% c("P21", "P22", "P23", "P24", "P26")) %>%
  dplyr::filter(!startsWith(sampleId, "RAS"))  # RAS samples

#save the raw driver gene file
write.table(datalist2 , file = file.path(WGS_plot_dir, "Driveroverview_raw_v2.txt"), sep="\t", col.names = TRUE,quote = F, row.names = FALSE)

# Update the overview creation
overview <- setNames(data.frame(matrix(ncol = length(unique(datalist2$gene)), 
                                       nrow = length(unique(datalist2$sampleId)))),
                     unique(datalist2$gene))
row.names(overview) <- unique(datalist2$sampleId)
overview <- overview %>% tibble::rownames_to_column(var = "sampleId")

# Make table with driver genes
for (geneID in unique(datalist2$gene)) {
  for (sampleID in unique(datalist2$sampleId)) {
    tt <- datalist2 %>% dplyr::filter(sampleId == sampleID) %>% dplyr::filter(gene == geneID)
    
    if (nrow(tt) > 1) {
      tt <- tt %>% dplyr::filter(isCanonical == "true")
    }
    
    if (nrow(tt) == 0) {
      next
    } else {
      simplemuts <- tt[c("missense", "nonsense", "splice", "inframe", "frameshift")]
      print(tt)
      
      # Use any() to check if any row meets the condition
      if (any(tt$likelihoodMethod == "DEL")) {
        overview[overview$sampleId == sampleID, geneID] <- "deep deletion"
      } else if (any(tt$likelihoodMethod == "AMP")) {
        if (sum(colSums(simplemuts)) == 0) {
          overview[overview$sampleId == sampleID, geneID] <- "amplification"
        } else if (sum(simplemuts$missense) > 0) {
          overview[overview$sampleId == sampleID, geneID] <- "amplification/missense"
        } else if (sum(simplemuts$nonsense) > 0) {
          overview[overview$sampleId == sampleID, geneID] <- "amplification/nonsense"
        } else if (sum(simplemuts$splice) > 0) {
          overview[overview$sampleId == sampleID, geneID] <- "amplification/splice"
        } else if (sum(simplemuts$inframe) > 0) {
          overview[overview$sampleId == sampleID, geneID] <- "amplification/inframe"
        } else if (sum(simplemuts$frameshift) > 0) {
          overview[overview$sampleId == sampleID, geneID] <- "amplification/frameshift"
        }
      } else if (any(tt$driver == "MUTATION")) {
        simplemuts_sum <- colSums(simplemuts)
        simplemuts_nonzero <- simplemuts_sum[simplemuts_sum > 0]
        if (length(simplemuts_nonzero) == 1) {
          overview[overview$sampleId == sampleID, geneID] <- names(simplemuts_nonzero)
        } else if (length(simplemuts_nonzero) >= 2) {
          overview[overview$sampleId == sampleID, geneID] <- paste(names(simplemuts_nonzero)[1], names(simplemuts_nonzero)[2], sep = "/")
        }
      }
    }
  }
}

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
                         

# When creating the plotting data
mat <- overview
write.table(mat , file = file.path(WGS_plot_dir, "Driveroverview_v2.txt"), sep="\t", col.names = TRUE,quote = F, row.names = FALSE)
row.names(mat) <- mat$sampleId
mat$sampleId <- NULL
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

# Remove sampleId columns and set row names
mat_filtered <- mat[, !colnames(mat) %in% c("sampleId", "sampleId.1")]
rownames(mat_filtered) <- mat$sampleId

# Calculate gene frequencies
gene_freq <- colMeans(!is.na(mat_filtered))

# Filter for genes that occur in more than 10% of samples
frequent_drivers <- names(gene_freq[gene_freq > 0.1])

# Subset the mat dataframe to include only frequent drivers
mat_filtered <- mat_filtered[, frequent_drivers]

# Order genes by frequency (high to low)
gene_order <- names(sort(colMeans(!is.na(mat_filtered)), decreasing = TRUE))
mat_filtered <- mat_filtered[, gene_order]

# Prepare data for plotting
plotting_data <- mat_filtered

# Melt the data for ggplot
m <- as.matrix(plotting_data)
dt <- reshape2::melt(m)

# Convert Var1 and Var2 to factors with correct levels
dt$Var1 <- factor(dt$Var1, levels = rownames(plotting_data))
dt$Var2 <- factor(dt$Var2, levels = colnames(plotting_data))

# Split the value column
dt <- dt %>% 
  mutate(value = as.character(value)) %>%
  separate_rows(value, sep = "/")

# Calculate positions for plotting
dt <- dt %>%
  group_by(Var1, Var2) %>%
  mutate(
    N = n(),
    shift = (seq_len(N) / N) - (1 / (2 * N)) - 0.5,
    height = 1 / N,
    test = as.numeric(Var2) + shift,
    test2 = round(test),
    height2 = 1
  ) %>%
  ungroup()

# Update gene names with percentages
gene_names <- sapply(colnames(plotting_data), function(gene) {
  freq <- mean(!is.na(plotting_data[, gene])) * 100
  sprintf("%s (%.1f%%)", gene, freq)
})

# Create the plot
plot <- ggplot() + 
  geom_tile(data = dt, aes(Var1, y=test, fill=value, height=height),
            color="grey", size=0.1) + 
  xlab('Sample ID') + 
  ylab('Gene') +
  scale_y_continuous(breaks = seq(1, ncol(plotting_data), by=1), 
                     labels = gene_names,
                     limits = c(0.5, ncol(plotting_data) + 0.5)) +
  scale_x_discrete(labels = levels(dt$Var1)) +
  scale_fill_manual(values=colorpalette_drivers, na.value = "white") +
  guides(fill=guide_legend(title="Effect")) +
  theme(axis.title.y=element_text(size=20,vjust=3),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=20),
        axis.text.x = element_text(colour="grey20",size=10,angle=90,hjust=1,vjust=0.5,face="plain"),
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
        legend.title=element_text(size=20)) +
  geom_tile(data = dt, aes(Var1, y=test2, height=height2),
            color="gray46", size=0.1, alpha=0)


plot
# Save the plot
ggsave(file=sprintf("%s/%s.pdf", WGS_plot_dir, "drivers_and_suppressors_filtered"), plot, width = 13, height = 8.5)




#seperate plot 

mat_N <- mat_filtered[grep("^N", rownames(mat_filtered)), ]
mat_P <- mat_filtered[grep("^P", rownames(mat_filtered)), ]

# Keep the same frequent drivers and order from the original analysis
frequent_drivers <- colnames(mat_filtered)
gene_order <- colnames(mat_filtered)

# Subsetting the matrices to keep only the frequent drivers
mat_N <- mat_N[, frequent_drivers]
mat_P <- mat_P[, frequent_drivers]

# Order the columns based on the gene_order
mat_N <- mat_N[, gene_order]
mat_P <- mat_P[, gene_order]

# filter niet gescreende PDOs
mat_N<- mat_N[rownames(mat_N) != "N12", ]
mat_N<- mat_N[rownames(mat_N) != "N15", ]
# Function to create a plot
create_plot <- function(data, gene_order) {
  # Prepare data for plotting
  plotting_data <- data
  
  # Melt the data for ggplot
  m <- as.matrix(plotting_data)
  dt <- reshape2::melt(m)
  
  # Convert Var1 and Var2 to factors with correct levels
  dt$Var1 <- factor(dt$Var1, levels = rownames(plotting_data))
  dt$Var2 <- factor(dt$Var2, levels = colnames(plotting_data))
  
  # Split the value column
  dt <- dt %>% 
    mutate(value = as.character(value)) %>%
    separate_rows(value, sep = "/")
  
  # Calculate positions for plotting
  dt <- dt %>%
    group_by(Var1, Var2) %>%
    mutate(
      N = n(),
      shift = (seq_len(N) / N) - (1 / (2 * N)) - 0.5,
      height = 1 / N,
      test = as.numeric(Var2) + shift,
      test2 = round(test),
      height2 = 1
    ) %>%
    ungroup()
  
  # Update gene names with percentages
  gene_names <- sapply(colnames(plotting_data), function(gene) {
    freq <- mean(!is.na(plotting_data[, gene])) * 100
    sprintf("%s (%.1f%%)", gene, freq)
  })
  
  # Create the plot
  plot <- ggplot() + 
    geom_tile(data = dt, aes(Var1, y=test, fill=value, height=height),
              color="grey", size=0.1) + 
    xlab('Samples') + 
    ylab('Gene (frequency)') +
    scale_y_continuous(breaks = seq(1, ncol(plotting_data), by=1), 
                       labels = gene_names,
                       limits = c(0.5, ncol(plotting_data) + 0.5)) +
    scale_x_discrete(labels = levels(dt$Var1)) +
    scale_fill_manual(values=colorpalette_drivers, na.value = "white") +
    guides(fill=guide_legend(title="Effect")) +
    theme(axis.title.y=element_text(size=20,vjust=3),
          axis.text.y=element_text(size=15),
          axis.title.x=element_text(size=20),
          axis.text.x = element_text(colour="grey20",size=10,angle=90,hjust=1,vjust=0.5,face="plain"),
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
          legend.title=element_text(size=20)) +
    geom_tile(data = dt, aes(Var1, y=test2, height=height2),
              color="gray46", size=0.1, alpha=0)
  
  return(plot)
}

library(patchwork)

plot_N <- create_plot(mat_N, gene_order)
plot_N <- plot_N + theme(legend.position = "none")
plot_P <- create_plot(mat_P, gene_order)

combined_plot <- plot_N + plot_P

print(combined_plot)

# Save the plot
ggsave(file=sprintf("%s/%s.pdf", WGS_plot_dir, "drivers_and_suppressors_filtered_seperate"), combined_plot, width = 15, height = 5)


