# Version control
# v14 added color scale

#library(colorspace) #diverging colors
#q6 <- diverging_hcl(10, palette = "Berlin") 
#"#E0E1E0" "#D8E7F6" "#7FBFF5" "#005682" "#241211" "#4B201D" "#7F3C37" "#F8A29E" "#F15D3F"
# Voor boxplots gebruik ik "#7FBFF5" en "#F8A29E"

#---- load libraries and data ----
# Load dplyr and readxl packages for data manipulation and reading Excel files
library(tidyverse)
library(gridExtra)
library(grid)
library(gtable)

rm(list=ls()) #Clear existing objects in the workspace

# Define necessary directories for accessing and saving files
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
home_dir <- dirname(dirname(script_dir))
WGS_dir <- file.path(home_dir, "WGS")
WGS_data_dir <- file.path(WGS_dir, "Data_HMF")
WGS_plot_dir <- file.path(home_dir, "Analyse 2.0/7_WGS")

#------import WGS data Hartwig (validation set)-----
driving <- read.csv(file.path(WGS_plot_dir, "Driveroverview_raw_v2.txt"), sep="\t")

str(driving)

# Create summary of mutated genes per sample
mutation_summary <- driving %>%
  # Filter for mutations only
  filter(driver == "MUTATION") %>%
  # Group by sample and create gene list
  group_by(sampleId) %>%
  summarise(
    mutated_genes = paste(sort(unique(gene)), collapse = ", "),
    mutation_count = n()
  ) %>%
  # Sort by sampleId
  arrange(sampleId) %>%
  filter(startsWith(sampleId, "N")) %>%
  # Rename columns for display
  rename(
    "Sample ID" = sampleId,
    "Mutated Genes" = mutated_genes,
    "Total Mutations" = mutation_count
  ) 
  

# Function to create table with grid
create_table_plot <- function(data) {
  # Convert data frame to table grob
  table <- tableGrob(data, rows = NULL, theme = ttheme_minimal(
    core = list(fg_params = list(hjust = 0, x = 0.1),
                padding = unit(c(4, 4), "mm")),
    colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.1))
  ))
  
  # Adjust column widths
  table$widths <- unit(c(0.15, 0.7, 0.15), "npc")
  
  return(table)
}

# Create and save plot
pdf(file.path(WGS_plot_dir, "mutation_table_chemonaive.pdf"), width = 11, height = 0.3 * nrow(mutation_summary))
grid.newpage()
grid.draw(create_table_plot(mutation_summary))
dev.off()

# Print preview of the table
print(mutation_summary)
