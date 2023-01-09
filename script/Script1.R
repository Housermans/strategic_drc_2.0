### SCRIPT 1: 
# Go from RAW SpectraMax files in the input folder ("1_raw_files") to a 
# file in a folder containing values on individual organoid-experiment combinations
# with drug concentrations and GR values of an individual plate. 
# The file is saved in the "2_organoid_data" folder. 
library(readxl)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(openxlsx)

rm(list=ls())

# Use code below to get the script directory if your working in Rstudio, if your 
# not working in Rstudio but from the Terminal, use alternative option commented 
# out directly below. 
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# script.dir <- dirname(sys.frame(1)$ofile) 
home_dir <- dirname(script_dir)
raw_dir <- file.path(home_dir, "1_raw_files")
organoid_dir <- file.path(home_dir, "2_organoid_data")
resource_dir <- file.path(home_dir, "resources")


