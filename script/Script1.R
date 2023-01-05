### SCRIPT 1: Go from RAW SpectraMax files in the input folder to a file containing values on drug concentrations and GR values of an individual plate. The file is saved in the output folder. 
library(readxl)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(openxlsx)

rm(list=ls())

# Onderstaande variabele script_dir is een manier om de locatie van het script in jouw context te krijgen. 
# werkt alleen in RStudio. Voor alternatief vanuit de terminal (https://stackoverflow.com/questions/47044068/get-the-path-of-current-script)
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)


list.files(sources)
