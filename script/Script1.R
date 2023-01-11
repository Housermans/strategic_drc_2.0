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
org_data_dir <- file.path(home_dir, "2_organoid_data")
resource_dir <- file.path(home_dir, "resources")

# screening_df moet een bestand inlezen waar elke organoid en elke screening
Screening_df <- read_excel(file.path(resource_dir, "Screening_overview.xlsx"))

# the list_raw_files function returns a list of all raw files in a experiment 
# folder
# TODO: use the list.files directly in the function to write an organoid file, 
# skip the 'function'. 
list_raw_files <- function(exp_id) {
  list.files(file.path(raw_dir, exp_id))
}

# TODO: make a overarching function that joins ctrl file and exp file. It should 
# use the org_name variable, which is either only the organoid name (org_id) or 
# the combination of org_id and the experimental condition (exp_cond), in case
# the differentiation was made

# This function reads the ctrl_file at the location and orientation listed in 
# the "Screening_overview.xlsx"-file located in the resources folder
# the file should have the 'ctrl' phrase in the filename. 
read_ctrl_file <- function(exp_id, organoid_name) {
  # TODO: plot quality of ctrl read out? Check if shouldn't have been inverted after all?
  # TODO: read out normally if position 9+, then invert automatically unless switch?
  # TODO: generate error if there is no ctrl file present
  ctrl_file_name <- list.files(file.path(raw_dir, exp_id), pattern="ctrl")
  ctrl_file <- read_excel(file.path(raw_dir, exp_id, ctrl_file_name), range="C1:Z17")
  ctrl_df <- read_excel(file.path(resource_dir, "standard_file.xlsx"))
  ctrl_df <- ctrl_df[-1,]
  
  # The orientation should be listed in the screening_overview file is pulled here 
  org_row <- Screening_df %>%
    filter(STR_ID==exp_id & org_name==organoid_name)
  inverted <-  pull(org_row, D0_inverted)
  position <-  pull(org_row, D0_rowstart)
  # a 2x24 matrix is created from the raw file at the location that the organoid is supposed to be
  # either the first or last column (depending on orientation) is left empty for the 'fluorescence' aka empty wells. 
  # the file is then converted to match the experiment file so that they can be joined in the read_organoid function
  if (inverted==1) {
    for (i in 0:1) {
      ctrl_df <- ctrl_df %>%
        add_row(Organoid = organoid_name, condition="Fluorescence", Timepoint="D0", Value=ctrl_file[[position+(i), 24]])
      for (j in 1:23) {
        ctrl_df <- ctrl_df %>%
          add_row(Organoid = organoid_name, condition="D0_ctrl", Timepoint="D0", Value=ctrl_file[[position+(i), j]])
      }
    }
  } else {
    for (i in 0:1) {
      ctrl_df <- ctrl_df %>%
        add_row(Organoid = organoid_name, condition="Fluorescence", Timepoint="D0", Value=ctrl_file[[position+(i), 1]])
      for (j in 1:23) {
        ctrl_df <- ctrl_df %>%
          add_row(Organoid = organoid_name, condition="D0_ctrl", Timepoint="D0", Value=ctrl_file[[position+(i), 1+j]])
      }
    }
  }
  return(ctrl_df)
}

read_screen <- function(exp_id, organoid_name) {
  exp_file_name <- list.files(file.path(raw_dir, exp_id), pattern=paste0("d5_", organoid_name))
  exp_file <- read_excel(file.path(raw_dir, exp_id, exp_file_name), range="B1:Z17", .name_repair = "unique_quiet")
  org_row <- Screening_df %>%
    filter(STR_ID==exp_id & org_name==organoid_name)
  setup_path <- file.path(resource_dir, pull(org_row, screen_setup))
  if (file.exists(setup_path)) {
    exp_df <- read_excel(setup_path)
  } else {
    return(paste("ERROR: Can't find", setup_path))
  }
  exp_df$Organoid <- rep(organoid_name, 384)
  i = 1
  for (row in 1:16) {
    for (column in 2:25) {
      exp_df[i, 'Value'] <- exp_file[[row, column]]
      i <- i+1
    }
  }
  exp_df
}

read_organoid <- function(exp_id, organoid_name) {
  exp <- read_screen(exp_id, organoid_name)
  if (typeof(exp) == "character") {
    return(exp)
  }
  ctrl <- read_ctrl_file(exp_id, organoid_name)
  org_df <- rbind(ctrl, exp)
  
  fluo_reads <- org_df %>%
    filter(condition == "Fluorescence") %>%
    group_by(Timepoint) %>%
    summarize(avg = mean(Value))
  # TODO: use fluorescence values to substract this from the respective timepoints
  # then use these corrected values to summarize control values
  # then calculate appropriate GR values
  # TODO: create workaround for GR values of Tween conditions
  # TODO: calculate Z score for plate
}







