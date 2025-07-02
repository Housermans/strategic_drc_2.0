### SCRIPT 1: 
# Go from RAW SpectraMax files in the input folder ("1_raw_files") to a 
# file in a folder containing values on individual organoid-experiment combinations
# with drug concentrations and GR values of an individual plate. 
# The file is saved in the "2_organoid_data" folder. 
library(readxl)
library(ggplot2)
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

# This function reads the ctrl_file at the location and orientation listed in 
# the "Screening_overview.xlsx"-file located in the resources folder
# the file should have the 'ctrl' phrase in the filename.
# This function is called in the read_organoid function to bind it together
# with the corresponding screening. 
read_ctrl_file <- function(exp_id, organoid_name) {
  # TODO: plot quality of ctrl read out? Check if shouldn't have been inverted after all?
  # TODO: read out normally if position 9+, then invert automatically unless switch?
  # TODO: generate error if there is no ctrl file present
  ctrl_file_name <- list.files(file.path(raw_dir, exp_id), pattern="ctrl")
  ctrl_file <- read_excel(file.path(raw_dir, exp_id, ctrl_file_name), range="C1:Z17")
  org_row <- Screening_df %>%
    filter(STR_ID==exp_id & org_name==organoid_name)
  setup_path <- file.path(resource_dir, pull(org_row, screen_setup))
  
  if (file.exists(setup_path)) {
    ctrl_df <- read_excel(setup_path)
  } else {
    return(paste("ERROR: Can't find", setup_path))
  }
  
  ctrl_df <- ctrl_df[1,]
  ctrl_df <- ctrl_df %>%
    mutate(Organoid = as.character(Organoid),
           Value = as.numeric(Value),
           value_corr = as.numeric(value_corr), 
           GR = as.numeric(GR))
  ctrl_df <- ctrl_df[-1,] 
  
  # The orientation should be listed in the screening_overview file is pulled here 
  org_row <- Screening_df %>%
    filter(STR_ID==exp_id & org_name==organoid_name)
  # inverted <-  pull(org_row, D0_inverted)
  position <-  pull(org_row, D0_rowstart)
  if (position > 8) {
    inverted <- 1
  } else {
    inverted <- 0
  }
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

# The read_screen function returns a datafame with the values from the screening
# plate pasted in the correct lane of the drugscreen setup. 
# This works if the screeningsetup is emptied out and the names of the condition
# have been correctly changed to something interpretable (i.e. binimetinib_lapatinib
# instead of '2 fluids' where the binimetinib and lapatinib column contain dose
# values.)
# Returns an error if the setup file can't be found. 
# This function is used in the read_organoid function. 
read_screen <- function(exp_id, organoid_name) {
  org_row <- Screening_df %>%
    filter(STR_ID==exp_id & org_name==organoid_name)
  setup_path <- file.path(resource_dir, pull(org_row, screen_setup))
  
  exp_file_name <- list.files(file.path(raw_dir, exp_id), pattern=paste0("d5_", organoid_name, ".xlsx"))
  if (length(exp_file_name) == 0) {
    return(paste("ERROR: Combination of", exp_id, "and", organoid_name, "does not exist in this folder!"))
  }
  if (startsWith(pull(org_row, screen_setup), "G12CScreen")) {
    halfscreen <- TRUE
    print(paste("Reading out", exp_id, organoid_name, "as a Halfscreen experiment"))
    exp_file <- read_excel(file.path(raw_dir, exp_id, exp_file_name), range="B1:Z9", .name_repair = "unique_quiet")
  } else {
    # Normal full 384 well plate
    exp_file <- read_excel(file.path(raw_dir, exp_id, exp_file_name), range="B1:Z17", .name_repair = "unique_quiet")
  }
  
  if (file.exists(setup_path)) {
    exp_df <- read_excel(setup_path)
  } else {
    return(paste("ERROR: Can't find", setup_path))
  }
  if (startsWith(pull(org_row, screen_setup), "G12Cscreen")) {
    # reading out a half screen experiment
    exp_df$Organoid <- rep(organoid_name, 192)
    i = 1
    for (row in 1:8) {
      for (column in 2:25) {
        exp_df[i, 'Value'] <- exp_file[[row, column]]
        i <- i+1
      }
    }
  } else {
    # reading out a normal screen (full plate) experiment
    exp_df$Organoid <- rep(organoid_name, 384)
    i = 1
    for (row in 1:16) {
      for (column in 2:25) {
        exp_df[i, 'Value'] <- exp_file[[row, column]]
        i <- i+1
      }
    }
  }
  exp_df
}

get_fluo <- function(org_df) {
  fluo_reads <- org_df %>%
    filter(condition == "Fluorescence") %>%
    group_by(Timepoint) %>%
    summarize(avg = mean(Value))
  fluo_vec <- pull(fluo_reads, avg)
  unlist(fluo_vec)
}

get_vehicle <- function(org_df) {
  vehicle_reads <- c(0,0,0)
  vehicle_reads[1] <- org_df %>%
    filter(condition == "Tween") %>%
    summarize(avg = mean(value_corr))
  vehicle_reads[2] <- org_df %>%
    filter(condition == "DMSO_1") %>%
    summarize(avg = mean(value_corr))
  vehicle_reads[3] <- org_df %>%
    filter(condition == "D0_ctrl") %>%
    summarize(avg = mean(value_corr))
  unlist(vehicle_reads)
}
  

apply_fluo <- function(Value, Timepoint, fluo_vec) {
  if (Timepoint == "D0") {
    value_corr2 <- Value - fluo_vec[1]
  } else {
    value_corr2 <- Value - fluo_vec[2]
  }
  if (value_corr2 < 0) {
    value_corr2 <- 0
  }
  value_corr2
}

apply_GR <- function(value_corr, DMSO_status, Tween, DMSO, Day0) {
  if (is.na(DMSO_status)) {
    return(NA)
  }
  if (value_corr <= 0) {
    GR <- -1
  } else if (DMSO_status == "Y") {
    GR <- 2^(log2(value_corr/Day0)/log2(DMSO/Day0))-1
  } else if (DMSO_status == "N") {
    GR <- 2^(log2(value_corr/Day0)/log2(Tween/Day0))-1
  }
}

# The read_organoid function returns a collated dataframe of the control and
# screening plate, unless the read_screen function returns an error message
# this usually means that the drug screen setup file wasn't found.
# it will then:
#   1. Subtract the fluorescence values of D0 and D5 from the ctrl and screening
#   values respectively
#   2. Use the negative control values (DMSO_1 and Tween) to calculate GR scores
#   3. Save the excel file to a file containing screening date, STR-ID, and
#   organoid_name variables. 
# TODO: use fluorescence values to substract this from the respective timepoints
# then use these corrected values to summarize control values
# then calculate appropriate GR values
# TODO: create workaround for GR values of Tween conditions
# TODO: calculate Z score for plate
read_organoid <- function(exp_id, organoid_name) {
  exp <- read_screen(exp_id, organoid_name)
  if (typeof(exp) == "character") {
    return(exp)
  }
  ctrl <- read_ctrl_file(exp_id, organoid_name)
  org_df <- rbind(ctrl, exp)
  
  org_df <- org_df %>%
    mutate(DMSO_status = case_when(
      DMSO_pct > 0.0001 ~ "Y", 
      DMSO_pct < 0.0001 ~ "N" 
    ))
  
  # distills fluorescence control values from the experiment and then distracts these from the 
  # appropriate value points in the experiment. If a value goes below 0, it is then set to 0
  # as negative values have no meaning here. 
  fluo_vec <- get_fluo(org_df)
  for (i in 1:nrow(org_df)) {
    # print(paste0("analyzing row ",i))
    org_df[i,]$value_corr <- apply_fluo(org_df[i,]$Value, org_df[i,]$Timepoint, fluo_vec)
    
  }
  
  # Summarizes the vehicle numbers and applies the appropriate vehicle control to each condition
  # based on the percentage of DMSO in the condition. If the percentage is above 0,1%, the DMSO status = "Y"
  # below, the DMSO status = "N".
  Tween <- unlist(org_df %>%
                    filter(condition == "Tween") %>%
                    summarize(avg = mean(value_corr))) 
  DMSO <- unlist(org_df %>%
                   filter(condition == "DMSO_1") %>%
                   summarize(avg = mean(value_corr)))
  Day0 <-unlist(org_df %>%
                  filter(condition == "D0_ctrl") %>%
                  summarize(avg = mean(value_corr)))
  
  for (i in 1:nrow(org_df)) {
    org_df[i,]$GR <- apply_GR(org_df[i,]$value_corr, org_df[i,]$DMSO_status, Tween, DMSO, Day0)
  }
  org_row <- Screening_df %>%
     filter(STR_ID==exp_id & org_name==organoid_name)
  date_long <-  pull(org_row, Date_D0)
  date_short <- substr(date_long, 1, 10)
  write.xlsx(org_df, file.path(org_data_dir, paste0(paste(date_short, exp_id, organoid_name, sep="_"), ".xlsx")))

  org_df
  
}

read_experiment <- function(exp_id) {
  experiment_orgs <-Screening_df %>%
    filter(STR_ID==exp_id)
  experiment_orgs <- unique(experiment_orgs$org_name)
  for (i in 1:length(experiment_orgs)) {
    print(paste("Reading organoid:", experiment_orgs[i]))
    org <- read_organoid(exp_id, experiment_orgs[i])
    if (typeof(org) == "character") {
      print(org)
    }
  }
}

 # read_screen(exp_id = "STR-D3", organoid_name = "OPT0024_Demi")
 # read_experiment("STR-D3")
 read_experiment("STR19")
 read_experiment("STRG12C1")
 # all_exp <- unique(Screening_df$STR_ID)
 # lapply(all_exp, read_experiment)
