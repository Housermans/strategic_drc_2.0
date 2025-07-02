#load packages
library(tidyverse)
library(boot) 
library(table1, quietly=TRUE)
library(rmarkdown)

rm(list=ls())

#TODO: OPTIC RASTRIC cohort


#Load files
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# db_clinical <- readRDS(file = file.path(script_dir,"db_clinic_only_table1.Rds"))
db_optic_original <- readRDS(file= file.path(script_dir, "OPTIC_clinical_data.rds"))
db_rastric_original <- read_delim(file.path(script_dir, "Gegevens_RASTRIC.csv"), delim=",")
db_rastric_biopsyloc <- read_delim(file.path(script_dir, "RASTRIC_biopsy_locations.csv"), delim=";")
db_strategic_original <- read_delim(file.path(script_dir, "STRATEGIC_overview2.csv"), delim=";")


db_strategic <- db_strategic_original %>% 
  filter(STRATEGIC_analysis == "yes") %>%
  mutate(organoid_name = str_trim(`Organoid line`), 
         MAPK_status = str_trim(`MAPK status`),
         Palliative_duration = as.numeric(Palliative_duration)) %>%
  mutate(
    Palliative_duration = ifelse(is.na(Palliative_duration), 0, Palliative_duration),
    Palliative_duration = case_when(
      Palliative_duration > 100 ~ floor(Palliative_duration / 30.5),
      TRUE ~ Palliative_duration
      ),
  ) %>%
  dplyr::select(organoid_name, `Sample type`, MAPK_status, Adjuvant_received, 
                Adjuvant_5FU, Adjuvant_Ox, Palliative_duration, STRATEGIC_analysis)
  
db_optic <- db_optic_original %>% mutate(
  creation_date = dmy_hms(Participant.Creation.Date),
  creation_year = year(creation_date),
  Age = creation_year - Birth_year, 
  Sex = Gender_Pt, 
  organoid_name = ID_OPTIC,
  origin_location = Biopsy_tissue
) %>%
  dplyr::select(organoid_name, Age, Sex, Sidedness, origin_location)

db_rastric <- db_rastric_original %>%
  mutate(
    organoid_name = str_replace(SUBJID, "RAS0", "RAS"),
    regisdat = dmy(regisdat),
    birthday= dmy(dob), 
    Age = as.integer(floor(interval(birthday, regisdat) / years(1))), 
    # PS = str_sub(whops_bl, -1, -1), 
    Sidedness = case_when(
      initmhloc %in% c("sigmoid", "descending colon") ~ "Left-sided (splenic flexure-sigmoid)",
      initmhloc == "cecum" ~ "Right-sided (coecum-transverse colon)",
      initmhloc %in% c("rectum", "Rectosigmoid") ~ "Rectum (rectosigmoid/rectal)",
      TRUE ~ NA_character_
    ),
    # RASMut = case_when(
    #   str_detect(krasmut, "G12") ~ "KRAS codon 12 mutation", 
    #   str_detect(krasmut, "Gly12") ~ "KRAS codon 12 mutation",
    #   str_detect(krasmut, "c.34G") ~ "KRAS codon 12 mutation", 
    #   str_detect(krasmut, "c.35G") ~ "KRAS codon 12 mutation", 
    #   str_detect(krasmut, "exon 2") ~ "KRAS codon 12 mutation", 
    #   str_detect(krasmut, "G13") ~ "KRAS codon 13 mutation", 
    #   str_detect(krasmut, "c.38G") ~ "KRAS codon 13 mutation",
    #   str_detect(krasmut, "436G") ~ "KRAS non-exon 2 mutation",
    #   str_detect(krasmut, "183A") ~ "KRAS non-exon 2 mutation",
    #   str_detect(krasmut, "873G") ~ "KRAS wildtype amplification",
    #   str_detect(nrasmut, "61") ~ "NRAS mutation"
    # ), 
    # nmeta_sites = case_when(
    #   nmeta > 3 ~ "≥ 4",
    #   .default = as.character(nmeta)
    # ),
    # nmeta_sites = factor(nmeta_sites, levels=c("1", "2", "3", "≥ 4")),
    Sex = sex
  ) %>% 
  left_join(db_rastric_biopsyloc, c("organoid_name" = "Study Number")) %>%
  mutate(
    origin_location = str_replace(Location, "metastasis", ""),
    origin_location = str_replace(origin_location, "l biopsy", ""),
    origin_location = case_when(
      str_detect(origin_location,"(skin only)") ~ NA_character_, 
      .default = origin_location
    )) %>%
  mutate(
    origin_location = str_replace(origin_location, "\\(.*?\\)", ""),
  ) %>%
  dplyr::select(organoid_name, Age, Sex, Sidedness, origin_location)

# Combine all data frames
db_combined <- db_strategic %>%
  left_join(db_optic, by = "organoid_name") %>%
  left_join(db_rastric, by = "organoid_name") %>%
  mutate(
    Age = coalesce(Age.x, Age.y),
    Sex = coalesce(Sex.x, Sex.y),
    Sidedness = coalesce(Sidedness.x, Sidedness.y),
    origin_location = coalesce(origin_location.x, origin_location.y),
    `Sample type` = str_to_title(`Sample type`)
  ) %>%
  dplyr::select(-Age.x, -Age.y, -Sex.x, -Sex.y, -Sidedness.x, -Sidedness.y, -origin_location.x, -origin_location.y) %>%
  mutate(
    Sidedness = ifelse(Sidedness == "", "Unknown", Sidedness), 
    MAPK_status = case_when(
      MAPK_status == "KRAS" ~ "RAS-mutant", 
      MAPK_status == "NRAS" ~ "RAS-mutant", 
      MAPK_status == "WT" ~ "RAS/BRAF-wildtype",
      MAPK_status == "WT KRAS amplification" ~ "KRAS wildtype amplification"),
    origin_location = str_trim(origin_location)
  ) %>% 
  filter(STRATEGIC_analysis == "yes") 

# # Identify and print non-matching names
# non_matching_strategic_optic <- db_strategic %>%
#   anti_join(db_optic, by = "organoid_name")
# 
# non_matching_strategic_rastric <- db_strategic %>%
#   anti_join(db_rastric, by = "organoid_name")
# 
# cat("Non-matching names between strategic and optic:\n")
# print(non_matching_strategic_optic)
# 
# cat("Non-matching names between strategic and rastric:\n")
# print(non_matching_strategic_rastric)

label(db_combined$Age) <- "Age" 
label(db_combined$Sex) <- "Sex" 
label(db_combined$Sidedness) <- "Primary tumour location"
label(db_combined$Adjuvant_received) <- "Adjuvant chemotherapy"
label(db_combined$MAPK_status) <- "Mutational status"
label(db_combined$Palliative_duration) <- "Months of palliative chemotherapy"
label(db_combined$origin_location) <- "Location of metastasis biopsy"


my.render.cont <- function(x) {
  with(stats.default(x), 
       c("",
         
         
         "Median (Min, Max)" = sprintf("%s (%s, %s)",
                                       round_pad(MEDIAN, 0), 
                                       round_pad(MIN, 0), 
                                       round_pad(MAX, 0)))
  )
}

table_combined <- table1(~ (Age) + (Sex)  + (Sidedness) + (MAPK_status) + (origin_location) + (Adjuvant_received) + (Palliative_duration) | `Sample type`, 
                         data=db_combined, 
                         render.continuous=my.render.cont)

table_combined

write.table(table_combined, file.path(script_dir,"table_1_v3.csv"), col.names = T, row.names=F, append= T)


rmarkdown::render(file.path(script_dir, "table1_v2.Rmd"))

#---- Biobank voorbeeld ----
# db_clinical$Prim_tum_sidedn [db_clinical$Prim_tum_sidedn == "Multiple primary tumours (with different sidedness)"] <- "Right-sided (coecum-transverse colon)"
# 
# #remove patient with 2 organoids
# baseline <- db_clinical[-c(11), ] #HUB-02-C2-154IV
# baseline <- baseline[-c(24), ] #M20-00016II
# 
# #liveronly
# baseline$metast_number1 <- ifelse(baseline$`met_location#Lung`!=0 | baseline$`met_location#Peritoneal`!=0 | baseline$`met_location#Other_abdominal_sites_than_liverperitoneal`!=0, ">1", "1")
# baseline$liveronly <- ifelse(baseline$metast_number1 == ">1","No","Yes")
#  
# #KRAS_BRAF
# baseline$mut_status2 <- ifelse(baseline$mut_status == "KRAS-m / BRAF-wt", "RAS-mutant",
#                                ifelse(baseline$mut_status == "KRAS-wt / BRAF-m","BRAF-mutant","RAS/BRAF-wildtype"))
# 
# #birth year to age
# baseline$age <- 2021-baseline$pat_birth_year 
# baseline$agecut <- cut(baseline$age, breaks=c(0, 56, 66, 76, 100), right = FALSE)
# baseline$agecutfact<- factor(baseline$agecut,
#                                   levels=c("[0,56)","[56,66)","[66,76)","[76,100)"),
#                                   labels=c("<56 years", "56-65 years","66-75 years",">75 years"))
# #combine pdo prim + meta
# baseline$meta_prim_capital <- factor(baseline$org_prim_met,
#                                           levels=c("metastatic","primary"), 
#                                           labels=c("Metastatic","Primary"))
# 
# baseline$prior <- ifelse(baseline$chemoexposed == "yes" & baseline$adjuv_given =="No", "Yes", "No")
# 
# #label and units
# label(baseline$meta_prim_capital) <- "PDO origin"
# label(baseline$pat_sex) <- "Sex"
# label(baseline$age) <- "Age"
# label(baseline$Prim_tum_sidedn) <- "Primary tumour location"
# label(baseline$prim_tum_res) <- "Primary tumour resection"
# label(baseline$adjuv_given) <- "Adjuvant chemotherapy"
# label(baseline$metast_number1) <- "Number of metastatic sites"
# label(baseline$liveronly) <- "Liver-only disease"
# label(baseline$mut_status2) <- "Mutational status"
# #label(baseline$PFS1) <- "Progression-free survival of 1st line treatment"
# label(baseline$prior) <- "Palliative treatment line(s) prior to PDO establishment"
# label(baseline$agecutfact) <- "Age"
# 
# my.render.cont <- function(x) {
#   with(stats.default(x), 
#        c("",
#          
#          
#          "Median (Min, Max)" = sprintf("%s (%s, %s)",
#                                        round_pad(MEDIAN, 0), 
#                                        round_pad(MIN, 0), 
#                                        round_pad(MAX, 0)))
#   )
# }
# 
# table_baseline <- table1(~ (age) + (Prim_tum_sidedn) + (prim_tum_res) + 
#                            (adjuv_given) + (liveronly) + (mut_status2) +(prior)+ (meta_prim_capital), 
#                          data=baseline, 
#                          render.continuous=my.render.cont)
# 
# table_baseline
# setwd("D:/SURFdrive/Lsmabers/PhD/Biobank/Scripts biobank Emerens/Analysis Screen 2021_2022_v2/Tables")
# write.table(table_baseline, "table_baseline.csv", col.names = T, row.names=F, append= T)
# 
# 
# 
# 
