#!/usr/bin/env Rscript
library(lme4)
library(lmerTest)
library(RMINC)
library(tidyverse)

setwd("/home/m/mchakrav/lanicupo/scratch/20230728_embryo_organ_lmer")
data <- read.csv("cropped_rel_jacs.csv")
dems <- read.csv("demographics_for_analysis.csv")
sex <- read.csv("transnetyx_all_results.csv")
qc <- read.csv("squish_qc.csv")
vols <- read.csv("label_vols.csv")


sex_dems <- merge(sex, dems, by = "Pup_ID")
sex_dems_vols <- merge(sex_dems, vols, by = c("merge"))

data1 <- merge(data,sex_dems_vols, by = c("merge"))
data1 <- merge(data1, qc, by = c("merge"))
data1 <- data1 %>%
	mutate(coil_cat = case_when(
		coil %in% c(1,4,13,16) ~"corner",
		coil %in% c(2,3,5,9,8,12,14,15) ~ "edge",
		coil %in% c(6,7,10,11) ~ "center",
		TRUE ~ NA_character_))

data1$treatment <- as.factor(data1$treatment)
data1$treatment = relevel(data1$treatment, ref = "Sal")

mask = "head_only_brain_mask_new_eroded3.mnc"

mod <-
  mincLmer(file.x ~ treatment * sex + volume+ (1|litter_size) + (1|coil_cat)
           , data1
           , mask = mask
           , parallel = c("snowfall",80))

save(mod, file = "iugr_mod_cropped_lmer.Rdata")

mincWriteVolume(mod,"iugr_cropped_main_treatmentTHC.mnc",column = 'tvalue-treatmentTHC')
mincWriteVolume(mod,"iugr_cropped_main_treatmentNull.mnc",column = 'tvalue-treatmentNull')
mincWriteVolume(mod,"iugr_cropped_main_sexMale.mnc",column = 'tvalue-sexmale')
mincWriteVolume(mod,"iugr_cropped_treatmentXsex_Null.mnc",column = 'tvalue-treatmentNull:sexmale')
mincWriteVolume(mod,"iugr_cropped_treatmentXsex_THC.mnc",column = 'tvalue-treatmentTHC:sexmale')
mincWriteVolume(mod,"iugr_cropped_volume.mnc",column = 'tvalue-volume')

