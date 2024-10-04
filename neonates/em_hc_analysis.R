#author: Lani Cupo
#written: 20230530
#this script will first run lmers to analyze my cell data, then shift function to compare groups.
suppressPackageStartupMessages({
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(ggpubr)
  library(rogme) # for shift function
  library(dplyr)
  library(stringr) #to extract pup ID
  library(car) #for levene test
})

#set working directory

#Set up the data --------------------------------------------------------

setwd("/data/chamal/projects/lani/neonate/electron_microscopy/")

data <- read.csv("neo_cell_counts_mm2.csv", header = T)
dems <- read.csv("../demographics.csv")
one_dem <- dems %>%
  distinct(dems$ID, .keep_all=TRUE)


data_dems <- merge(data, one_dem, by = c("ID"))
one_dem2 <- data_dems %>%
  distinct(data_dems$ID, .keep_all=T)
#run models------------------------------------------------------------------
neuron_mod <- lmer(total_dark_neuron ~ condition + sex  + (1|ID), data= data_dems )
summary(neuron_mod)

dglia_mod <- lmer(total_dark_micro ~ condition + sex  + (1|ID), data= data_dems )
summary(dglia_mod)

mglia_mod <- lmer(total_microglia ~ condition + sex  + (1|ID), data = data_dems)
summary(mglia_mod)

apop_mod <- lmer(total_apop ~ condition + sex  + (1|ID), data= data_dems)
summary(apop_mod)

#div_mod <- lmer(total_dividing ~ condition + sex  + (1|ID), data= data_dems)
#summary(div_mod)

#------------------------------------------------------------------------
data_dems <- subset(data_dems, !is.na(data_dems$total_dark_micro))

filtered.data.tibble <- as_tibble(data_dems)
filtered.data.tibble$condition <- as.factor(filtered.data.tibble$condition)
sf <- shifthd_pbci(data = filtered.data.tibble, formula = total_apop ~ condition)

plot_sf(sf)
p = plot_sf(sf)
p = add_sf_lab(p, sf=sf)

# Plots------------------------------------------------------------------
tot_micro <-ggplot(data = data_dems, aes(x=condition, y=total_microglia, color=condition, group=condition)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3, height = 0))+
  #geom_boxplot(data = fitdata_jacob, aes(y=fit),size=1, alpha=0.8) +
  #geom_ribbon(data = fitdata_jacob, aes(y=fit, ymin=lower, ymax=upper, color=NULL), alpha=0.1) +
  scale_color_manual('Treatment', values = c('sal' = '#420085', 'thc' = '#018700'))+
  labs(x = "Treatment", y= bquote("Cell/"~mm^2), title = "Microglia in the Neonate HC")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=20), legend.title=element_text(size=18)) +
  theme_classic()
tot_micro

apop <-ggplot(data = data_dems, aes(x=condition, y=total_apop, color=condition, group=condition)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3, height = 0))+
  #geom_boxplot(data = fitdata_jacob, aes(y=fit),size=1, alpha=0.8) +
  #geom_ribbon(data = fitdata_jacob, aes(y=fit, ymin=lower, ymax=upper, color=NULL), alpha=0.1) +
  scale_color_manual('Treatment', values = c('sal' = '#420085', 'thc' = '#018700'))+
  labs(x = "Treatment", y= bquote("Cell/"~mm^2), title = "Apoptosis in the Neonate Hippocampus")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=20), legend.title=element_text(size=18)) +
  theme_classic()
apop

d_micro <-ggplot(data = data_dems, aes(x=condition, y=total_dark_micro, color=condition, group=condition)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3, height = 0))+
  #geom_boxplot(data = fitdata_jacob, aes(y=fit),size=1, alpha=0.8) +
  #geom_ribbon(data = fitdata_jacob, aes(y=fit, ymin=lower, ymax=upper, color=NULL), alpha=0.1) +
  scale_color_manual('Treatment', values = c('sal' = '#420085', 'thc' = '#018700'))+
  labs(x = "Treatment", y= bquote("Cell/"~mm^2), title = "Dark Microglia in the Neonate Hippocampus")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=20), legend.title=element_text(size=18)) +
  theme_classic()
d_micro

d_neur <-ggplot(data = data_dems, aes(x=condition, y=total_dark_neuron, color=condition, group=condition)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3, height = 0))+
  #geom_boxplot(data = fitdata_jacob, aes(y=fit),size=1, alpha=0.8) +
  #geom_ribbon(data = fitdata_jacob, aes(y=fit, ymin=lower, ymax=upper, color=NULL), alpha=0.1) +
  scale_color_manual('Treatment', values = c('sal' = '#420085', 'thc' = '#018700'))+
  labs(x = "Treatment", y= bquote("Cell/"~mm^2), title = "Dark Neurons in the Neonate Hippocampus")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=20), legend.title=element_text(size=18)) +
  theme_classic()
d_neur

ggarrange(tot_micro, apop, d_micro, d_neur)

png("../../pte_paper/figures/neonate_microglia.png", units="in", width=3, height=3, res=300)
tot_micro
dev.off()
