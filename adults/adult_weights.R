# Author: Lani Cupo Last updated: 20220923
# This script analyzes the weights from thc exposure during pregnancy to see if it has an impact on weight gain in adults. 

library(lme4)
library(lmerTest)
library(dplyr)
library(viridis)
library(ggplot2)
library(ggpubr)
library(effects)
library(splines)
library(tidyverse)
library(sjPlot)

#set the directory where the files are
setwd("/data/chamal/projects/lani/")

#Set theme for ggplots------------------------------

myTheme <- theme_classic()+ theme(axis.text=element_text(size=16), axis.title=element_text(size=16), 
                                  legend.text=element_text(size=16), legend.key.size = unit(1, 'cm'),
                                  legend.title=element_text(size=18), plot.title = element_text(size = 16))

#data <- read.csv("demographics.csv", encoding = "UTF-8", dec = ".", header = T)
data <- read.csv("adult_thc/demographics.csv")
data <- subset(data, data$scanned == "yes")
unique <- data %>%
  distinct(data$ID, .keep_all = T)

data <- subset(data,!is.na(data$age))
mod <- lmer(weight ~ treatment*poly(age,2)* sex +(1|ID) + (1|mom_id), data)
summary(mod)
tab_model(mod, file = "pte_paper/stats_tables/adult_weight_lmer.doc")


fitdata_mod = as.data.frame(Effect(c("sex", "treatment","age"),mod, xlevels=list(age=seq(25,90,1))))

weight_condition <- ggplot(data = data, aes(y=weight,x=age, color=treatment, group=treatment)) +
  geom_line(data = fitdata_mod, aes(y=fit),size=1, alpha=0.8) +
  geom_ribbon(data = fitdata_mod, aes(y=fit, ymin=lower, ymax=upper), alpha=0.4, color = NA) +
  geom_jitter(alpha=0.7, size = 3) +
  scale_color_manual('Treatment', values = c('SAL' = '#420085', 'THC' = '#018700'))+
  xlab("Age (days)") +
  ylab("Adult Weight (g)") + 
  ggtitle("Adult Weights") +
  scale_x_continuous(breaks = c(25,35,60,90))+
  facet_wrap(~ sex) +
  myTheme +
  theme(strip.text = element_text(
    size = 20))
weight_condition
data <- subset(data, !is.na(data$weight ))
data$weight_resid <- resid(mod)

write.csv(x = data, file = "adult_thc/dems_with_weight_resids.csv")



