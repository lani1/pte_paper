# Author: Lani Cupo Last updated: 20210712
# This script analyzes the weights from thc exposure during pregnancy to see if it has an impact on weight gain in pups.

library(lme4)
library(lmerTest)
library(dplyr)
library(viridis)
library(ggplot2)
library(ggpubr)
library(effects)
library(splines)
library(tidyverse)

#set the directory where the files are
setwd("/data/chamal/projects/lani/neonate")

# load data
data <- read.csv("demographics.csv", encoding = "UTF-8", dec = ".", header = T)
data$scanned <-relevel(data$scanned, ref = "no")
scan_qc <- read.csv("./qc/20211018_qc.csv")
data_qc <- merge(data,scan_qc, by= c("ID","age"))

#Set theme for ggplots------------------------------
myTheme <- theme_classic()+ theme(axis.text=element_text(size=16), axis.title=element_text(size=16), 
                                  legend.text=element_text(size=16), legend.key.size = unit(1, 'cm'),
                                  legend.title=element_text(size=18), plot.title = element_text(size = 16))

#mark those that were "incomplete" as "scanned" 
data <- data %>%
  mutate(total_scanned = ifelse(scanned == "no", "no", "yes"))
data$total_scanned <- as.factor(data$total_scanned)

data <- subset(data, data$MnCl2 == "manganese")

#calculate total, complete sample size
unique_data <- data %>%
  distinct(data$ID, .keep_all = T)
unique_data_qc <- data_qc %>%
  distinct(data_qc$ID, .keep_all = T)

incomplete_unique <- subset(unique_data, unique_data$scanned %in% c("yes", "incomplete"))
complete_unique <- subset(unique_data, unique_data$scanned %in% c("yes"))

#run model
data <- subset(data, !is.na(data$weight))
mod <- lmer(weight ~ condition*age *sex + total_scanned + (1|ID) + (1|mom_id), data)
summary(mod)
tab_model(mod, file = "../pte_paper/stats_tables/neonate_weight_lmer.doc")
  
data$weight_resid <- resid(mod)

write.csv(x = data, file = "data_with_weight_resids.csv")

mod_scan <- lmer(weight ~ condition*age + total_scanned + sex + (1|ID) + (1|mom_id), data)
summary(mod_scan)
(res.table <- as.data.frame(coef(summary(mod_scan))))
xtable(res.table, type = "latex")

plot <- ggplot(data, aes(x=condition, y=weight, color=condition))+
  geom_boxplot()+
  geom_jitter()+
  labs(x="Treatment", y="Weight", title="Daily Weights") +
  #scale_color_manual('Sex', values = c('F' = '#da91ac', 'M' = '#3c5839'))+
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title=element_blank()) +
  facet_wrap(~ sex) +
  theme_classic()

fitdata_mod = as.data.frame(Effect(c("age", "condition"),mod, xlevels=list(Day=seq(0,15,1))))
fitdata_mod_scan = as.data.frame(Effect(c("age", "total_scanned"),mod, xlevels=list(Day=seq(0,15,1))))

condition <- ggplot(data = data, aes(y=weight,x=age, color=condition, group=condition)) +
  geom_line(data = fitdata_mod, aes(y=fit),size=1, alpha=0.8) +
  geom_ribbon(data = fitdata_mod, aes(y=fit, ymin=lower, ymax=upper), alpha=0.4, color = NA) +
  geom_jitter(alpha=0.7, size = 3, aes(shape = factor(sex))) +
  scale_color_manual('Condition', values = c('sal' = '#420085', 'thc' = '#018700'))+
  theme(axis.text=element_text(size=14), 
        legend.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5,size = 16, face = "bold"))+
  xlab("Age (days)") +
  ylab("Neonate Weight (g)") + 
  ggtitle("Pup Weights by Condition") +
  scale_x_continuous(breaks = c(seq(3, 14, 2)))+
  myTheme
condition

scan <- ggplot(data = data, aes(y=weight,x=age, color=total_scanned, group=total_scanned)) +
  geom_point(alpha=0.7, size = 2) +
  geom_line(data = fitdata_mod_scan, aes(y=fit),size=1, alpha=0.8) +
  geom_ribbon(data = fitdata_mod_scan, aes(y=fit, ymin=lower, ymax=upper), alpha=0.4) +
  geom_jitter(alpha=0.7, size = 3, aes(shape = factor(sex))) +
  xlab("Age (days)") +
  ylab("Weight (g)") + 
  ggtitle("Pup Weights if Scanned") +
  scale_x_continuous(breaks = c(seq(3, 14, 2))) +
  labs(colour = "Scanned?")+
  theme(axis.text=element_text(size=26), axis.title=element_text(size=26,face="bold"), legend.text=element_text(size=26),
        legend.title=element_text(size=28), plot.title = element_text(size = 26, face = "bold"))+
  scale_x_continuous(breaks = c(seq(3, 14, 2))) +
  theme_classic()+ theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=16),
                         legend.title=element_text(size=18), plot.title = element_text(size = 16, face = "bold"))
scan

png("visuals/neonate_weight.png", units="in", width=7, height=5, res=300)
ggarrange(condition,scan)
dev.off()

