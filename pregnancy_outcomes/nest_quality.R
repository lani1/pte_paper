#Author: Lani Cupo
#Last updated: 20220617
# This script was written to analyze pilot data for the nest quality from my moms

#Load Libraries
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(ggplot2)
  library(ggpubr)
  library(ggeffects)
  library(grid)
  library(tidyverse)
  library(effects)
  library(vcd)
})

setwd("/data/chamal/projects/lani/adult_thc//behavior/raw_data/")

#Load data
data <- read.csv("nest_quality.csv", header = T)

data$Session = relevel(data$Session, ref = "Pre-Injection")

nest_test <- lmer(nest_qual ~ treatment*Session + (1|Mom), data = data)
summary(nest_test)

data <- data %>% 
  group_by(Session) %>% 
  mutate(average = mean(nest_qual, na.rm = T), .keep = "all") %>%
  ungroup()

fitdata = as.data.frame(Effect(c("treatment", "Session"),nest_test, xlevels=list(Session=seq(1,2,1))))
nest_quality_plot <-  ggplot(data = data, aes(y=nest_qual,color=treatment,x=Session)) +
  geom_boxplot(outlier.shape = NA)+
  #geom_line(position = position_dodge(width = 0.75), aes(group = Mom), alpha = 0.3) +
  geom_jitter(position=position_dodge(width=0.75), aes(group=treatment))+
  xlab("Session") +
  ylab(bquote('Nest Quality')) + 
  labs(fill = "Treatment", color="Treatment", title = "Nest Quality")+
  scale_color_manual('Treatment', values = c('SAL' = '#420085', 'THC' = '#018700'))+
  scale_size_manual(values=c(2,2)) +
  theme_classic() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=16),
        legend.title=element_text(size=18), plot.title = element_text(size = 16, face = "bold"))
nest_quality_plot

#eof

