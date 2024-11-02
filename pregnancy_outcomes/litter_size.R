#Author: Lani Cupo
#Written: 20230713
# Analyze weight data from all dams togethers

library(lme4)
library(lmerTest)
library(dplyr)
library(viridis)
library(ggplot2)
library(ggpubr)
library(effects)
library(splines)
library(scales)

#set up dataframe
em_sex <- read.csv("/data/chamal/projects/lani/embryos/transnetyx/transnetyx_all_results.csv")
em_dems <- read.csv("/data/chamal/projects/lani/embryos/demographics_for_analysis.csv")
em_data <- merge(em_sex, em_dems, by = c("Pup_ID"))
em_data <- em_data[, c('litter_size', 'Mom_ID','treatment')] 

neo_data <- read.csv("/data/chamal/projects/lani/neonate/demographics.csv")
neo_data <- neo_data[, c('litter_size', 'mom_id','condition')] 
names(neo_data) <- c("litter_size","Mom_ID","treatment")

adult_data <-read.csv("/data/chamal/projects/lani/adult_thc/behavior/p1_weights.csv")
adult_data <- adult_data[, c('litter_size', 'Mom_ID','Treatment')] 
names(adult_data) <- c("litter_size","Mom_ID","treatment")

em_neo <- rbind(em_data, neo_data)
data <- rbind(em_neo, adult_data)
data_lit <- data %>% distinct(Mom_ID, .keep_all = TRUE)
data_lit$treatment <- gsub("sal","Sal",data_lit$treatment)
data_lit$treatment <- gsub("thc","THC",data_lit$treatment)
data_lit <- subset(data_lit, !is.na(data_lit$treatment))
data_lit$treatment <- relevel(as.factor(data_lit$treatment), ref = "Null")

data_lit$treatment <- relevel(data_lit$treatment, ref = "Sal")
litt_size <- lm(litter_size ~treatment, data = data_lit)
summary(litt_size)

lit_plot <- ggplot(data = data_lit, aes(x=treatment, y=litter_size, color=treatment, group=treatment)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.1, height = 0))+
  scale_color_manual('Treatment', values = c('Sal' = '#420085', 'THC' = '#018700', "Null" = "#0F52BA"))+
  labs(x = "Treatment", y= "Litter Size", title = "THC dams have smaller litters")+
  theme(plot.title = element_text(hjust = 0.5,size = 16, face = "bold"))+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=16), 
        legend.title=element_text(size=18)) +
  geom_signif(comparisons = list(c("Null","Sal"),c("Sal","THC")),
              annotations = c(paste("p=",round(summary(litt_size)$coefficients[2,4],digits=3)), paste("p=",signif(summary(litt_size)$coefficients[3,4],digits=3))),
              y_position = c(max(data_lit$litter_size)+0.5, max(data_lit$litter_size) + 1),textsize=4, colour = "black") +
  theme_classic()
lit_plot

#eof
