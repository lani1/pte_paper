#Author: Lani Cupo
#Written: Apr 16, 2023
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

#Maternal weight--------------------------------------

#set up dataframe
data <-read.csv("/data/chamal/projects/lani/neonate/milestones_weight/dam_weights.csv")
data <- subset(data, data$Pregnant == "Yes")
data$Treatment <- relevel(as.factor(data$Treatment), ref = "Sal")
data <- subset(data, data$Day <15)

#Create model
mod <- lmer(Weight ~ Treatment*Day + (1|unique_ID), data)
summary(mod)

fitdata_mod = as.data.frame(Effect(c("Day", "Treatment"),mod, xlevels=list(Day=seq(0,10,1))))

treatment <- ggplot(data = data, aes(y=Weight,x=Day, color=Treatment, group=Treatment, fill = Treatment))+
  geom_point() +
  geom_line(data = fitdata_mod, aes(y=fit)) +
  geom_ribbon(data = fitdata_mod, aes(y=fit, ymin=lower, ymax=upper), alpha=0.3,linetype=0) +
  scale_color_manual('Treatment', values = c('Sal' = '#420085', 'THC' = '#018700', "Null" = "#0F52BA"))+
  scale_fill_manual('Treatment', values = c('Sal' = '#420085', 'THC' = '#018700', "Null" = "#0F52BA"))+
  xlab("Gestational Day") +
  ylab("Weight (g)") + 
  #facet_grid(.~Sex)+
  labs(fill = "Treatment", group="Treatment", color="Treatment")+
  ggtitle("Treatment impacts dam weight gain") +
  scale_x_continuous(labels = number_format(accuracy = 1))+
  theme_classic()
treatment

data0 <- subset(data, data$Day == "0")
data0 <- subset(data0, !is.na(data0$Weight))
pre_treat <- lm(Weight ~Treatment, data0)
summary(pre_treat)

pre_treat_plot <-ggplot(data = data0, aes(x=Treatment, y=Weight, color=Treatment, group=Treatment)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.1, height = 0))+
  scale_color_manual('Treatment', values = c('Sal' = '#420085', 'THC' = '#018700', "Null" = "#0F52BA"))+
  labs(x = "Treatment", y= "Dam Weight", title = "Initial weight difference at GD0")+
  theme(plot.title = element_text(hjust = 0.5,size = 16, face = "bold"))+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=16), 
        legend.title=element_text(size=18)) +
  geom_signif(comparisons = list(c("Null","Sal"),c("Sal","THC")),
              annotations = c(paste("p=",round(summary(pre_treat)$coefficients[2,4],digits=3)), paste("p=",signif(summary(pre_treat)$coefficients[3,4],digits=3))),
              y_position = c(max(data0$Weight)+.5, max(data0$Weight) +2),textsize=4, colour = "black") +
  theme_classic()
pre_treat_plot

ggarrange(treatment,pre_treat_plot,common.legend = T, labels = c("A","B"))
