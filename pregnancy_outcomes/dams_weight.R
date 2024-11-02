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

data <- read.csv("/data/chamal/projects/lani/adult_thc/behavior/dam_weights.csv")
data <- subset(data, data$pregnant == "yes")

mod <- lmer(weight ~ treatment*GD + (1|mom_id), data)
summary(mod)
(res.table <- as.data.frame(coef(summary(mod))))
xtable(res.table, type = "latex")

fitdata_mod = as.data.frame(Effect(c("GD", "treatment"),mod, xlevels=list(Day=seq(0,10,1))))

treatment <- ggplot(data = data, aes(y=weight,x=GD, color=treatment, group=treatment, fill = treatment))+
  geom_point() +
  geom_line(data = fitdata_mod, aes(y=fit)) +
  geom_ribbon(data = fitdata_mod, aes(y=fit, ymin=lower, ymax=upper), alpha=0.3,linetype=0) +
  scale_color_manual('Treatment', values = c('SAL' = '#420085', 'THC' = '#018700'))+
  scale_fill_manual('Treatment', values = c('SAL' = '#420085', 'THC' = '#018700'))+
  xlab("Gestational Day") +
  ylab("Weight (g)") + 
  #facet_grid(.~Sex)+
  labs(fill = "Treatment", group="Treatment", color="Treatment")+
  ggtitle("Treatment impacts dam weight gain") +
  scale_x_continuous(labels = number_format(accuracy = 1))+
  theme_classic()
treatment

data0 <- subset(data, data$GD == "0")
pre_treat <- lm(weight ~treatment, data)
summary(pre_treat)

pre_treat_plot <-ggplot(data = data0, aes(x=treatment, y=weight, color=treatment, group=treatment)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.1, height = 0))+
  scale_color_manual('Treatment', values = c('SAL' = '#420085', 'THC' = '#018700'))+
  labs(x = "Treatment", y= "Dam Weight", title = "Initial weight difference at GD0")+
  theme(plot.title = element_text(hjust = 0.5,size = 16, face = "bold"))+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=16), 
        legend.title=element_text(size=18)) +
  geom_signif(comparisons = list(c("THC","SAL")),
              annotations = c(paste("p=",round(summary(pre_treat)$coefficients[2,4],digits=3))),
              y_position = c(max(data0$weight)+.5),textsize=4, colour = "black") +
  theme_classic()
pre_treat_plot

png("/data/chamal/projects/lani/adult_thc/visuals/dam_weights.png", units="in", width=7, height=5, res=300)

ggarrange(treatment,pre_treat_plot,common.legend = T, labels = c("A","B"))
dev.off()
#eof
