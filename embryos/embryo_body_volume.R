#Author: Lani Cupo
#Date: 20220426
#This script was written to analyze the difference between body volume for my embryos 

#run with runRneonate.sh

#Load Libraries
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(ggplot2)
  library(ggpubr)
  library(ggsignif)
  library(EnvStats)
  library(dplyr)
  library(texreg)
  library(sjPlot)
})
library(clipr)
#Set working directory locally
setwd("/data/chamal/projects/lani/embryos/")

#Load in data --------------------------------------
dems <- read.csv("demographics_for_analysis.csv")
sex <- read.csv("transnetyx_all_results.csv")
sex_dems <- merge(sex, dems, by = "Pup_ID")
vols <- read.csv("label_vols.csv")
data <- merge(sex_dems, vols, by = c("merge"))
qc <- read.csv("qc_csvs/squish_qc.csv")
data <- merge(data,qc, by = c("merge") )
no_squish_data <- subset(data, data$Final_qc <1)
no_squish_data$treatment <- as.factor(no_squish_data$treatment)
  
no_squish_data$treatment <- relevel(no_squish_data$treatment, ref = "Sal")

#Set theme for ggplots------------------------------

myTheme <- theme_classic()+ theme(axis.text=element_text(size=16), axis.title=element_text(size=16), 
                                  legend.text=element_text(size=16), legend.key.size = unit(1, 'cm'),
                         legend.title=element_text(size=18), plot.title = element_text(size = 16))
#Run model--------------------------------------
mod <- lmer(scale(volume) ~ treatment * sex + (1|litter_size) + (1|coil), data = no_squish_data)
summary(mod)
tab_model(mod, file = "../pte_paper/stats_tables/embryo_volume_lmer.doc")
write_clip(res.table)
no_squish_data$volume_resid <- resid(mod)

write.csv(x = no_squish_data, file = "data_with_volume_resids.csv")
write.csv(x = res.table, file = "../stats_tables/embryo_volume_lmer.csv")
# make plot--------------------------------------
no_squish_data$treatment <- factor(no_squish_data$treatment , levels=c("Sal", "Null", "THC"))
png("visuals/embryo_volume.png", units="in", width=5, height=5, res=300)

em_plot <- ggplot(data = no_squish_data, aes(x=treatment, y=volume, color=treatment, group=treatment)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(shape = sex, size = 3), position = position_jitter(width = 0.3, height = 0))+ 
  # geom_signif(comparisons = list(c("Null","Sal"),c("Sal","THC")),
  #             annotations = c(paste("p=",round(summary(mod)$coefficients[2,5],digits=3)), 
  #                             paste("p=",signif(summary(mod)$coefficients[3,5],digits=3))),
  #             y_position = c(1100, 1170),textsize=6, colour = "black") +
  scale_color_manual('Treatment', values = c('Sal' = '#420085', 'THC' = '#018700', "Null" = "#0F52BA"))+
  labs(x = "Treatment", title = "Treatment alters embryo volume")+ ylab(bquote('Embryo Volume '(mm^3))) +
  myTheme
em_plot
# insert ggplot code
dev.off()

#Examine percent diff-------------------------------
avg_txt <- no_squish_data %>% group_by(treatment) %>% 
  summarise(mean_vol=mean(volume),
            .groups = 'drop')
avg_txt = as.data.frame(avg_txt)

#percent diff null and sal
null_sal_bot <- (avg_txt$mean_vol[1] + avg_txt$mean_vol[2])/2
null_sal_pct <- (avg_txt$mean_vol[2] - avg_txt$mean_vol[1])/null_sal_bot *100
#percent diff thc and sal
thc_sal_bot <- (avg_txt$mean_vol[1] + avg_txt$mean_vol[3])/2
thc_sal_pct <- (avg_txt$mean_vol[3] - avg_txt$mean_vol[1])/thc_sal_bot *100
null_sal_pct
thc_sal_pct
