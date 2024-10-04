#Authors: Lani Cupo and Katerina Bradshaw
#Date: 20230705
#this script was written to analyze OFT data from the adult experiment

# Run with /data/chamal/projects/lani/runRneonate.sh

#Load necessary packages
library(lme4)
library(lmerTest)
library(dplyr)
library(viridis)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(rogme)

#setting up the environment
setwd("/data/chamal/projects/lani/adult_thc/")
dems <- read.csv("demographics.csv")
oft <- read.csv("behavior/derived_data/oft/master_oft.csv")

data <- merge(oft, dems, by= c("ID"))
data <- subset(data, data$timepoint == "OFT")

#Analyze total distance moved-----------------------------------------------

#Zone 5 is the center

data$total_distance_moved <- ave(data$Distance_moved, data$ID, FUN = sum)

data <- as.data.frame(data)
head(data)
data$treatment[data$treatment == "SAL"] <- "Sal"

distance_unique <- distinct(data, ID, total_distance_moved, .keep_all = T)

head(distance_unique)

#statistically compare the groups
mod_distance <- lmer(total_distance_moved ~treatment * sex + (1|mom_id), data = distance_unique)
summary(mod_distance)


# Fixed effects:
#   Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)           5596.48     401.53    25.06  13.938 2.61e-13 ***
#   treatmentTHC         -1442.81     550.21    22.98  -2.622   0.0152 *  
#   sexmale              -1018.64     497.21    35.39  -2.049   0.0480 *  
#   treatmentTHC:sexmale  1593.93     695.51    34.99   2.292   0.0281 * 

fig_dist <- ggplot(data = distance_unique, aes(x=treatment, y=total_distance_moved, color=treatment, group=treatment))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(alpha=0.7, size = 2, width = 0.3) +
  labs(x="Treatment", y="Distance moved (cm)", title="Total distance moved") +
  scale_color_manual('Treatment', values = c('Sal' = '#420085', 'THC' = '#018700'))+
  theme_classic()+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=16),
        legend.title=element_text(size=18), plot.title = element_text(size = 16, face = "bold"))+
  facet_wrap(~ sex) 
fig_dist

time_center <- subset(data, data$Zone=="5")

#statistically compare the groups
mod_center <- lmer(Duration ~treatment * sex + (1|mom_id), data = time_center)
summary(mod_center)

(res.table <- as.data.frame(coef(summary(mod_center))))
xtable(res.table, type = "latex")

# Fixed effects:
#   Estimate Std. Error     df t value Pr(>|t|)    
# (Intercept)            158.45      21.39  30.22   7.408 2.84e-08 ***
#   treatmentTHC           -93.27      29.11  27.80  -3.204  0.00339 ** 
#   sexmale                -35.41      29.05  36.34  -1.219  0.23074    
# treatmentTHC:sexmale    81.16      40.69  35.84   1.995  0.05374 .  

#fitdata_center = as.data.frame(Effect(c("treatment","sex"),mod_center))
fig_center <- ggplot(data = time_center, aes(x=treatment, y=Duration, color=treatment, group=treatment))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(alpha=0.7, size = 2, width = 0.3) +
  theme(axis.text=element_text(size=14), 
        legend.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5,size = 16, face = "bold"))+
  theme_classic()+
  labs(x="Treatment", y="Time spent in center (sec)", title="Time spent in center") +
  scale_color_manual('Treatment', values = c('Sal' = '#420085', 'THC' = '#018700'))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=16),
        legend.title=element_text(size=18), plot.title = element_text(size = 16, face = "bold"))+
  facet_wrap(~ sex)
fig_center

#statistically compare the groups
mod_center_freq <- lmer(Frequency ~treatment * sex + (1|mom_id), data = time_center)
tab_model(mod_distance, file = "../pte_paper/stats_tables/oft_dist_lmer.doc")
xtable(res.table, type = "latex")
# Fixed effects:
#   Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)            58.334      5.743  23.586  10.158 4.36e-10 ***
#   treatmentTHC          -23.609      7.886  21.657  -2.994  0.00676 ** 
#   sexmale               -11.256      6.889  35.139  -1.634  0.11118    
# treatmentTHC:sexmale   22.946      9.633  34.759   2.382  0.02282 * 

#fitdata_center = as.data.frame(Effect(c("treatment","sex"),mod_center))
fig_center_freq <- ggplot(data = time_center, aes(x=treatment, y=Frequency, color=treatment, group=treatment))+
  geom_boxplot(outlier.alpha = NA)+
  geom_jitter(alpha=0.7, size = 2, width = 0.3) + theme_classic()+
  labs(x="Treatment", y="Frequency of passes", title="Passes through center by condition") +
  scale_color_manual('Treatment', values = c('Sal' = '#420085', 'THC' = '#018700'))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=16),
        legend.title=element_text(size=18), plot.title = element_text(size = 16, face = "bold"))+
  facet_wrap(~ sex) 
fig_center_freq

figure2<-ggarrange(fig_dist + rremove("xlab")  , fig_center + rremove("xlab"),
                   #fig_center_freq + rremove("xlab"),
                   #labels = c("A", "B", "C"),
                   ncol = 2,# nrow = 2,
                   common.legend = TRUE, legend = "bottom",
                   font.label = list(size = 10, color = "black", face = "bold",  position = "top"))
figure2_annotated<-annotate_figure(figure2,
                                   bottom = text_grob("Treatment", size = 20, face = "bold"))
png("../pte_paper/figures/oft_only.png", units = "in", width=8, height=4, res=300)

figure2_annotated
dev.off()

#Plot with centered third plot
blank <- ggplot() + theme_void()

three_plot <- ggarrange(
  ggarrange(fig_dist+ rremove("xlab"), fig_center+ rremove("xlab"), blank, nrow = 1,
            widths = c(1, 1, 0), legend = "none"),
  ggarrange(blank, fig_center_freq+ rremove("xlab"), blank, nrow = 1,
            widths = c(0.5, 1, 0.5), legend = "none"),
  nrow = 2,
  legend = "none"
)
figure_annotated<-annotate_figure(three_plot,
                                   bottom = text_grob("Treatment", size = 20, face = "bold"))
figure_annotated

#eof

