#Author: Lani Cupo
#Written: 20230731
#This script is used to analyze the data from maternal observations. 
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lme4)
library(lmerTest)

setwd("/data/chamal/projects/lani/adult_thc/behavior/derived_data/maternal_observations/")

#set up data frames
data <- read.csv("master_maternal_observations_old_thc.csv")
dems <- read.csv("../../raw_data/nest_quality.csv")
dems <- dems %>% distinct(Mom, .keep_all = T)
data_dems <- merge(data,dems, by.x = c("ID"), by.y = c("Mom"))

#First, check correlation between ethovision and manual duration and frequency -------------------------
etho_dat <- subset(data_dems, !is.na(data_dems$etho_duration_nest))

dur_check <- lm(manual_duration_nest ~ etho_duration_nest, etho_dat)
summary(dur_check)

freq_check <- lm(manual_freq_nest ~ etho_freq_nest, etho_dat)
summary(freq_check)

dur_check_plot <- ggplot(etho_dat, aes(x=manual_duration_nest, y=etho_duration_nest)) + geom_point() + geom_abline(slope = 1, intercept = 0)
freq_check_plot <- ggplot(etho_dat, aes(x=manual_freq_nest, y = etho_freq_nest)) + geom_point() + geom_abline(slope = 1, intercept = 0)
ggarrange(dur_check_plot, freq_check_plot)

#Prepare duration and frequency for analysis----------------------------------------------------------
data_dems$total_time <- data_dems$manual_duration_nest + data_dems$manual_duration_of
data_dems$percent_dur <- (data_dems$manual_duration_nest / data_dems$total_time) * 100

#Next, Analyze duration with manual accounts------------------------------------------------------------
duration_mod <- lmer(percent_dur ~ treatment*Day*Time + (1|ID), data = data_dems)
summary(duration_mod)

#plot duration
duration_plot <- ggplot(data = data_dems, aes(y=percent_dur,x=Day, color=treatment)) + facet_grid(. ~ Time)+
  geom_boxplot(aes(group=Day), outlier.shape = NA) +
  geom_line(aes(group = ID), alpha = 0.3) + geom_smooth(se = F) +
  xlab("Day") +  scale_x_continuous(limits = c(1, 9), breaks = c( 2, 4, 6, 8)) +
  ylab(bquote('Percent duration on nest')) + 
  labs(fill = "Treatment", group="Treatment", color="Treatment", title = "Maternal observation: duration")+
  scale_color_manual('Treatment', values = c('SAL' = '#420085', 'THC' = '#018700'))+
  scale_size_manual(values=c(2,2)) +
  theme_classic() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=16),
        legend.title=element_text(size=18), plot.title = element_text(size = 16, face = "bold"))
duration_plot

#Analyze frequency--------------------------------------------------------------------------------------
freq_mod <- lmer(manual_freq_nest ~ treatment*Day*Time + (1|ID), data = data_dems)
summary(freq_mod)

freq_plot <- ggplot(data = data_dems, aes(y=manual_freq_nest,x=Day, color=treatment)) + facet_grid(. ~ Time)+
  geom_boxplot(aes(group=Day), outlier.shape = NA) +
  geom_line(aes(group = ID), alpha = 0.3) + geom_smooth(se = F) +
  xlab("Day") +  scale_x_continuous(limits = c(1, 9), breaks = c( 2, 4, 6, 8)) +
  ylab(bquote('Frequency of passes on nest')) + 
  labs(fill = "Treatment", group="Treatment", color="Treatment", title = "Maternal observation: frequency")+
  scale_color_manual('Treatment', values = c('SAL' = '#420085', 'THC' = '#018700'))+
  scale_size_manual(values=c(2,2)) +
  theme_classic() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=16),
        legend.title=element_text(size=18), plot.title = element_text(size = 16, face = "bold"))
freq_plot

figure2<-ggarrange(duration_plot + rremove("xlab")  , freq_plot + rremove("xlab"),
                   labels = c("A", "B"),
                   ncol = 2, nrow = 1,
                   common.legend = TRUE, legend = "bottom",
                   font.label = list(size = 10, color = "black", face = "bold",  position = "top"))
figure2_annotated<-annotate_figure(figure2,
                                   bottom = text_grob("PND", size = 20, face = "bold"))

png("../../../../pte_paper/figures/mat_obs.png", units="in", width=9, height=6, res=300)
figure2_annotated
dev.off()
#eof