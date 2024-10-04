#Author: Lani Cupo
#Last Updated: 20210727
# This script was written to analyze the first batch of data. 
#should be run with runRnew.sh

setwd('/data/chamal/projects/lani/neonate/usvs/')
suppressPackageStartupMessages({
  library(ggplot2)  
  library(lme4)
  library(lmerTest)
  library(graphics)
  library(rogme)
  library(summarytools)
  library(tidyverse)
  library(ggpubr)
  library(gridExtra)
  library(dplyr)
})

calls <- read.csv("master_usv.csv")
dems <- read.csv("../demographics.csv")

data <- merge(calls, dems, by = "ID")
data$sex <- as.factor(data$sex)
data$sex <- relevel(data$sex, ref = "female")
data$condition <- as.factor(data$condition)
data$condition <- relevel(data$condition, ref = "sal")
data$scanned <- as.factor(data$scanned)
data$scanned <- relevel(data$scanned, ref = "no")

# remove calls longer than 300ms. 
filtered.data=subset(data, data$Duration <300)

#for sample size
#calculate total, complete sample size
unique_data <- data %>%
  distinct(data$ID, .keep_all = T)

ftable(unique_data$scanned, unique_data$condition, unique_data$sex)


#Do stats on the USVs
#drop "incomplete" from scanned variable
filtered.data<-droplevels(filtered.data)
sal_only <- subset(filtered.data, filtered.data$condition == "sal")
sal_only <- droplevels(sal_only)
#filtered.data$tx=interaction(filtered.data$sex)
filtered.data.tibble <- as_tibble(filtered.data)
sal.tibble <- as_tibble(sal_only)

unique_data_sal <- sal_only %>%
  distinct(sal_only$ID, .keep_all = T)

table(unique_data_sal$scanned)

#Right now, filtered.data.tibble contains repeated calls because there is an entry for each timepoint. We want to only keep
#one row per ID per call. WE can do this by choosing unique ID + begin times.
filtered.data.tibble$ID_begin.time <- paste(filtered.data.tibble$ID, filtered.data.tibble$Begin.time, sep="_")
one_call <- filtered.data.tibble %>%
  distinct(filtered.data.tibble$ID_begin.time, .keep_all = T)

#uncomment to check the results in the mncl2 exposed subset
# sub_data <- subset(one_call, one_call$MnCl2 == "manganese")
# sub_data <- subset(sub_data, sub_data$scanned == "yes")
# one_call <- sub_data

##As an example why we shouldn't just examine the means with a t-test, examine the means with a t-test
t.test(Duration~condition, one_call)
lmer_mod <- lmer(Duration ~ condition*sex + (1|mom_id) + (1|ID), data = one_call)
summary(lmer_mod)
png("../visuals/mean_calls.png", units="in", width=4, height=3, res=300)

ggplot(data = one_call, aes(x=condition, y=Duration, color=condition, group=condition)) +
  geom_boxplot(outlier.shape = NA)+
  #geom_point(alpha=0.2, size = 2) +
  #geom_jitter(position = position_jitter(width = 0.3, height = 0),aes(alpha=0.1))+
  scale_color_manual('Treatment', values = c('sal' = '#420085', 'thc' = '#018700'))+
  labs(x = "Treatment", y= "Duration", title = "Mean Duration of Calls")+
  theme(plot.title = element_text(hjust = 0.5,size = 16, face = "bold"))+
  theme_classic()+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"), legend.text=element_text(size=16), 
        legend.title=element_text(size=18))
dev.off()

one_call$sex_treatment <- as.factor(paste0(one_call$sex, one_call$condition))
predictors <- c("scanned")
one_call$scanned <- sub ("yes" , "scanned" , one_call$scanned)
one_call$scanned <- sub ("no" , "unscanned" , one_call$scanned)
one_call$scanned <- as.factor(one_call$scanned)
#predictors <- c("sex", "scanned", "condition")
for (each in predictors){
  sf <- shifthd_pbci(data = one_call, formula = formula(paste0("Duration ~ ",each)))
  p = plot_sf(sf)
  p = add_sf_lab(p, sf=sf)
  p <- plot_scat2(data = one_call, formula = formula(paste0("Duration ~ ",each)),
                  xlabel = "",
                  ylabel = "Call length (ms)",
                  alpha = 1,
                  shape = 21,
                  colour = "grey20",
                  fill = "grey90",
                  size = 3) +
    #scale_x_discrete(breaks=c("female", "male"),
     #                labels=c("female", "male")) +
    theme(axis.text.y = element_text(angle = 90, hjust = .5))
 p <- plot_hd_bars(p,
                   col = "black",
                   q_size = 0.5,
                   md_size = 1.5,
                   alpha = 1)
  p <- p + coord_flip() #> flip axes
  pscat <- p
  pscat
}
png("../visuals/neo_usv_distro_scan.png", units="in", width=8, height=5, res=300)
pscat
dev.off()

sf_cond <- shifthd_pbci(data = one_call, formula = Duration ~ sex_treatment, doall = T)
sf_scan <- shifthd_pbci(data = one_call, formula = Duration ~ scanned, todo = list(c("scanned","unscanned")))
#set this to change between looking at condition or scan!
sf <- sf_cond

png("../../pte_paper/figures/neonate_usvs_cond_20240627.png", units="in", width=4, height=9, res=300)

plot_sf(sf)
p = plot_sf(sf)
#p = add_sf_lab(p, sf=sf)
sub = p[c(1,2,6)]
do.call("grid.arrange", c(sub, ncol=1))

dev.off()



#make scaterplot colored by ID
ggplot(data = one_call,aes(x=Duration, y = sex_treatment, color = ID)) +
  geom_jitter()
#eof

