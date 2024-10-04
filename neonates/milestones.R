#Author: Lani Cupo 
#Last updated 20220311
#This script was written to analyse the milestone data from my neonates

#First I will just see if there is a difference for the "first appearance" of each milestone. Second I will 
#analyze time days since first checked, since if a milestone was checked for one group on PND 10 and the other on PND11
#it may appear that group 2 is delayed compared to group 1, but this way they would both be assigned 0. Third, I will use survival 
#curves to measure the development of milestones.

suppressPackageStartupMessages({
library(Matrix)
library(lme4)
library(lmerTest)
library(dplyr)
library(viridisLite)
library(viridis)
library(ggplot2)
library(magrittr)
library(ggpubr)
library(splines)})

setwd("/data/chamal/projects/lani/neonate/")
data <- read.csv("./milestones_weight/milestones_999_as_NA.csv", encoding = "UTF-8", dec = ".", header = T)

data <- subset(data,scanned %in% c("yes","no"))
#Loop through the relevant variables
pvalue = matrix(, nrow = 0, ncol = 4) 
index = seq(8,13)
for (varname in names(data)[index]) {
  print(varname)
  mod = lmer(as.formula(paste0(varname,"~ condition + scanned + sex + (1|mom_id)")), data=data)
  summary_mod = summary(mod)
  pvalue=rbind(pvalue, summary_mod$coefficients[,"Pr(>|t|)"]) 
}
  
rownames(pvalue) = names(data)[index]
pvalue
#Plot milestones-------------------------------------------------------------------
plot_milestone <- function(variable,title,term){
  if(term[1] %in% c("thc","sal")){
  p <- ggplot(data, aes(x=term, y=variable, color=term)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter()+ labs(x="Treatment", y="Day of Appearance", title=title)+
    scale_color_manual('Treatment', values = c('sal' = '#420085', 'thc' = '#018700'))+
    theme(axis.text = element_text(size= 20)) + theme(axis.title = element_text(size = 18)) + 
    theme(legend.text = element_text(size = 18)) + theme(legend.title = element_blank())
  }else if(term[1] %in% c("male","female")){
    p <- ggplot(data, aes(x=term, y=variable, color=term)) +
      geom_boxplot(outlier.shape = NA) + 
      geom_jitter()+ labs(x="Sex", y="Day of Appearance", title=title)+
      scale_color_manual('Sex', values = c('female' = '#ea88ac', 'male' = '#88eab2'))+
      theme(axis.text = element_text(size= 20)) + theme(axis.title = element_text(size = 18)) + 
      theme(legend.text = element_text(size = 18)) + theme(legend.title = element_blank())
  }else if(term[1] %in% c("yes","no","incomplete")){
    p <- ggplot(data, aes(x=term, y=variable, color=term)) +
      geom_boxplot(outlier.shape = NA) + 
      geom_jitter()+ labs(x="Sex", y="Day of Appearance", title=title)+
      theme(axis.text = element_text(size= 20)) + theme(axis.title = element_text(size = 18)) + 
      theme(legend.text = element_text(size = 18)) + theme(legend.title = element_blank())
  }
}

#to change the term, change this to data$sex, data$condition, or data$scanned
input_term <- data$sex
incisor_eruption <- plot_milestone(data$incisor_eruption,"Incisor Eruption",input_term)
fur_development <- plot_milestone(data$fur_development, "Fur Development", input_term)
pinnae_detachment <- plot_milestone(data$pinnae_detachment, "Pinnae Detachment",input_term)
vib_plac_resp <- plot_milestone(data$vibrissa_placing_response, "Vibrissa Placing Response",input_term)
ear_twitch <- plot_milestone(data$ear_twitch,"Ear Twitch",input_term)
aud_startle <- plot_milestone(data$auditory_startle_response,"Auditory Startle",input_term)
surf_right <- plot_milestone(data$surface_righting_reflex, "Surface Righting",input_term)
grasp_reflex <- plot_milestone(data$grasp_reflex, "Grasp Reflex",input_term)

figure<-ggarrange(incisor_eruption + rremove("ylab") + rremove("xlab"), fur_development + rremove("ylab") + rremove("xlab"),  pinnae_detachment+ rremove("ylab") + rremove("xlab"), vib_plac_resp+ rremove("ylab") + rremove("xlab"), # remove axis labels from plots
                  ear_twitch + rremove("ylab") + rremove("xlab"), aud_startle + rremove("ylab") + rremove("xlab"), surf_right + rremove("ylab") + rremove("xlab"), grasp_reflex + rremove("ylab") + rremove("xlab"),
                  ncol = 4, nrow = 2,
                  common.legend = TRUE, legend = "bottom",
                  font.label = list(size = 10, color = "black", face = "bold",  position = "top"))
figure_annotated<-annotate_figure(figure, left = text_grob("Day of Appearance", rot = 90, vjust = 1, size = 20, face = "bold"),
                                  bottom = text_grob("Sex", size = 20, face = "bold"))
figure_annotated

# check proportion of 999s per group---------------------------------------------
#Function to calculate z-score 
#change these three variables if you want to look at condition, scanned, or sex
data <- read.csv("./milestones_weight/milestones.csv", encoding = "UTF-8", dec = ".", header = T)
control = "no"
compare = "yes"
term2 = data$scanned

find_z <- function(milestone, term_name){
  con999 <- length(which(milestone=="999" & term2==control))
  contot <- length(which(term2==control))
  com999 <- length(which(milestone=="999" & term2==compare))
  comtot <- length(which(term2==compare))
  print(c(con999,contot,com999,comtot))
  res <- prop.test(x = c(con999,com999), n = c(contot,comtot))
  return(print(paste0(term_name[1]," is ",res$p.value)))
}

#incisor_eruption_z <- find_z(data$incisor_eruption) #no 999
fur_development_z <- find_z(data$fur_development, term2)
#pinnae_detachment_z <- find_z(data$pinnae_detachment) #too few 999
#vib_plac_resp_z <- find_z(data$vibrissa_placing_response) #too few
ear_twitch_z <- find_z(data$ear_twitch, term2)
aud_startle_z <- find_z(data$auditory_startle_response, term2)
surf_right_z <- find_z(data$surface_righting_reflex, term2)
#grasp_reflex_z <- find_z(data$grasp_reflex) #no 999

res <- prop.test(x = c(3, 11), n = c(24, 51))
# Printing the results
res 
#eof