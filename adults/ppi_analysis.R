#Author: Lani Cupo Date: 20230628
#This script was written to analyze my behavioral data: ppI

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(dplyr)
  library(ggpubr)
  library(ggplot2)
  library(multcomp)
  library(reshape)
  library(splines)
  library(viridis)
  library(car) # necessary for levene's test
  library(effects)
  library(ggpubr)
})

#set up data-----------------------------------
setwd("/data/chamal/projects/lani/adult_thc/")
dems <- read.csv("demographics.csv")
#read in data
avg_trials <- read.csv("behavior/derived_data/ppi/avg_trials_ppi.csv")
avg_means <- read.csv("behavior/derived_data/ppi/avg_means_ppi.csv")
max_trials <- read.csv("behavior/derived_data/ppi/max_trials_ppi.csv")
max_means <- read.csv("behavior/derived_data/ppi/max_means_ppi.csv")
#merge data
dems_avg_trials <- merge(dems,avg_trials,by.x = c("ID"), by.y = "Subject")
dems_avg_means <- merge(dems, avg_means, by.x = c("ID"), by.y = "Subject")
dems_max_trials <- merge(dems, max_trials, by.x = c("ID"), by.y = "Subject")
dems_max_means <- merge(dems, max_means, by.x = c("ID"), by.y = "Subject")

#Calculate PPI--------------------------------------------
#Calculate Percent PPI and store as new column
ppi_calc_trials <- function(data){ data <- data %>%
  mutate(ppi03 = (((Start_mid - pp3) / Start_mid) * 100),
         ppi06 = (((Start_mid - pp6) / Start_mid) * 100),
         ppi09 = (((Start_mid - pp9) / Start_mid) * 100),
         ppi12 = (((Start_mid - pp12) / Start_mid) * 100),
         ppi15 = (((Start_mid - pp15) / Start_mid) * 100))
  #Melt the data frame to convert PPI into a variable with levels
  names = colnames(data)
  m.data <- melt(data, id.vars = names[1:19],
                 measure.vars = names[20:24])
 
   m.data <- m.data %>%
      dplyr::rename(PPI = variable, PPI.value = value)
  return(m.data)
}

ppi_calc_means <- function(data){ data <- data %>%
  mutate(ppi03 = (((STAR120 - PP3) / STAR120) * 100),
         ppi06 = (((STAR120 - PP6) / STAR120) * 100),
         ppi09 = (((STAR120 - PP9) / STAR120) * 100),
         ppi12 = (((STAR120 - PP12) / STAR120) * 100),
         ppi15 = (((STAR120 - PP15) / STAR120) * 100))
#Melt the data frame to convert PPI into a variable with levels
names = colnames(data)
m.data <- melt(data, id.vars = names[1:17],
               measure.vars = names[18:22])

m.data <- m.data %>%
  dplyr::rename(PPI = variable, PPI.value = value)
return(m.data)
}

clean_ppi <- function(data){
  #Specify unique columns for ID and PPI
  unique_rows <- distinct(data, ID, PPI, .keep_all = T)
  
  #Change values to characters to rename
  unique_rows$PPI <- as.character((unique_rows$PPI))
  
  #remove "ppi" from each label
  unique_rows$PPI[unique_rows$PPI == "ppi03"] <- "03"
  unique_rows$PPI[unique_rows$PPI == "ppi06"] <- "06"
  unique_rows$PPI[unique_rows$PPI == "ppi09"] <- "09"
  unique_rows$PPI[unique_rows$PPI == "ppi12"] <- "12"
  unique_rows$PPI[unique_rows$PPI == "ppi15"] <- "15"
  #Change PPI to a numeric column
  unique_rows$PPI <- as.numeric((unique_rows$PPI))
  
  #Include only positive %PPI
  unique_rows_pos <- subset(unique_rows, unique_rows$PPI.value >=0)
  
  return(unique_rows_pos)
}

dems_avg_ppi <- ppi_calc_trials(dems_avg_trials)
dems_max_ppi <- ppi_calc_trials(dems_max_trials)
avg_tri_pos <- clean_ppi(dems_avg_ppi)
max_tri_pos <- clean_ppi(dems_max_ppi)

dems_avg_means_ppi <- ppi_calc_means(dems_avg_means)
avg_mean_pos <- clean_ppi(dems_avg_means_ppi)

#Compare means produced by collate PPI vs calculated manually based on Elisa's work----
simp_tri_pos <- avg_tri_pos[, c('ID', 'PPI',"PPI.value")]
simp_mean_pos <- avg_mean_pos[,c("ID","PPI","PPI.value")]
merged <- merge(simp_tri_pos, simp_mean_pos, by = c("ID","PPI"))
names(merged)[3] = "PPI.value.tri"
names(merged)[4] = "PPI.value.mean"

ggplot(data = merged, aes(x = PPI.value.tri, y = PPI.value.mean))+
  geom_point()+
  geom_abline(slope = 1)
#there doesn't seem to be a systematic bias using just the output means, but I'll use Elisa's method for now.
#Analyze differences between groups---------------------------

#Run model 
avg_mod <- lmer(PPI.value ~treatment * PPI * sex  + (1|ID) + (1|mom_id), avg_tri_pos)
summary(avg_mod)
(res.table <- as.data.frame(coef(summary(avg_mod))))
xtable(res.table, type = "latex")

max_mod <- lmer(PPI.value ~treatment * PPI * sex  + (1|ID) + (1|mom_id), max_tri_pos)
summary(max_mod)
#Make Plots--------------------------------------------------------

fitdata_mod = as.data.frame(Effect(c("PPI", "treatment","sex"),avg_mod, xlevels=list(PPI_Level=seq(3,15,3))))
fig_ppi <- ggplot(data = avg_tri_pos, aes(x=PPI, y=PPI.value, color=treatment, group=treatment))+
  geom_line(data = fitdata_mod, aes(y=fit),size=1, alpha=0.8) +
  geom_ribbon(data = fitdata_mod, aes(y=fit, ymin=lower, ymax=upper), alpha=0.4, color = NA) +
  geom_jitter(alpha=0.7, size = 2, width = 0.3) +
  theme_classic()+
  labs(x="PPI Level (dB)", y="% PPI", title="Average PPI") +
  scale_color_manual('Treatment', values = c('SAL' = '#420085', 'THC' = '#018700'))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=16),
        legend.title=element_text(size=18), plot.title = element_text(size = 16, face = "bold"))+
  facet_wrap(~ sex)
fig_ppi

fitdata_mod_max = as.data.frame(Effect(c("PPI", "treatment","sex"),max_mod, xlevels=list(PPI_Level=seq(3,15,3))))
fig_ppi_max <- ggplot(data = max_tri_pos, aes(x=PPI, y=PPI.value, color=treatment, group=treatment))+
  geom_line(data = fitdata_mod_max, aes(y=fit),size=1, alpha=0.8) +
  geom_ribbon(data = fitdata_mod_max, aes(y=fit, ymin=lower, ymax=upper), alpha=0.4, color = NA) +
  geom_jitter(alpha=0.7, size = 2, width = 0.3) +
  theme_classic()+
  labs(x="PPI Level (dB)", y="% PPI", title="Maximum PPI") +
  scale_color_manual('Treatment', values = c('SAL' = '#420085', 'THC' = '#018700'))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=16),
        legend.title=element_text(size=18), plot.title = element_text(size = 16, face = "bold"))+
  facet_wrap(~ sex)
fig_ppi_max

#check overall startle by condition/sex------------------------------------------------
#Eliminate duplicates
unique_avg_tri <- avg_tri_pos %>%
  distinct(ID, .keep_all = TRUE)

mod_startle <- lmer(Start_mid ~treatment * sex + Chamber + (1|mom_id), unique_avg_tri)
summary(mod_startle)
 
fig_startle <- ggplot(data = unique_avg_tri, aes(x=treatment, y=Start_mid, color=treatment, group=treatment))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(alpha=0.7, size = 2, width = 0.3) +
  theme(axis.text=element_text(size=14), 
        legend.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5,size = 16, face = "bold"))+
  theme_classic()+
  labs(x="Treatment", y="Startle", title="Startle") +
  scale_color_manual('Treatment', values = c('SAL' = '#420085', 'THC' = '#018700'))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=16),
        legend.title=element_text(size=18), plot.title = element_text(size = 16, face = "bold"))+
  facet_wrap(~ sex)
fig_startle

#eof

