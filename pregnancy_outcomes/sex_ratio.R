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

#set up dataframe------------------------
#embryos
em_sex <- read.csv("/data/chamal/projects/lani/embryos/transnetyx/transnetyx_all_results.csv")
em_dems <- read.csv("/data/chamal/projects/lani/embryos/demographics_for_analysis.csv")
em_data <- merge(em_sex, em_dems, by = c("Pup_ID"))
head(em_data)

#Use sex because data is just for scanned mice
em_sex$Mom_ID <- as.numeric(substr(em_sex$Pup_ID, start = 1, stop = 2))

#create a dataframe with the number of females
em_sex_ratio <- em_sex %>% 
  group_by(Mom_ID)%>%
  summarise(n_female = sum(sex == "female", na.rm = T))
#merge sex and sex ratio by mom id
em_sex_ratio <- merge(em_sex, em_sex_ratio, by = c("Mom_ID"))
#calculate the proportion of females from total litter size
em_sex_ratio$prop_female <- em_sex_ratio$n_female / em_sex_ratio$litter_size
#merge with dems for treatment
em_ratio_treat <- merge(em_sex_ratio, em_dems, by = c("Mom_ID"))
#reduce to a single value per mom again
em_ratio_treat <- em_ratio_treat %>% distinct(Mom_ID, .keep_all = TRUE)
em_ratio_treat$Mom_ID <- paste0("e",em_ratio_treat$Mom_ID)
em_ratio_treat <- em_ratio_treat[, c('prop_female', 'Mom_ID','treatment')] 
names(em_ratio_treat)[1] <- c("ratio_female")

#neonates
neo_dems <- read.csv("/data/chamal/projects/lani/neonate/demographics.csv")
neo_data_sex <- neo_dems %>% distinct(ID, .keep_all = TRUE)
neo_mom_count_male <-neo_data_sex %>% group_by(mom_id) %>% summarise(count_male = sum(sex == "male"))
neo_mom_count_female <- neo_data_sex %>% group_by(mom_id) %>% summarise(count_female = sum(sex == "female"))
neo_ratio_female <- neo_mom_count_female$count_female/(neo_mom_count_male$count_male+neo_mom_count_female$count_female)
neo_mom_count_male$ratio_female <- neo_ratio_female
neo_data_merge <- merge(neo_data_sex, neo_mom_count_male, by = c("mom_id"))
neo_data_mom <- neo_data_merge %>% distinct(mom_id, .keep_all = TRUE)
neo_data_mom$mom_id <- paste0("n", neo_data_mom$mom_id)
neo_data_mom <- neo_data_mom[, c('ratio_female', 'mom_id','condition')] 
names(neo_data_mom) <- c("ratio_female","Mom_ID","treatment")

#adults
ad_dat <- read.csv("/data/chamal/projects/lani/adult_thc/behavior/sex_ratio.csv")
ad_dat$mom_id <- paste0("a", ad_dat$mom_id)
ad_dat <- ad_dat[, c('n_prop_female', 'mom_id','condition')] 
names(ad_dat) <- c("ratio_female","Mom_ID","treatment")

em_neo <- rbind(em_ratio_treat, neo_data_mom)
data <- rbind(em_neo, ad_dat)
data$treatment <- gsub("sal","Sal",data$treatment)
data$treatment <- gsub("thc","THC",data$treatment)

data$treatment <-relevel(as.factor(data$treatment), ref = "Sal")
sex_rat <- lm(ratio_female ~treatment, data = data)
summary(sex_rat)

sex_plot <- ggplot(data = data, aes(x=treatment, y=ratio_female, color=treatment, group=treatment)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.1, height = 0))+
  scale_color_manual('Treatment', values = c('Sal' = '#420085', 'THC' = '#018700', "Null" = "#0F52BA"))+
  labs(x = "Treatment", y= "Ratio of females", title = "Treatment does not alter ratio of females")+
  theme(plot.title = element_text(hjust = 0.5,size = 16, face = "bold"))+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=16), 
        legend.title=element_text(size=18)) +
  geom_signif(comparisons = list(c("Null","Sal"),c("Sal","THC")),
              annotations = c(paste("p=",round(summary(sex_rat)$coefficients[2,4],digits=3)), paste("p=",signif(summary(sex_rat)$coefficients[3,4],digits=3))),
              y_position = c(max(data$ratio_female)+0.15, max(data$ratio_female) + 0.3),textsize=4, colour = "black") +
  theme_classic()
sex_plot
