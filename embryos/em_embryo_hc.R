#author: Lani Cupo
#written: 20230530
#this script will first run lmers to analyze my cell data, then shift function to compare groups.
suppressPackageStartupMessages({
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(ggpubr)
  library(rogme) # for shift function
  library(dplyr)
  library(stringr) #to extract pup ID
  library(car) #for levene test
})

#set working directory

#Set up the data --------------------------------------------------------

setwd("/data/chamal/projects/lani/embryos/electron_microscopy/data/")

data <- read.csv("em_embryo_hc_mm2csv.csv", header = T)
dems <- read.csv("../../demographics_for_analysis.csv")
sex <- read.csv("../../transnetyx/transnetyx_all_results.csv")
# 
# #extract ID and number of areas per animal from the filename
# reg_ID <- "[0-9]{2}_[0-9]{2}"
# data <- data %>% mutate(Pup_ID = str_extract(Original_name, reg_ID))
# #number the samples within each animal
# data <- data %>%                              # Create numbering variable
#   group_by(Pup_ID) %>%
#   mutate(sample = row_number())
# 
# #turn each numbered sample into a letter for modelling
# data$sample <- sapply(data$sample, function(i) letters[i])
# 
# #Calculate the cells per mm2----------------------------------------------------
# data$dark_neuron_per_mm2 <- data$dark_neuron/data$area_mm2
# data$dark_glia_per_mm2 <- data$dark_glia/data$area_mm2
# data$total_dark_per_mm2 <- data$total_dark_cell/data$area_mm2
# 
# data$dark_neuron_per_um2 <- data$dark_neuron/data$area_um2
# data$dark_glia_per_um2 <- data$dark_glia/data$area_um2
# data$total_dark_per_um2 <- data$total_dark_cell/data$area_um2

#merge with demographics----------------------------------------------
data_dems <- merge(data, dems, by.x="ID", by.y = "Pup_ID")
data_dems <- merge(data_dems,sex, by.x="ID", by.y = "Pup_ID")

data_dems$apoptotic <- data_dems$Apoptotic.Cells.Condensed+data_dems$Apoptotic.Advanced.Necrosis+data_dems$Apoptotic.Necklas
data_dems$treatment <- as.factor(data_dems$treatment)
data_dems$treatment <- relevel(data_dems$treatment, ref = "Sal")
#run models------------------------------------------------------------------
neuron_mod <- lmer(Dark.Neuron ~ treatment + sex  + (1|ID), data= data_dems )
summary(neuron_mod)

glia_mod <- lmer(Dark.Glia ~ treatment + sex  + (1|ID), data= data_dems )
summary(glia_mod)

apop_mod <- lmer(apoptotic ~ treatment + sex  + (1|ID), data= data_dems)
summary(apop_mod)

div_mod <- lmer(Potentially.Dividing ~ treatment + sex  + (1|ID), data= data_dems)
summary(div_mod)
#create plots-------------------------------------------------------------------

div_plot <-ggplot(data = data_dems, aes(x=treatment, y=Potentially.Dividing, color=treatment, group=treatment)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3, height = 0))+
  #geom_boxplot(data = fitdata_jacob, aes(y=fit),size=1, alpha=0.8) +
  #geom_ribbon(data = fitdata_jacob, aes(y=fit, ymin=lower, ymax=upper, color=NULL), alpha=0.1) +
  scale_color_manual('Treatment', values = c('Sal' = '#420085', 'THC' = '#018700', "Null" = "#0F52BA"))+
  labs(x = "Treatment", y= bquote("Cell/"~mm^2), title = "Dividing Cells in the Embryo HC")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=20), legend.title=element_text(size=18)) +
  theme_classic()
png("../../../pte_paper/figures/embryo_dividing.png", units="in", width=3.4, height=3, res=300)
div_plot
dev.off()

apop <-ggplot(data = data_dems, aes(x=treatment, y=apoptotic, color=treatment, group=treatment)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3, height = 0))+
  #geom_boxplot(data = fitdata_jacob, aes(y=fit),size=1, alpha=0.8) +
  #geom_ribbon(data = fitdata_jacob, aes(y=fit, ymin=lower, ymax=upper, color=NULL), alpha=0.1) +
  scale_color_manual('Treatment', values = c('Sal' = '#420085', 'THC' = '#018700', "Null" = "#0F52BA"))+
  labs(x = "Treatment", y= bquote("Cell/"~mm^2), title = "Apoptotic Cells in the Embryo Hippocampus")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=20), legend.title=element_text(size=18)) +
  theme_classic()
apop

neurons <-ggplot(data = data_dems, aes(x=treatment, y=Dark.Neuron, color=treatment, group=treatment)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3, height = 0))+
  #geom_boxplot(data = fitdata_jacob, aes(y=fit),size=1, alpha=0.8) +
  #geom_ribbon(data = fitdata_jacob, aes(y=fit, ymin=lower, ymax=upper, color=NULL), alpha=0.1) +
  scale_color_manual('Treatment', values = c('Sal' = '#420085', 'THC' = '#018700', "Null" = "#0F52BA"))+
  labs(x = "Treatment", y= bquote("Cell/"~mm^2), title = "Dark neurons in the Embryo Hippocampus")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=20), legend.title=element_text(size=18)) +
  theme_classic()
neurons

glia <-ggplot(data = data_dems, aes(x=treatment, y=Dark.Glia, color=treatment, group=treatment)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3, height = 0))+
  #geom_boxplot(data = fitdata_jacob, aes(y=fit),size=1, alpha=0.8) +
  #geom_ribbon(data = fitdata_jacob, aes(y=fit, ymin=lower, ymax=upper, color=NULL), alpha=0.1) +
  scale_color_manual('Treatment', values = c('Sal' = '#420085', 'THC' = '#018700', "Null" = "#0F52BA"))+
  labs(x = "Treatment", y= bquote("Cell/"~mm^2), title = "Dark glia in the Embryo Hippocampus")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=20), legend.title=element_text(size=18)) +
  theme_classic()
glia

ggarrange(neurons, glia, apop, div_plot)
# Color by pup_ID
div_pupID <-ggplot(data = data_dems, aes(x=treatment, y=Potentially.Dividing, color=, group=treatment)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3, height = 0)) + 
  labs(x = "Treatment", y= "Cell/mm2", title = "Dividing Cells in the Embryo Hippocampus")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=20), legend.title=element_text(size=18)) +
  theme_classic()

div_pupID

#use the shift function------------------------------------------------------------------
data_dems <- subset(data_dems, !is.na(data_dems$Potentially.Dividing))

filtered.data.tibble <- as_tibble(data_dems)

sf <- shifthd_pbci(data = filtered.data.tibble, formula = Potentially.Dividing ~ treatment, todo = list(c("Null","Sal"),c("Null","THC"),c("Sal","THC")))
sf

plot_sf(sf)
p = plot_sf(sf)
#p = add_sf_lab(p, sf=sf)
sub = p[c(3)]
png("../../../pte_paper/figures/embryo_div_pairwise.png", units="in", width=5, height=3, res=300)
sub
dev.off()

p <- plot_scat2(data = filtered.data.tibble, formula = Potentially.Dividing ~ treatment,
                xlabel = "",
                ylabel = "Cell density per mm^2",
                alpha = 1,
                shape = 21,
                colour = "grey20",
                fill = "grey90",
                size = 3) +
  scale_x_discrete(breaks=c("Null","Sal", "THC"),
                   labels=c("Null","Sal", "THC")) +
  theme(axis.text.y = element_text(angle = 90, hjust = .5))
p <- plot_hd_bars(p,
                  col = "black",
                  q_size = 0.5,
                  md_size = 1.5,
                  alpha = 1)
p <- p + coord_flip() #> flip axes
pscat <- p
png("../../../pte_paper/figures/embryo_div_distro.png", units="in", width=5, height=3, res=300)
pscat
dev.off()
#check variance--------------------------------------------------
levene_neuron = leveneTest(dark_neuron_per_mm2 ~ treatment, data_dems)
levene_glia = leveneTest(dark_glia_per_mm2 ~ treatment, data_dems)
levene_total = leveneTest(total_dark_per_mm2 ~ treatment, data_dems)
print(levene_neuron)
print(levene_glia)
print(levene_total)

#do an anova on residuals manually to get group differences-----------------
#follows example here: https://stackoverflow.com/questions/43646987/multiple-comparison-post-hoc-test-for-levenes-test
data_dems <- data_dems %>% 
  mutate(med = median(total_dark_per_mm2, na.rm=TRUE))
data_dems$med.res<-abs(data_dems$total_dark_per_mm2-data_dems$med)
levene.total.aov<-aov(med.res~treatment,data_dems)
summary(levene.total.aov)
TukeyHSD(levene.total.aov)

#Nonparametric stats----------------------------------------------------
kruskal.test(total_dark_per_mm2 ~ treatment, data = data_dems)

res<- pairwise.wilcox.test(data_dems$total_dark_per_mm2, data_dems$treatment, p.adjust.method = 'BH')

#eof
