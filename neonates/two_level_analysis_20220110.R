#Author: Lani Cupo
#Date: 20210718
#This script was written to analyze my neonates.

#run with runRnew.sh  (run in R version 3.5.1)

#Load Libraries
suppressPackageStartupMessages({
library(lme4)
library(lmerTest)
library(ggplot2)
library(RMINC)
library(dplyr)
library(viridis)
library(ggeffects)
library(grid)
library(tidyverse)
library(MRIcrotome)
library(magrittr)
library(effects)
library(ggpubr)
  library(xtable)
})

#---------------------------------------Set Up Environment-----------------------------------------------------
#Set working directory locally
setwd("/data/scratch2/cuplan/neonate/neo_thc/20230508_gabe_dbm_outputs/")
#setwd("/data/scratch/cuplan/neonate/neo_thc/alternative_dbm/")

#Load data
jac <- read.csv(file = "full_smooth/full_jacs.csv", header = F)
#jac <- read.csv(file = "relative_smooth/rel_jacs.csv", header = F)
dems <- read.csv(file = "/data/chamal/projects/lani/neonate/data_with_weight_resids.csv") #, header = TRUE)
colnames(dems)[colnames(dems) == "timepoint"] <- "age"
colnames(dems)[colnames(dems) == "id"] <- "ID"
#--------------------------------Set Up spreadsheet----------------------------
#first parse the text to include a variable that is the ID and PND
colnames(jac) <- c("file")
jac$tempID <- sub(".*MCH_NEO_THC_", "", jac$file)
jac$ID <- substring(jac$tempID, 1, 5)
jac$ID

#rassign age automatically from filename
jac$tempPND <- sub(".*MCH_NEO_THC_[0-9]{2}_[0-9]{2}_", "", jac$file)
jac$age <- substring(jac$tempPND,1,2)
jac$age <- as.numeric(jac$age)
jac$age

qc <- read.csv("/data/chamal/projects/lani/neonate/qc/agreement_qc.csv", header = T)
data <- merge(jac,dems,by = c("ID", "age"))
data <- merge(data, qc, by.x = c("ID","age"), by.y = c("ID","timepoint"))

#cluster set up
#options(RMINC_BATCH_CONF = "/data/chamal/projects/lani/MIA_lani/sge_resources.R")

#set mask
mask = "/data/chamal/projects/lani/neonate/make_mask/eroded_extracted_mask.mnc"
#Set Female and control as default
data$sex = relevel(data$sex, ref = "female")
data$condition = relevel(data$condition, ref = "sal")

#Drop those with only 1 age
two_times <- data %>% group_by(ID) %>% filter(n()>1)
two_times <- droplevels(two_times)
#two_times <- two_times %>% ungroup(two_times)
ftable(two_times$age, two_times$condition, two_times$sex)

#I already dropped poor qc (3 or worse), but I should run the model twice, once with only those who got a 1 or 1.5
high_qc <- subset(two_times, two_times$final_qc <= 1.5)
#-------------------------------Compare models--------------------------------------
mlm <-
  mincLmer(file.x ~ condition * I(age-3) * sex + weight_resid + (1|ID) + (1|litter_size)
           , two_times#two_timess
           , mask = mask)
           #, parallel = c("sge", 50))
saveRDS(mlm, "neonate_mod_full.rds")
mlm = mincLmerEstimateDF(mlm)   
FDR = mincFDR(mlm, mask= mask)
FDR

#absolute jacobians
# FDR Thresholds:
#   tvalue-(Intercept) tvalue-conditionthc tvalue-I(age - 3) tvalue-sexmale tvalue-weight_resid tvalue-conditionthc:I(age - 3) tvalue-conditionthc:sexmale
# 0.01           2.901372                  NA          2.613756             NA            3.306791                       3.722925                          NA
# 0.05           2.111322            3.956427          1.980399             NA            2.368218                       2.524727                          NA
# 0.1            1.740884            3.316921          1.658517             NA            1.910166                       1.962680                    4.108715
# 0.15           1.510320            2.868056          1.450803             NA            1.638579                       1.668639                    3.715055
# 0.2            1.336338            1.723101          1.290433             NA            1.438034                       1.460460                    2.285121
# tvalue-I(age - 3):sexmale tvalue-conditionthc:I(age - 3):sexmale
# 0.01                        NA                                     NA
# 0.05                        NA                                     NA
# 0.1                         NA                                     NA
# 0.15                        NA                                     NA
# 0.2                         NA                                     NA

mincWriteVolume(mlm, "/data/chamal/projects/lani/neonate/visuals/absolute_ageXcondition.mnc", column = 'tvalue-conditionthc:I(age - 3)')
mincWriteVolume(mlm, "/data/chamal/projects/lani/neonate/visuals/absolute_age.mnc", column = 'tvalue-I(age - 3)')
mincWriteVolume(mlm, "/data/chamal/projects/lani/neonate/visuals/absolute_condition.mnc", column = 'tvalue-conditionthc')
mincWriteVolume(mlm, "/data/chamal/projects/lani/neonate/visuals/absolute_weight_resid.mnc", column = 'tvalue-weight_resid')
mincWriteVolume(mlm, "/data/chamal/projects/lani/neonate/visuals/converged.mnc", column = "converged")
#relative jacobians

# > FDR
# Multidimensional MINC volume
# Columns:       qvalue-(Intercept) qvalue-conditionthc qvalue-I(age - 3) qvalue-sexmale qvalue-weight_resid qvalue-conditionthc:I(age - 3) qvalue-conditionthc:sexmale qvalue-I(age - 3):sexmale qvalue-conditionthc:I(age - 3):sexmale 
# [1] "/data/scratch2/cuplan/neonate/neo_thc/20230508_gabe_dbm_outputs/relative_smooth/MCH_NEO_THC_01_01_10_1-2-1_fix_lsq6_manual_extracted_fwhm_4vox.mnc"
# Degrees of Freedom: 155.999999996382 155.999999996304 155.999999996202 155.999999996341 155.999999995989 155.999999995927 155.999999996179 155.999999996047 155.999999996203 
# FDR Thresholds:
#   tvalue-(Intercept) tvalue-conditionthc tvalue-I(age - 3) tvalue-sexmale tvalue-weight_resid tvalue-conditionthc:I(age - 3) tvalue-conditionthc:sexmale
# 0.01           2.789325                  NA          2.733594             NA                  NA                       4.080589                          NA
# 0.05           2.133394            3.994845          2.087102             NA            3.410460                       3.323977                          NA
# 0.1            1.800196            3.495641          1.758135             NA            2.812698                       2.886403                          NA
# 0.15           1.583826            3.206469          1.543652             NA            2.457065                       2.564825                          NA
# 0.2            1.416863            2.984159          1.378279             NA            2.196660                       2.282591                          NA
# tvalue-I(age - 3):sexmale tvalue-conditionthc:I(age - 3):sexmale
# 0.01                        NA                                     NA
# 0.05                        NA                                     NA
# 0.1                         NA                                     NA
# 0.15                        NA                                     NA
# 0.2                         NA                                     NA

mincWriteVolume(mlm, "/data/chamal/projects/lani/neonate/visuals/relative_ageXcondition.mnc", column = 'tvalue-conditionthc:I(age - 3)')
mincWriteVolume(mlm, "/data/chamal/projects/lani/neonate/visuals/relative_age.mnc", column = 'tvalue-I(age - 3)')
mincWriteVolume(mlm, "/data/chamal/projects/lani/neonate/visuals/relative_condition.mnc", column = 'tvalue-conditionthc')
mincWriteVolume(mlm, "/data/chamal/projects/lani/neonate/visuals/relative_weight_resid.mnc", column = 'tvalue-weight_resid')
mincWriteVolume(mlm, "/data/chamal/projects/lani/neonate/visuals/converged.mnc", column = "converged")
# Visualize brain results----------------------------
#anatVol <- mincArray(mincGetVolume("/20210716_niftis_for_dbm/output/secondlevel_template0.mnc"))
#DBM

averagemask <- mincArray(mincGetVolume("/data/chamal/projects/lani/neonate/make_mask/eroded_extracted_mask.mnc"))

model =mlm
thresholds = attr(mincFDR(model), "thresholds")
#Here we figure out the extent of the mask, so we can use it to limit
# the FOV
# You may want to adjust these slighty buy adding/subtracting, so you don't
# quite go to the edge of the mask
#sagittal
dim1_begin <- min(which(averagemask == 1, arr.ind=TRUE)[,"dim1"])# +20
dim1_end <-max(which(averagemask == 1, arr.ind=TRUE)[,"dim1"]) #-15
# coronal
dim2_begin <- min(which(averagemask == 1, arr.ind=TRUE)[,"dim2"])# +10
dim2_end <- max(which(averagemask == 1, arr.ind=TRUE)[,"dim2"]) #-15
#Axial
dim3_begin <- min(which(averagemask == 1, arr.ind=TRUE)[,"dim3"])# +25
dim3_end <- max(which(averagemask == 1, arr.ind=TRUE)[,"dim3"]) #-15

anatVol <- mincArray(mincGetVolume("template_sharpen_shapeupdate.mnc"))[dim1_begin:dim1_end,dim2_begin:dim2_end, dim3_begin:dim3_end]

jacobian = "absolute"

predictors = "tvalue-conditionthc:I(age - 3)"
for (predictor in predictors){#dimnames(thresholds)[[2]][c(-1)]) {
 png(paste0("head_only_",jacobian,"_",predictor,".png"), units = "in", height = 5, width = 3.5, res = 300)
  variable_to_plot <- mincArray(model, predictor)[dim1_begin:dim1_end,dim2_begin:dim2_end, dim3_begin:dim3_end]
  begin1 = 20
  end1 = dim1_end - dim1_begin - 15
  begin2 = 20
  end2 = dim2_end - dim2_begin -15
  begin3 = 30
  end3 = dim3_end - dim3_begin -15
  tryCatch({
    lowerthreshold = round(thresholds["0.05",predictor],digits=2)
    #Sometimes, there isn't an 0.01 threshold, when thats the case, we use the max instead, be careful to read the threshold array printed above
    upperthreshold = round(ifelse(is.na(thresholds["0.01",predictor]), max(c(max(mincArray(model, predictor)),abs(min(mincArray(model, predictor))))), thresholds["0.01",predictor]),digits=2)
    # Here is the plotting code, we do all three slice directions in one figure. This was optimized for a human brain, you may need to
    # adjust the nrow and ncol to get exactly the figure you want
    # You will also need to adjust the anatVol low and high values to correspond to good thresholds for your template
    sliceSeries(nrow = 4, ncol = 2, dimension = 2, begin = begin2, end = end2) %>%
      anatomy(anatVol, low=1, high=30) %>%
      addtitle("Neonate") %>%
      overlay(variable_to_plot,
              low=lowerthreshold,
              high=upperthreshold,
              symmetric = T, alpha=0.6) %>%
      contourSliceIndicator(anatVol, c(1, 20)) %>%
      # sliceSeries(nrow = 7, ncol= 2, dimension = 1, begin = begin1, end = end1) %>%
      # anatomy(anatVol, low=1, high=30) %>%
      # addtitle("Sagittal") %>%
      # overlay(variable_to_plot,
      #         low=lowerthreshold,
      #         high=upperthreshold,
      #         symmetric = T, alpha=0.6) %>%
      # sliceSeries(nrow = 5, ncol= 2, dimension = 3, begin = begin3, end = end3) %>%
      # anatomy(anatVol, low=1, high=30) %>%
      # addtitle("Axial") %>%
      # overlay(variable_to_plot,
      #         low=lowerthreshold,
      #         high=upperthreshold,
      #         symmetric = T, alpha=0.6) %>%
      legend(predictor) %>%
      draw()}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  #If you are saving to file, also uncomment this
  dev.off()
}

# ------------------------------------------------ give it an atlas and extract the names of peak voxels.----------------------------------
# Load string with path to atlas definitions
defs = read.csv("/data/scratch2/cuplan/neonate/neo_thc/volume_from_dbm/atlas/P7_ECox_mapping_of_labels.csv")
# Load atlas volume
#Labels.
atlasLabel <- mincArray(mincGetVolume("/data/chamal/projects/lani/neonate/make_mask/labels_on_average.mnc"))
condition_age_linear<-mincArray(mincGetVolume("/data/chamal/projects/lani/neonate/visuals/absolute_ageXcondition.mnc"))
predictor = "tvalue-conditionthc:I(age - 3)"
lowerthreshold = round(thresholds["0.01",predictor],digits=2)
peaks <- mincFindPeaks(condition_age_linear,direction="both", threshold = lowerthreshold, minDistance = 0.5)

# Function to label peaks
label_peaks <- function(defs){
  #get information from the atlas to associate with the peaks (peaks need to be associated with a value to match with the Structure)
  peaks$label <- atlasLabel[as.matrix(peaks[,1:3])]
  peaks$label <-as.integer(peaks$label)
  #Recast data frame so the label name (side) become one column (variable) and the value becomes the other (value)
  mdefs <- defs[,1:3] %>% gather("variable","value", c("left.label", "right.label"))
  #make variable and structure characters
  mdefs$variable <- as.character(mdefs$variable)
  mdefs$Structure <- as.character(mdefs$Structure)
  #remove left and right
  mdefs$variable <- sub('.label','', mdefs$variable, fixed=T)
  #add left and right to the structure name. 
  mdefs$Structure <- paste(mdefs[,2], mdefs[,1])
  #add a 0 label for unlabelled regions
  mdefs <- rbind(mdefs, c("unlabelled","both",0))
  #turn value (the label) into a numeric
  mdefs$value <- as.numeric(mdefs$value)
  #go throught peaks, add a column labelling the peak. 
  for (i in 1:nrow(peaks)){
    if(any(peaks$label[i]==mdefs$value)){
      ars<-which(peaks$label[i]==mdefs$value)
      peaks$label[i]<-mdefs$Structure[ars]
    }
  }
  #define two functions that do...?
  numbers_only <- function(x) !grepl("\\D", x)
  numbers_only_neg <- function(x) grepl("^-.*", x)
  peaks$labels_missing_neg<-numbers_only_neg(peaks$label)
  peaks$labels_missing<-numbers_only(peaks$label)
  peaks$label2<-as.numeric(peaks$label)
  for (i in 1:nrow(peaks)){
    if (peaks$labels_missing[i]||peaks$labels_missing_neg[i]){
      x<-mdefs$value
      target <- which(abs(x - peaks$label2[i]) == min(abs(x - peaks$label2[i])))
      y<-x[max(target)]
      ars<-which(mdefs$value==max(y))
      peaks$label[i]<-mdefs$Structure[ars]
    }   
  }
  peaks<-select(peaks,-9:-11)
  return(peaks)
}
# Call function and output dataframe with labeled peaks
peaks<-label_peaks(defs)

# Save peak labels for future use as csv
write.csv(peaks, "/data/chamal/projects/lani/neonate/visuals/absolute_conditionXage_labelled_peaks.csv")

# ----------------------------- Plot Peak Voxels -------------------------------------------------
peaks <- read.csv("/data/chamal/projects/lani/neonate/visuals/absolute_conditionXage_labelled_peaks.csv")
two_times$voxel <-mincGetWorldVoxel(filenames = two_times$file.x, 1.3, 1.4, 3.2)

plot_peak_voxel <- function(Xw, Yw, Zw, region){
  two_times$voxel <-mincGetWorldVoxel(filenames = two_times$file.x , v1 = Xw, v2 = Yw, v3 = Zw)
  model_voxel = lmer(voxel ~ condition * I(age-3) +sex + (1|ID) + (1|litter_size), data = two_times)
  fitdata_voxel = as.data.frame(Effect(c("condition", "age"),model_voxel, xlevels=list(age=seq(3,12,1))))
  p <-  ggplot(data = two_times, aes(y=voxel,x=age, color=condition, group=condition)) +
    geom_point() +
    geom_line(data = fitdata_voxel, aes(y=fit),size = 2) +
    geom_ribbon(data = fitdata_voxel, aes(y=fit, ymin=lower, ymax=upper), alpha=0.2) +
    xlab("Age") +
    ylab(bquote('Log-Jacobians, Absolute')) + 
    labs(fill = "Condition", group="Condition", color="Condition", title = paste0(region))+
    scale_color_manual('Treatment', values = c('sal' = '#420085', 'thc' = '#018700'))+
    scale_size_manual(values=c(2,2))+
    theme(plot.title = element_text(hjust = 0.5,size = 16, face = "bold"))+
    theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), 
          legend.text=element_text(size=16), legend.title=element_text(size=18))+
    scale_x_continuous(breaks=seq(3,11,by=2)) +
    theme_classic()
  print(p)
  return(p)
}

#absolute jacobians, interaction
all_plots <- peaks %>% 
  rowwise() %>% 
  mutate(all_plots = list(plot_peak_voxel(x,y,z,label))) %>% 
  ungroup()

figure<-ggarrange(all_plots$all_plots[[2]] + rremove("ylab") + rremove("xlab"),all_plots$all_plots[[5]] + rremove("ylab") + rremove("xlab"),
                  all_plots$all_plots[[7]] + rremove("ylab") + rremove("xlab"),all_plots$all_plots[[17]] + rremove("ylab") + rremove("xlab"), # remove axis labels from plots
                 # all_plots$all_plots[[13]] + rremove("ylab") + rremove("xlab"),all_plots$all_plots[[23]] + rremove("ylab") + rremove("xlab"),
                #all_plots$all_plots[[29]] + rremove("ylab") + rremove("xlab"),all_plots$all_plots[[32]] + rremove("ylab") + rremove("xlab"),
                  ncol = 2, nrow = 2,
                  common.legend = TRUE, legend = "bottom",
                  font.label = list(size = 10, color = "black", face = "bold",  position = "top"))

figure_annotated<-annotate_figure(figure, left = text_grob("Log-Jacobian, Absolute", rot = 90, vjust = 1, size = 20, face = "bold"),
                                  bottom = text_grob("Age (days)", size = 20, face = "bold"))
png("/data/chamal/projects/lani/pte_paper/figures/neonate_peak_vox.png", units = "in", height = 5, width = 7.5, res = 300)
figure_annotated
dev.off()


#eof-----------------------------------------------------------

