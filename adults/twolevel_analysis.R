#Author: Lani Cupo
#Date: 20230918
#This script was written to analyze my adult data.

#run with runRneonates.sh

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
})

#---------------------------------------Set Up Environment-----------------------------------------------------
#Set working directory locally
#setting the wd determines whether you get full or relative jacs! 
setwd("/data/chamal/projects/lani/adult_thc/dbm_outputs/full_smooth/")
full_jac <- read.csv("../abs_jac.csv")
dems <- read.csv(file = "/data/chamal/projects/lani/adult_thc/dems_with_weight_resids.csv") #, header = TRUE)

data <- merge(full_jac,dems,by.x = c("ID", "timepoint"), by.y = c("ID", "scan_tp"))

#set mask
mask = "/data/chamal/projects/lani/adult_thc/make_mask/eroded_adult_mask.mnc"

#al data represented at least twice!
#two_times <- data %>% group_by(ID) %>% filter(n()>1)
#two_times <- droplevels(two_times)

sample_size <- data %>%
  distinct(data$ID, .keep_all=TRUE)

table(sample_size$treatment,sample_size$sex)

#choose whether or not to do a sex-specific models-------------------------
AIC_summary <-
  function(mmod)
    c(RMINC:::fixef_summary(mmod), AIC = extractAIC(mmod)[2])

mlm2 <-
  mincLmer(file ~ treatment * poly(scale(I(age-24)),2) * sex + scale(weight_resid)+ (1|ID) + (1|mom_id)
           , data
           , mask = mask)
           #, parallel=c("slurm",8))
           #, REML = F
           #, summary_type = AIC_summary)
saveRDS(mlm2, "adult_mod_full.rds")
mlm2 = mincLmerEstimateDF(mlm2)   
FDR = mincFDR(mlm2, mask= mask)
FDR

mlm <-
  mincLmer(file ~ treatment * poly(I(age-24),2) + sex + weight_resid+ (1|ID) + (1|mom_id)
           , data
           , mask = mask
           , parallel = c("sge", 200)
           , REML = F
           , summary_type = AIC_summary)

mincWriteVolume(mlm, "../../visuals/model1_AIC.mnc", column = 'AIC')
mincWriteVolume(mlm2, "../../visuals/model2_AIC.mnc", column = 'AIC')
mincWriteVolume(mlm2, "../../visuals/treat_sex_beta_scaled.mnc", column = "beta-treatmentTHC:sexmale")
mincWriteVolume(mlm2, "../../visuals/adult_thc.mnc", column = "tvalue-treatmentTHC")
mincWriteVolume(mlm2, "../../visuals/adult_thcXsex.mnc", column = "tvalue-treatmentTHC:sexmale")

mlm_den <- density(mlm[,"AIC"])
plot(mlm_den)
mlm2_den <- density(mlm2[,"AIC"])
plot(mlm2_den)

# Multidimensional MINC volume
# Columns:       qvalue-(Intercept) qvalue-treatmentTHC qvalue-poly(I(age - 24), 2)1 qvalue-poly(I(age - 24), 2)2 qvalue-sexmale qvalue-weight_resid qvalue-treatmentTHC:poly(I(age - 24), 2)1 qvalue-treatmentTHC:poly(I(age - 24), 2)2 qvalue-treatmentTHC:sexmale qvalue-poly(I(age - 24), 2)1:sexmale qvalue-poly(I(age - 24), 2)2:sexmale qvalue-treatmentTHC:poly(I(age - 24), 2)1:sexmale qvalue-treatmentTHC:poly(I(age - 24), 2)2:sexmale 
# [1] "sub-MCH_OLD_THC_32_01_ses-1_FLASH_lsq6_fwhm_4vox.mnc"
# Degrees of Freedom: 22.4055163956464 22.8226326709432 121.335195382056 125.999736076909 38.6069619440316 129.939340260602 123.446480899006 131.226054903755 38.9388809224157 122.959649461468 127.18521797812 122.099080585887 135.000027287197 
# FDR Thresholds:
#   tvalue-(Intercept) tvalue-treatmentTHC tvalue-poly(I(age - 24), 2)1 tvalue-poly(I(age - 24), 2)2 tvalue-sexmale tvalue-weight_resid
# 0.01           3.301942                  NA                           NA                     2.794593       4.002696            2.900143
# 0.05           2.430223            3.473147                           NA                     2.128835       3.005686            2.190246
# 0.1            2.021846            2.648180                           NA                     1.794475       2.534275            1.838021
# 0.15           1.765684            2.211156                           NA                     1.577356       2.225762            1.610599
# 0.2            1.571857            1.900103                           NA                     1.410732       1.987362            1.436536
# tvalue-treatmentTHC:poly(I(age - 24), 2)1 tvalue-treatmentTHC:poly(I(age - 24), 2)2 tvalue-treatmentTHC:sexmale tvalue-poly(I(age - 24), 2)1:sexmale
# 0.01                                        NA                                        NA                          NA                                   NA
# 0.05                                        NA                                        NA                    3.811104                                   NA
# 0.1                                         NA                                        NA                    3.305266                                   NA
# 0.15                                        NA                                        NA                    2.943021                                   NA
# 0.2                                         NA                                        NA                    2.647863                                   NA
# tvalue-poly(I(age - 24), 2)2:sexmale tvalue-treatmentTHC:poly(I(age - 24), 2)1:sexmale tvalue-treatmentTHC:poly(I(age - 24), 2)2:sexmale
# 0.01                             3.725654                                                NA                                                NA
# 0.05                             2.877506                                                NA                                                NA
# 0.1                              2.355118                                                NA                                                NA
# 0.15                             1.998775                                                NA                                                NA
# 0.2                              1.743136                                                NA                                                NA

#Relative Jacobians
# > FDR
# Multidimensional MINC volume
# Columns:       qvalue-(Intercept) qvalue-treatmentTHC qvalue-poly(I(age - 24), 2)1 qvalue-poly(I(age - 24), 2)2 qvalue-sexmale qvalue-weight_resid qvalue-treatmentTHC:poly(I(age - 24), 2)1 qvalue-treatmentTHC:poly(I(age - 24), 2)2 qvalue-treatmentTHC:sexmale qvalue-poly(I(age - 24), 2)1:sexmale qvalue-poly(I(age - 24), 2)2:sexmale qvalue-treatmentTHC:poly(I(age - 24), 2)1:sexmale qvalue-treatmentTHC:poly(I(age - 24), 2)2:sexmale 
# [1] "sub-MCH_OLD_THC_32_01_ses-1_FLASH_lsq6_fwhm_4vox.mnc"
# Degrees of Freedom: 23.7604497449454 23.9788809137044 122.417624913042 127.221058912389 39.0549808437112 132.707886917293 124.480456190313 133.269349741805 39.8175407810997 123.550572171365 128.230188287337 122.400136038421 135.651465792006 
# FDR Thresholds:
#   tvalue-(Intercept) tvalue-treatmentTHC tvalue-poly(I(age - 24), 2)1 tvalue-poly(I(age - 24), 2)2 tvalue-sexmale tvalue-weight_resid
# 0.01           3.465859                  NA                           NA                     2.862301       3.894094            2.958705
# 0.05           2.585685                  NA                           NA                     2.195914       2.826220            2.225936
# 0.1            2.167743                  NA                           NA                     1.857128       2.361969            1.866236
# 0.15           1.903797                  NA                     3.518083                     1.636704       2.077065            1.636208
# 0.2            1.704317                  NA                     3.175943                     1.465983       1.865829            1.461220
# tvalue-treatmentTHC:poly(I(age - 24), 2)1 tvalue-treatmentTHC:poly(I(age - 24), 2)2 tvalue-treatmentTHC:sexmale tvalue-poly(I(age - 24), 2)1:sexmale
# 0.01                                        NA                                        NA                          NA                                   NA
# 0.05                                        NA                                        NA                    3.886807                                   NA
# 0.1                                         NA                                        NA                    3.323360                                   NA
# 0.15                                        NA                                        NA                    2.976541                                   NA
# 0.2                                         NA                                        NA                    2.715981                                   NA
# tvalue-poly(I(age - 24), 2)2:sexmale tvalue-treatmentTHC:poly(I(age - 24), 2)1:sexmale tvalue-treatmentTHC:poly(I(age - 24), 2)2:sexmale
# 0.01                             3.994582                                                NA                                                NA
# 0.05                             3.416264                                                NA                                                NA
# 0.1                              3.101497                                                NA                                                NA
# 0.15                             2.873804                                                NA                                                NA
# 0.2                              2.689967                                                NA                                                NA

#Visualize results-----------------------------
anatVol <- mincArray(mincGetVolume("../template_sharpen_shapeupdate.mnc"))
averagemask <- mincArray(mincGetVolume("../../make_mask/eroded_adult_mask.mnc"))

model =mlm2
thresholds = attr(mincFDR(model), "thresholds")
#Here we figure out the extent of the mask, so we can use it to limit
# the FOV
# You may want to adjust these slighty buy adding/subtracting, so you don't
# quite go to the edge of the mask
#sagittal
dim1_begin <- min(which(averagemask == 1, arr.ind=TRUE)[,"dim1"]) +15
dim1_end <-max(which(averagemask == 1, arr.ind=TRUE)[,"dim1"]) -15
# coronal
dim2_begin <- min(which(averagemask == 1, arr.ind=TRUE)[,"dim2"]) +20
dim2_end <- max(which(averagemask == 1, arr.ind=TRUE)[,"dim2"]) -50
#Axial
dim3_begin <- min(which(averagemask == 1, arr.ind=TRUE)[,"dim3"]) +15
dim3_end <- max(which(averagemask == 1, arr.ind=TRUE)[,"dim3"]) -15

jacobian = "absolute"

for (predictor in c("tvalue-treatmentTHC")){#){# dimnames(thresholds)[[2]][c(-1)]) { #
  tryCatch({
   png(paste0("adult_",jacobian,"_",predictor,".png"), units = "in", height = 5, width = 3.5, res = 300)
    lowerthreshold = round(thresholds["0.1",predictor],digits=2)
    #lowerthreshold = 2.5
    #Sometimes, there isn't an 0.01 threshold, when thats the case, we use the max instead, be careful to read the threshold array printed above
    upperthreshold = round(ifelse(is.na(thresholds["0.01",predictor]), max(c(max(mincArray(model, predictor)),abs(min(mincArray(model, predictor))))), thresholds["0.01",predictor]),digits=2)
    # Here is the plotting code, we do all three slice directions in one figure. This was optimized for a human brain, you may need to
    # adjust the nrow and ncol to get exactly the figure you want
    # You will also need to adjust the anatVol low and high values to correspond to good thresholds for your template
    sliceSeries(nrow = 5, ncol = 2, dimension = 2, begin = dim2_begin, end = dim2_end) %>%  
      anatomy(anatVol, low=1, high=7) %>%
      addtitle("Adults") %>%
      overlay(mincArray(model, predictor),
              low=lowerthreshold, 
              high=upperthreshold, 
              symmetric = T, alpha=0.6) %>%
      contourSliceIndicator(anatVol, c(1, 20)) %>%
      #  sliceSeries(nrow = 7, ncol= 2, dimension = 1, begin = dim1_begin, end = dim1_end) %>%
      #  anatomy(anatVol, low=1, high=7) %>%
      #  addtitle("Sagittal") %>%
      # overlay(mincArray(model, predictor),
      #         low=lowerthreshold,
      #         high=upperthreshold,
      #         symmetric = T, alpha=0.6) %>%
      # sliceSeries(nrow = 5, ncol= 2, dimension = 3, begin = dim3_begin, end = dim3_end) %>%
      # anatomy(anatVol, low=1, high=7) %>%
      # addtitle("Axial") %>%
      # overlay(mincArray(model, predictor),
      #        low=lowerthreshold,
      #        high=upperthreshold,
      #        symmetric = T, alpha=0.6) %>%
      legend(predictor) %>%
      draw()}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  dev.off()
}

#from file, read in
model = mincArray(mincGetVolume("../../visuals/adult_thc.mnc"))
tFDR5 =2.67
tmax = 6.03

#png(paste0("/data/chamal/projects/lani/pte_paper/figures/adult_main_condition.png"), units = "in", height = 5, width = 3.5, res = 300)

sliceSeries(nrow = 8, ncol = 2, dimension = 2, begin = dim2_begin, end = dim2_end) %>%  
  anatomy(anatVol, low=1, high=7) %>%
  addtitle("Adults") %>%
  overlay(mincArray(model, predictor),
          low=tFDR5, 
          high=tmax, 
          symmetric = T, alpha=0.6) %>%
  sliceSeries(nrow = 7, ncol= 2, dimension = 1, begin = dim1_begin, end = dim1_end) %>%
  anatomy(anatVol, low=1, high=7) %>%
  addtitle("Sagittal") %>%
  overlay(mincArray(model, predictor),
          low=tFDR5,
          high=tmax,
          symmetric = T, alpha=0.6) %>%
  sliceSeries(nrow = 5, ncol= 2, dimension = 3, begin = dim3_begin, end = dim3_end) %>%
  anatomy(anatVol, low=1, high=7) %>%
  addtitle("Axial") %>%
  overlay(mincArray(model, predictor),
          low=tFDR5,
          high=tmax,
          symmetric = T, alpha=0.6) %>%
  legend("tvalue-treatmentTHC") %>%
  draw()
#dev.off()
# ------------------------------------------------ give it an atlas and extract the names of peak voxels.----------------------------------

model <-readRDS("adult_mod_full.rds")
# Load string with path to atlas definitions
defs = read.csv("../../make_mask/DSURQE_40micron_R_mapping.csv")
#Labels.
atlasLabel <- mincArray(mincGetVolume("../../make_mask/labels_on_average.mnc"))
find_peaks <- function(file_name, term){
  file <- mincArray((mincGetVolume(paste0("../../visuals/",file_name,".mnc"))))
  #file <- mlm2[,'tvalue-sexmale']
  #lowerthreshold = round(thresholds["0.05",term],digits=2)
  lowerthreshold = 3.47
  #peaks <- mincFindPeaks(inputStats = mlm2,column = term,direction="both", threshold = lowerthreshold, minDistance = 1)
  peaks <- mincFindPeaks(inputStats = file, direction = "both",threshold = lowerthreshold, minDistance = 1 )
  # Function to label peaks
  label_peaks <- function(defs){
    #get information from the atlas to associate with the peaks (peaks need to be associated with a value to match with the Structure)
    peaks$label <- as.integer(atlasLabel[as.matrix(peaks[,1:3])])
    #Recast data frame so the label name (side) become one column (variable) and the value becomes the other (value)
    mdefs <- defs[,1:3] %>% gather("variable","value", c("left.label", "right.label"))
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
    return(peaks)
  }
  # Call function and output dataframe with labeled peaks
  peaks<-label_peaks(defs)
  # write.csv(peaks, paste0("/data/chamal/projects/lani/neonate/visuals/",file_name,"_labelled_peaks.csv"))
  return(peaks)
}
con_peaks <- find_peaks("adult_thc","tvalue-treatmentTHC")
con_peaks_regions <- con_peaks %>% distinct(label, .keep_all = TRUE)
con_peaks_regions
write.csv(x = con_peaks_regions, file = "../../visuals/sexXage_quad_regions_0p05.csv")

#Plot peak voxel------------------------------------------------------------
data$voxel <-mincGetWorldVoxel(filenames = data$file , 1.3, 1.4, 3.2)
plot_peak_voxel <- function(Xw, Yw, Zw, region){
  data$voxel <-mincGetWorldVoxel(filenames = data$file , v1 = Xw, v2 = Yw, v3 = Zw)
  model_voxel = lmer(voxel ~ treatment * poly(I(age-24),2)*sex + (1|ID) + (1|mom_id), data = data)
  fitdata_voxel = as.data.frame(Effect(c("treatment", "age","sex"),model_voxel, xlevels=list(age=seq(20,95,5))))
  p <-  ggplot(data = data, aes(y=voxel,x=age, color=treatment, group=treatment)) +
    geom_point() +
    geom_line(data = fitdata_voxel, aes(y=fit),size = 2) +
    geom_ribbon(data = fitdata_voxel, aes(y=fit, ymin=lower, ymax=upper), alpha=0.2) +
    xlab("Age") +
    ylab(bquote('Log-Jacobians, Absolute')) + 
    theme_classic()+
    labs(fill = "Treatment", group="Treatment", color="Treatment", title = paste0( region))+
    scale_color_manual('Treatment', values = c('SAL' = '#420085', 'THC' = '#018700'))+
    scale_size_manual(values=c(2,2))+
    theme(plot.title = element_text(hjust = 0.5,size = 16))+
    theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), 
          legend.text=element_text(size=16), legend.title=element_text(size=18))+
    scale_x_continuous(breaks=c(20,40,60,90))+
    facet_wrap(~sex)
}

#absolute jacobians, interaction
all_plots <- con_peaks %>% 
  rowwise() %>% 
  mutate(all_plots = list(plot_peak_voxel(x,y,z,label))) %>% 
  ungroup()
figure<-ggarrange(all_plots$all_plots[[11]] + rremove("ylab") + rremove("xlab") + theme(axis.title.x=element_blank(),
                                                                                       axis.text.x=element_blank(),
                                                                                       axis.ticks.x=element_blank()),
                  all_plots$all_plots[[16]] + rremove("ylab") + rremove("xlab") + theme(axis.title.x=element_blank(),
                                                                                        axis.text.x=element_blank(),
                                                                                        axis.ticks.x=element_blank()),
                  all_plots$all_plots[[15]] + rremove("ylab") + rremove("xlab")+ theme(
                                                                                      axis.ticks.x=element_blank(),
                                                                                      strip.background = element_blank(),
                                                                                      strip.text.x = element_blank()),
                  all_plots$all_plots[[7]] + rremove("ylab") + rremove("xlab") + theme(
                                                                                       strip.background = element_blank(),
                                                                                        strip.text.x = element_blank()), # remove axis labels from plots
                  ncol = 2, nrow = 2,
                  common.legend = TRUE, legend = "bottom",
                  font.label = list(size = 10, color = "black", face = "bold",  position = "top"))

figure_annotated<-annotate_figure(figure, left = text_grob("Log-Jacobian, Absolute", rot = 90, vjust = 1, size = 20, face = "bold"),
                                  bottom = text_grob("Age (days)", size = 20, face = "bold"))
png("/data/chamal/projects/lani/pte_paper/figures/adult_peak_vox.png", units = "in", height = 5, width = 7.5, res = 300)
figure_annotated
dev.off()
#Quick check if mom 52's pups are weird---------------------
plot_peak_voxel1 <- function(Xw, Yw, Zw, region){
  data$voxel <-mincGetWorldVoxel(filenames = data$file , v1 = Xw, v2 = Yw, v3 = Zw)
  p <-  ggplot(data = data, aes(y=voxel,x=age, color=factor(mom_id), group=factor(mom_id))) +
    geom_point() +
    theme(plot.title = element_text(hjust = 0.5,size = 16, face = "bold"))+
    theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), 
          legend.text=element_text(size=16), legend.title=element_text(size=18))+
    scale_x_continuous(breaks=c(20,40,60,90))+
    facet_wrap(~sex)
}
a <- plot_peak_voxel1(-0.47,  7.61, 1.5,"left Olfactory bulb")

#Maybe output volumes for PLS? ----------------------

pls_ppi <- read.csv("../../behavior/derived_data/ppi/avg_trials_ppi.csv")
pls_oft <-read.csv("../../behavior/derived_data/oft/master_oft.csv")
pls_oft <- pls_oft %>%
  group_by(ID) %>%
  mutate(total_distance_moved = sum(Distance_moved)) %>%
  ungroup()

pls_oft <- subset(pls_oft, pls_oft$Zone == "4")
pls_behav <- merge(pls_ppi, pls_oft, by.x = "Subject", by.y = "ID")
pls_dems <- merge(timepoint2, pls_behav, by.x = "ID", by.y = "Subject")
pls_dems <- pls_dems[c("ID","weight","Start_mid","pp3","pp6","pp12","pp9","pp15","Frequency","Duration","total_distance_moved")]

timepoint2 <- subset(data, data$timepoint == "2")
timepoint2 <- subset(timepoint2, timepoint2$ID %in% pls_dems$ID)
t2_volumes <- anatGetAll(timepoint2$file, atlas = "../../make_mask/labels_on_average.mnc", defs = "../../make_mask/DSURQE_40micron_R_mapping.csv")

pls_dems <- pls_dems[c("weight","Start_mid","pp3","pp6","pp12","pp9","pp15","Frequency","Duration","total_distance_moved")]
write.csv(t2_volumes, "../../pls/t2_volumes.csv")
write.csv(pls_dems, "../../pls/t2_demographics.csv")
write.csv(names(pls_dems), "../../pls/voxelwise_pls/behavior_vars.csv")

#eof