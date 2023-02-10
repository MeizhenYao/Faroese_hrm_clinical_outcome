library(car)
library(readr)
library(lattice)
library(nlme)
library(ggplot2)
library(GGally)
library(nnet)
library(foreign)
library(biotools)
library(glmmML)
library(MASS)
library(lme4)
library(multcomp)
library(dplyr)
library(qwraps2)
library(knitr)
library(xtable)
library(kableExtra)
library(DT)
library(glmnet)
library(corrplot)
library(ggpubr)
library("EnvStats")
library(lmerTest)
library("merTools")
library(reshape2)
library(ggplot2)
library(GGally)
library(mgcv)
library(gplots)
library(tidyr)
library(cluster)
library(factoextra)
library(psych)
library(xgboost)
library(SHAPforxgboost)
library(bkmr)
library(mi)
library(stargazer)
library(factoextra) 
library(spatstat)
library(Hmisc)
library(robustHD)
library(gWQS)
library(gridExtra)
library(ggcorrplot)
library(BWQS)
library(qwraps2)
library(mgcv)
library(ggcorrplot)
library("MatchIt")
library(cobalt)
library(WeightIt)
library(ggalluvial)
library(mice)
library(data.table)
library(omu)
library(renv)
require(foreign)
require(sandwich)


setwd("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/git/Faroese_PFAS_vs_Metabolites/HILIC")
d <- as.data.frame(read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/New_faroese/HILIC/17.12.19_December2020.csv"))
d$id <- as.character(d$id)
d1 <- as.data.frame(read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/New_faroese/HILIC/eligible_ids.csv"))
d1$id <- as.character(as.numeric(as.character(d1$id)))
y <- rep(NA_real_, length(d1$id))
for(i in 1:length(d1$id)){
  y[i] = sum(d$id %in% d1$id[i])
}
dt <- d[which(d$id %in% d1$id),]
d <- dt[,c("id",'pfos_0', 'pfos_7', 'pfos_14', 'pfos_22', 'pfos_28', 'pfoa_0', 'pfoa_7', 'pfoa_14', 'pfoa_22', 'pfoa_28', 
           'pfhxs_0', 'pfhxs_7', 'pfhxs_14', 'pfhxs_22', 'pfhxs_28',  'pfna_0', 'pfna_7', 'pfna_14', 'pfna_22', 'pfna_28', 
           'pfda_0', 'pfda_7', 'pfda_14', 'pfda_22', 'pfda_28', 'bmi_28', 'waistcirc_28', 'diaBP_dxt_28', 'sysBP_dxt_28', 
           'uhdl_28y', 'trig_28y', 'dldl_28y', 'chol_28y','cir_28', 'insulinauc_28y', 'glucoseauc_28y', 'matsuda_28', 
           'homair_28', 'igi_28', "sex","mage","mbmi","parity","smokepreg_2","matfishpreg_cat2","breastfed_tot",
           "age_7","age_14","age22","age28","diaBP_dxt_28","IDF_metS")]

d$id <- as.character(d$id)



# Impute Data using MICE

init = mice(d, maxit=0) 
meth = init$method
predM = init$predictorMatrix
set.seed(1234)
imputed = mice(d, method="norm", predictorMatrix=predM, m=3, ntree = 100)
imputed <- mice::complete(imputed, action = 2) # this is your final imputed dataset


d_subset <- imputed[,c("id",'pfos_0', 'pfos_7', 'pfos_14', 'pfos_22', 'pfos_28', 'pfoa_0', 'pfoa_7', 'pfoa_14', 'pfoa_22', 
                       'pfoa_28', 'pfhxs_0', 'pfhxs_7', 'pfhxs_14', 'pfhxs_22', 'pfhxs_28',  'pfna_0', 'pfna_7', 'pfna_14', 
                       'pfna_22', 'pfna_28', 'pfda_0', 'pfda_7', 'pfda_14', 'pfda_22', 'pfda_28', 'bmi_28', 'waistcirc_28', 
                       'diaBP_dxt_28', 'sysBP_dxt_28', 'uhdl_28y', 'trig_28y', 'dldl_28y', 'chol_28y','cir_28', 'insulinauc_28y', 
                       'glucoseauc_28y', 'matsuda_28', 'homair_28', 'igi_28', "sex","mage","mbmi","parity","smokepreg_2",
                       "matfishpreg_cat2","breastfed_tot",
                       "age_7","age_14","age22","age28")]

d_long <- reshape(d_subset, idvar = "id", direction = "long", timevar = "Year", 
                  varying = list(c('pfos_0', 'pfos_7', 'pfos_14', 'pfos_22', 'pfos_28'), 
                                 c('pfoa_0', 'pfoa_7', 'pfoa_14', 'pfoa_22', 'pfoa_28'),
                                 c('pfhxs_0', 'pfhxs_7', 'pfhxs_14', 'pfhxs_22', 'pfhxs_28'),
                                 c('pfna_0', 'pfna_7', 'pfna_14', 'pfna_22', 'pfna_28'),
                                 c('pfda_0', 'pfda_7', 'pfda_14', 'pfda_22', 'pfda_28')),
                  v.names = c("PFOS","PFOA","PFHxS","PFNA","PFDA"))

d_long$Year <- (d_long$Year - 1)*7
d_long$Year <- ifelse(d_long$Year == 21, 22, d_long$Year)
d_long$Year <- factor(d_long$Year, levels = c("0","7","14","22","28"))



## BMI_28

## waistcirc_28: (male:40, female:35)
merged_omics_HILIC$waist <- merged_omics_HILIC$waistcirc_28
for (i in 1:nrow(merged_omics_HILIC)){
  if (merged_omics_HILIC$sex[i]==0 & merged_omics_HILIC$waistcirc_28[i] > 35*2.54){merged_omics_HILIC$waist[i]=1}
  else if(merged_omics_HILIC$sex[i]==0 & merged_omics_HILIC$waistcirc_28[i] <= 35*2.54){merged_omics_HILIC$waist[i]=0}
  else if (merged_omics_HILIC$sex[i]==1 & merged_omics_HILIC$waistcirc_28[i] > 40*2.54){merged_omics_HILIC$waist[i]=1}
  else if(merged_omics_HILIC$sex[i]==1 & merged_omics_HILIC$waistcirc_28[i] <= 40*2.54){merged_omics_HILIC$waist[i]=0}
}


## (merged_omics$hypertension <- as.numeric(merged_omics$sysBP_dxt_28 > 130 & merged_omics$diaBP_dxt_28 > 85) 
bmi
bmi_bi
# waist/height ratio

##--------------------- metabolic syndrome:yes or no
metabolic syndrome 
## HLD
## LDL
## TRI
## CHOL

##--------------------- maybe pre-hypertension
merged_omics_HILIC$pre_hypertension <- as.numeric(merged_omics_HILIC$sysBP_dxt_28 > 120 & merged_omics_HILIC$diaBP_dxt_28 > 80) 
# 
# # https://www.nhlbi.nih.gov/health/high-blood-pressure
# 
# merged_omics$metabolic_disorder <- 

##--------------------- continuous
"insulinauc_28y"   
"glucoseauc_28y"  
"matsuda_28"       
"homair_28"        
"igi_28"


##--------------------- metabolic syndrome:yes or no
merged_omics_HILIC_age7$mets <- merged_omics_HILIC_age7$trig_28y
for (i in 1:nrow(merged_omics_HILIC_age7)){
  if (merged_omics_HILIC_age7$sex[i]==0 & merged_omics_HILIC_age7$uhdl_28y[i] < 1.293 & merged_omics_HILIC_age7$trig_28y[i] >= 1.6935){merged_omics_HILIC_age7$mets[i]="mixed"}
  else if (merged_omics_HILIC_age7$sex[i]==0 & merged_omics_HILIC_age7$uhdl_28y[i] < 1.293 & merged_omics_HILIC_age7$trig_28y[i] < 1.6935){merged_omics_HILIC_age7$mets[i]="hypo"}
  else if (merged_omics_HILIC_age7$sex[i]==0 & merged_omics_HILIC_age7$uhdl_28y[i] > 1.293 & merged_omics_HILIC_age7$trig_28y[i] >= 1.6935){merged_omics_HILIC_age7$mets[i]="hyper"}
  else if (merged_omics_HILIC_age7$sex[i]==0 & merged_omics_HILIC_age7$uhdl_28y[i] > 1.293 & merged_omics_HILIC_age7$trig_28y[i] < 1.6935){merged_omics_HILIC_age7$mets[i]="no"}
  else if (merged_omics_HILIC_age7$sex[i]==1 & merged_omics_HILIC_age7$uhdl_28y[i] < 1.0344 & merged_omics_HILIC_age7$trig_28y[i] >= 1.6935){merged_omics_HILIC_age7$mets[i]="mixed"}
  else if (merged_omics_HILIC_age7$sex[i]==1 & merged_omics_HILIC_age7$uhdl_28y[i] < 1.0344 & merged_omics_HILIC_age7$trig_28y[i] < 1.6935){merged_omics_HILIC_age7$mets[i]="hypo"}
  else if (merged_omics_HILIC_age7$sex[i]==1 & merged_omics_HILIC_age7$uhdl_28y[i] > 1.0344 & merged_omics_HILIC_age7$trig_28y[i] >= 1.6935){merged_omics_HILIC_age7$mets[i]="hyper"}
  else if (merged_omics_HILIC_age7$sex[i]==1 & merged_omics_HILIC_age7$uhdl_28y[i] > 1.0344 & merged_omics_HILIC_age7$trig_28y[i] < 1.6935){merged_omics_HILIC_age7$mets[i]="no"}
}





