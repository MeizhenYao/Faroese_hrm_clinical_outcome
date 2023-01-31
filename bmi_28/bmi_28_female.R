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
library(ggrepel)
library(wesanderson)


###################################################################################################################
###################################################################################################################
###################################################################################################################


new_quantile <- function(x, cuts ){
  
  y <- x[x!= 0 & !is.na(x)]
  qi <- unique(quantile(y, probs = seq(0, 1, by = 1/cuts), na.rm = TRUE))
  
  if(length(qi) == 1){ 
    qi = c(-Inf, qi)
  } else{ 
    qi[1] <- -Inf
    qi[length(qi)] <- Inf
  }
  
  x[which(x!= 0 & !is.na(x))] = cut(x[x!= 0 & !is.na(x)], breaks = qi, labels = FALSE, include.lowest = TRUE)
  
  return(x)
  
}

###################################################################################################################
###################################################################################################################
###################################################################################################################

data_HILIC <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/data_hilic.csv", check.names = F)
data_C18 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/data_c18.csv", check.names = F)

merged_omics_HILIC <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/merged_omics_hilic.csv", check.names = F)
merged_omics_C18 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese//hrm_clinical_outcome/merged_omics_c18.csv", check.names = F)

###################################################################################################################
###################################################################################################################
###################################################################################################################
merged_omics_HILIC$bmi_28_female<- ifelse(merged_omics_HILIC$bmi_28 < 25, 0, 1)
merged_omics_C18$bmi_28_female<- ifelse(merged_omics_C18$bmi_28 < 25, 0, 1)
# restrict to sex=0
merged_omics_HILIC<- merged_omics_HILIC %>% 
                     filter(sex==0)
merged_omics_C18<- merged_omics_C18 %>% 
                     filter(sex==0)
# bmi_28 and met at age 7

## HILIC

merged_omics_HILIC_age7 <- merged_omics_HILIC[merged_omics_HILIC$Year == "7",]
qqt <- as.data.frame(apply(merged_omics_HILIC_age7[,colnames(merged_omics_HILIC_age7) %in% paste0("Met",seq(1:nrow(data_HILIC)))], 2, function(x) new_quantile(x, cuts = 10)))
data.bmi_28_female.met_age7 = cbind(qqt[,c(paste0("Met",seq(1:nrow(data_HILIC))))], merged_omics_HILIC_age7[,c( 'mage',  'mbmi', 'age7', "bmi_28_female")])
data.bmi_28_female.met_age7$beta <- rep(NA_real_, nrow(data.bmi_28_female.met_age7))
data.bmi_28_female.met_age7$model_pval <- rep(NA_real_, nrow(data.bmi_28_female.met_age7))
d_lm_status <- data.frame(Met_id = NA_character_, Beta = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)


## fit model between clinical outcomes and metabolites
for(i in 1:1991){
  s_lm <- (glm(bmi_28_female ~ data.bmi_28_female.met_age7[,i] + mage + age7, family = binomial(), data = data.bmi_28_female.met_age7))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = exp(coef(s_lm))
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  d_lm_status <- rbind(d_lm_status,c(colnames(data.bmi_28_female.met_age7)[i], as.numeric(r.est[2,c(1,2,3,4)])))
}


### create data for Metaboanalyst and plot
d_lm_status <- d_lm_status[-1,]
d_lm_status$z.value <- as.numeric(d_lm_status$z.value)
d_lm_status$p.value <- as.numeric(d_lm_status$p.value)
data.bmi_28_female.met_age7 <- cbind(d_lm_status, data_HILIC[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Class")])

mumm_bmi_28_female.met_age7 <- data.bmi_28_female.met_age7[,c("mz","p.value","z.value","time")]
colnames(mumm_bmi_28_female.met_age7) <- c("m.z", "p.value", "t.score", "rt")
mumm_bmi_28_female.met_age7 <- mumm_bmi_28_female.met_age7[order(mumm_bmi_28_female.met_age7$p.value),]
write.table(mumm_bmi_28_female.met_age7, sep = "\t",
            "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/mummichog_HILIC_bmi_28_female_age7.txt",
            row.names = F, col.names = T)


write.csv(data.bmi_28_female.met_age7,
          "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_HILIC_bmi_28_female_age7.csv",
          row.names = F)

### calculate p.value based on number of effect test
exwas_HILIC_bmi_28_female_age7 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_HILIC_bmi_28_female_age7.csv")
exwas_HILIC_bmi_28_female_age7$Mode <- rep("HILIC",nrow(exwas_HILIC_bmi_28_female_age7))

ght <- cor(cbind(merged_omics_HILIC_age7[,paste0("Met",seq(1:nrow(data_HILIC)))]))
et <- eigen(ght)
p_HILIC_bmi_28_female_age7<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
p_HILIC_bmi_28_female_age7
# 0.0005197505

exwas_HILIC_bmi_28_female_age7$Met_id[exwas_HILIC_bmi_28_female_age7$p.value < p_HILIC_bmi_28_female_age7]
# character(0)


###################################################################################################################

## C18

merged_omics_C18_age7 <- merged_omics_C18[merged_omics_C18$Year == "7",]
qqt <- as.data.frame(apply(merged_omics_C18_age7[,colnames(merged_omics_C18_age7) %in% paste0("Met",seq(1:nrow(data_C18)))], 2, function(x) new_quantile(x, cuts = 10)))
data.bmi_28_female.met_age7 = cbind(qqt[,c(paste0("Met",seq(1:nrow(data_C18))))], merged_omics_C18_age7[,c( 'mage',  'mbmi', 'age7', "bmi_28_female")])
data.bmi_28_female.met_age7$beta <- rep(NA_real_, nrow(data.bmi_28_female.met_age7))
data.bmi_28_female.met_age7$model_pval <- rep(NA_real_, nrow(data.bmi_28_female.met_age7))
d_lm_status <- data.frame(Met_id = NA_character_, Beta = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:787){
  s_lm <- (glm(bmi_28_female ~ data.bmi_28_female.met_age7[,i] +  mage + age7 ,family = binomial(),  data = data.bmi_28_female.met_age7))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = exp(coef(s_lm))
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  d_lm_status <- rbind(d_lm_status,c(colnames(data.bmi_28_female.met_age7)[i], as.numeric(r.est[2,c(1,2,3,4)])))
}

d_lm_status <- d_lm_status[-1,]
d_lm_status$z.value <- as.numeric(d_lm_status$z.value)
d_lm_status$p.value <- as.numeric(d_lm_status$p.value)
data.bmi_28_female.met_age7 <- cbind(d_lm_status, data_C18[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Class")])

mumm_bmi_28_female.met_age7 <- data.bmi_28_female.met_age7[,c("mz","p.value","z.value","time")]
colnames(mumm_bmi_28_female.met_age7) <- c("m.z", "p.value", "t.score", "rt")
mumm_bmi_28_female.met_age7 <- mumm_bmi_28_female.met_age7[order(mumm_bmi_28_female.met_age7$p.value),]
write.table(mumm_bmi_28_female.met_age7, sep = "\t",
            "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/mummichog_C18_bmi_28_female_age7.txt",
            row.names = F, col.names = T)


write.csv(data.bmi_28_female.met_age7,
          "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_C18_bmi_28_female_age7.csv",
          row.names = F)

exwas_C18_bmi_28_female_age7 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_C18_bmi_28_female_age7.csv")
exwas_C18_bmi_28_female_age7$Mode <- rep("C18",nrow(exwas_C18_bmi_28_female_age7))

ght <- cor(cbind(merged_omics_C18_age7[,paste0("Met",seq(1:nrow(data_C18)))]))
et <- eigen(ght)
p_C18_bmi_28_female_age7<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# 0.001388889

exwas_C18_bmi_28_female_age7$Met_id[exwas_C18_bmi_28_female_age7$p.value < p_C18_bmi_28_female_age7]
# "(0)


###################################################################################################################
###################################################################################################################
###################################################################################################################

# bmi_28_female and met at age 14

## HILIC

merged_omics_HILIC_age14 <- merged_omics_HILIC[merged_omics_HILIC$Year == "14",]
qqt <- as.data.frame(apply(merged_omics_HILIC_age14[,colnames(merged_omics_HILIC_age14) %in% paste0("Met",seq(1:nrow(data_HILIC)))], 2, function(x) new_quantile(x, cuts = 10)))
data.bmi_28_female.met_age14 = cbind(qqt[,c(paste0("Met",seq(1:nrow(data_HILIC))))], merged_omics_HILIC_age14[,c( 'mage',  'mbmi', 'age14', "bmi_28_female")])
data.bmi_28_female.met_age14$beta <- rep(NA_real_, nrow(data.bmi_28_female.met_age14))
data.bmi_28_female.met_age14$model_pval <- rep(NA_real_, nrow(data.bmi_28_female.met_age14))
d_lm_status <- data.frame(Met_id = NA_character_, Beta = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:1991){
  s_lm <- (glm(bmi_28_female ~ data.bmi_28_female.met_age14[,i] +  mage + age14 ,family = binomial(),  data = data.bmi_28_female.met_age14))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = exp(coef(s_lm))
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  d_lm_status <- rbind(d_lm_status,c(colnames(data.bmi_28_female.met_age14)[i], as.numeric(r.est[2,c(1,2,3,4)])))
}

d_lm_status <- d_lm_status[-1,]
d_lm_status$z.value <- as.numeric(d_lm_status$z.value)
d_lm_status$p.value <- as.numeric(d_lm_status$p.value)
data.bmi_28_female.met_age14 <- cbind(d_lm_status, data_HILIC[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Class")])

mumm_bmi_28_female.met_age14 <- data.bmi_28_female.met_age14[,c("mz","p.value","z.value","time")]
colnames(mumm_bmi_28_female.met_age14) <- c("m.z", "p.value", "t.score", "rt")
mumm_bmi_28_female.met_age14 <- mumm_bmi_28_female.met_age14[order(mumm_bmi_28_female.met_age14$p.value),]
write.table(mumm_bmi_28_female.met_age14, sep = "\t",
            "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/mummichog_HILIC_bmi_28_female_age14.txt",
            row.names = F, col.names = T)


write.csv(data.bmi_28_female.met_age14,
          "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_HILIC_bmi_28_female_age14.csv",
          row.names = F)

exwas_HILIC_bmi_28_female_age14 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_HILIC_bmi_28_female_age14.csv")
exwas_HILIC_bmi_28_female_age14$Mode <- rep("HILIC",nrow(exwas_HILIC_bmi_28_female_age14))

ght <- cor(cbind(merged_omics_HILIC_age14[,paste0("Met",seq(1:nrow(data_HILIC)))]))
et <- eigen(ght)
p_HILIC_bmi_28_female_age14<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# 0.0005197505

exwas_HILIC_bmi_28_female_age14$Met_id[exwas_HILIC_bmi_28_female_age14$p.value < p_HILIC_bmi_28_female_age14]
# "Met640"


###################################################################################################################

## C18

merged_omics_C18_age14 <- merged_omics_C18[merged_omics_C18$Year == "14",]
qqt <- as.data.frame(apply(merged_omics_C18_age14[,colnames(merged_omics_C18_age14) %in% paste0("Met",seq(1:nrow(data_C18)))], 2, function(x) new_quantile(x, cuts = 10)))
data.bmi_28_female.met_age14 = cbind(qqt[,c(paste0("Met",seq(1:nrow(data_C18))))], merged_omics_C18_age14[,c( 'mage',  'mbmi', 'age14', "bmi_28_female")])
data.bmi_28_female.met_age14$beta <- rep(NA_real_, nrow(data.bmi_28_female.met_age14))
data.bmi_28_female.met_age14$model_pval <- rep(NA_real_, nrow(data.bmi_28_female.met_age14))
d_lm_status <- data.frame(Met_id = NA_character_, Beta = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:787){
  s_lm <- (glm(bmi_28_female ~ data.bmi_28_female.met_age14[,i] +  mage + age14 , family = binomial(), data = data.bmi_28_female.met_age14))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = exp(coef(s_lm))
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  d_lm_status <- rbind(d_lm_status,c(colnames(data.bmi_28_female.met_age14)[i], as.numeric(r.est[2,c(1,2,3,4)])))
}

d_lm_status <- d_lm_status[-1,]
d_lm_status$z.value <- as.numeric(d_lm_status$z.value)
d_lm_status$p.value <- as.numeric(d_lm_status$p.value)
data.bmi_28_female.met_age14 <- cbind(d_lm_status, data_C18[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Class")])

mumm_bmi_28_female.met_age14 <- data.bmi_28_female.met_age14[,c("mz","p.value","z.value","time")]
colnames(mumm_bmi_28_female.met_age14) <- c("m.z", "p.value", "t.score", "rt")
mumm_bmi_28_female.met_age14 <- mumm_bmi_28_female.met_age14[order(mumm_bmi_28_female.met_age14$p.value),]
write.table(mumm_bmi_28_female.met_age14, sep = "\t",
            "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/mummichog_C18_bmi_28_female_age14.txt",
            row.names = F, col.names = T)


write.csv(data.bmi_28_female.met_age14,
          "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_C18_bmi_28_female_age14.csv",
          row.names = F)

exwas_C18_bmi_28_female_age14 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_C18_bmi_28_female_age14.csv")
exwas_C18_bmi_28_female_age14$Mode <- rep("C18",nrow(exwas_C18_bmi_28_female_age14))

ght <- cor(cbind(merged_omics_C18_age14[,paste0("Met",seq(1:nrow(data_C18)))]))
et <- eigen(ght)
p_C18_bmi_28_female_age14<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
#  0.001388889

exwas_C18_bmi_28_female_age14$Met_id[exwas_C18_bmi_28_female_age14$p.value < p_C18_bmi_28_female_age14]
#  "Met366"


###################################################################################################################
###################################################################################################################
###################################################################################################################

# bmi_28_female and met at age 22

## HILIC

merged_omics_HILIC_age22 <- merged_omics_HILIC[merged_omics_HILIC$Year == "22",]
qqt <- as.data.frame(apply(merged_omics_HILIC_age22[,colnames(merged_omics_HILIC_age22) %in% paste0("Met",seq(1:nrow(data_HILIC)))], 2, function(x) new_quantile(x, cuts = 10)))
data.bmi_28_female.met_age22 = cbind(qqt[,c(paste0("Met",seq(1:nrow(data_HILIC))))], merged_omics_HILIC_age22[,c( 'mage',  'mbmi', 'age22', "bmi_28_female")])
data.bmi_28_female.met_age22$beta <- rep(NA_real_, nrow(data.bmi_28_female.met_age22))
data.bmi_28_female.met_age22$model_pval <- rep(NA_real_, nrow(data.bmi_28_female.met_age22))
d_lm_status <- data.frame(Met_id = NA_character_, Beta = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:1991){
  s_lm <- (glm(bmi_28_female ~ data.bmi_28_female.met_age22[,i] +  mage + age22 ,family = binomial(),  data = data.bmi_28_female.met_age22))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = exp(coef(s_lm))
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  d_lm_status <- rbind(d_lm_status,c(colnames(data.bmi_28_female.met_age22)[i], as.numeric(r.est[2,c(1,2,3,4)])))
}

d_lm_status <- d_lm_status[-1,]
d_lm_status$z.value <- as.numeric(d_lm_status$z.value)
d_lm_status$p.value <- as.numeric(d_lm_status$p.value)
data.bmi_28_female.met_age22 <- cbind(d_lm_status, data_HILIC[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Class")])

mumm_bmi_28_female.met_age22 <- data.bmi_28_female.met_age22[,c("mz","p.value","z.value","time")]
colnames(mumm_bmi_28_female.met_age22) <- c("m.z", "p.value", "t.score", "rt")
mumm_bmi_28_female.met_age22 <- mumm_bmi_28_female.met_age22[order(mumm_bmi_28_female.met_age22$p.value),]
write.table(mumm_bmi_28_female.met_age22, sep = "\t",
            "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/mummichog_HILIC_bmi_28_female_age22.txt",
            row.names = F, col.names = T)


write.csv(data.bmi_28_female.met_age22,
          "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_HILIC_bmi_28_female_age22.csv",
          row.names = F)

exwas_HILIC_bmi_28_female_age22 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_HILIC_bmi_28_female_age22.csv")
exwas_HILIC_bmi_28_female_age22$Mode <- rep("HILIC",nrow(exwas_HILIC_bmi_28_female_age22))

ght <- cor(cbind(merged_omics_HILIC_age22[,paste0("Met",seq(1:nrow(data_HILIC)))]))
et <- eigen(ght)
p_HILIC_bmi_28_female_age22<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# 0.0005197505

exwas_HILIC_bmi_28_female_age22$Met_id[exwas_HILIC_bmi_28_female_age22$p.value < p_HILIC_bmi_28_female_age22]
# character(0)


###################################################################################################################

## C18

merged_omics_C18_age22 <- merged_omics_C18[merged_omics_C18$Year == "22",]
qqt <- as.data.frame(apply(merged_omics_C18_age22[,colnames(merged_omics_C18_age22) %in% paste0("Met",seq(1:nrow(data_C18)))], 2, function(x) new_quantile(x, cuts = 10)))
data.bmi_28_female.met_age22 = cbind(qqt[,c(paste0("Met",seq(1:nrow(data_C18))))], merged_omics_C18_age22[,c( 'mage',  'mbmi', 'age22', "bmi_28_female")])
data.bmi_28_female.met_age22$beta <- rep(NA_real_, nrow(data.bmi_28_female.met_age22))
data.bmi_28_female.met_age22$model_pval <- rep(NA_real_, nrow(data.bmi_28_female.met_age22))
d_lm_status <- data.frame(Met_id = NA_character_, Beta = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:787){
  s_lm <- (glm(bmi_28_female ~ data.bmi_28_female.met_age22[,i] +  mage + age22 ,family = binomial(),  data = data.bmi_28_female.met_age22))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = exp(coef(s_lm))
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  d_lm_status <- rbind(d_lm_status,c(colnames(data.bmi_28_female.met_age22)[i], as.numeric(r.est[2,c(1,2,3,4)])))
}

d_lm_status <- d_lm_status[-1,]
d_lm_status$z.value <- as.numeric(d_lm_status$z.value)
d_lm_status$p.value <- as.numeric(d_lm_status$p.value)
data.bmi_28_female.met_age22 <- cbind(d_lm_status, data_C18[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Class")])

mumm_bmi_28_female.met_age22 <- data.bmi_28_female.met_age22[,c("mz","p.value","z.value","time")]
colnames(mumm_bmi_28_female.met_age22) <- c("m.z", "p.value", "t.score", "rt")
mumm_bmi_28_female.met_age22 <- mumm_bmi_28_female.met_age22[order(mumm_bmi_28_female.met_age22$p.value),]
write.table(mumm_bmi_28_female.met_age22, sep = "\t",
            "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/mummichog_C18_bmi_28_female_age22.txt",
            row.names = F, col.names = T)


write.csv(data.bmi_28_female.met_age22,
          "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_C18_bmi_28_female_age22.csv",
          row.names = F)

exwas_C18_bmi_28_female_age22 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_C18_bmi_28_female_age22.csv")
exwas_C18_bmi_28_female_age22$Mode <- rep("C18",nrow(exwas_C18_bmi_28_female_age22))

ght <- cor(cbind(merged_omics_C18_age22[,paste0("Met",seq(1:nrow(data_C18)))]))
et <- eigen(ght)
p_C18_bmi_28_female_age22<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# 0.001388889

exwas_C18_bmi_28_female_age22$Met_id[exwas_C18_bmi_28_female_age22$p.value < p_C18_bmi_28_female_age22]
# "Met128"


###################################################################################################################
###################################################################################################################
###################################################################################################################

# bmi_28_female and met at age 28

## HILIC

merged_omics_HILIC_age28 <- merged_omics_HILIC[merged_omics_HILIC$Year == "28",]
qqt <- as.data.frame(apply(merged_omics_HILIC_age28[,colnames(merged_omics_HILIC_age28) %in% paste0("Met",seq(1:nrow(data_HILIC)))], 2, function(x) new_quantile(x, cuts = 10)))
data.bmi_28_female.met_age28 = cbind(qqt[,c(paste0("Met",seq(1:nrow(data_HILIC))))], merged_omics_HILIC_age28[,c( 'mage',  'mbmi', 'age28', "bmi_28_female")])
data.bmi_28_female.met_age28$beta <- rep(NA_real_, nrow(data.bmi_28_female.met_age28))
data.bmi_28_female.met_age28$model_pval <- rep(NA_real_, nrow(data.bmi_28_female.met_age28))
d_lm_status <- data.frame(Met_id = NA_character_, Beta = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:1991){
  s_lm <- (glm(bmi_28_female ~ data.bmi_28_female.met_age28[,i] +  mage + age28 , family = binomial(), data = data.bmi_28_female.met_age28))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = exp(coef(s_lm))
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  d_lm_status <- rbind(d_lm_status,c(colnames(data.bmi_28_female.met_age28)[i], as.numeric(r.est[2,c(1,2,3,4)])))
}

d_lm_status <- d_lm_status[-1,]
d_lm_status$z.value <- as.numeric(d_lm_status$z.value)
d_lm_status$p.value <- as.numeric(d_lm_status$p.value)
data.bmi_28_female.met_age28 <- cbind(d_lm_status, data_HILIC[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Class")])

mumm_bmi_28_female.met_age28 <- data.bmi_28_female.met_age28[,c("mz","p.value","z.value","time")]
colnames(mumm_bmi_28_female.met_age28) <- c("m.z", "p.value", "t.score", "rt")
mumm_bmi_28_female.met_age28 <- mumm_bmi_28_female.met_age28[order(mumm_bmi_28_female.met_age28$p.value),]
write.table(mumm_bmi_28_female.met_age28, sep = "\t",
            "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/mummichog_HILIC_bmi_28_female_age28.txt",
            row.names = F, col.names = T)


write.csv(data.bmi_28_female.met_age28,
          "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_HILIC_bmi_28_female_age28.csv",
          row.names = F)

exwas_HILIC_bmi_28_female_age28 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_HILIC_bmi_28_female_age28.csv")
exwas_HILIC_bmi_28_female_age28$Mode <- rep("HILIC",nrow(exwas_HILIC_bmi_28_female_age28))

ght <- cor(cbind(merged_omics_HILIC_age28[,paste0("Met",seq(1:nrow(data_HILIC)))]))
et <- eigen(ght)
p_HILIC_bmi_28_female_age28<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# 0.0005197505

exwas_HILIC_bmi_28_female_age28$Met_id[exwas_HILIC_bmi_28_female_age28$p.value < p_HILIC_bmi_28_female_age28]
# character(0)


###################################################################################################################

## C18

merged_omics_C18_age28 <- merged_omics_C18[merged_omics_C18$Year == "28",]
qqt <- as.data.frame(apply(merged_omics_C18_age28[,colnames(merged_omics_C18_age28) %in% paste0("Met",seq(1:nrow(data_C18)))], 2, function(x) new_quantile(x, cuts = 10)))
data.bmi_28_female.met_age28 = cbind(qqt[,c(paste0("Met",seq(1:nrow(data_C18))))], merged_omics_C18_age28[,c( 'sex','mage',  'mbmi', 'age28', "bmi_28_female")])
data.bmi_28_female.met_age28$beta <- rep(NA_real_, nrow(data.bmi_28_female.met_age28))
data.bmi_28_female.met_age28$model_pval <- rep(NA_real_, nrow(data.bmi_28_female.met_age28))
d_lm_status <- data.frame(Met_id = NA_character_, Beta = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:787){
  s_lm <- (glm(bmi_28_female ~ data.bmi_28_female.met_age28[,i] +  mage + age28 ,family = binomial(),  data = data.bmi_28_female.met_age28))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = exp(coef(s_lm))
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  d_lm_status <- rbind(d_lm_status,c(colnames(data.bmi_28_female.met_age28)[i], as.numeric(r.est[2,c(1,2,3,4)])))
}

d_lm_status <- d_lm_status[-1,]
d_lm_status$z.value <- as.numeric(d_lm_status$z.value)
d_lm_status$p.value <- as.numeric(d_lm_status$p.value)
data.bmi_28_female.met_age28 <- cbind(d_lm_status, data_C18[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Class")])

mumm_bmi_28_female.met_age28 <- data.bmi_28_female.met_age28[,c("mz","p.value","z.value","time")]
colnames(mumm_bmi_28_female.met_age28) <- c("m.z", "p.value", "t.score", "rt")
mumm_bmi_28_female.met_age28 <- mumm_bmi_28_female.met_age28[order(mumm_bmi_28_female.met_age28$p.value),]
write.table(mumm_bmi_28_female.met_age28, sep = "\t",
            "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/mummichog_C18_bmi_28_female_age28.txt",
            row.names = F, col.names = T)


write.csv(data.bmi_28_female.met_age28,
          "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_C18_bmi_28_female_age28.csv",
          row.names = F)

exwas_C18_bmi_28_female_age28 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_C18_bmi_28_female_age28.csv")
exwas_C18_bmi_28_female_age28$Mode <- rep("C18",nrow(exwas_C18_bmi_28_female_age28))

ght <- cor(cbind(merged_omics_C18_age28[,paste0("Met",seq(1:nrow(data_C18)))]))
et <- eigen(ght)
p_C18_bmi_28_female_age28<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# 0.001388889

exwas_C18_bmi_28_female_age28$Met_id[exwas_C18_bmi_28_female_age28$p.value < p_C18_bmi_28_female_age28]
#  character(0)


################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

# Plots
## At age7

exwas_C18_bmi_28_female_age7 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_C18_bmi_28_female_age7.csv")
exwas_C18_bmi_28_female_age7$Mode <- rep("C18",nrow(exwas_C18_bmi_28_female_age7))
exwas_C18_bmi_28_female_age7 <- exwas_C18_bmi_28_female_age7[exwas_C18_bmi_28_female_age7$p.value < p_C18_bmi_28_female_age7,]
exwas_C18_bmi_28_female_age7$age <- rep("Age 7", nrow(exwas_C18_bmi_28_female_age7))

exwas_HILIC_bmi_28_female_age7 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_HILIC_bmi_28_female_age7.csv")
exwas_HILIC_bmi_28_female_age7$Mode <- rep("HILIC",nrow(exwas_HILIC_bmi_28_female_age7))
exwas_HILIC_bmi_28_female_age7 <- exwas_HILIC_bmi_28_female_age7[exwas_HILIC_bmi_28_female_age7$p.value < p_HILIC_bmi_28_female_age7,]
exwas_HILIC_bmi_28_female_age7$age <- rep("Age 7", nrow(exwas_HILIC_bmi_28_female_age7))

## At age14

exwas_C18_bmi_28_female_age14 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_C18_bmi_28_female_age14.csv")
exwas_C18_bmi_28_female_age14$Mode <- rep("C18",nrow(exwas_C18_bmi_28_female_age14))
exwas_C18_bmi_28_female_age14 <- exwas_C18_bmi_28_female_age14[exwas_C18_bmi_28_female_age14$p.value < p_C18_bmi_28_female_age14,]
exwas_C18_bmi_28_female_age14$age <- rep("Age 14", nrow(exwas_C18_bmi_28_female_age14))


exwas_HILIC_bmi_28_female_age14 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_HILIC_bmi_28_female_age14.csv")
exwas_HILIC_bmi_28_female_age14$Mode <- rep("HILIC",nrow(exwas_HILIC_bmi_28_female_age14))
exwas_HILIC_bmi_28_female_age14 <- exwas_HILIC_bmi_28_female_age14[exwas_HILIC_bmi_28_female_age14$p.value < p_HILIC_bmi_28_female_age14,]
exwas_HILIC_bmi_28_female_age14$age <- rep("Age 14", nrow(exwas_HILIC_bmi_28_female_age14))

## At age22

exwas_C18_bmi_28_female_age22 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_C18_bmi_28_female_age22.csv")
exwas_C18_bmi_28_female_age22$Mode <- rep("C18",nrow(exwas_C18_bmi_28_female_age22))
exwas_C18_bmi_28_female_age22 <- exwas_C18_bmi_28_female_age22[exwas_C18_bmi_28_female_age22$p.value < p_C18_bmi_28_female_age22,]
exwas_C18_bmi_28_female_age22$age <- rep("Age 22", nrow(exwas_C18_bmi_28_female_age22))


exwas_HILIC_bmi_28_female_age22 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_HILIC_bmi_28_female_age22.csv")
exwas_HILIC_bmi_28_female_age22$Mode <- rep("HILIC",nrow(exwas_HILIC_bmi_28_female_age22))
exwas_HILIC_bmi_28_female_age22 <- exwas_HILIC_bmi_28_female_age22[exwas_HILIC_bmi_28_female_age22$p.value < p_HILIC_bmi_28_female_age22,]
exwas_HILIC_bmi_28_female_age22$age <- rep("Age 22", nrow(exwas_HILIC_bmi_28_female_age22))

## At age28

exwas_C18_bmi_28_female_age28 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_C18_bmi_28_female_age28.csv")
exwas_C18_bmi_28_female_age28$Mode <- rep("C18",nrow(exwas_C18_bmi_28_female_age28))
exwas_C18_bmi_28_female_age28 <- exwas_C18_bmi_28_female_age28[exwas_C18_bmi_28_female_age28$p.value < p_C18_bmi_28_female_age28,]
exwas_C18_bmi_28_female_age28$age <- rep("Age 28", nrow(exwas_C18_bmi_28_female_age28))


exwas_HILIC_bmi_28_female_age28 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_HILIC_bmi_28_female_age28.csv")
exwas_HILIC_bmi_28_female_age28$Mode <- rep("HILIC",nrow(exwas_HILIC_bmi_28_female_age28))
exwas_HILIC_bmi_28_female_age28 <- exwas_HILIC_bmi_28_female_age28[exwas_HILIC_bmi_28_female_age28$p.value < p_HILIC_bmi_28_female_age28,]
exwas_HILIC_bmi_28_female_age28$age <- rep("Age 28", nrow(exwas_HILIC_bmi_28_female_age28))


##########################################################################################
##########################################################################################
##########################################################################################

# Naming the metabolites and Corresponding adducts

exwas_bmi_28_female <- rbind(exwas_HILIC_bmi_28_female_age7,exwas_C18_bmi_28_female_age7,
                         exwas_HILIC_bmi_28_female_age14,exwas_C18_bmi_28_female_age14,
                         exwas_HILIC_bmi_28_female_age22,exwas_C18_bmi_28_female_age22,
                         exwas_HILIC_bmi_28_female_age28,exwas_C18_bmi_28_female_age28)


conf_met <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/confirmed_metabolites.csv",fileEncoding = "Latin1")
conf_met_HILIC <- conf_met[conf_met$Method.RT == "HILIC+",]
conf_met_C18 <- conf_met[conf_met$Method.RT == "C18-",]

stage4_hilic <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/Stage4.csv")
stage4_c18 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/Stage4.csv")

all_sig_hits <- exwas_bmi_28_female
all_sig_hits$chemical_ID <- rep(NA_character_, dim(all_sig_hits)[1])
all_sig_hits$Annotation.confidence.score <- rep(NA_character_, dim(all_sig_hits)[1])
all_sig_hits$Name <- rep(NA_character_, dim(all_sig_hits)[1])
all_sig_hits$Adduct <- rep(NA_character_, dim(all_sig_hits)[1])

for(i in 1:nrow(all_sig_hits)){
  
  if(all_sig_hits$Mode[i] == "HILIC"){
    
    rrt <- which(abs(all_sig_hits$mz[i] - stage4_hilic$mz) < 5 & abs(all_sig_hits$time[i] - stage4_hilic$time) < 30)
    dft <- stage4_hilic[unique(rrt[!is.na(rrt)]),]
    b <- dft[dft$Formula %in% intersect(conf_met_HILIC$Elemental.Composition, unique(dft$Formula)),]
    
    if(nrow(b)!= 0){
      
      score <- max(stage4_hilic$Annotation.confidence.score[stage4_hilic$chemical_ID %in% b$chemical_ID])
      dft2 <- stage4_hilic[stage4_hilic$chemical_ID %in% b$chemical_ID,c("chemical_ID","Annotation.confidence.score","Name","Adduct")]
      all_sig_hits$chemical_ID[i] = knitr::combine_words(unique(dft2$chemical_ID[dft2$Annotation.confidence.score == score]), and = "", sep = "/")
      all_sig_hits$Annotation.confidence.score[i] = knitr::combine_words(unique(dft2$Annotation.confidence.score[dft2$Annotation.confidence.score == score]), and = "", sep = "/")
      all_sig_hits$Name[i] = knitr::combine_words(unique(dft2$Name[dft2$Annotation.confidence.score == score]), and = "", sep = "/")
      
      names <- unlist(strsplit(all_sig_hits$Name[i],"/"))
      tp <- NA_character_
      for(j in 1:length(names)){
        
        tp <- c(tp,unique(dft2$Adduct[dft2$Name == names[j]]))
        
      }
      tp <- tp[-1]
      all_sig_hits$Adduct[i] = knitr::combine_words(tp, and = "", sep = "/")
      
    }
  }
  
  else if(all_sig_hits$Mode[i] == "C18"){
    
    rrt <- which(abs(all_sig_hits$mz[i] - stage4_c18$mz) < 5 & abs(all_sig_hits$time[i] - stage4_c18$time) < 30)
    dft <- stage4_c18[unique(rrt[!is.na(rrt)]),]
    b <- dft[dft$Formula %in% intersect(conf_met_C18$Elemental.Composition, unique(dft$Formula)),]
    
    if(nrow(b)!= 0){
      
      score <- max(stage4_c18$Annotation.confidence.score[stage4_c18$chemical_ID %in% b$chemical_ID])
      dft2 <- stage4_c18[stage4_c18$chemical_ID %in% b$chemical_ID,c("chemical_ID","Annotation.confidence.score","Name","Adduct")]
      all_sig_hits$chemical_ID[i] = knitr::combine_words(unique(dft2$chemical_ID[dft2$Annotation.confidence.score == score]), and = "", sep = "/")
      all_sig_hits$Annotation.confidence.score[i] = knitr::combine_words(unique(dft2$Annotation.confidence.score[dft2$Annotation.confidence.score == score]), and = "", sep = "/")
      all_sig_hits$Name[i] = knitr::combine_words(unique(dft2$Name[dft2$Annotation.confidence.score == score]), and = "", sep = "/")
      
      names <- unlist(strsplit(all_sig_hits$Name[i],"/"))
      tp <- NA_character_
      for(j in 1:length(names)){
        
        tp <- c(tp,unique(dft2$Adduct[dft2$Name == names[j]]))
        
      }
      tp <- tp[-1]
      all_sig_hits$Adduct[i] = knitr::combine_words(tp, and = "", sep = "/")
      
      
    }
  }
  
}



all_sig_hit <- all_sig_hits[, c("Met_id", "Beta" , "Std.Error", "z.value" , "p.value" ,"mz" ,
                                "time" , "Annotation.confidence.score",
                                "Class" , "Mode", 
                                "age" , "chemical_ID" ,"Name",
                                "Adduct")]
write.csv(all_sig_hits, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_allsig_metabolites_bmi_28_female.csv", row.names = F )

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
# significant Metabolites plot
# Without the ones at age 28

exwas_bmi_28_female <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/exwas_allsig_metabolites_bmi_28_female.csv")
exwas_bmi_28_female <- exwas_bmi_28_female[exwas_bmi_28_female$age!= "Age 28",]

exwas_bmi_28_female$log10.p.value <- -log(exwas_bmi_28_female$p.value,10)
exwas_bmi_28_female$age <- factor(exwas_bmi_28_female$age, levels = c('Age 14', 'Age 22'))

exwas_bmi_28_female$Metabolite <- rep(NA_character_,nrow(exwas_bmi_28_female))

for(i in 1:nrow(exwas_bmi_28_female)) exwas_bmi_28_female$Metabolite[i] <- strsplit(exwas_bmi_28_female$Name,"/")[[i]][1]

exwas_bmi_28_female<- exwas_bmi_28_female %>% # if there are duplicates, only keep one metabolite with lowest p.value at each timepoint
  arrange(age, p.value)%>%
  group_by(age) %>% 
  distinct(Metabolite, .keep_all = TRUE) 


pal <- wes_palette("Zissou1", 100, type = "continuous")
vol <- (ggplot(exwas_bmi_28_female, aes(x=age, y=Beta, color=log10.p.value, label=Metabolite, shape = Mode)) +# Show all points
          geom_point(size = 3)+  
          scale_shape_manual(values = c(15,17)) + # change shape
          geom_hline(yintercept= 1, color = "black", size = 1 )   +
          scale_color_gradientn(colours = pal) + # change color
          labs(x = "Age at metabolomic assessment", y = "OR", 
               title = "(A) Prospective metabolic exposures and Overweight at age 28 (Female)",
               col = "-log10(p.value)") +
          theme_bw()+
          geom_label_repel(size = 6, family = 'serif',
                           fontface = 'bold',
                           box.padding = unit(0.5, "lines"),
                           max.iter = 2e4,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                           force = 2, force_pull = 2, show.legend = F)) # add label for each point 

vol <-  (vol + theme(plot.title=element_text(size=16,face="bold"),
                     axis.title=element_text(size=14,face="bold"),
                     plot.tag = element_text(size = 14,face = "bold"),
                     axis.text.y = element_text(size=14,face="bold"),
                     axis.text.x = element_text(size=14,face="bold"),
                     strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
                     legend.title = element_text(face = "bold"), legend.position = 'bottom',
                     legend.text = element_text(size = 12),
                     legend.background = element_rect(fill="white", 
                                                      linewidth=0.5, linetype="solid",  colour ="darkblue"),
                     plot.margin=unit(c(0.5,0.5,1,1.5), "cm")) +
           ylim(c(-0.5,2.5)) +
           annotate("rect", xmin = c(0.75, 1.75), 
                    xmax = c(1.25, 2.25), 
                    ymin = c(-0.5, -0.5), 
                    ymax = c(2.5, 2.5),
                    alpha = 0.05))


jpeg("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/plot_allsig_metabolites_bmi_28_female.jpeg",
     units="in", width=15, height=10, res=600)

vol


dev.off()


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

# Pathway Plot

C18_bmi_28_female_age7  <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/metaboanalyst/mummichog_pathway_enrichment_C18_bmi_28_female_age7.csv")
C18_bmi_28_female_age7 <- C18_bmi_28_female_age7[C18_bmi_28_female_age7$Hits.sig >=3 ,]
C18_bmi_28_female_age7$age <- rep(7, nrow(C18_bmi_28_female_age7))

C18_bmi_28_female_age14  <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/metaboanalyst/mummichog_pathway_enrichment_C18_bmi_28_female_age14.csv")
C18_bmi_28_female_age14 <- C18_bmi_28_female_age14[C18_bmi_28_female_age14$Hits.sig >=3 ,]
C18_bmi_28_female_age14$age <- rep(14, nrow(C18_bmi_28_female_age14))

C18_bmi_28_female_age22  <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/metaboanalyst/mummichog_pathway_enrichment_C18_bmi_28_female_age22.csv")
C18_bmi_28_female_age22 <- C18_bmi_28_female_age22[C18_bmi_28_female_age22$Hits.sig >=3 ,]
C18_bmi_28_female_age22$age <- rep(22, nrow(C18_bmi_28_female_age22))

C18_bmi_28_female <- rbind(C18_bmi_28_female_age7,C18_bmi_28_female_age14, C18_bmi_28_female_age22)
C18_bmi_28_female$Mode <- rep("C18-", nrow(C18_bmi_28_female))


HILIC_bmi_28_female_age7  <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/metaboanalyst/mummichog_pathway_enrichment_HILIC_bmi_28_female_age7.csv")
HILIC_bmi_28_female_age7 <- HILIC_bmi_28_female_age7[HILIC_bmi_28_female_age7$Hits.sig >=3 ,]
HILIC_bmi_28_female_age7$age <- rep(7, nrow(HILIC_bmi_28_female_age7))

HILIC_bmi_28_female_age14  <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/metaboanalyst/mummichog_pathway_enrichment_HILIC_bmi_28_female_age14.csv")
HILIC_bmi_28_female_age14 <- HILIC_bmi_28_female_age14[HILIC_bmi_28_female_age14$Hits.sig >=3 ,]
HILIC_bmi_28_female_age14$age <- rep(14, nrow(HILIC_bmi_28_female_age14))

HILIC_bmi_28_female_age22  <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_female/metaboanalyst/mummichog_pathway_enrichment_HILIC_bmi_28_female_age22.csv")
HILIC_bmi_28_female_age22 <- HILIC_bmi_28_female_age22[HILIC_bmi_28_female_age22$Hits.sig >=3 ,]
HILIC_bmi_28_female_age22$age <- rep(22, nrow(HILIC_bmi_28_female_age22))

HILIC_bmi_28_female <- rbind(HILIC_bmi_28_female_age7,HILIC_bmi_28_female_age14, HILIC_bmi_28_female_age22)
HILIC_bmi_28_female$Mode <- rep("HILIC+", nrow(HILIC_bmi_28_female))

pathways_bmi_28_female <- rbind(C18_bmi_28_female, HILIC_bmi_28_female)
colnames(pathways_bmi_28_female)[1] <- "Pathways"
pathways_bmi_28_female$age <- factor(pathways_bmi_28_female$age, levels = c("7","14"))

pfas_pathway_subset <- pathways_bmi_28_female[pathways_bmi_28_female$FET < 0.05,]

all_pathways <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/all_pathways.csv")
pfas_pathway_subset <- merge(pfas_pathway_subset, all_pathways, by.x = "Pathways")


pfas_pathway_subset$Group <- factor(pfas_pathway_subset$Group, levels = c("Metabolism of amino acids and derivatives", 
                                                                          "Xenobiotics biodegradation and metabolism"
))

groups <- c("Metabolism of amino acids and derivatives", 
            "Xenobiotics biodegradation and metabolism")

pfas_pathway_subset <- pfas_pathway_subset[order(pfas_pathway_subset$Group),]
uni_paths <- pfas_pathway_subset$Pathways
uni_paths <- uni_paths[!duplicated(uni_paths)]
pfas_pathway_subset$Pathways <- factor(pfas_pathway_subset$Pathways, levels = uni_paths)

ret <- c(length(unique(pfas_pathway_subset$Pathways[pfas_pathway_subset$Group ==  unique(pfas_pathway_subset$Group)[1]])),
         length(unique(pfas_pathway_subset$Pathways[pfas_pathway_subset$Group ==  unique(pfas_pathway_subset$Group)[2]]))
)

gp <- (ggplot(pfas_pathway_subset, aes(x = age,  y = Pathways)) + 
         geom_point(aes(size = -log(FET , base = 10), shape = Mode), 
                    alpha  = 5, 
                    position = position_dodge2(width = 0.5, preserve  = "total", padding = 0.5, reverse = T)) + 
         theme_bw()   + scale_shape_manual(values = c(15,17)) +
         ylab(NULL) + xlab("Age at Metabolomic Assessment") +
         ggtitle("(B) Enriched pathways in association between\nserum Metabolites and Overweight at age 28 (Female)") + 
         labs(size = "-log10(p.value)") +
         theme(plot.title=element_text(size=16,face="bold"), 
               axis.title=element_text(size=14,face="bold"),
               plot.tag = element_text(size = 13,face = "bold"),
               axis.text.y = element_text(size=15,face="bold"),
               axis.text.x = element_text(size=12,face="bold"),
               strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
               legend.title = element_text(face = "bold"), legend.position = 'right',
               legend.text = element_text(size = 12),
               legend.background = element_rect(fill="white", 
                                                size=0.5, linetype="solid",  colour ="darkblue"),
               plot.margin=unit(c(0.5,0.5,1,0.5), "cm")) +
         annotate("rect", xmin = c(0.75,1.75), # change x groups
                   xmax = c(1.25,2.25),
                   ymin = rep(0,2), ymax =rep(4.5,2),
                   alpha = .1)+  
        guides(color = guide_legend(override.aes = list(size = 4)),
                shape = guide_legend(override.aes = list(size = 4))))

gp <- (gp +  
        geom_hline(yintercept = c(cumsum(ret))+0.5)+
        geom_hline(yintercept = 4.5, size = 2)+
        annotate("text", x = rep(4.9, length(unique(pfas_pathway_subset$Group))), 
                  y = c(1.8,  4), # adjust position of red label
                  label = c("Metabolism of amino acids\nand derivatives", 
                            "Xenobiotics biodegradation\nand metabolism"), color = "red",
                  size = rep(5, length(groups)), fontface  = 'bold')+
        coord_cartesian(xlim = c(1.2, 6), ylim = c(0.5, 4), clip = "off") )

jpeg("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese//hrm_clinical_outcome/bmi/bmi_28_female/pathways_bmi_28_female.jpeg",
     units="in", width=12, height=10, res=600)

gp

dev.off()

##################################################

# combine plots

jpeg("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese//hrm_clinical_outcome/bmi/bmi_28_female/plot_metabolites_pathways_bmi_28_female.jpeg",
     units="in", width=28, height=15, res=800)

ggpubr::ggarrange(vol, gp, nrow = 1, ncol = 2, common.legend = F, 
                  legend = "bottom") + theme_bw()


dev.off()
