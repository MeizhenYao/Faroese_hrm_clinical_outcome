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

## HILIC

merged_omics_HILIC_age7 <- merged_omics_HILIC[merged_omics_HILIC$Year == "7",]
qqt <- as.data.frame(apply(merged_omics_HILIC_age7[,colnames(merged_omics_HILIC_age7) %in% paste0("Met",seq(1:nrow(data_HILIC)))], 2, function(x) new_quantile(x, cuts = 10)))
data.homair_28.met_age7 = cbind(qqt[,c(paste0("Met",seq(1:nrow(data_HILIC))))], merged_omics_HILIC_age7[,c('sex', 'mage',  'mbmi', 'age7', "homair_28")])
data.homair_28.met_age7$beta <- rep(NA_real_, nrow(data.homair_28.met_age7))
data.homair_28.met_age7$model_pval <- rep(NA_real_, nrow(data.homair_28.met_age7))
d_lm_status <- data.frame(Met_id = NA_character_, Beta = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)


## fit model between clinical outcomes and metabolites
for(i in 1:1991){
  s_lm <- (lm(homair_28 ~ data.homair_28.met_age7[,i] + sex + age7 , data = data.homair_28.met_age7))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  d_lm_status <- rbind(d_lm_status,c(colnames(data.homair_28.met_age7)[i], as.numeric(r.est[2,c(1,2,3,4)])))
}


### create data for Metaboanalyst and plot
d_lm_status <- d_lm_status[-1,]
d_lm_status$z.value <- as.numeric(d_lm_status$z.value)
d_lm_status$p.value <- as.numeric(d_lm_status$p.value)
data.homair_28.met_age7 <- cbind(d_lm_status, data_HILIC[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Class")])

mumm_homair_28.met_age7 <- data.homair_28.met_age7[,c("mz","p.value","z.value","time")]
colnames(mumm_homair_28.met_age7) <- c("m.z", "p.value", "t.score", "rt")
mumm_homair_28.met_age7 <- mumm_homair_28.met_age7[order(mumm_homair_28.met_age7$p.value),]
write.table(mumm_homair_28.met_age7, sep = "\t",
            "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/mummichog_HILIC_homair_28_age7.txt",
            row.names = F, col.names = T)


write.csv(data.homair_28.met_age7,
          "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_HILIC_homair_28_age7.csv",
          row.names = F)

### calculate p.value based on number of effect test
exwas_HILIC_homair_28_age7 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_HILIC_homair_28_age7.csv")
exwas_HILIC_homair_28_age7$Mode <- rep("HILIC",nrow(exwas_HILIC_homair_28_age7))

ght <- cor(cbind(merged_omics_HILIC_age7[,paste0("Met",seq(1:nrow(data_HILIC)))]))
et <- eigen(ght)
1/ (sum((et$values>1 + 0)* (et$values - 1)))
# 0.0005356186

exwas_HILIC_homair_28_age7$Met_id[exwas_HILIC_homair_28_age7$p.value < 0.0005356186]
# character(0)


###################################################################################################################

## C18

merged_omics_C18_age7 <- merged_omics_C18[merged_omics_C18$Year == "7",]
qqt <- as.data.frame(apply(merged_omics_C18_age7[,colnames(merged_omics_C18_age7) %in% paste0("Met",seq(1:nrow(data_C18)))], 2, function(x) new_quantile(x, cuts = 10)))
data.homair_28.met_age7 = cbind(qqt[,c(paste0("Met",seq(1:nrow(data_C18))))], merged_omics_C18_age7[,c('sex', 'mage',  'mbmi', 'age7', "homair_28")])
data.homair_28.met_age7$beta <- rep(NA_real_, nrow(data.homair_28.met_age7))
data.homair_28.met_age7$model_pval <- rep(NA_real_, nrow(data.homair_28.met_age7))
d_lm_status <- data.frame(Met_id = NA_character_, Beta = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:787){
  s_lm <- (lm(homair_28 ~ data.homair_28.met_age7[,i] + sex + age7 , data = data.homair_28.met_age7))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  d_lm_status <- rbind(d_lm_status,c(colnames(data.homair_28.met_age7)[i], as.numeric(r.est[2,c(1,2,3,4)])))
}

d_lm_status <- d_lm_status[-1,]
d_lm_status$z.value <- as.numeric(d_lm_status$z.value)
d_lm_status$p.value <- as.numeric(d_lm_status$p.value)
data.homair_28.met_age7 <- cbind(d_lm_status, data_C18[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Class")])

mumm_homair_28.met_age7 <- data.homair_28.met_age7[,c("mz","p.value","z.value","time")]
colnames(mumm_homair_28.met_age7) <- c("m.z", "p.value", "t.score", "rt")
mumm_homair_28.met_age7 <- mumm_homair_28.met_age7[order(mumm_homair_28.met_age7$p.value),]
write.table(mumm_homair_28.met_age7, sep = "\t",
            "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/mummichog_C18_homair_28_age7.txt",
            row.names = F, col.names = T)


write.csv(data.homair_28.met_age7,
          "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_C18_homair_28_age7.csv",
          row.names = F)

exwas_C18_homair_28_age7 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_C18_homair_28_age7.csv")
exwas_C18_homair_28_age7$Mode <- rep("C18",nrow(exwas_C18_homair_28_age7))

ght <- cor(cbind(merged_omics_C18_age7[,paste0("Met",seq(1:nrow(data_C18)))]))
et <- eigen(ght)
1/ (sum((et$values>1 + 0)* (et$values - 1)))
# 0.001507876

exwas_C18_homair_28_age7$Met_id[exwas_C18_homair_28_age7$p.value < 0.001507876]
# "Met183"


###################################################################################################################
###################################################################################################################
###################################################################################################################

# homair_28 and met at age 14

## HILIC

merged_omics_HILIC_age14 <- merged_omics_HILIC[merged_omics_HILIC$Year == "14",]
qqt <- as.data.frame(apply(merged_omics_HILIC_age14[,colnames(merged_omics_HILIC_age14) %in% paste0("Met",seq(1:nrow(data_HILIC)))], 2, function(x) new_quantile(x, cuts = 10)))
data.homair_28.met_age14 = cbind(qqt[,c(paste0("Met",seq(1:nrow(data_HILIC))))], merged_omics_HILIC_age14[,c('sex', 'mage',  'mbmi', 'age14', "homair_28")])
data.homair_28.met_age14$beta <- rep(NA_real_, nrow(data.homair_28.met_age14))
data.homair_28.met_age14$model_pval <- rep(NA_real_, nrow(data.homair_28.met_age14))
d_lm_status <- data.frame(Met_id = NA_character_, Beta = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:1991){
  s_lm <- (lm(homair_28 ~ data.homair_28.met_age14[,i] + sex + age14 , data = data.homair_28.met_age14))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  d_lm_status <- rbind(d_lm_status,c(colnames(data.homair_28.met_age14)[i], as.numeric(r.est[2,c(1,2,3,4)])))
}

d_lm_status <- d_lm_status[-1,]
d_lm_status$z.value <- as.numeric(d_lm_status$z.value)
d_lm_status$p.value <- as.numeric(d_lm_status$p.value)
data.homair_28.met_age14 <- cbind(d_lm_status, data_HILIC[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Class")])

mumm_homair_28.met_age14 <- data.homair_28.met_age14[,c("mz","p.value","z.value","time")]
colnames(mumm_homair_28.met_age14) <- c("m.z", "p.value", "t.score", "rt")
mumm_homair_28.met_age14 <- mumm_homair_28.met_age14[order(mumm_homair_28.met_age14$p.value),]
write.table(mumm_homair_28.met_age14, sep = "\t",
            "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/mummichog_HILIC_homair_28_age14.txt",
            row.names = F, col.names = T)


write.csv(data.homair_28.met_age14,
          "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_HILIC_homair_28_age14.csv",
          row.names = F)

exwas_HILIC_homair_28_age14 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_HILIC_homair_28_age14.csv")
exwas_HILIC_homair_28_age14$Mode <- rep("HILIC",nrow(exwas_HILIC_homair_28_age14))

ght <- cor(cbind(merged_omics_HILIC_age14[,paste0("Met",seq(1:nrow(data_HILIC)))]))
et <- eigen(ght)
1/ (sum((et$values>1 + 0)* (et$values - 1)))
# 0.0005356186

exwas_HILIC_homair_28_age14$Met_id[exwas_HILIC_homair_28_age14$p.value < 0.0005356186]
# "Met147" "Met193"


###################################################################################################################

## C18

merged_omics_C18_age14 <- merged_omics_C18[merged_omics_C18$Year == "14",]
qqt <- as.data.frame(apply(merged_omics_C18_age14[,colnames(merged_omics_C18_age14) %in% paste0("Met",seq(1:nrow(data_C18)))], 2, function(x) new_quantile(x, cuts = 10)))
data.homair_28.met_age14 = cbind(qqt[,c(paste0("Met",seq(1:nrow(data_C18))))], merged_omics_C18_age14[,c('sex', 'mage',  'mbmi', 'age14', "homair_28")])
data.homair_28.met_age14$beta <- rep(NA_real_, nrow(data.homair_28.met_age14))
data.homair_28.met_age14$model_pval <- rep(NA_real_, nrow(data.homair_28.met_age14))
d_lm_status <- data.frame(Met_id = NA_character_, Beta = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:787){
  s_lm <- (lm(homair_28 ~ data.homair_28.met_age14[,i] + sex + age14 , data = data.homair_28.met_age14))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  d_lm_status <- rbind(d_lm_status,c(colnames(data.homair_28.met_age14)[i], as.numeric(r.est[2,c(1,2,3,4)])))
}

d_lm_status <- d_lm_status[-1,]
d_lm_status$z.value <- as.numeric(d_lm_status$z.value)
d_lm_status$p.value <- as.numeric(d_lm_status$p.value)
data.homair_28.met_age14 <- cbind(d_lm_status, data_C18[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Class")])

mumm_homair_28.met_age14 <- data.homair_28.met_age14[,c("mz","p.value","z.value","time")]
colnames(mumm_homair_28.met_age14) <- c("m.z", "p.value", "t.score", "rt")
mumm_homair_28.met_age14 <- mumm_homair_28.met_age14[order(mumm_homair_28.met_age14$p.value),]
write.table(mumm_homair_28.met_age14, sep = "\t",
            "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/mummichog_C18_homair_28_age14.txt",
            row.names = F, col.names = T)


write.csv(data.homair_28.met_age14,
          "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_C18_homair_28_age14.csv",
          row.names = F)

exwas_C18_homair_28_age14 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_C18_homair_28_age14.csv")
exwas_C18_homair_28_age14$Mode <- rep("C18",nrow(exwas_C18_homair_28_age14))

ght <- cor(cbind(merged_omics_C18_age14[,paste0("Met",seq(1:nrow(data_C18)))]))
et <- eigen(ght)
1/ (sum((et$values>1 + 0)* (et$values - 1)))
# 0.001507908

exwas_C18_homair_28_age14$Met_id[exwas_C18_homair_28_age14$p.value < 0.001507908]
#  "Met8"   "Met12"  "Met134" "Met217" "Met384"


###################################################################################################################
###################################################################################################################
###################################################################################################################

# homair_28 and met at age 22

## HILIC

merged_omics_HILIC_age22 <- merged_omics_HILIC[merged_omics_HILIC$Year == "22",]
qqt <- as.data.frame(apply(merged_omics_HILIC_age22[,colnames(merged_omics_HILIC_age22) %in% paste0("Met",seq(1:nrow(data_HILIC)))], 2, function(x) new_quantile(x, cuts = 10)))
data.homair_28.met_age22 = cbind(qqt[,c(paste0("Met",seq(1:nrow(data_HILIC))))], merged_omics_HILIC_age22[,c('sex', 'mage',  'mbmi', 'age22', "homair_28")])
data.homair_28.met_age22$beta <- rep(NA_real_, nrow(data.homair_28.met_age22))
data.homair_28.met_age22$model_pval <- rep(NA_real_, nrow(data.homair_28.met_age22))
d_lm_status <- data.frame(Met_id = NA_character_, Beta = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:1991){
  s_lm <- (lm(homair_28 ~ data.homair_28.met_age22[,i] + sex + age22 , data = data.homair_28.met_age22))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  d_lm_status <- rbind(d_lm_status,c(colnames(data.homair_28.met_age22)[i], as.numeric(r.est[2,c(1,2,3,4)])))
}

d_lm_status <- d_lm_status[-1,]
d_lm_status$z.value <- as.numeric(d_lm_status$z.value)
d_lm_status$p.value <- as.numeric(d_lm_status$p.value)
data.homair_28.met_age22 <- cbind(d_lm_status, data_HILIC[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Class")])

mumm_homair_28.met_age22 <- data.homair_28.met_age22[,c("mz","p.value","z.value","time")]
colnames(mumm_homair_28.met_age22) <- c("m.z", "p.value", "t.score", "rt")
mumm_homair_28.met_age22 <- mumm_homair_28.met_age22[order(mumm_homair_28.met_age22$p.value),]
write.table(mumm_homair_28.met_age22, sep = "\t",
            "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/mummichog_HILIC_homair_28_age22.txt",
            row.names = F, col.names = T)


write.csv(data.homair_28.met_age22,
          "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_HILIC_homair_28_age22.csv",
          row.names = F)

exwas_HILIC_homair_28_age22 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_HILIC_homair_28_age22.csv")
exwas_HILIC_homair_28_age22$Mode <- rep("HILIC",nrow(exwas_HILIC_homair_28_age22))

ght <- cor(cbind(merged_omics_HILIC_age22[,paste0("Met",seq(1:nrow(data_HILIC)))]))
et <- eigen(ght)
1/ (sum((et$values>1 + 0)* (et$values - 1)))
# 0.0005356186

exwas_HILIC_homair_28_age22$Met_id[exwas_HILIC_homair_28_age22$p.value < 0.0005356186]
# character(0)


###################################################################################################################

## C18

merged_omics_C18_age22 <- merged_omics_C18[merged_omics_C18$Year == "22",]
qqt <- as.data.frame(apply(merged_omics_C18_age22[,colnames(merged_omics_C18_age22) %in% paste0("Met",seq(1:nrow(data_C18)))], 2, function(x) new_quantile(x, cuts = 10)))
data.homair_28.met_age22 = cbind(qqt[,c(paste0("Met",seq(1:nrow(data_C18))))], merged_omics_C18_age22[,c('sex', 'mage',  'mbmi', 'age22', "homair_28")])
data.homair_28.met_age22$beta <- rep(NA_real_, nrow(data.homair_28.met_age22))
data.homair_28.met_age22$model_pval <- rep(NA_real_, nrow(data.homair_28.met_age22))
d_lm_status <- data.frame(Met_id = NA_character_, Beta = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:787){
  s_lm <- (lm(homair_28 ~ data.homair_28.met_age22[,i] + sex + age22 , data = data.homair_28.met_age22))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  d_lm_status <- rbind(d_lm_status,c(colnames(data.homair_28.met_age22)[i], as.numeric(r.est[2,c(1,2,3,4)])))
}

d_lm_status <- d_lm_status[-1,]
d_lm_status$z.value <- as.numeric(d_lm_status$z.value)
d_lm_status$p.value <- as.numeric(d_lm_status$p.value)
data.homair_28.met_age22 <- cbind(d_lm_status, data_C18[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Class")])

mumm_homair_28.met_age22 <- data.homair_28.met_age22[,c("mz","p.value","z.value","time")]
colnames(mumm_homair_28.met_age22) <- c("m.z", "p.value", "t.score", "rt")
mumm_homair_28.met_age22 <- mumm_homair_28.met_age22[order(mumm_homair_28.met_age22$p.value),]
write.table(mumm_homair_28.met_age22, sep = "\t",
            "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/mummichog_C18_homair_28_age22.txt",
            row.names = F, col.names = T)


write.csv(data.homair_28.met_age22,
          "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_C18_homair_28_age22.csv",
          row.names = F)

exwas_C18_homair_28_age22 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_C18_homair_28_age22.csv")
exwas_C18_homair_28_age22$Mode <- rep("C18",nrow(exwas_C18_homair_28_age22))

ght <- cor(cbind(merged_omics_C18_age22[,paste0("Met",seq(1:nrow(data_C18)))]))
et <- eigen(ght)
1/ (sum((et$values>1 + 0)* (et$values - 1)))
# 0.001508217

exwas_C18_homair_28_age22$Met_id[exwas_C18_homair_28_age22$p.value < 0.001508217]
#  "Met252" "Met339" "Met673" "Met727"

###################################################################################################################
###################################################################################################################
###################################################################################################################

# homair_28 and met at age 28

## HILIC

merged_omics_HILIC_age28 <- merged_omics_HILIC[merged_omics_HILIC$Year == "28",]
qqt <- as.data.frame(apply(merged_omics_HILIC_age28[,colnames(merged_omics_HILIC_age28) %in% paste0("Met",seq(1:nrow(data_HILIC)))], 2, function(x) new_quantile(x, cuts = 10)))
data.homair_28.met_age28 = cbind(qqt[,c(paste0("Met",seq(1:nrow(data_HILIC))))], merged_omics_HILIC_age28[,c('sex', 'mage',  'mbmi', 'age28', "homair_28")])
data.homair_28.met_age28$beta <- rep(NA_real_, nrow(data.homair_28.met_age28))
data.homair_28.met_age28$model_pval <- rep(NA_real_, nrow(data.homair_28.met_age28))
d_lm_status <- data.frame(Met_id = NA_character_, Beta = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:1991){
  s_lm <- (lm(homair_28 ~ data.homair_28.met_age28[,i] + sex + age28 , data = data.homair_28.met_age28))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  d_lm_status <- rbind(d_lm_status,c(colnames(data.homair_28.met_age28)[i], as.numeric(r.est[2,c(1,2,3,4)])))
}

d_lm_status <- d_lm_status[-1,]
d_lm_status$z.value <- as.numeric(d_lm_status$z.value)
d_lm_status$p.value <- as.numeric(d_lm_status$p.value)
data.homair_28.met_age28 <- cbind(d_lm_status, data_HILIC[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Class")])

mumm_homair_28.met_age28 <- data.homair_28.met_age28[,c("mz","p.value","z.value","time")]
colnames(mumm_homair_28.met_age28) <- c("m.z", "p.value", "t.score", "rt")
mumm_homair_28.met_age28 <- mumm_homair_28.met_age28[order(mumm_homair_28.met_age28$p.value),]
write.table(mumm_homair_28.met_age28, sep = "\t",
            "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/mummichog_HILIC_homair_28_age28.txt",
            row.names = F, col.names = T)


write.csv(data.homair_28.met_age28,
          "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_HILIC_homair_28_age28.csv",
          row.names = F)

exwas_HILIC_homair_28_age28 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_HILIC_homair_28_age28.csv")
exwas_HILIC_homair_28_age28$Mode <- rep("HILIC",nrow(exwas_HILIC_homair_28_age28))

ght <- cor(cbind(merged_omics_HILIC_age28[,paste0("Met",seq(1:nrow(data_HILIC)))]))
et <- eigen(ght)
1/ (sum((et$values>1 + 0)* (et$values - 1)))
# 0.0005356186

exwas_HILIC_homair_28_age28$Met_id[exwas_HILIC_homair_28_age28$p.value < 0.0005356186]
# "Met5"    "Met12"   "Met15"   "Met23"   "Met31"   "Met57"   "Met99"   "Met377"  "Met412"  "Met424"  "Met433" 
# [12] "Met496"  "Met509"  "Met510"  "Met606"  "Met712"  "Met723"  "Met728"  "Met1017" "Met1022" "Met1131" "Met1242"
# [23] "Met1261" "Met1277" "Met1311" "Met1317" "Met1431" "Met1441" "Met1470" "Met1597" "Met1741" "Met1772" "Met1883"
# [34] "Met1948"


###################################################################################################################

## C18

merged_omics_C18_age28 <- merged_omics_C18[merged_omics_C18$Year == "28",]
qqt <- as.data.frame(apply(merged_omics_C18_age28[,colnames(merged_omics_C18_age28) %in% paste0("Met",seq(1:nrow(data_C18)))], 2, function(x) new_quantile(x, cuts = 10)))
data.homair_28.met_age28 = cbind(qqt[,c(paste0("Met",seq(1:nrow(data_C18))))], merged_omics_C18_age28[,c('sex', 'mage',  'mbmi', 'age28', "homair_28")])
data.homair_28.met_age28$beta <- rep(NA_real_, nrow(data.homair_28.met_age28))
data.homair_28.met_age28$model_pval <- rep(NA_real_, nrow(data.homair_28.met_age28))
d_lm_status <- data.frame(Met_id = NA_character_, Beta = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:787){
  s_lm <- (lm(homair_28 ~ data.homair_28.met_age28[,i] + sex + age28 , data = data.homair_28.met_age28))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  d_lm_status <- rbind(d_lm_status,c(colnames(data.homair_28.met_age28)[i], as.numeric(r.est[2,c(1,2,3,4)])))
}

d_lm_status <- d_lm_status[-1,]
d_lm_status$z.value <- as.numeric(d_lm_status$z.value)
d_lm_status$p.value <- as.numeric(d_lm_status$p.value)
data.homair_28.met_age28 <- cbind(d_lm_status, data_C18[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Class")])

mumm_homair_28.met_age28 <- data.homair_28.met_age28[,c("mz","p.value","z.value","time")]
colnames(mumm_homair_28.met_age28) <- c("m.z", "p.value", "t.score", "rt")
mumm_homair_28.met_age28 <- mumm_homair_28.met_age28[order(mumm_homair_28.met_age28$p.value),]
write.table(mumm_homair_28.met_age28, sep = "\t",
            "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/mummichog_C18_homair_28_age28.txt",
            row.names = F, col.names = T)


write.csv(data.homair_28.met_age28,
          "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_C18_homair_28_age28.csv",
          row.names = F)

exwas_C18_homair_28_age28 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_C18_homair_28_age28.csv")
exwas_C18_homair_28_age28$Mode <- rep("C18",nrow(exwas_C18_homair_28_age28))

ght <- cor(cbind(merged_omics_C18_age28[,paste0("Met",seq(1:nrow(data_C18)))]))
et <- eigen(ght)
1/ (sum((et$values>1 + 0)* (et$values - 1)))
# 0.001508217

exwas_C18_homair_28_age28$Met_id[exwas_C18_homair_28_age28$p.value < 0.001508293]
# [1] "Met5"   "Met45"  "Met69"  "Met86"  "Met96"  "Met196" "Met216" "Met219" "Met273" "Met425" "Met475" "Met500"
# [13] "Met554" "Met708"


################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

# Plots
## At age7

exwas_C18_homair_28_age7 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_C18_homair_28_age7.csv")
exwas_C18_homair_28_age7$Mode <- rep("C18",nrow(exwas_C18_homair_28_age7))
exwas_C18_homair_28_age7 <- exwas_C18_homair_28_age7[exwas_C18_homair_28_age7$p.value < 0.001507876,]
exwas_C18_homair_28_age7$age <- rep("Age 7", nrow(exwas_C18_homair_28_age7))

exwas_HILIC_homair_28_age7 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_HILIC_homair_28_age7.csv")
exwas_HILIC_homair_28_age7$Mode <- rep("HILIC",nrow(exwas_HILIC_homair_28_age7))
exwas_HILIC_homair_28_age7 <- exwas_HILIC_homair_28_age7[exwas_HILIC_homair_28_age7$p.value < 0.0005356186,]
exwas_HILIC_homair_28_age7$age <- rep("Age 7", nrow(exwas_HILIC_homair_28_age7))

## At age14

exwas_C18_homair_28_age14 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_C18_homair_28_age14.csv")
exwas_C18_homair_28_age14$Mode <- rep("C18",nrow(exwas_C18_homair_28_age14))
exwas_C18_homair_28_age14 <- exwas_C18_homair_28_age14[exwas_C18_homair_28_age14$p.value < 0.001507908,]
exwas_C18_homair_28_age14$age <- rep("Age 14", nrow(exwas_C18_homair_28_age14))


exwas_HILIC_homair_28_age14 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_HILIC_homair_28_age14.csv")
exwas_HILIC_homair_28_age14$Mode <- rep("HILIC",nrow(exwas_HILIC_homair_28_age14))
exwas_HILIC_homair_28_age14 <- exwas_HILIC_homair_28_age14[exwas_HILIC_homair_28_age14$p.value < 0.0005356186,]
exwas_HILIC_homair_28_age14$age <- rep("Age 14", nrow(exwas_HILIC_homair_28_age14))

## At age22

exwas_C18_homair_28_age22 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_C18_homair_28_age22.csv")
exwas_C18_homair_28_age22$Mode <- rep("C18",nrow(exwas_C18_homair_28_age22))
exwas_C18_homair_28_age22 <- exwas_C18_homair_28_age22[exwas_C18_homair_28_age22$p.value < 0.001508217,]
exwas_C18_homair_28_age22$age <- rep("Age 22", nrow(exwas_C18_homair_28_age22))


exwas_HILIC_homair_28_age22 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_HILIC_homair_28_age22.csv")
exwas_HILIC_homair_28_age22$Mode <- rep("HILIC",nrow(exwas_HILIC_homair_28_age22))
exwas_HILIC_homair_28_age22 <- exwas_HILIC_homair_28_age22[exwas_HILIC_homair_28_age22$p.value < 0.0005356186,]
exwas_HILIC_homair_28_age22$age <- rep("Age 22", nrow(exwas_HILIC_homair_28_age22))

## At age28

exwas_C18_homair_28_age28 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_C18_homair_28_age28.csv")
exwas_C18_homair_28_age28$Mode <- rep("C18",nrow(exwas_C18_homair_28_age28))
exwas_C18_homair_28_age28 <- exwas_C18_homair_28_age28[exwas_C18_homair_28_age28$p.value < 0.001508293,]
exwas_C18_homair_28_age28$age <- rep("Age 28", nrow(exwas_C18_homair_28_age28))


exwas_HILIC_homair_28_age28 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_HILIC_homair_28_age28.csv")
exwas_HILIC_homair_28_age28$Mode <- rep("HILIC",nrow(exwas_HILIC_homair_28_age28))
exwas_HILIC_homair_28_age28 <- exwas_HILIC_homair_28_age28[exwas_HILIC_homair_28_age28$p.value < 0.0005356186,]
exwas_HILIC_homair_28_age28$age <- rep("Age 28", nrow(exwas_HILIC_homair_28_age28))


##########################################################################################
##########################################################################################
##########################################################################################

# Naming the metabolites and Corresponding adducts

exwas_homair_28 <- rbind(exwas_HILIC_homair_28_age7,exwas_C18_homair_28_age7,
                      exwas_HILIC_homair_28_age14,exwas_C18_homair_28_age14,
                      exwas_HILIC_homair_28_age22,exwas_C18_homair_28_age22,
                      exwas_HILIC_homair_28_age28,exwas_C18_homair_28_age28)


conf_met <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/confirmed_metabolites.csv",fileEncoding = "Latin1")
conf_met_HILIC <- conf_met[conf_met$Method.RT == "HILIC+",]
conf_met_C18 <- conf_met[conf_met$Method.RT == "C18-",]

stage4_hilic <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/Stage4.csv")
stage4_c18 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/Stage4.csv")

all_sig_hits <- exwas_homair_28
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
write.csv(all_sig_hits, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_allsig_metabolites_homair_28.csv", row.names = F )

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
# significant Metabolites plot
# Without the ones at age 28

exwas_homair_28 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_allsig_metabolites_homair_28.csv")
exwas_homair_28 <- exwas_homair_28[exwas_homair_28$age!= "Age 28",]
  
exwas_homair_28$log10.p.value <- -log(exwas_homair_28$p.value,10)
exwas_homair_28$age <- factor(exwas_homair_28$age, levels = c('Age 7', 'Age 14', 'Age 22'))

exwas_homair_28$Metabolite <- rep(NA_character_,nrow(exwas_homair_28))

for(i in 1:nrow(exwas_homair_28)) exwas_homair_28$Metabolite[i] <- strsplit(exwas_homair_28$Name,"/")[[i]][1]

exwas_homair_28<- exwas_homair_28 %>% # if there are duplicates, only keep one metabolite with lowest p.value at each timepoint
  arrange(age, p.value)%>%
  group_by(age) %>% 
  distinct(Metabolite, .keep_all = TRUE) 


pal <- wes_palette("Zissou1", 100, type = "continuous")
vol <- (ggplot(exwas_homair_28, aes(x=age, y=Beta, color=log10.p.value, label=Metabolite, shape = Mode)) +# Show all points
          geom_point(size = 3)+  
          scale_shape_manual(values = c(15,17)) + # change shape
          geom_hline(yintercept= 0, color = "black", size = 1 )   +
          scale_color_gradientn(colours = pal) + # change color
          labs(x = "Age at metabolomic assessment", y = "Beta Coefficients", 
               title = "(A) Prospective metabolic exposures and HOMA-IR at age 28",
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
           ylim(c(-0.3,0.3)) +
           annotate("rect", xmin = c(0.75, 1.75, 2.75), 
                    xmax = c(1.25, 2.25, 3.25), 
                    ymin = c(-0.3, -0.3, -0.3), 
                    ymax = c(0.3, 0.3, 0.3),
                    alpha = 0.05))



jpeg("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/plot_allsig_metabolites_homair_28.jpeg",
     units="in", width=15, height=10, res=600)

vol

dev.off()


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

# Pathway Plot

C18_homair_28_age7  <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/metaboanalyst/mummichog_pathway_enrichment_C18_homair_28_age7.csv")
C18_homair_28_age7 <- C18_homair_28_age7[C18_homair_28_age7$Hits.sig >=3 ,]
C18_homair_28_age7$age <- rep(7, nrow(C18_homair_28_age7))

C18_homair_28_age14  <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/metaboanalyst/mummichog_pathway_enrichment_C18_homair_28_age14.csv")
C18_homair_28_age14 <- C18_homair_28_age14[C18_homair_28_age14$Hits.sig >=3 ,]
C18_homair_28_age14$age <- rep(14, nrow(C18_homair_28_age14))

C18_homair_28_age22  <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/metaboanalyst/mummichog_pathway_enrichment_C18_homair_28_age22.csv")
C18_homair_28_age22 <- C18_homair_28_age22[C18_homair_28_age22$Hits.sig >=3 ,]
C18_homair_28_age22$age <- rep(22, nrow(C18_homair_28_age22))

C18_homair_28 <- rbind(C18_homair_28_age7,C18_homair_28_age14, C18_homair_28_age22)
C18_homair_28$Mode <- rep("C18-", nrow(C18_homair_28))


HILIC_homair_28_age7  <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/metaboanalyst/mummichog_pathway_enrichment_HILIC_homair_28_age7.csv")
HILIC_homair_28_age7 <- HILIC_homair_28_age7[HILIC_homair_28_age7$Hits.sig >=3 ,]
HILIC_homair_28_age7$age <- rep(7, nrow(HILIC_homair_28_age7))

HILIC_homair_28_age14  <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/metaboanalyst/mummichog_pathway_enrichment_HILIC_homair_28_age14.csv")
HILIC_homair_28_age14 <- HILIC_homair_28_age14[HILIC_homair_28_age14$Hits.sig >=3 ,]
HILIC_homair_28_age14$age <- rep(14, nrow(HILIC_homair_28_age14))

HILIC_homair_28_age22  <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/metaboanalyst/mummichog_pathway_enrichment_HILIC_homair_28_age22.csv")
HILIC_homair_28_age22 <- HILIC_homair_28_age22[HILIC_homair_28_age22$Hits.sig >=3 ,]
HILIC_homair_28_age22$age <- rep(22, nrow(HILIC_homair_28_age22))

HILIC_homair_28 <- rbind(HILIC_homair_28_age7,HILIC_homair_28_age14, HILIC_homair_28_age22)
HILIC_homair_28$Mode <- rep("HILIC+", nrow(HILIC_homair_28))

pathways_homair_28 <- rbind(C18_homair_28, HILIC_homair_28)
colnames(pathways_homair_28)[1] <- "Pathways"
pathways_homair_28$age <- factor(pathways_homair_28$age, levels = c("14","22"))

pfas_pathway_subset <- pathways_homair_28[pathways_homair_28$FET < 0.05,]

all_pathways <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/all_pathways.csv")
pfas_pathway_subset <- merge(pfas_pathway_subset, all_pathways, by.x = "Pathways")
write.csv(pfas_pathway_subset, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/allsig_pathway_homair_28.csv", row.names = F )





pfas_pathway_subset<- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/allsig_pathway_homair_28.csv")

pfas_pathway_subset$Group <- factor(pfas_pathway_subset$Group, levels = c("Lipid metabolism",
                                                                          "Metabolism of amino acids and derivatives", 
                                                                          "Nucleotide metabolism"
                                                                          ))

groups <- c("Lipid metabolism",
            "Metabolism of amino acids and derivatives", 
            "Nucleotide metabolism")
            
pfas_pathway_subset <- pfas_pathway_subset[order(pfas_pathway_subset$Group),]
uni_paths <- pfas_pathway_subset$Pathways
uni_paths <- uni_paths[!duplicated(uni_paths)]
pfas_pathway_subset$Pathways <- factor(pfas_pathway_subset$Pathways, levels = uni_paths)

ret <- c(length(unique(pfas_pathway_subset$Pathways[pfas_pathway_subset$Group ==  unique(pfas_pathway_subset$Group)[1]])),
         length(unique(pfas_pathway_subset$Pathways[pfas_pathway_subset$Group ==  unique(pfas_pathway_subset$Group)[2]])),
         length(unique(pfas_pathway_subset$Pathways[pfas_pathway_subset$Group ==  unique(pfas_pathway_subset$Group)[3]])))

gp <- (ggplot(pfas_pathway_subset, aes(x = age,  y = Pathways)) + 
         geom_point(aes(size = -log(FET , base = 10), shape = Mode), 
                    alpha  = 5, 
                    position = position_dodge2(width = 0.5, preserve  = "total", padding = 0.5, reverse = T)) + 
         theme_bw()   + scale_shape_manual(values = c(15,17)) +
         ylab(NULL) + xlab("Age at Metabolomic Assessment") +
         ggtitle("(B) Enriched pathways in association between\nserum Metabolites and HOMA-IR at age 28") + labs(size = "-log10(p.value)") 
       + theme(plot.title=element_text(size=16,face="bold"), axis.title=element_text(size=14,face="bold"),
               plot.tag = element_text(size = 13,face = "bold"),
               axis.text.y = element_text(size=15,face="bold"),
               axis.text.x = element_text(size=12,face="bold"),
               strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
               legend.title = element_text(face = "bold"), legend.position = 'right',
               legend.text = element_text(size = 12),
               legend.background = element_rect(fill="white", 
                                                size=0.5, linetype="solid",  colour ="darkblue"),
               plot.margin=unit(c(0.5,0.5,1,0.5), "cm"))
       +  annotate("rect", xmin = c(0.75,1.75), 
                   xmax = c(1.25,2.25), 
                   ymin = rep(0,2), ymax =rep(6,2),
                   alpha = .1)  
       + guides(color = guide_legend(override.aes = list(size = 4)),
                shape = guide_legend(override.aes = list(size = 4))))

gp <- (gp +  geom_hline(yintercept = c(cumsum(ret))+0.5)
       + geom_hline(yintercept = 6, size = 2)
       + annotate("text", x = rep(4.9, length(unique(pfas_pathway_subset$Group))), 
                  y = c(0.8, 3, 5), 
                  label = c("Lipid metabolism",
                            "Metabolism of amino acids\nand derivatives", 
                            "Nucleotide metabolism"), color = "red",
                  size = rep(5, length(groups)), fontface  = 'bold')
       + coord_cartesian(xlim = c(1.2, 6), ylim = c(0.5, 5.5), clip = "off") )

jpeg("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese//hrm_clinical_outcome/homair/homair_28/pathways_homair_28.jpeg",
     units="in", width=12, height=10, res=600)

gp

dev.off()

##################################################

# combine plots

jpeg("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese//hrm_clinical_outcome/homair/homair_28/plot_metabolites_pathways_homair_28.jpeg",
     units="in", width=32, height=15, res=800)

ggpubr::ggarrange(vol, gp, nrow = 1, ncol = 2, common.legend = F, 
                  legend = "bottom") + theme_bw()


dev.off()
