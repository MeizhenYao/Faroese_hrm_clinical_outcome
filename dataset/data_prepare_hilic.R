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

# Set working directory in local repository

# DDD

# setwd("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/git/Faroese_PFAS_vs_Metabolites/HILIC")


# Always start with renv to capture generate workflow

# renv::init()

# Epi data
t<-as.data.frame(read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/17.12.19_December2020.csv"))
d <- as.data.frame(read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/17.12.19_December2020.csv"))
d$id <- as.character(d$id)

# eligible id: 500 ids

d1 <- as.data.frame(read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/eligible_ids.csv"))
d1$id <- as.character(as.numeric(as.character(d1$id)))

# Match the eligible ids with epi data

y <- rep(NA_real_, length(d1$id))
for(i in 1:length(d1$id)){
  y[i] = sum(d$id %in% d1$id[i])
}

dt <- d[which(d$id %in% d1$id),]

# Select the PFAS exposure columns and covariates

d <- dt[,c("id",'pfos_0', 'pfos_7', 'pfos_14', 'pfos_22', 'pfos_28', 'pfoa_0', 'pfoa_7', 'pfoa_14', 'pfoa_22', 'pfoa_28', 
          'pfhxs_0', 'pfhxs_7', 'pfhxs_14', 'pfhxs_22', 'pfhxs_28',  'pfna_0', 'pfna_7', 'pfna_14', 'pfna_22', 'pfna_28', 
          'pfda_0', 'pfda_7', 'pfda_14', 'pfda_22', 'pfda_28', 'bmi_28', 'waistcirc_28', 'diaBP_dxt_28', 'sysBP_dxt_28', 
          'uhdl_28y', 'trig_28y', 'dldl_28y', 'chol_28y','cir_28', 'insulinauc_28y', 'glucoseauc_28y', 'matsuda_28', 
          'homair_28', 'igi_28', "sex","mage","mbmi","parity","smokepreg_2","matfishpreg_cat2","breastfed_tot",
          "age_7","age_14","age22","age28")]

d$id <- as.character(d$id)

# 493 ids

# Impute Data using MICE

init = mice(d, maxit=0) 
meth = init$method
predM = init$predictorMatrix
set.seed(1234)
imputed = mice(d, method="norm", predictorMatrix=predM, m=3, ntree = 100)
imputed <- complete(imputed, action = 2) # this is your final imputed dataset


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


# Read Metabolomics data - HILIC+

mapping_hilic <- read.delim("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/Perfluoroalkyl_ALL_mapping_hilicpos.txt")

# read median summarized HILIC+ table
median_sum_feature_hilic <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/MetabComBat_HILICpos_Merged_mediansummarized_featuretable_final.csv")

#7109

### remove VT_201114_M456_133_A, VT_201114_M456_139_A, VT_201114_M456_133_B, VT_201114_M456_139_B (Not patient sample - duplicated)

median_sum_feature_hilic <- median_sum_feature_hilic[,!(colnames(median_sum_feature_hilic) %in% c("VT_201114_M456_133_A","VT_201114_M456_139_A",
                                                                                                 "VT_201114_M456_133_B","VT_201114_M456_139_B"))]

# Extract File name ID_YEAR
x <- colnames(median_sum_feature_hilic)[3:ncol(median_sum_feature_hilic)]
y <- rep(NA_character_, length(x))
for(i in 1:length(x)){
  y[i] <- sum(mapping_hilic$File.Name %in% x[i])
}

# median_sum_feature_hilic <- median_sum_feature_hilic[,!(colnames(median_sum_feature_hilic) %in% x[which(y == "0")])]
x <- colnames(median_sum_feature_hilic)[3:ncol(median_sum_feature_hilic)]
y <- rep(NA_character_, length(x))
for(i in 1:length(x)){
  y[i] = mapping_hilic$Sample.ID[which(mapping_hilic$File.Name %in% x[i])] 
}
colnames(median_sum_feature_hilic)[3:ncol(median_sum_feature_hilic)] <- y


# valid columns - columns with year_id
valid_colnames <- NA_character_
for(i in 1:length(colnames(median_sum_feature_hilic))){
  if(is.na(strsplit(colnames(median_sum_feature_hilic)[i],"_")[[1]][2])==F & (is.na((strsplit(strsplit(colnames(median_sum_feature_hilic)[i],"_")[[1]][2],"")[[1]][3] == "y")) == F & (strsplit(strsplit(colnames(median_sum_feature_hilic)[i],"_")[[1]][2],"")[[1]][3] == "y")) ||
     (is.na((strsplit(strsplit(colnames(median_sum_feature_hilic)[i],"_")[[1]][2],"")[[1]][2] == "y")) == F & (strsplit(strsplit(colnames(median_sum_feature_hilic)[i],"_")[[1]][2],"")[[1]][2] == "y"))){
    valid_colnames <- c(valid_colnames, colnames(median_sum_feature_hilic)[i])
  }
}

valid_colnames <- unique(valid_colnames[-1])
length(valid_colnames)/4      # 309  - year_id

# Data filtering
# for each metabolite, non-zero value > 75%
median_sum_feature_hilic$seq <- paste0("chem_",seq(1, nrow(median_sum_feature_hilic)))
median_sum_feature_hilic<- median_sum_feature_hilic[,c(valid_colnames,"mz","time","seq")]
median_sum_feature_hilic$zero<- rowSums(median_sum_feature_hilic[,valid_colnames]==0) #!!!!!!!!!!!!! by timepoint?
median_sum_feature_hilic$zero_per<- median_sum_feature_hilic$zero/length(valid_colnames)*100

median_sum_feature_hilic<- median_sum_feature_hilic[median_sum_feature_hilic$zero_per<25,]
#5067

# for each metabolite, CV% > 25%
median_sum_feature_hilic<- median_sum_feature_hilic %>% 
                           rowwise() %>% 
                           mutate(
                             mean = mean(c_across(valid_colnames[1]:valid_colnames[1236])),
                             sd = sd(c_across(valid_colnames[1]:valid_colnames[1236])),
                             cv = 100*(sd/mean))
median_sum_feature_hilic<- median_sum_feature_hilic[median_sum_feature_hilic$cv>30,]
#4205




# export data
data_met_clinical_hilic <- median_sum_feature_hilic[,c(valid_colnames,"mz","time","seq")]
data_met_clinical_hilic$Met_id <- paste0("Met",seq(1:nrow(data_met_clinical_hilic)))
dim(data_met_clinical_hilic) 
# [1] 4205 1240

write.csv(data_met_clinical_hilic, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/data_met_clinical_hilic.csv", row.names = F)


###################################################################################################################

data_met_clinical_hilic <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/data_met_clinical_hilic.csv", check.names = F)

# transpose met data
t_data_met_clinical_hilic <- data.table::transpose(data_met_clinical_hilic)

# get row and colnames in order
Met_name<- data_met_clinical_hilic$Met_id
colnames(t_data_met_clinical_hilic) <- Met_name
t_data_met_clinical_hilic$id_Year <- colnames(data_met_clinical_hilic)
dim(t_data_met_clinical_hilic)
# [1] 1238 4206


t_data_met_clinical_hilic$id <- rep(NA_character_,nrow(t_data_met_clinical_hilic))
t_data_met_clinical_hilic$Year <- rep(NA_character_,nrow(t_data_met_clinical_hilic))

for(i in 1:nrow(t_data_met_clinical_hilic)){
  t_data_met_clinical_hilic$id[i] <- strsplit(t_data_met_clinical_hilic$id_Year[i],"_")[[1]][1]
  t_data_met_clinical_hilic$Year[i] <- as.numeric(knitr::combine_words(strsplit(strsplit(t_data_met_clinical_hilic$id_Year[i],"_")[[1]][2],"y")[[1]][1],and = "", sep =""))
}


# combine with Epi data
merged_met_clinical_hilic <- merge(d_long,t_data_met_clinical_hilic, by = c("id","Year"))
merged_met_clinical_hilic$id <- as.character(merged_met_clinical_hilic$id)
dim(merged_met_clinical_hilic)


# check available sample size at different timepoints
nrow(merged_met_clinical_hilic %>% filter(Year == "7"))
nrow(merged_met_clinical_hilic %>% filter(Year == "14"))
nrow(merged_met_clinical_hilic %>% filter(Year == "22"))
nrow(merged_met_clinical_hilic %>% filter(Year == "28"))

## Subset for only 125 ids
merged_met_clinical_hilic <- merged_met_clinical_hilic[which(merged_met_clinical_hilic$id %in% as.vector(t_data_met_clinical_hilic[t_data_met_clinical_hilic$Year == "7","id"])),]
dim(merged_met_clinical_hilic)
# [1] 500 4238



# Met data standardization
merged_met_clinical_hilic<- merged_met_clinical_hilic %>% 
                      mutate_at(Met_name, ~(as.numeric(.) %>% as.vector))  %>% 
                      mutate_at(Met_name, ~(log(. + 1, base = 2) %>% as.vector))%>% 
                      mutate_at(Met_name,  ~(scale(.) %>% as.vector))


# combine with dichotomous PFAS
d_subset$cpfoa0 <- ifelse(d_subset$pfoa_0 > median(d_subset$pfoa_0), 1, 0)
d_subset$cpfos0 <- ifelse(d_subset$pfos_0 > median(d_subset$pfos_0), 1, 0)
d_subset$cpfhxs0 <- ifelse(d_subset$pfhxs_0 > median(d_subset$pfhxs_0), 1, 0)
d_subset$cpfna0 <- ifelse(d_subset$pfna_0 > median(d_subset$pfna_0), 1, 0)
d_subset$cpfda0 <- ifelse(d_subset$pfda_0 > median(d_subset$pfda_0), 1, 0)


d_subset$cpfoa7 <- ifelse(d_subset$pfoa_7 > median(d_subset$pfoa_7), 1, 0)
d_subset$cpfos7 <- ifelse(d_subset$pfos_7 > median(d_subset$pfos_7), 1, 0)
d_subset$cpfhxs7 <- ifelse(d_subset$pfhxs_7 > median(d_subset$pfhxs_7), 1, 0)
d_subset$cpfna7 <- ifelse(d_subset$pfna_7 > median(d_subset$pfna_7), 1, 0)
d_subset$cpfda7 <- ifelse(d_subset$pfda_7 > median(d_subset$pfda_7), 1, 0)

d_subset$cpfoa14 <- ifelse(d_subset$pfoa_14 > median(d_subset$pfoa_14), 1, 0)
d_subset$cpfos14 <- ifelse(d_subset$pfos_14 > median(d_subset$pfos_14), 1, 0)
d_subset$cpfhxs14 <- ifelse(d_subset$pfhxs_14 > median(d_subset$pfhxs_14), 1, 0)
d_subset$cpfna14 <- ifelse(d_subset$pfna_14 > median(d_subset$pfna_14), 1, 0)
d_subset$cpfda14 <- ifelse(d_subset$pfda_14 > median(d_subset$pfda_14), 1, 0)

d_subset$cpfoa22 <- ifelse(d_subset$pfoa_22 > median(d_subset$pfoa_22), 1, 0)
d_subset$cpfos22 <- ifelse(d_subset$pfos_22 > median(d_subset$pfos_22), 1, 0)
d_subset$cpfhxs22 <- ifelse(d_subset$pfhxs_22 > median(d_subset$pfhxs_22), 1, 0)
d_subset$cpfna22 <- ifelse(d_subset$pfna_22 > median(d_subset$pfna_22), 1, 0)
d_subset$cpfda22 <- ifelse(d_subset$pfda_22 > median(d_subset$pfda_22), 1, 0)

d_subset$cpfoa28 <- ifelse(d_subset$pfoa_28 > median(d_subset$pfoa_28), 1, 0)
d_subset$cpfos28 <- ifelse(d_subset$pfos_28 > median(d_subset$pfos_28), 1, 0)
d_subset$cpfhxs28 <- ifelse(d_subset$pfhxs_28 > median(d_subset$pfhxs_28), 1, 0)
d_subset$cpfna28 <- ifelse(d_subset$pfna_28 > median(d_subset$pfna_28), 1, 0)
d_subset$cpfda28 <- ifelse(d_subset$pfda_28 > median(d_subset$pfda_28), 1, 0)


merged_met_clinical_hilic <- merge(merged_met_clinical_hilic, d_subset[,c("pfoa_0","pfos_0","pfhxs_0","pfna_0","pfda_0","id",
                                                              "cpfoa0","cpfos0","cpfhxs0","cpfna0","cpfda0",
                                                              "cpfoa7","cpfos7","cpfhxs7","cpfna7","cpfda7",
                                                              "cpfoa14","cpfos14","cpfhxs14","cpfna14","cpfda14",
                                                              "cpfoa22","cpfos22","cpfhxs22","cpfna22","cpfda22",
                                                              "cpfoa28","cpfos28","cpfhxs28","cpfna28","cpfda28")],by = "id")



# reformate variables
merged_met_clinical_hilic$Year <- as.character(merged_met_clinical_hilic$Year)
merged_met_clinical_hilic$Year <- factor(merged_met_clinical_hilic$Year, levels  = c("7","14","22","28"))

merged_met_clinical_hilic$cmbmi <- ifelse(merged_met_clinical_hilic$mbmi >= 25, "High","Low")
merged_met_clinical_hilic$cmage <- ifelse(merged_met_clinical_hilic$mage >= as.numeric(quantile(merged_met_clinical_hilic$mage)[4]), "High","Low")
merged_met_clinical_hilic$cparity <- ifelse(merged_met_clinical_hilic$parity >= 1, "1","0")
merged_met_clinical_hilic$csmokepreg <- ifelse(merged_met_clinical_hilic$smokepreg_2 == 1, "1","0")
merged_met_clinical_hilic$cbreastfed_tot <- ifelse(merged_met_clinical_hilic$breastfed_tot >= as.numeric(quantile(merged_met_clinical_hilic$breastfed_tot)[2]), "High","Low")
merged_met_clinical_hilic$cmatfishpreg <- as.numeric(ifelse(merged_met_clinical_hilic$matfishpreg_cat2 >= 0.5, "1","0"))
merged_met_clinical_hilic[,(colnames(merged_met_clinical_hilic) %in% c(data_c18$Met_id[order(data_c18$Class,data_c18$Met_id)]))] <- (apply(merged_met_clinical_hilic[,(colnames(merged_met_clinical_hilic) %in% c(data_c18$Met_id[order(data_c18$Class,data_c18$Met_id)]))],2,as.numeric))

dim(merged_met_clinical_hilic[merged_met_clinical_hilic$Year == 7,(colnames(merged_met_clinical_hilic) %in% Met_name)])
dim(merged_met_clinical_hilic[merged_met_clinical_hilic$Year == 14,(colnames(merged_met_clinical_hilic) %in% Met_name)])
dim(merged_met_clinical_hilic[merged_met_clinical_hilic$Year == 22,(colnames(merged_met_clinical_hilic) %in% Met_name)])
dim(merged_met_clinical_hilic[merged_met_clinical_hilic$Year == 28,(colnames(merged_met_clinical_hilic) %in% Met_name)])
#[1]  125 4205

write.csv(merged_met_clinical_hilic, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/merged_met_clinical_hilic.csv", row.names = F)

