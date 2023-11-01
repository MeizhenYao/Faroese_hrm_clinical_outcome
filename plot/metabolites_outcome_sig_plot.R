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

# hypertension
exwas_hypertension_28_bi <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/hypertension/hypertension_28_bi/exwas_allsig_metabolites_hypertension_28_bi.csv")
exwas_hypertension_28_bi <- exwas_hypertension_28_bi[exwas_hypertension_28_bi$age!= "Age 28",]

exwas_hypertension_28_bi$log10.p.value <- -log(exwas_hypertension_28_bi$p.value,10)
exwas_hypertension_28_bi$age <- factor(exwas_hypertension_28_bi$age, levels = c('Age 7'))

exwas_hypertension_28_bi$Metabolite <- rep(NA_character_,nrow(exwas_hypertension_28_bi))

for(i in 1:nrow(exwas_hypertension_28_bi)) exwas_hypertension_28_bi$Metabolite[i] <- strsplit(exwas_hypertension_28_bi$Name,"/")[[i]][1]


exwas_hypertension_28_bi$outcome <- rep("Hypertension",nrow(exwas_hypertension_28_bi))

# matsuda
exwas_matsuda_28 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/matsuda/matsuda_28/exwas_allsig_metabolites_matsuda_28.csv")
exwas_matsuda_28 <- exwas_matsuda_28[exwas_matsuda_28$age!= "Age 28",]

exwas_matsuda_28$log10.p.value <- -log(exwas_matsuda_28$p.value,10)
exwas_matsuda_28$age <- factor(exwas_matsuda_28$age, levels = c('Age 7', 'Age 14', 'Age 22'))

exwas_matsuda_28$Metabolite <- rep(NA_character_,nrow(exwas_matsuda_28))

for(i in 1:nrow(exwas_matsuda_28)) exwas_matsuda_28$Metabolite[i] <- strsplit(exwas_matsuda_28$Name,"/")[[i]][1]

exwas_matsuda_28$outcome <- rep("Matsuda",nrow(exwas_matsuda_28))

# bmi
exwas_bmi_28 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28/exwas_allsig_metabolites_bmi_28.csv")
exwas_bmi_28 <- exwas_bmi_28[exwas_bmi_28$age!= "Age 28",]

exwas_bmi_28$log10.p.value <- -log(exwas_bmi_28$p.value,10)
exwas_bmi_28$age <- factor(exwas_bmi_28$age, levels = c('Age 14', 'Age 22'))

exwas_bmi_28$Metabolite <- rep(NA_character_,nrow(exwas_bmi_28))

for(i in 1:nrow(exwas_bmi_28)) exwas_bmi_28$Metabolite[i] <- strsplit(exwas_bmi_28$Name,"/")[[i]][1]

exwas_bmi_28$outcome <- rep("BMI",nrow(exwas_bmi_28))

# overweight
exwas_bmi_28_bi <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_bi/exwas_allsig_metabolites_bmi_28_bi.csv")
exwas_bmi_28_bi <- exwas_bmi_28_bi[exwas_bmi_28_bi$age!= "Age 28",]

exwas_bmi_28_bi$log10.p.value <- -log(exwas_bmi_28_bi$p.value,10)
exwas_bmi_28_bi$age <- factor(exwas_bmi_28_bi$age, levels = c('Age 7','Age 14', 'Age 22'))

exwas_bmi_28_bi$Metabolite <- rep(NA_character_,nrow(exwas_bmi_28_bi))

for(i in 1:nrow(exwas_bmi_28_bi)) exwas_bmi_28_bi$Metabolite[i] <- strsplit(exwas_bmi_28_bi$Name,"/")[[i]][1]

exwas_bmi_28_bi$outcome <- rep("Overweight",nrow(exwas_bmi_28_bi))

# GlucoseAUC
exwas_glucoseauc_28y <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/glucoseauc/glucoseauc_28y/exwas_allsig_metabolites_glucoseauc_28y.csv")
exwas_glucoseauc_28y <- exwas_glucoseauc_28y[exwas_glucoseauc_28y$age!= "Age 28",]

exwas_glucoseauc_28y$log10.p.value <- -log(exwas_glucoseauc_28y$p.value,10)
exwas_glucoseauc_28y$age <- factor(exwas_glucoseauc_28y$age, levels = c('Age 7','Age 14', 'Age 22'))

exwas_glucoseauc_28y$Metabolite <- rep(NA_character_,nrow(exwas_glucoseauc_28y))

for(i in 1:nrow(exwas_glucoseauc_28y)) exwas_glucoseauc_28y$Metabolite[i] <- strsplit(exwas_glucoseauc_28y$Name,"/")[[i]][1]

exwas_glucoseauc_28y$outcome <- rep("GlucoseAUC",nrow(exwas_glucoseauc_28y))

# HOMA-IR
exwas_homair_28 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/exwas_allsig_metabolites_homair_28.csv")
exwas_homair_28 <- exwas_homair_28[exwas_homair_28$age!= "Age 28",]

exwas_homair_28$log10.p.value <- -log(exwas_homair_28$p.value,10)
exwas_homair_28$age <- factor(exwas_homair_28$age, levels = c('Age 7', 'Age 14', 'Age 22'))

exwas_homair_28$Metabolite <- rep(NA_character_,nrow(exwas_homair_28))

for(i in 1:nrow(exwas_homair_28)) exwas_homair_28$Metabolite[i] <- strsplit(exwas_homair_28$Name,"/")[[i]][1]

exwas_homair_28$outcome <- rep("HOMA-IR",nrow(exwas_homair_28))

# InsulinAUC
exwas_insulinauc_28y <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/insulinauc/insulinauc_28y/exwas_allsig_metabolites_insulinauc_28y.csv")
exwas_insulinauc_28y <- exwas_insulinauc_28y[exwas_insulinauc_28y$age!= "Age 28",]

exwas_insulinauc_28y$log10.p.value <- -log(exwas_insulinauc_28y$p.value,10)
exwas_insulinauc_28y$age <- factor(exwas_insulinauc_28y$age, levels = c('Age 14'))

exwas_insulinauc_28y$Metabolite <- rep(NA_character_,nrow(exwas_insulinauc_28y))

for(i in 1:nrow(exwas_insulinauc_28y)) exwas_insulinauc_28y$Metabolite[i] <- strsplit(exwas_insulinauc_28y$Name,"/")[[i]][1]

exwas_insulinauc_28y$outcome <- rep("InsulinAUC",nrow(exwas_insulinauc_28y))

# Mets
exwas_mets_28_bi <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/mets/mets_28_bi/exwas_allsig_metabolites_mets_28_bi.csv")
exwas_mets_28_bi <- exwas_mets_28_bi[exwas_mets_28_bi$age!= "Age 28",]

exwas_mets_28_bi$log10.p.value <- -log(exwas_mets_28_bi$p.value,10)
exwas_mets_28_bi$age <- factor(exwas_mets_28_bi$age, levels = c('Age 7','Age 14', 'Age 22'))

exwas_mets_28_bi$Metabolite <- rep(NA_character_,nrow(exwas_mets_28_bi))

for(i in 1:nrow(exwas_mets_28_bi)) exwas_mets_28_bi$Metabolite[i] <- strsplit(exwas_mets_28_bi$Name,"/")[[i]][1]

exwas_mets_28_bi$outcome <- rep("Dyslipidemia",nrow(exwas_mets_28_bi))

# combine metabolites from all outcomes
all_outcomes_sig_metabolites<- rbind(exwas_bmi_28,
                                     exwas_bmi_28_bi,
                                     exwas_hypertension_28_bi,
                                     exwas_glucoseauc_28y,
                                     exwas_insulinauc_28y,
                                     exwas_homair_28,
                                     exwas_matsuda_28,
                                     exwas_mets_28_bi) %>% 
                               select(Metabolite, z.value, age, outcome, log10.p.value)

all_outcomes_sig_metabolites$age<- factor(all_outcomes_sig_metabolites$age,
                                          levels = c("Age 7", "Age 14", "Age 22"))

all_outcomes_sig_metabolites$outcome<- factor(all_outcomes_sig_metabolites$outcome,
                                          levels = c("BMI", "Overweight", "Hypertension", "GlucoseAUC",
                                                     "InsulinAUC", "HOMA-IR", "Matsuda", "Dyslipidemia"))


# plot
pal <- wes_palette("Zissou1", 100, type = "continuous")
vol <- (ggplot(all_outcomes_sig_metabolites, aes(x=outcome, y=Metabolite, color=z.value, size=log10.p.value, shape=age)) +# Show all points
          geom_point(position = position_dodge2(width = 0.8, preserve  = "total", padding = 0.6, reverse = T))+  
          scale_color_gradientn(colours = pal) + # change color
          scale_shape_manual(values = c(15,17,19))+ # change shape
          labs(x = "Metabolic Outcomes", y = "Metabolites", 
               title = "Prospective metabolic exposures and Metabolic outcomes at age 28",
               col = "Beta z.value") +
          theme_bw())
vol<-   (vol + 
         theme(plot.title=element_text(size=16,face="bold"),
               axis.title=element_text(size=14,face="bold"),
               plot.tag = element_text(size = 14,face = "bold"),
               axis.text.y = element_text(size=10,face="bold"),
               axis.text.x = element_text(size=10,face="bold"),
               strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
               legend.title = element_text(face = "bold"), legend.position = 'bottom',
               legend.text = element_text(size = 12),
               legend.background = element_rect(fill="white", 
                                                linewidth=0.5, linetype="solid",  colour ="darkblue"),
               plot.margin=unit(c(0.5,0.5,1,1.5), "cm"))+ 
          annotate("rect", xmin = seq(from=0.6,to=7.6,by=1), 
                    xmax = seq(from=1.4,to=8.4,by=1), 
                    ymin = rep(0,8), ymax =rep(36),
                    alpha = .1)+
          guides(shape = guide_legend(override.aes = list(size = 4))))

jpeg("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese//hrm_clinical_outcome/plot/metabolites_outcome_sig_plot.jpeg",
     units="in", width=18, height=10, res=600)
vol

dev.off()








