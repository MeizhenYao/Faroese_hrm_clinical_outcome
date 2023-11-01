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

# import data
hypertension<- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/hypertension/hypertension_28_bi/allsig_pathway_hypertension_28_bi.csv")
insulinauc<- read.csv( "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/insulinauc/insulinauc_28y/allsig_pathway_insulinauc_28y.csv")
matsuda<- read.csv( "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/matsuda/matsuda_28/allsig_pathway_matsuda_28.csv")
mets<- read.csv( "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/mets/mets_28_bi/allsig_pathway_mets_28_bi.csv")
bmi<- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28/allsig_pathway_bmi_28.csv")
glucoseauc<- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/glucoseauc/glucoseauc_28y/allsig_pathway_glucoseauc_28y.csv")
homair<- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/homair/homair_28/allsig_pathway_homair_28.csv")
bmi_bi<- read.csv( "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/hrm_clinical_outcome/bmi/bmi_28_bi/allsig_pathway_bmi_28_bi.csv")


hypertension$outcome <- rep("Hypertension",nrow(hypertension))
insulinauc$outcome <- rep("InsulinAUC",nrow(insulinauc))
matsuda$outcome <- rep("Matsuda",nrow(matsuda))
mets$outcome <- rep("Dyslipidemia",nrow(mets))
bmi$outcome <- rep("BMI",nrow(bmi))
glucoseauc$outcome <- rep("GlucoseAUC",nrow(glucoseauc))
homair$outcome <- rep("HOMA-IR",nrow(homair))
bmi_bi$outcome <- rep("Overweight",nrow(bmi_bi))


# combine data
all_outcomes_sig_pathways<- rbind(hypertension,
                                  insulinauc,
                                  matsuda,
                                  mets,
                                  bmi,
                                  glucoseauc,
                                  homair,
                                  bmi_bi) %>% 
                            mutate(log_FET=-log(FET , base = 10)) %>% 
                            select(Pathways, outcome, log_FET, age, Group)

# format variables
all_outcomes_sig_pathways$age<- factor(all_outcomes_sig_pathways$age,
                                       levels = c("7", "14", "22"),
                                       labels = c("Age 7", "Age 14", "Age 22"))

all_outcomes_sig_pathways$outcome<- factor(all_outcomes_sig_pathways$outcome,
                                        levels = c("BMI", "Overweight", "Hypertension", "GlucoseAUC",
                                                   "InsulinAUC", "HOMA-IR", "Matsuda", "Dyslipidemia"))

all_outcomes_sig_pathways$Group <- factor(all_outcomes_sig_pathways$Group, levels = c("Carbohydrate metabolism",
                                                                                      "Glycan biosynthesis and metabolism",
                                                                                      "Lipid metabolism",
                                                                                      "Metabolism of amino acids and derivatives",
                                                                                      "Metabolism of vitamins and cofactors",
                                                                                      "Nucleotide metabolism",
                                                                                      "Steroid hormone biosynthesis",
                                                                                      "Xenobiotics biodegradation and metabolism"))

groups <- c("Carbohydrate metabolism",
            "Glycan biosynthesis and metabolism",
            "Lipid metabolism",
            "Metabolism of amino acids and derivatives",
            "Metabolism of vitamins and cofactors",
            "Nucleotide metabolism",
            "Steroid hormone biosynthesis",
            "Xenobiotics biodegradation and metabolism")

all_outcomes_sig_pathways <- all_outcomes_sig_pathways[order(all_outcomes_sig_pathways$Group),]
uni_paths <- all_outcomes_sig_pathways$Pathways
uni_paths <- uni_paths[!duplicated(uni_paths)]
all_outcomes_sig_pathways$Pathways <- factor(all_outcomes_sig_pathways$Pathways, levels = uni_paths)

ret <- c(length(unique(all_outcomes_sig_pathways$Pathways[all_outcomes_sig_pathways$Group ==  unique(all_outcomes_sig_pathways$Group)[1]])),
         length(unique(all_outcomes_sig_pathways$Pathways[all_outcomes_sig_pathways$Group ==  unique(all_outcomes_sig_pathways$Group)[2]])),
         length(unique(all_outcomes_sig_pathways$Pathways[all_outcomes_sig_pathways$Group ==  unique(all_outcomes_sig_pathways$Group)[3]])),
         length(unique(all_outcomes_sig_pathways$Pathways[all_outcomes_sig_pathways$Group ==  unique(all_outcomes_sig_pathways$Group)[4]])),
         length(unique(all_outcomes_sig_pathways$Pathways[all_outcomes_sig_pathways$Group ==  unique(all_outcomes_sig_pathways$Group)[5]])),
         length(unique(all_outcomes_sig_pathways$Pathways[all_outcomes_sig_pathways$Group ==  unique(all_outcomes_sig_pathways$Group)[6]])),
         length(unique(all_outcomes_sig_pathways$Pathways[all_outcomes_sig_pathways$Group ==  unique(all_outcomes_sig_pathways$Group)[7]])),
         length(unique(all_outcomes_sig_pathways$Pathways[all_outcomes_sig_pathways$Group ==  unique(all_outcomes_sig_pathways$Group)[8]])))




# plot
pal <- wes_palette("Zissou1", 100, type = "continuous")
vol <- (ggplot(all_outcomes_sig_pathways, aes(x=outcome, y=Pathways, color=log_FET, shape=age)) +# Show all points
          geom_point(size=3, position = position_dodge2(width = 0.8, preserve  = "total", padding = 0.6, reverse = T))+  
          scale_color_gradientn(colours = pal) + # change color
          scale_shape_manual(values = c(15,17,19))+ # change shape
          labs(x = "Metabolic Outcomes",
               title = "Enriched pathways in association between\nprospective serum Metabolites and Metabolic outcomes at age 28",
               col = "-log10(p.value)") +
          theme_bw()+
          theme(plot.title=element_text(size=16,face="bold"), axis.title=element_text(size=14,face="bold"),
                plot.tag = element_text(size = 13,face = "bold"),
                axis.text.y = element_text(size=10,face="bold"),
                axis.text.x = element_text(size=10,face="bold"),
                strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
                legend.title = element_text(face = "bold"), legend.position = 'bottom',
                legend.text = element_text(size = 12),
                legend.background = element_rect(fill="white", 
                                                 size=0.5, linetype="solid",  colour ="darkblue"),
                plot.margin=unit(c(0.5,0.5,1,0.5), "cm"))+ 
          annotate("rect", xmin = seq(from=0.6,to=7.6,by=1), 
                   xmax = seq(from=1.4,to=8.4,by=1), 
                   ymin = rep(0,8), ymax =rep(31),
                   alpha = .1))

vol<- (vol +
         geom_hline(yintercept = c(cumsum(ret))+0.5)+
         geom_hline(yintercept = 31, size = 1.5)+
         annotate("text", x = rep(9.5, length(unique(all_outcomes_sig_pathways$Group))), 
                  y = c(4, 9, 11.5, 17.8, 24, 27.5, 29, 30), 
                  label = c("Carbohydrate metabolism",
                            "Glycan biosynthesis and metabolism",
                            "Lipid metabolism",
                            "Metabolism of amino acids and derivatives",
                            "Metabolism of vitamins and cofactors",
                            "Nucleotide metabolism",
                            "Steroid hormone biosynthesis",
                            "Xenobiotics biodegradation and metabolism"), 
                  color = "red",
                  size = rep(3, length(groups)), 
                  fontface  = 'bold')+
         coord_cartesian(xlim = c(1.2, 10), ylim = c(0.5, 30.5), clip = "off"))

jpeg("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese//hrm_clinical_outcome/plot/pathways_outcome_sig_plot.jpeg",
     units="in", width=20, height=10, res=600)
vol

dev.off()







