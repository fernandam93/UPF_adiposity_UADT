###############################################################################
#                UPF - HEAD & NECK AND OESOPHAGEAL CANCERS PROJECT            #
###############################################################################
#last modified: 14 March 2023

#-----------------------------------------------------------------------------#
#                                  Housekeeping                               #
#-----------------------------------------------------------------------------#

rm(list=ls()) #Remove any existing objects in R 

setwd("DIRECTORY") 

if (!require("pacman")) install.packages("pacman")
pacman::p_load("haven", "tidyverse", "dplyr", "data.table", "expss", "table1",
               "dplyr", "flextable", "magrittr", "officer", "survival", "survminer",
               "broom", "rms", "stats", "gtsummary", "plotrix", "rstatix", "ggsurvfit", "gridExtra", "corrplot", "splines", "lattice")

#-----------------------------------------------------------------------------#
#                                  Read dataset                               #
#-----------------------------------------------------------------------------#

data <- readRDS(file = "FILE.RDS")

#-----------------------------------------------------------------------------#
#            Create categorical UPF consumption variable                      #
#-----------------------------------------------------------------------------#

data$QE_P_M_N4_halfs <- ntile(data$QE_P_M_N4, 2) 

#-----------------------------------------------------------------------------#
#                        Investigate assumptions                              #
#-----------------------------------------------------------------------------#

########################## Proportionality assumption - Schoenfeld residuals #########################

proportionality_model3_HNC <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, HN_Mal_INHANCE_2007_SiteCode) ~ QE_P_M_N4 + 
                                            strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat_C) + 
                                            factor(L_School_C) + Height_C + factor(Pa_Index_C) + Alc_Re,
                                          data = data, method = "breslow") 
plot_zph_HNC <- ggcoxzph(cox.zph(proportionality_model3_HNC, transform = "km", global=TRUE))
ggsave("FILE.png", arrangeGrob(grobs = plot_zph_HNC), width = 50, height= 40, units = "cm")


proportionality_model3_adeno <- survival::coxph(Surv(Age_Recr, Agexit_Stom, Esoph_Mal_Adeno) ~ QE_P_M_N4 + 
                                                strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat_C) + 
                                                factor(L_School_C) + Height_C + factor(Pa_Index_C) + Alc_Re,
                                              data = data, method = "breslow") 
plot_zph_adeno <- ggcoxzph(cox.zph(proportionality_model3_adeno, transform = "km", global=TRUE))
ggsave("FILE.png", arrangeGrob(grobs = plot_zph_adeno), width = 50, height= 40, units = "cm")




########################## Proportionality assumption - loglog plot #########################

proportionality_model3_cat_HNC <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, HN_Mal_INHANCE_2007_SiteCode) ~ strata(QE_P_M_N4_halfs) + Sex + Cntr_C + Age_Recr_1y + 
                                                Smoke_Stat_C + L_School_C + Height_C + Pa_Index_C  + Alc_Re,
                                              data = data, method = "breslow") 
png("FILE.png", width = 2500, height = 1500, res = 300)
plot(survfit(proportionality_model3_cat_HNC), ylab="log(-log(Survival Probability)", xlab="Age in Years", 
     col=c("orange", "skyblue"), fun="cloglog", main="Head and Neck Cancer log-log Survival Plot") 
legend("bottomright", inset = .01, legend=c("Lower UPF consumption", "Higher UPF consumption"),
       col=c("orange", "skyblue"), lty=1:1, horiz=F, box.lty=0 )
dev.off()

proportionality_model3_cat_Adeno <- survival::coxph(Surv(Age_Recr, Agexit_Stom, Esoph_Mal_Adeno) ~ strata(QE_P_M_N4_halfs) + Sex + Cntr_C + Age_Recr_1y + 
                                                    Smoke_Stat_C + L_School_C + Height_C + Pa_Index_C + Alc_Re,
                                                  data = data, method = "breslow") 
png("FILE.png", width = 2500, height = 1500, res = 300)
plot(survfit(proportionality_model3_cat_Adeno), ylab="log(-log(Survival Probability)", xlab="Age in Years", 
             col=c("orange", "skyblue"), fun="cloglog", main="Oesophageal Adenocarcinoma log-log Survival Plot") 
  legend("bottomright", inset = .01, legend=c("Lower UPF consumption", "Higher UPF consumption"),
         col=c("orange", "skyblue"), lty=1:1, horiz=F, box.lty=0 )
dev.off()


########################## Correlation matrix #########################

#save correlation matrix
new_data <- data %>%  ungroup()  %>% select("UPF intake (%g/d)"=QE_P_M_N4, "Smoking status"=Smoke_Stat_C, "Education level"=L_School_C, "Height"=Height_C, "Physical activity"=Pa_Index_C, "Alcohol intake"=Alc_Re)
correlation <- cor(x = new_data)  #no categorical vars in matrix
png("FILE.png", width = 2000, height = 2000, res = 300)
corrplot(correlation, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
dev.off()

########################## Multicollinearity #########################

#save vif
proportionality_model3_continuous_HN <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, HN_Mal_INHANCE_2007_SiteCode) ~ QE_P_M_N4 + 
                                                       strata(Sex, Cntr_C, Age_Recr_1y) + Smoke_Stat_C + 
                                                       L_School_C + Height_C + Pa_Index_C + Alc_Re,
                                                     data = data, method = "breslow") 
vif_HN <- rms::vif(proportionality_model3_continuous_HN) %>% as.data.frame() %>% add_rownames()   #no evidence of multicollinearity because the Variance Inflation Factors (VIF) are close to 1
colnames(vif_HN)[2] <- "VIF for head and neck cancer"
vif_HN



proportionality_model3_continuous_adeno <- survival::coxph(Surv(Age_Recr, Agexit_Stom, Esoph_Mal_Adeno) ~ QE_P_M_N4 + 
                                                           strata(Sex, Cntr_C, Age_Recr_1y) + Smoke_Stat_C + 
                                                           L_School_C + Height_C + Pa_Index_C + Alc_Re,
                                                         data = data, method = "breslow") 
vif_adeno <- rms::vif(proportionality_model3_continuous_adeno) %>% as.data.frame() %>% add_rownames()    #no evidence of multicollinearity because the Variance Inflation Factors (VIF) are close to 1
colnames(vif_adeno)[2] <- "VIF for oesophageal adenocarcinoma"
vif_adeno

set_flextable_defaults(table.layout = "autofit")
sect_properties <- prop_section(page_size = page_size(orient = "landscape"))

vif_table <- cbind(vif_HN, vif_adeno)
vif_table[, c(3,5)] <- NULL
vif_table %>% flextable() %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)



########################## Non-linearity #########################

#check linearity HNC
base_model0_HNC <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, HN_Mal_INHANCE_2007_SiteCode) ~ QE_P_M_N4 ,
                               data = data, method = "breslow") 

splines_model0_HNC <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, HN_Mal_INHANCE_2007_SiteCode) ~ ns(QE_P_M_N4, 4) ,
                                  data = data, method = "breslow") 

anova(base_model0_HNC, splines_model0_HNC) #p=0.5418, no evidence against linearity

#create a plot to depict relationship
ND_HNC <- with(data, data.frame(QE_P_M_N4=seq(min(QE_P_M_N4),
                                          max(QE_P_M_N4), length.out=500)))

prs<- predict(splines_model0_HNC, newdata=ND_HNC, type = "lp", se.fit = T)
ND_HNC$pred <- prs[[1]]
ND_HNC$se <- prs[[2]]
ND_HNC$lo <- ND_HNC$pred - 1.96*ND_HNC$se
ND_HNC$up <- ND_HNC$pred + 1.96*ND_HNC$se

png("FILE.png", width = 2000, height = 2000, res = 300)
xyplot(pred + lo + up ~ QE_P_M_N4, data= ND_HNC,
       type="l", col= "black", lwd = 2, lty=c(1, 2, 2),
       abline=list(h=0, lty=2, lwd=2, col="red"),
       xlab = "UPF consumption (in %g/d)", ylab = "logHR")
dev.off()





















