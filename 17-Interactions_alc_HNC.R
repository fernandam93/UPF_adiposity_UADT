###############################################################################
#                UPF - HEAD & NECK AND OESOPHAGEAL CANCERS PROJECT            #
###############################################################################

#EXPLORING INTERACTIONS 

#last modified: 29 April 2023

#-----------------------------------------------------------------------------#
#                                  Housekeeping                               #
#-----------------------------------------------------------------------------#

rm(list=ls()) #Remove any existing objects in R 

setwd("DIRECTORY") 

if (!require("pacman")) install.packages("pacman")
pacman::p_load("haven", "tidyverse", "dplyr", "data.table", "expss", "table1",
               "dplyr", "flextable", "magrittr", "officer", "survival", "survminer",
               "broom", "rms", "stats", "gtsummary", "plotrix", "rstatix")

#-----------------------------------------------------------------------------#
#                                  Read dataset                               #
#-----------------------------------------------------------------------------#

data <- readRDS(file = "FILE.RDS")

#-----------------------------------------------------------------------------#
#                            Create alcohol variables                          #
#-----------------------------------------------------------------------------#
table(data$Alc_Re_C)

#categories as in https://onlinelibrary.wiley.com/doi/full/10.1002/ijc.29559
#1=non-drinkers, very light or occasional drinkers [0.1–6 (men); 0.1-3 (women) g/d]. 
#2=moderate low drinkers [6.1–12 (men); 3.1-12 (women) g/d]. moderate high drinkers [12.1-24 g/d]. 
#3=heavy drinkers [24.1-60 g/d]. very heavy drinkers [>60 g/d].
data$Alc_categories <- data$Alc_Re_C
table(data$Alc_categories)
data$Alc_categories[data$Alc_categories<2.1] <- 1 
data$Alc_categories[data$Alc_categories==3 | data$Alc_categories==4] <- 2 
data$Alc_categories[data$Alc_categories==5 | data$Alc_categories==6 | data$Alc_categories==7] <- 3 
table(data$Alc_categories)


#-----------------------------------------------------------------------------#
#                      Interaction models for alcohol                         #
#-----------------------------------------------------------------------------#

#write functions to run cox regression model 2 with interactions for alcohol intake = to obtain p value for interaction
#interaction model
cox2_int_alc_func <- function(ageexit, outcome, exposure) {
  x <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + exposure*factor(Alc_categories) + strata(Sex, Cntr_C, Age_Recr_1y, Alc_categories) + 
                         factor (Smoke_Stat_C) + factor(L_School_C) + factor(Pa_Index_C) + Height_C, data = data)
  #summary(x)
}
#reference model - no interaction
cox2_ref_alc_func <- function(ageexit, outcome, exposure) {
  x <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y, Alc_categories) + 
                         factor (Smoke_Stat_C) + factor(L_School_C) + factor(Pa_Index_C) + Height_C, data = data)
  #summary(x)
}


#apply functions to HN models
cox2_int_alc_HN_Mal_INHANCE_2007_SiteCode_N4 <- cox2_int_alc_func(data$Agexit_Uadt, data$HN_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4)
cox2_ref_alc_HN_Mal_INHANCE_2007_SiteCode_N4 <- cox2_ref_alc_func(data$Agexit_Uadt, data$HN_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4)

#run test
anova(cox2_ref_alc_HN_Mal_INHANCE_2007_SiteCode_N4, cox2_int_alc_HN_Mal_INHANCE_2007_SiteCode_N4) #no evidence of interaction based on p value


#stratified analysis
  #among non-drinkers and light drinkers - alc_categories variable
cox2_int_HN_Mal_INHANCE_2007_SiteCode_alc_0_light_N4 <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, HN_Mal_INHANCE_2007_SiteCode) ~ QE_P_M_N4 + strata(Sex, Cntr_C, Age_Recr_1y) + 
                                                                  factor (Smoke_Stat_C) + factor(L_School_C) + factor(Pa_Index_C) + Height_C, data = data[data$Alc_categories == 1, ]) #1.02 per 1% increase in UPF
  #among moderate drinkers - alc_categories variable
cox2_int_HN_Mal_INHANCE_2007_SiteCode_alc_moderate_N4 <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, HN_Mal_INHANCE_2007_SiteCode) ~ QE_P_M_N4 + strata(Sex, Cntr_C, Age_Recr_1y) + 
                                                                   factor (Smoke_Stat_C) + factor(L_School_C) + factor(Pa_Index_C) + Height_C, data = data[data$Alc_categories == 2, ]) #1.03 per 1% increase in UPF
  #among heavy drinkers - alc_categories variable
cox2_int_HN_Mal_INHANCE_2007_SiteCode_alc_heavy_N4 <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, HN_Mal_INHANCE_2007_SiteCode) ~ QE_P_M_N4 + strata(Sex, Cntr_C, Age_Recr_1y) + 
                                                                factor (Smoke_Stat_C) + factor(L_School_C) + factor(Pa_Index_C) + Height_C, data = data[data$Alc_categories == 3, ]) # 1.00 per 1% increase in UPF

#-----------------------------------------------------------------------------#
#                             Create tables  %g/d                             #
#-----------------------------------------------------------------------------#
set_flextable_defaults(table.layout = "autofit")

#create dataframe with interaction results
int_func <- function(cox_object, interaction_level, exposure_var, outcome_var) {
  x <- tidy(cox_object) %>% 
    filter(., row_number()==1) %>%
    mutate(.,"estimate_10" = estimate*10, "std.error_10" = std.error*10) %>% 
    mutate(., "HR_10" = exp(estimate_10), "conf.low_10" = exp(estimate_10-1.96*std.error_10), "conf.high_10" = exp(estimate_10+1.96*std.error_10)) %>%
    add_column(., "model" = interaction_level, "exposure" = exposure_var, "outcome" = outcome_var, "N" = cox_object$n, "Nevent" = cox_object$nevent) %>%
    mutate(conf.low_10=format(round(conf.low_10, 2), nsmall = 2), conf.high_10=format(round(conf.high_10, 2), nsmall = 2)) %>%
    mutate("95% CI" = paste(conf.low_10, "-", conf.high_10, sep = "")) %>%
    select(., model, exposure, outcome, N, Nevent, HR_10, "95% CI", p.value) %>%
    mutate(HR_10=round(HR_10,2), p.value=p_format(p_round(p.value), digits=3, accuracy = 0.001)) %>%
    rename(., "Interaction level" = model, "Exposure" = exposure, "Outcome" = outcome, "N total" = N, 
           "N events" = Nevent, "HR" = HR_10, "P-value" = p.value) 
  print(x)
}

#apply function to head and neck cancer models
int_alc_light_HN_Mal_table_N4 <- int_func(cox2_int_HN_Mal_INHANCE_2007_SiteCode_alc_0_light_N4, "no/light alcohol intake", "10% increase in UPF in %g/d", "Head and neck cancer")
int_alc_moderate_HN_Mal_table_N4 <- int_func(cox2_int_HN_Mal_INHANCE_2007_SiteCode_alc_moderate_N4, "moderate alcohol intake", "10% increase in UPF in %g/d", "Head and neck cancer")
int_alc_heavy_HN_Mal_table_N4 <- int_func(cox2_int_HN_Mal_INHANCE_2007_SiteCode_alc_heavy_N4, "heavy alcohol intake", "10% increase in UPF in %g/d", "Head and neck cancer")



int_alc_HN_Mal_table_all <- rbind(int_alc_light_HN_Mal_table_N4, int_alc_moderate_HN_Mal_table_N4, int_alc_heavy_HN_Mal_table_N4) %>% flextable(.)


#-----------------------------------------------------------------------------#
#                                 Save tables                                 #
#-----------------------------------------------------------------------------#
#export tables to word
sect_properties <- prop_section(page_size = page_size(orient = "landscape"))

#save %g/d tables
int_alc_HN_Mal_table_all %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)












