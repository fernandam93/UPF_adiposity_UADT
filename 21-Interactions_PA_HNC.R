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
#                      Interaction models for Pa_Index_C                         #
#-----------------------------------------------------------------------------#

#write functions to run cox regression model 3 with interactions for Pa_Index_C intake = to obtain p value for interaction
#interaction model
cox3_int_Pa_Index_C_func <- function(ageexit, outcome, exposure) {
  x <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + exposure*factor(Pa_Index_C) + strata(Sex, Cntr_C, Age_Recr_1y, Pa_Index_C) + 
                         factor(Smoke_Stat_C) + factor(L_School_C) + Height_C + Alc_Re, data = data)
  #summary(x)
}
#reference model - no interaction
cox3_ref_Pa_Index_C_func <- function(ageexit, outcome, exposure) {
  x <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y, Pa_Index_C) + 
                         factor(Smoke_Stat_C) + factor(L_School_C) + Height_C + Alc_Re, data = data)
  #summary(x)
}


#apply functions to HN models
cox3_int_Pa_Index_C_HN_Mal_INHANCE_2007_SiteCode_N4 <- cox3_int_Pa_Index_C_func(data$Agexit_Uadt, data$HN_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4)
cox3_ref_Pa_Index_C_HN_Mal_INHANCE_2007_SiteCode_N4 <- cox3_ref_Pa_Index_C_func(data$Agexit_Uadt, data$HN_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4)

#run test
anova(cox3_ref_Pa_Index_C_HN_Mal_INHANCE_2007_SiteCode_N4, cox3_int_Pa_Index_C_HN_Mal_INHANCE_2007_SiteCode_N4) 


#stratified analysis
  #among inactive - Pa_Index_C variable
cox3_int_HN_Mal_INHANCE_2007_SiteCode_Pa_Index_C_1_N4 <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, HN_Mal_INHANCE_2007_SiteCode) ~ QE_P_M_N4 + strata(Sex, Cntr_C, Age_Recr_1y) + 
                                                                   factor(Smoke_Stat_C) + factor(L_School_C) + Height_C + Alc_Re, data = data[data$Pa_Index_C == 1, ]) 
  #among moderately inactive - Pa_Index_C variable
cox3_int_HN_Mal_INHANCE_2007_SiteCode_Pa_Index_C_2_N4 <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, HN_Mal_INHANCE_2007_SiteCode) ~ QE_P_M_N4 + strata(Sex, Cntr_C, Age_Recr_1y) + 
                                                                   factor(Smoke_Stat_C) + factor(L_School_C) + Height_C + Alc_Re, data = data[data$Pa_Index_C == 2, ]) 
#among moderately active - Pa_Index_C variable
cox3_int_HN_Mal_INHANCE_2007_SiteCode_Pa_Index_C_3_N4 <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, HN_Mal_INHANCE_2007_SiteCode) ~ QE_P_M_N4 + strata(Sex, Cntr_C, Age_Recr_1y) + 
                                                                   factor(Smoke_Stat_C) + factor(L_School_C) + Height_C + Alc_Re, data = data[data$Pa_Index_C == 3, ]) 
#among active - Pa_Index_C variable
cox3_int_HN_Mal_INHANCE_2007_SiteCode_Pa_Index_C_4_N4 <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, HN_Mal_INHANCE_2007_SiteCode) ~ QE_P_M_N4 + strata(Sex, Cntr_C, Age_Recr_1y) + 
                                                                   factor(Smoke_Stat_C) + factor(L_School_C) + Height_C + Alc_Re, data = data[data$Pa_Index_C == 4, ]) 


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
int_Pa_Index_C_1_HN_Mal_table_N4 <- int_func(cox3_int_HN_Mal_INHANCE_2007_SiteCode_Pa_Index_C_1_N4, "inactive", "10% increase in UPF in %g/d", "Head and neck cancer")
int_Pa_Index_C_2_HN_Mal_table_N4 <- int_func(cox3_int_HN_Mal_INHANCE_2007_SiteCode_Pa_Index_C_2_N4, "moderately inactive", "10% increase in UPF in %g/d", "Head and neck cancer")
int_Pa_Index_C_3_HN_Mal_table_N4 <- int_func(cox3_int_HN_Mal_INHANCE_2007_SiteCode_Pa_Index_C_3_N4, "moderately active", "10% increase in UPF in %g/d", "Head and neck cancer")
int_Pa_Index_C_4_HN_Mal_table_N4 <- int_func(cox3_int_HN_Mal_INHANCE_2007_SiteCode_Pa_Index_C_4_N4, "active", "10% increase in UPF in %g/d", "Head and neck cancer")


int_Pa_Index_C_HN_Mal_table_all <- rbind(int_Pa_Index_C_1_HN_Mal_table_N4, int_Pa_Index_C_2_HN_Mal_table_N4, int_Pa_Index_C_3_HN_Mal_table_N4, int_Pa_Index_C_4_HN_Mal_table_N4) %>% flextable(.)


#-----------------------------------------------------------------------------#
#                                 Save tables                                 #
#-----------------------------------------------------------------------------#
#export tables to word
sect_properties <- prop_section(page_size = page_size(orient = "landscape"))

#save %g/d tables
int_Pa_Index_C_HN_Mal_table_all %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)








