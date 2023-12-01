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

data <- readRDS(file = "FILE.rds")


#-----------------------------------------------------------------------------#
#                      Interaction models for Smoke_Stat_C                         #
#-----------------------------------------------------------------------------#

#write functions to run cox regression model 3 with interactions for Smoke_Stat_C intake = to obtain p value for interaction
#interaction model
cox3_int_Smoke_Stat_C_func <- function(ageexit, outcome, exposure) {
  x <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + exposure*factor(Smoke_Stat_C) + strata(Sex, Cntr_C, Age_Recr_1y, Smoke_Stat_C) + 
                         factor(L_School_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data)
  #summary(x)
}
#reference model - no interaction
cox3_ref_Smoke_Stat_C_func <- function(ageexit, outcome, exposure) {
  x <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y, Smoke_Stat_C) + 
                         factor(L_School_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data)
  #summary(x)
}


#apply functions to adeno models
cox3_int_Smoke_Stat_C_Esoph_Mal_Adeno_N4 <- cox3_int_Smoke_Stat_C_func(data$Agexit_Uadt, data$Esoph_Mal_Adeno, data$QE_P_M_N4)
cox3_ref_Smoke_Stat_C_Esoph_Mal_Adeno_N4 <- cox3_ref_Smoke_Stat_C_func(data$Agexit_Uadt, data$Esoph_Mal_Adeno, data$QE_P_M_N4)

#run test
anova(cox3_ref_Smoke_Stat_C_Esoph_Mal_Adeno_N4, cox3_int_Smoke_Stat_C_Esoph_Mal_Adeno_N4) 


#stratified analysis
  #among never - Smoke_Stat_C variable
cox3_int_Esoph_Mal_Adeno_Smoke_Stat_C_1_N4 <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, Esoph_Mal_Adeno) ~ QE_P_M_N4 + strata(Sex, Cntr_C, Age_Recr_1y) + 
                                                               factor(L_School_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data[data$Smoke_Stat_C == 1, ]) 
  #among former - Smoke_Stat_C variable
cox3_int_Esoph_Mal_Adeno_Smoke_Stat_C_2_N4 <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, Esoph_Mal_Adeno) ~ QE_P_M_N4 + strata(Sex, Cntr_C, Age_Recr_1y) + 
                                                                factor(L_School_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data[data$Smoke_Stat_C == 2, ]) 
#among smoker - Smoke_Stat_C variable
cox3_int_Esoph_Mal_Adeno_Smoke_Stat_C_3_N4 <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, Esoph_Mal_Adeno) ~ QE_P_M_N4 + strata(Sex, Cntr_C, Age_Recr_1y) + 
                                                                     factor(L_School_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data[data$Smoke_Stat_C == 3, ]) 



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

#apply function to Oesophageal adenocarcinoma models
int_Smoke_Stat_C_1_adeno_Mal_table_N4 <- int_func(cox3_int_Esoph_Mal_Adeno_Smoke_Stat_C_1_N4, "Never", "10% increase in UPF in %g/d", "Oesophageal adenocarcinoma")
int_Smoke_Stat_C_2_adeno_Mal_table_N4 <- int_func(cox3_int_Esoph_Mal_Adeno_Smoke_Stat_C_2_N4, "Former", "10% increase in UPF in %g/d", "Oesophageal adenocarcinoma")
int_Smoke_Stat_C_3_adeno_Mal_table_N4 <- int_func(cox3_int_Esoph_Mal_Adeno_Smoke_Stat_C_3_N4, "Smoker", "10% increase in UPF in %g/d", "Oesophageal adenocarcinoma")



int_Smoke_Stat_C_adeno_Mal_table_all <- rbind(int_Smoke_Stat_C_1_adeno_Mal_table_N4, int_Smoke_Stat_C_2_adeno_Mal_table_N4,int_Smoke_Stat_C_3_adeno_Mal_table_N4) %>% flextable(.)


#-----------------------------------------------------------------------------#
#                                 Save tables                                 #
#-----------------------------------------------------------------------------#
#export tables to word
sect_properties <- prop_section(page_size = page_size(orient = "landscape"))

#save %g/d tables
int_Smoke_Stat_C_adeno_Mal_table_all %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)







