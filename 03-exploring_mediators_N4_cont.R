###############################################################################
#                UPF - HEAD & NECK AND OESOPHAGEAL CANCERS PROJECT            #
###############################################################################
#last modified: 13 Mar 2023

#-----------------------------------------------------------------------------#
#                                  Housekeeping                               #
#-----------------------------------------------------------------------------#

rm(list=ls()) #Remove any existing objects in R 

setwd("WORKING_DIR") 

if (!require("pacman")) install.packages("pacman")
pacman::p_load("haven", "tidyverse", "dplyr", "data.table", "expss", "table1",
               "dplyr", "flextable", "magrittr", "officer", "survival", "survminer",
               "broom", "rms", "stats", "gtsummary", "plotrix", "rstatix")

#-----------------------------------------------------------------------------#
#                                  Read dataset                               #
#-----------------------------------------------------------------------------#

data <- readRDS(file = "FILE.RDS")

#-----------------------------------------------------------------------------#
#                      Exposure-mediator associations                         #
#-----------------------------------------------------------------------------#

#linear regression exp- continuous mediator
lin_expmed_func<- function(exposure, mediator) {
  x <- glm(mediator ~ exposure + Sex + factor(Cntr_C) + Age_Recr_1y + factor(Smoke_Stat_C) +
                                 factor(L_School_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data)
}

#WHR
lin_expmed_percentgrams_whr <- lin_expmed_func(data$QE_P_M_N4, data$Whr_Adj*100)
lin_expmed_percentgrams_whr$coefficients[[2]]

#BMI
lin_expmed_percentgrams_bmi <- lin_expmed_func(data$QE_P_M_N4, data$Bmi_Adj)
lin_expmed_percentgrams_bmi$coefficients[[2]]


#-----------------------------------------------------------------------------#
#                         Mediator-outcome associations                       #
#-----------------------------------------------------------------------------#

#write function to run basic cox regression model
cox_medout_func <- function(ageexit, outcome, exposure, mediator) {
  x <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ mediator + strata(Sex, Cntr_C, Age_Recr_1y) + exposure + 
                         factor(Smoke_Stat_C) + factor(L_School_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data)
  #summary(x)
}

#WHR
cox_medout_whr_HN_Mal_INHANCE_2007_SiteCode <- cox_medout_func(data$Agexit_Uadt, data$HN_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4, data$Whr_Adj*100)
cox_medout_whr_Esoph_Mal_Adeno <- cox_medout_func(data$Agexit_Stom, data$Esoph_Mal_Adeno, data$QE_P_M_N4, data$Whr_Adj*100)

#BMI
cox_medout_bmi_HN_Mal_INHANCE_2007_SiteCode <- cox_medout_func(data$Agexit_Uadt, data$HN_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4, data$Bmi_Adj)
cox_medout_bmi_Esoph_Mal_Adeno <- cox_medout_func(data$Agexit_Stom, data$Esoph_Mal_Adeno, data$QE_P_M_N4, data$Bmi_Adj)

#-----------------------------------------------------------------------------#
#                     Direct effects (adjusting for mediators)                #
#-----------------------------------------------------------------------------#

cox_expout_medadjusted_func <- function(ageexit, outcome, exposure, mediator) {
  x <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y) + mediator + 
                         factor(Smoke_Stat_C) + factor(L_School_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data)
  #summary(x)
}

#WHR
cox_expout_medadjusted_whr_HN_Mal_INHANCE_2007_SiteCode <- cox_expout_medadjusted_func(data$Agexit_Uadt, data$HN_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4, data$Whr_Adj*100)
cox_expout_medadjusted_whr_Esoph_Mal_Adeno <- cox_expout_medadjusted_func(data$Agexit_Stom, data$Esoph_Mal_Adeno, data$QE_P_M_N4, data$Whr_Adj*100)

#BMI
cox_expout_medadjusted_bmi_HN_Mal_INHANCE_2007_SiteCode <- cox_expout_medadjusted_func(data$Agexit_Uadt, data$HN_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4, data$Bmi_Adj)
cox_expout_medadjusted_bmi_Esoph_Mal_Adeno <- cox_expout_medadjusted_func(data$Agexit_Stom, data$Esoph_Mal_Adeno, data$QE_P_M_N4, data$Bmi_Adj)



#-----------------------------------------------------------------------------#
#                                    Create tables                            #
#-----------------------------------------------------------------------------#

set_flextable_defaults(table.layout = "autofit")

#create dataframe with cox regression results
cox_func_mediator <- function(cox_object, model_type, mediator_var, outcome_var) {
  x <- tidy(cox_object) %>% 
    filter(., row_number()==1) %>%
    mutate(., "HR" = exp(estimate), "conf.low" = exp(estimate-1.96*std.error), "conf.high" = exp(estimate+1.96*std.error)) %>%
    add_column(., "model" = model_type, "mediator" = mediator_var, "outcome" = outcome_var, "N" = cox_object$n, "Nevent" = cox_object$nevent) %>%
    mutate(conf.low=format(round(conf.low, 2), nsmall = 2), conf.high=format(round(conf.high, 2), nsmall = 2)) %>%
    mutate("95% CI" = paste(conf.low, "-", conf.high, sep = "")) %>%
    select(., model, mediator, outcome, N, Nevent, HR, "95% CI", p.value) %>%
    mutate(HR=round(HR,2), p.value=p_format(p_round(p.value), digits=3, accuracy = 0.001)) %>%
    rename(., "Model" = model, "Mediator" = mediator, "Outcome" = outcome, "N total" = N, 
           "N events" = Nevent, "P-value" = p.value) 
  print(x)
}

#apply function to whr models
cox_medout_whr_HN_Mal_INHANCE_2007_SiteCode_table <- cox_func_mediator(cox_medout_whr_HN_Mal_INHANCE_2007_SiteCode, "Cox regression w/continuous mediator", "Waist-to-hip ratio", "Head and neck cancer")
cox_medout_whr_Esoph_Mal_Adeno_table <- cox_func_mediator(cox_medout_whr_Esoph_Mal_Adeno, "Cox regression w/continuous mediator", "Waist-to-hip ratio", "Oesophageal adenocarcinoma")

#apply function to bmi models
cox_medout_bmi_HN_Mal_INHANCE_2007_SiteCode_table <- cox_func_mediator(cox_medout_bmi_HN_Mal_INHANCE_2007_SiteCode, "Cox regression w/continuous mediator", "Body mass index", "Head and neck cancer")
cox_medout_bmi_Esoph_Mal_Adeno_table <- cox_func_mediator(cox_medout_bmi_Esoph_Mal_Adeno, "Cox regression w/continuous mediator", "Body mass index", "Oesophageal adenocarcinoma")

cox_table_med<- rbind(cox_medout_whr_HN_Mal_INHANCE_2007_SiteCode_table, cox_medout_whr_Esoph_Mal_Adeno_table,
                      cox_medout_bmi_HN_Mal_INHANCE_2007_SiteCode_table, cox_medout_bmi_Esoph_Mal_Adeno_table)  %>% flextable(.)


#create dataframe with linear regression results
lin_func_mediator <- function(cox_object, model_type, exposure_var, mediator_var) {
  x <- tidy(cox_object) %>% 
    filter(., row_number()==2) %>%
    mutate(.,"estimate" = estimate*10, "std.error" = std.error*10) %>%
    mutate(., "coefficient" = estimate, "conf.low" = (estimate-1.96*std.error), "conf.high" = (estimate+1.96*std.error)) %>%
    add_column(., "model" = model_type, "mediator" = mediator_var, "exposure" = exposure_var) %>%
    mutate(conf.low=format(round(conf.low, 2), nsmall = 2), conf.high=format(round(conf.high, 2), nsmall = 2)) %>%
    mutate("95% CI" = paste(conf.low, "-", conf.high, sep = "")) %>%
    select(., model, exposure, mediator, coefficient, "95% CI", p.value) %>%
    mutate(coefficient=round(coefficient,2), p.value=p_format(p_round(p.value), digits=3, accuracy = 0.001)) %>%
    rename(., "Model" = model, "Mediator" = mediator, "Exposure" = exposure, 
           "P-value" = p.value, "Coefficient" = coefficient) 
  print(x)
}

lin_expmed_percentgrams_whr_table <- lin_func_mediator(lin_expmed_percentgrams_whr, "Linear regression", "UPF intake in %g/d", "Waist-to-hip ratio")
lin_expmed_percentgrams_bmi_table <- lin_func_mediator(lin_expmed_percentgrams_bmi, "Linear regression", "UPF intake in %g/d", "Body mass index")
 
lin_table <- rbind(lin_expmed_percentgrams_whr_table,
                   lin_expmed_percentgrams_bmi_table) %>% flextable(.)



cox_func_direct <- function(cox_object, model_type, exposure_var, mediator_var, outcome_var) {
  x <- tidy(cox_object) %>% 
    filter(., row_number()==1) %>%
    mutate(.,"estimate" = estimate*10, "std.error" = std.error*10) %>%
    mutate(., "HR" = exp(estimate), "conf.low" = exp(estimate-1.96*std.error), "conf.high" = exp(estimate+1.96*std.error)) %>%
    add_column(., "model" = model_type, "exposure" = exposure_var, "mediator" = mediator_var, "outcome" = outcome_var, "N" = cox_object$n, "Nevent" = cox_object$nevent) %>%
    mutate(conf.low=format(round(conf.low, 2), nsmall = 2), conf.high=format(round(conf.high, 2), nsmall = 2)) %>%
    mutate("95% CI" = paste(conf.low, "-", conf.high, sep = "")) %>%
    select(., model, exposure, mediator, outcome, N, Nevent, HR, "95% CI", p.value) %>%
    mutate(HR=round(HR,2), p.value=p_format(p_round(p.value), digits=3, accuracy = 0.001)) %>%
    rename(., "Model" = model, "Exposure" = exposure, "Mediator" = mediator, "Outcome" = outcome, "N total" = N, 
           "N events" = Nevent, "P-value" = p.value) 
  print(x)
}


#apply function to whr models
cox_expout_medadjusted_whr_HN_Mal_INHANCE_2007_SiteCode_table <- cox_func_direct(cox_expout_medadjusted_whr_HN_Mal_INHANCE_2007_SiteCode, "Cox regression adjusted for mediator", "UPF intake in %g/d", "Waist-to-hip ratio", "Head and neck cancer")
cox_expout_medadjusted_whr_Esoph_Mal_Adeno_table <- cox_func_direct(cox_expout_medadjusted_whr_Esoph_Mal_Adeno, "Cox regression adjusted for mediator", "UPF intake in %g/d", "Waist-to-hip ratio", "Oesophageal adenocarcinoma")
 
#apply function to bmi models
cox_expout_medadjusted_bmi_HN_Mal_INHANCE_2007_SiteCode_table <- cox_func_direct(cox_expout_medadjusted_bmi_HN_Mal_INHANCE_2007_SiteCode, "Cox regression adjusted for mediator", "UPF intake in %g/d", "Body mass index", "Head and neck cancer")
cox_expout_medadjusted_bmi_Esoph_Mal_Adeno_table <- cox_func_direct(cox_expout_medadjusted_bmi_Esoph_Mal_Adeno, "Cox regression adjusted for mediator", "UPF intake in %g/d", "Body mass index", "Oesophageal adenocarcinoma")

cox_table_direct <- rbind(cox_expout_medadjusted_whr_HN_Mal_INHANCE_2007_SiteCode_table, cox_expout_medadjusted_whr_Esoph_Mal_Adeno_table,
                          cox_expout_medadjusted_bmi_HN_Mal_INHANCE_2007_SiteCode_table, cox_expout_medadjusted_bmi_Esoph_Mal_Adeno_table) %>% flextable(.)

#-----------------------------------------------------------------------------#
#                                 Save tables                                 #
#-----------------------------------------------------------------------------#
#export tables to word
sect_properties <- prop_section(page_size = page_size(orient = "landscape"))

#save %g/d tables
lin_table %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)
cox_table_med %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)
cox_table_direct %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)





