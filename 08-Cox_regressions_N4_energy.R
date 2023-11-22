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
               "broom", "rms", "stats", "gtsummary", "plotrix", "rstatix", "openxlsx")

#-----------------------------------------------------------------------------#
#                                  Read dataset                               #
#-----------------------------------------------------------------------------#

data <- readRDS(file = "FILE.RDS")

#-----------------------------------------------------------------------------#
#                               Cox regressions                               #
#-----------------------------------------------------------------------------#
#NOTE: R will only run the analysis using complete cases (if missing it will drop the row)
# we use the breslow method for handling ties (default)

#-----------------------------------------------------------------------------#
#                             Create tables  %kcal/d                             #
#-----------------------------------------------------------------------------#

set_flextable_defaults(table.layout = "autofit")

#create dataframe with cox regression results
cox_table_func <- function(cox_object, model_number, exposure_var, outcome_var) {
          x_table <- tidy(cox_object) %>% 
            filter(., row_number()==1) %>%
            mutate(.,"estimate_10" = estimate*10, "std.error_10" = std.error*10) %>% 
            mutate(., "HR_10" = exp(estimate_10), "conf.low_10" = exp(estimate_10-1.96*std.error_10), "conf.high_10" = exp(estimate_10+1.96*std.error_10)) %>%
            add_column(., "model" = model_number, "exposure" = exposure_var, "outcome" = outcome_var, "N" = cox_object$n, "Nevent" = cox_object$nevent) %>%
            mutate(conf.low_10=format(round(conf.low_10, 2), nsmall = 2), conf.high_10=format(round(conf.high_10, 2), nsmall = 2)) %>%
            mutate("95% CI" = paste(conf.low_10, "-", conf.high_10, sep = "")) %>%
            select(., model, exposure, outcome, N, Nevent, HR_10, "95% CI", p.value) %>%
            mutate(HR_10=round(HR_10,2), p.value=p_format(p_round(p.value), digits=3, accuracy = 0.001)) %>%
            rename(., "Model" = model, "Exposure" = exposure, "Outcome" = outcome, "N total" = N, 
                   "N events" = Nevent, "HR" = HR_10, "P-value" = p.value) 
          print(x_table)
}


cox_func <- function(ageexit, outcome, exposure, exposure_var, outcome_var) {
  #write function to run cox regression model 1
  x1 <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y), data = data)
  #write function to run cox regression model 2
  x2 <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat_C) +
                          factor(L_School_C) + factor(Pa_Index_C) + Height_C, data = data)
  #write function to run cox regression model 3
  x3 <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat_C) +
                          factor(L_School_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data)
  
  cox_table1 <- cox_table_func(x1, 1, exposure_var, outcome_var)
  cox_table2 <- cox_table_func(x2, 2, exposure_var, outcome_var)
  cox_table3 <- cox_table_func(x3, 3, exposure_var, outcome_var)
  
  cox_table_all <- rbind(cox_table1, cox_table2, cox_table3)
  
  print(cox_table_all)
  
}

cox_Esoph_Mal_Adeno_table <- cox_func(data$Agexit_Stom, data$Esoph_Mal_Adeno, data$QE_P_E_M_N4, "UPF intake in %kcal/d", "Oesophageal adenocarcinoma")

#create table for oesophageal cancer cox results 
Esoph_Mal_morp_table <- cox_Esoph_Mal_Adeno_table %>% flextable(.)


cox_HN_Mal_INHANCE_2007_SiteCode_table <- cox_func(data$Agexit_Uadt, data$HN_Mal_INHANCE_2007_SiteCode, data$QE_P_E_M_N4, "UPF intake in %kcal/d", "Head and neck cancer (INHANCE_2007)")

#create table for head and neck cancer cox results 
HN_Mal_table <- cox_HN_Mal_INHANCE_2007_SiteCode_table %>% flextable(.)


cox_Larynx_Mal_INHANCE_2007_table <- cox_func(data$Agexit_Uadt, data$Larynx_Mal_INHANCE_2007_SiteCode, data$QE_P_E_M_N4, "UPF intake in %kcal/d", "Larynx cancer (INHANCE_2007)")
cox_Oral_Mal_INHANCE_2007_table <- cox_func(data$Agexit_Uadt, data$Oral_cavity_Mal_INHANCE_2007_SiteCode, data$QE_P_E_M_N4, "UPF intake in %kcal/d", "Oral cancer (INHANCE_2007)")
cox_Oropharynx_Mal_INHANCE_2007_SiteCode_table <- cox_func(data$Agexit_Uadt, data$Oropharynx_Mal_INHANCE_2007_SiteCode, data$QE_P_E_M_N4, "UPF intake in %kcal/d", "Oropharynx cancer (INHANCE_2007)")
cox_Hypopharynx_Mal_INHANCE_2007_SiteCode_table <- cox_func(data$Agexit_Uadt, data$Hypopharynx_Mal_INHANCE_2007_SiteCode, data$QE_P_E_M_N4, "UPF intake in %kcal/d", "Hypopharynx cancer (INHANCE_2007)")
cox_Unspecified_Mal_INHANCE_2007_SiteCode_table <- cox_func(data$Agexit_Uadt, data$Unspecified_Mal_INHANCE_2007_SiteCode, data$QE_P_E_M_N4, "UPF intake in %kcal/d", "Unspecified/overlapping cancer (INHANCE_2007)")

#create table for head and neck cancer cox results 
HN_subtypes_Mal_table <- rbind(cox_Larynx_Mal_INHANCE_2007_table, cox_Oral_Mal_INHANCE_2007_table, cox_Oropharynx_Mal_INHANCE_2007_SiteCode_table,
                               cox_Hypopharynx_Mal_INHANCE_2007_SiteCode_table, cox_Unspecified_Mal_INHANCE_2007_SiteCode_table) %>% flextable(.)



#-----------------------------------------------------------------------------#
#                             Create tables  kcal/d                              #
#-----------------------------------------------------------------------------#
set_flextable_defaults(table.layout = "autofit")

#create dataframe with cox regression results
cox_table_func_kcal <- function(cox_object, model_num, exposure_var, outcome_var) {
  x <- tidy(cox_object) %>% 
    filter(., row_number()==1) %>%
    mutate(.,"estimate_100" = estimate*100, "std.error_100" = std.error*100) %>% 
    mutate(., "HR_100" = exp(estimate_100), "conf.low_100" = exp(estimate_100-1.96*std.error_100), "conf.high_100" = exp(estimate_100+1.96*std.error_100)) %>%
    add_column(., "model" = model_num, "exposure" = exposure_var, "outcome" = outcome_var, "N" = cox_object$n, "Nevent" = cox_object$nevent) %>%
    mutate(conf.low_100=format(round(conf.low_100, 2), nsmall = 2), conf.high_100=format(round(conf.high_100, 2), nsmall = 2)) %>%
    mutate("95% CI" = paste(conf.low_100, "-", conf.high_100, sep = "")) %>%
    select(., model, exposure, outcome, N, Nevent, HR_100, "95% CI", p.value) %>%
    mutate(HR_100=round(HR_100,2), p.value=p_format(p_round(p.value), digits=3, accuracy = 0.001)) %>%
    rename(., "Model" = model, "Exposure" = exposure, "Outcome" = outcome, "N total" = N, 
           "N events" = Nevent, "HR" = HR_100, "P-value" = p.value) 
  print(x)
}

cox_func_kcal <- function(ageexit, outcome, exposure, exposure_var, outcome_var) {
  #write function to run cox regression model 1
  x1 <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y), data = data)
  #write function to run cox regression model 2
  x2 <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat_C) +
                          factor(L_School_C) + factor(Pa_Index_C) + Height_C, data = data)
  #write function to run cox regression model 3
  x3 <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat_C) +
                          factor(L_School_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data)
  
  cox_table1 <- cox_table_func_kcal(x1, 1, exposure_var, outcome_var)
  cox_table2 <- cox_table_func_kcal(x2, 2, exposure_var, outcome_var)
  cox_table3 <- cox_table_func_kcal(x3, 3, exposure_var, outcome_var)
  
  cox_table_all <- rbind(cox_table1, cox_table2, cox_table3)
  
  print(cox_table_all)
  
}


cox_Esoph_Mal_Adeno_table_kcal <- cox_func_kcal(data$Agexit_Stom, data$Esoph_Mal_Adeno, data$QE_E_M_N4, "UPF intake in kcal/d", "Oesophageal adenocarcinoma")

#create table for oesophageal cancer cox results 
Esoph_Mal_morp_table_kcal <- cox_Esoph_Mal_Adeno_table_kcal %>% flextable(.)


cox_HN_Mal_INHANCE_2007_SiteCode_table_kcal <- cox_func_kcal(data$Agexit_Uadt, data$HN_Mal_INHANCE_2007_SiteCode, data$QE_E_M_N4, "UPF intake in kcal/d", "Head and neck cancer (INHANCE_2007)")

#create table for head and neck cancer cox results 
HN_Mal_table_kcal <-cox_HN_Mal_INHANCE_2007_SiteCode_table_kcal %>% flextable(.)


cox_Larynx_Mal_INHANCE_2007_table_kcal <- cox_func_kcal(data$Agexit_Uadt, data$Larynx_Mal_INHANCE_2007_SiteCode, data$QE_E_M_N4, "UPF intake in kcal/d", "Larynx cancer (INHANCE_2007)")
cox_Oral_Mal_INHANCE_2007_table_kcal <- cox_func_kcal(data$Agexit_Uadt, data$Oral_cavity_Mal_INHANCE_2007_SiteCode, data$QE_E_M_N4, "UPF intake in kcal/d", "Oral cancer (INHANCE_2007)")
cox_Oropharynx_Mal_INHANCE_2007_SiteCode_table_kcal <- cox_func_kcal(data$Agexit_Uadt, data$Oropharynx_Mal_INHANCE_2007_SiteCode, data$QE_E_M_N4, "UPF intake in kcal/d", "Oropharynx cancer (INHANCE_2007)")
cox_Hypopharynx_Mal_INHANCE_2007_SiteCode_table_kcal <- cox_func_kcal(data$Agexit_Uadt, data$Hypopharynx_Mal_INHANCE_2007_SiteCode, data$QE_E_M_N4, "UPF intake in kcal/d", "Hypopharynx cancer (INHANCE_2007)")
cox_Unspecified_Mal_INHANCE_2007_SiteCode_table_kcal <- cox_func_kcal(data$Agexit_Uadt, data$Unspecified_Mal_INHANCE_2007_SiteCode, data$QE_E_M_N4, "UPF intake in kcal/d", "Unspecified/overlapping cancer (INHANCE_2007)")

#create table for head and neck cancer cox results 
HN_subtypes_Mal_table_kcal <- rbind(cox_Larynx_Mal_INHANCE_2007_table_kcal, cox_Oral_Mal_INHANCE_2007_table_kcal, cox_Oropharynx_Mal_INHANCE_2007_SiteCode_table_kcal,
                               cox_Hypopharynx_Mal_INHANCE_2007_SiteCode_table_kcal, cox_Unspecified_Mal_INHANCE_2007_SiteCode_table_kcal) %>% flextable(.)


#-----------------------------------------------------------------------------#
#                                 Save tables                                 #
#-----------------------------------------------------------------------------#
#export tables to word
sect_properties <- prop_section(page_size = page_size(orient = "landscape"))

#save %kcal/d tables
Esoph_Mal_morp_table %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)

HN_Mal_table %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)

HN_subtypes_Mal_table %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)


#save kcal/d tables
Esoph_Mal_morp_table_kcal %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)

HN_Mal_table_kcal %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)

HN_subtypes_Mal_table_kcal %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)













############ export to excel file %kcal/d
cox_func_excel <- function(cox_object, model_num, exposure_var, outcome_var) {
  x <- tidy(cox_object) %>% 
    filter(., row_number()==1) %>%
    mutate(.,"b" = estimate*10, "se" = std.error*10) %>% 
    add_column(., "model" = model_num, "exposure" = exposure_var, "outcome" = outcome_var, "N" = cox_object$n, "Nevent" = cox_object$nevent) %>%
    select(., model, exposure, outcome, N, Nevent, b, "se", p.value)
  print(x)
}

excel_func <- function(ageexit, outcome, exposure, exposure_var, outcome_var) {
  #write function to run cox regression model 1
  x1 <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y), data = data)
  #write function to run cox regression model 2
  x2 <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat_C) +
                          factor(L_School_C) + factor(Pa_Index_C) + Height_C, data = data)
  #write function to run cox regression model 3
  x3 <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat_C) +
                          factor(L_School_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data)
  
  cox_table1 <- cox_func_excel(x1, 1, exposure_var, outcome_var)
  cox_table2 <- cox_func_excel(x2, 2, exposure_var, outcome_var)
  cox_table3 <- cox_func_excel(x3, 3, exposure_var, outcome_var)
  
  cox_table_all <- rbind(cox_table1, cox_table2, cox_table3)
  
  print(cox_table_all)
  
}

cox_HN_Mal_excel <- excel_func(data$Agexit_Uadt, data$HN_Mal_INHANCE_2007_SiteCode, data$QE_P_E_M_N4, "Ultra-processed foods", "Head and neck cancer")
cox_Esoph_Mal_Adeno_excel <- excel_func(data$Agexit_Stom, data$Esoph_Mal_Adeno, data$QE_P_E_M_N4, "Ultra-processed foods", "Oesophageal adenocarcinoma")


HN_O_Mal_excel <- rbind(cox_HN_Mal_excel,
                        cox_Esoph_Mal_Adeno_excel)

write.xlsx(HN_O_Mal_excel, file = "FILE.xlsx")






############ export to excel file kcal
cox_func_excel_kcal <- function(cox_object, model_num, exposure_var, outcome_var) {
  x <- tidy(cox_object) %>% 
    filter(., row_number()==1) %>%
    mutate(.,"b" = estimate*100, "se" = std.error*100) %>% 
    add_column(., "model" = model_num, "exposure" = exposure_var, "outcome" = outcome_var, "N" = cox_object$n, "Nevent" = cox_object$nevent) %>%
    select(., model, exposure, outcome, N, Nevent, b, "se", p.value)
  print(x)
}

excel_func_kcal <- function(ageexit, outcome, exposure, exposure_var, outcome_var) {
  #write function to run cox regression model 1
  x1 <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y), data = data)
  #write function to run cox regression model 2
  x2 <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat_C) +
                          factor(L_School_C) + factor(Pa_Index_C) + Height_C, data = data)
  #write function to run cox regression model 3
  x3 <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat_C) +
                          factor(L_School_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data)
  
  cox_table1 <- cox_func_excel_kcal(x1, 1, exposure_var, outcome_var)
  cox_table2 <- cox_func_excel_kcal(x2, 2, exposure_var, outcome_var)
  cox_table3 <- cox_func_excel_kcal(x3, 3, exposure_var, outcome_var)
  
  cox_table_all <- rbind(cox_table1, cox_table2, cox_table3)
  
  print(cox_table_all)
  
}

cox_HN_Mal_excel_kcal <- excel_func_kcal(data$Agexit_Uadt, data$HN_Mal_INHANCE_2007_SiteCode, data$QE_E_M_N4, "Ultra-processed foods", "Head and neck cancer")
cox_Esoph_Mal_Adeno_excel_kcal <- excel_func_kcal(data$Agexit_Stom, data$Esoph_Mal_Adeno, data$QE_E_M_N4, "Ultra-processed foods", "Oesophageal adenocarcinoma")


HN_O_Mal_kcal_excel <- rbind(cox_HN_Mal_excel_kcal, 
                              cox_Esoph_Mal_Adeno_excel_kcal)

write.xlsx(HN_O_Mal_kcal_excel, file = "FILE.xlsx")






