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
#                            Create education variables                       #
#-----------------------------------------------------------------------------#

#recode education levels
table(data$L_School_C)
data$Educ_categories <- data$L_School_C
table(data$Educ_categories)
data$Educ_categories[data$Educ_categories<1.1] <- 1 
data$Educ_categories[data$Educ_categories==2 | data$Educ_categories==3] <- 2 
data$Educ_categories[data$Educ_categories==4] <- 3 
table(data$Educ_categories)

#-----------------------------------------------------------------------------#
#                      Interaction models for Educ_categories                         #
#-----------------------------------------------------------------------------#

#write functions to run cox regression model 3 with interactions for Educ_categories intake = to obtain p value for interaction
#interaction model
cox3_int_Educ_categories_func <- function(ageexit, outcome, exposure) {
  x <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + exposure*factor(Educ_categories) + strata(Sex, Cntr_C, Age_Recr_1y, Educ_categories) + 
                         factor (Smoke_Stat_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data)
  #summary(x)
}
#reference model - no interaction
cox3_ref_Educ_categories_func <- function(ageexit, outcome, exposure) {
  x <- survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y, Educ_categories) + 
                         factor (Smoke_Stat_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data)
  #summary(x)
}


#apply functions to HN models
cox3_int_Educ_categories_HN_Mal_INHANCE_2007_SiteCode_N4 <- cox3_int_Educ_categories_func(data$Agexit_Uadt, data$HN_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4)
cox3_ref_Educ_categories_HN_Mal_INHANCE_2007_SiteCode_N4 <- cox3_ref_Educ_categories_func(data$Agexit_Uadt, data$HN_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4)

#run test
anova(cox3_ref_Educ_categories_HN_Mal_INHANCE_2007_SiteCode_N4, cox3_int_Educ_categories_HN_Mal_INHANCE_2007_SiteCode_N4) 


#stratified analysis
  #among primary or less - Educ_categories variable
cox3_int_HN_Mal_INHANCE_2007_SiteCode_Educ_categories_1_N4 <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, HN_Mal_INHANCE_2007_SiteCode) ~ QE_P_M_N4 + strata(Sex, Cntr_C, Age_Recr_1y) + 
                                                                        factor (Smoke_Stat_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data[data$Educ_categories == 1, ]) 
  #among secondary or technical school - Educ_categories variable
cox3_int_HN_Mal_INHANCE_2007_SiteCode_Educ_categories_2_N4 <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, HN_Mal_INHANCE_2007_SiteCode) ~ QE_P_M_N4 + strata(Sex, Cntr_C, Age_Recr_1y) + 
                                                                        factor (Smoke_Stat_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data[data$Educ_categories == 2, ]) 
#among higher education - Educ_categories variable
cox3_int_HN_Mal_INHANCE_2007_SiteCode_Educ_categories_3_N4 <- survival::coxph(Surv(Age_Recr, Agexit_Uadt, HN_Mal_INHANCE_2007_SiteCode) ~ QE_P_M_N4 + strata(Sex, Cntr_C, Age_Recr_1y) + 
                                                                        factor (Smoke_Stat_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data[data$Educ_categories == 3, ]) 



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
int_Educ_categories_1_HN_Mal_table_N4 <- int_func(cox3_int_HN_Mal_INHANCE_2007_SiteCode_Educ_categories_1_N4, "primary education or less", "10% increase in UPF in %g/d", "Head and neck cancer")
int_Educ_categories_2_HN_Mal_table_N4 <- int_func(cox3_int_HN_Mal_INHANCE_2007_SiteCode_Educ_categories_2_N4, "secondary or technical/professional education", "10% increase in UPF in %g/d", "Head and neck cancer")
int_Educ_categories_3_HN_Mal_table_N4 <- int_func(cox3_int_HN_Mal_INHANCE_2007_SiteCode_Educ_categories_3_N4, "higher education", "10% increase in UPF in %g/d", "Head and neck cancer")

int_Educ_categories_HN_Mal_table_all <- rbind(int_Educ_categories_1_HN_Mal_table_N4, int_Educ_categories_2_HN_Mal_table_N4, int_Educ_categories_3_HN_Mal_table_N4) %>% flextable(.)


#-----------------------------------------------------------------------------#
#                                 Save tables                                 #
#-----------------------------------------------------------------------------#
#export tables to word
sect_properties <- prop_section(page_size = page_size(orient = "landscape"))

#save %g/d tables
int_Educ_categories_HN_Mal_table_all %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)









