###############################################################################
#                UPF - HEAD & NECK AND OESOPHAGEAL CANCERS PROJECT            #
###############################################################################
#last modified: 10 Mar 2023

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
#                      Read datasets we need to merge                         #
#-----------------------------------------------------------------------------#

dataset_death <- read_sas("dq_upf_nova_mort_2022.sas7bdat", catalog_file=('formats.sas7bcat'))

dataset_cancer <- readRDS(file = "FILE.RDS")

#-----------------------------------------------------------------------------#
#                                Subset dataset                               #
#-----------------------------------------------------------------------------#

vars_keep_death <- c("Idepic", "Agexit", "Age_Recr", "C_Death_O", "Length", "Py")


vars_keep_cancer <- c("Idepic", "Sex", "Cntr_C","Height_C", "Alc_Re", "Smoke_Stat_C", "L_School_C", "Pa_Index_C", 
                      "QE_P_M_N1", "QE_P_M_N2", "QE_P_M_N3", "QE_P_M_N4", "Age_Recr_1y")


dataset_death <- dataset_death[vars_keep_death]
dataset_cancer <- dataset_cancer[vars_keep_cancer]

#-----------------------------------------------------------------------------#
#                                  Merge datasets                             #
#-----------------------------------------------------------------------------#

both_data <- merge(dataset_death, dataset_cancer, by = "Idepic", all.y = T)

both_data$Mortality_acc <- ifelse(grepl("^V|^W|^X0|^X1|^X2|^X3|^X4|^X5", both_data$C_Death_O), 1, 0)

#-----------------------------------------------------------------------------#
#                                  Save dataset                               #
#-----------------------------------------------------------------------------#

saveRDS(both_data, file = "FILE_mortality.rds") # death dataset

#-----------------------------------------------------------------------------#
#                           Read merged dataset                               ####
#-----------------------------------------------------------------------------#

data <- readRDS(file = "FILE_mortality.RDS")

#-----------------------------------------------------------------------------#
#                               Cox regressions                               #
#-----------------------------------------------------------------------------#
#NOTE: R will only run the analysis using complete cases (if missing it will drop the row)
# we use the breslow method for handling ties (default)

#-----------------------------------------------------------------------------#
#                             Create tables  %g/d                             #
#-----------------------------------------------------------------------------#
set_flextable_defaults(table.layout = "autofit")

#create dataframe with cox regression results
cox_table_func <- function(cox_object, model_num, exposure_var) {
          x <- tidy(cox_object) %>% 
            filter(., row_number()==1) %>%
            mutate(.,"estimate_10" = estimate*10, "std.error_10" = std.error*10) %>% 
            mutate(., "HR_10" = exp(estimate_10), "conf.low_10" = exp(estimate_10-1.96*std.error_10), "conf.high_10" = exp(estimate_10+1.96*std.error_10)) %>%
            add_column(., "model" = model_num, "exposure" = exposure_var, "Death_acc" = "Accidental death", "N" = cox_object$n, "Nevent" = cox_object$nevent) %>%
            mutate(conf.low_10=format(round(conf.low_10, 2), nsmall = 2), conf.high_10=format(round(conf.high_10, 2), nsmall = 2)) %>%
            mutate("95% CI" = paste(conf.low_10, "-", conf.high_10, sep = "")) %>%
            select(., model, exposure, Death_acc, N, Nevent, HR_10, "95% CI", p.value) %>%
            mutate(HR_10=round(HR_10,2), p.value=p_format(p_round(p.value), digits=3, accuracy = 0.001)) %>%
            rename(., "Model" = model, "Exposure" = exposure, "Death_acc" = Death_acc, "N total" = N, 
                   "N events" = Nevent, "HR" = HR_10, "P-value" = p.value) 
          print(x)
}

cox_func <- function(exposure, exposure_var) {
  #write function to run cox regression model 1
  x1 <- survival::coxph(Surv(Age_Recr, Agexit, Mortality_acc) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y), data = data, control = coxph.control(timefix = FALSE))
  #write function to run cox regression model 2
  x2 <- survival::coxph(Surv(Age_Recr, Agexit, Mortality_acc) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat_C) +
                         factor(L_School_C) + factor(Pa_Index_C) + Height_C, data = data, control = coxph.control(timefix = FALSE))
  #write function to run cox regression model 3
  x3 <- survival::coxph(Surv(Age_Recr, Agexit, Mortality_acc) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat_C) +
                         factor(L_School_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = data, control = coxph.control(timefix = FALSE))
  
  cox_table1 <- cox_table_func(x1, 1, exposure_var)
  cox_table2 <- cox_table_func(x2, 2, exposure_var)
  cox_table3 <- cox_table_func(x3, 3, exposure_var)
  
  cox_table_all <- rbind(cox_table1, cox_table2, cox_table3)
  
  print(cox_table_all)
  
}

#apply function to  death models
Death_table_N4 <- cox_func(data$QE_P_M_N4, "UPF intake in %g/d") %>% flextable(.)
Death_table_N3 <- cox_func(data$QE_P_M_N3, "Processed food intake in %g/d") %>% flextable(.)
Death_table_N1 <- cox_func(data$QE_P_M_N1, "Unprocessed/minimally processed food intake in %g/d") %>% flextable(.)

#-----------------------------------------------------------------------------#
#                                 Save tables                                 #
#-----------------------------------------------------------------------------#
#export tables to word
sect_properties <- prop_section(page_size = page_size(orient = "landscape"))

#save %g/d tables
Death_table_N4 %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)

Death_table_N3 %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)

Death_table_N1 %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)





###excel
cox_func_excel <- function(cox_object, model_num, exposure_var, outcome_var) {
  x <- tidy(cox_object) %>% 
    filter(., row_number()==1) %>%
    mutate(.,"b" = estimate*10, "se" = std.error*10) %>% 
    add_column(., "model" = model_num, "exposure" = exposure_var, "outcome" = outcome_var, "N" = cox_object$n, "Nevent" = cox_object$nevent) %>%
    select(., model, exposure, outcome, N, Nevent, b, "se", p.value)
  print(x)
}

cox1_Death_N4_excel <- cox_func_excel(cox1_Death_N4, 1, "Ultra-processed foods", "Accidental death")
cox2_Death_N4_excel <- cox_func_excel(cox2_Death_N4, 2, "Ultra-processed foods", "Accidental death")
cox3_Death_N4_excel <- cox_func_excel(cox3_Death_N4, 3, "Ultra-processed foods", "Accidental death")
cox1_Death_N3_excel <- cox_func_excel(cox1_Death_N3, 1, "Processed foods", "Accidental death")
cox2_Death_N3_excel <- cox_func_excel(cox2_Death_N3, 2, "Processed foods", "Accidental death")
cox3_Death_N3_excel <- cox_func_excel(cox3_Death_N3, 3, "Processed foods", "Accidental death")
cox1_Death_N1_excel <- cox_func_excel(cox1_Death_N1, 1, "Unprocessed/minimally processed foods", "Accidental death")
cox2_Death_N1_excel <- cox_func_excel(cox2_Death_N1, 2, "Unprocessed/minimally processed foods", "Accidental death")
cox3_Death_N1_excel <- cox_func_excel(cox3_Death_N1, 3, "Unprocessed/minimally processed foods", "Accidental death")

death_excel <- rbind(cox1_Death_N4_excel, cox2_Death_N4_excel, cox3_Death_N4_excel,
                        cox1_Death_N3_excel, cox2_Death_N3_excel, cox3_Death_N3_excel,
                        cox1_Death_N1_excel, cox2_Death_N1_excel, cox3_Death_N1_excel)

write.xlsx(death_excel, file = "FILE.xlsx")












################### checking if adiposity/diabetes is associated with accidental deaths in EPIC



data_adiposity <- readRDS(file = "FILE.RDS")
vars_keep <- c("Idepic", "Whr_Adj", "Bmi_Adj", "Bmi_C", "Diabet_C")
data_adiposity <- data_adiposity[vars_keep]



both_data <- merge(data, data_adiposity, by = "Idepic", all.x = T)




cox_whr <- survival::coxph(Surv(Age_Recr, Agexit, Mortality_acc) ~ Whr_Adj + strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat_C) +
                       factor(L_School_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = both_data, control = coxph.control(timefix = FALSE))
summary(cox_whr)

cox_bmi <- survival::coxph(Surv(Age_Recr, Agexit, Mortality_acc) ~ Bmi_Adj + strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat_C) +
                             factor(L_School_C) + factor(Pa_Index_C) + Height_C + Alc_Re, data = both_data, control = coxph.control(timefix = FALSE))
summary(cox_bmi)

cox_UPF_adjBMIDiabet <- survival::coxph(Surv(Age_Recr, Agexit, Mortality_acc) ~ QE_P_M_N4 + strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat_C) +
                             factor(L_School_C) + factor(Pa_Index_C) + Height_C + Alc_Re + Bmi_C + Diabet_C, data = both_data, control = coxph.control(timefix = FALSE))
summary(cox_UPF_adjBMIDiabet)
