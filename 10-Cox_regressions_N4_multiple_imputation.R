###############################################################################
#                UPF - HEAD & NECK AND OESOPHAGEAL CANCERS PROJECT            #
###############################################################################
#last modified: 12 July 2023

#-----------------------------------------------------------------------------#
#                                  Housekeeping                               #
#-----------------------------------------------------------------------------#

rm(list=ls()) #Remove any existing objects in R 

setwd("WORKING_DIR") 

if (!require("pacman")) install.packages("pacman")
pacman::p_load("haven", "tidyverse", "dplyr", "data.table", "expss", "table1",
               "dplyr", "flextable", "magrittr", "officer", "survival", "survminer",
               "broom", "rms", "stats", "gtsummary", "plotrix", "rstatix", "openxlsx", "mice")

#-----------------------------------------------------------------------------#
#                                  Read dataset                               #
#-----------------------------------------------------------------------------#

data <- readRDS(file = "FILE.RDS")

#-----------------------------------------------------------------------------#
#                               Multiple imputation                               #
#-----------------------------------------------------------------------------#

#select variables of interest and make a new dataset
data.preimp <- select(data, Idepic, Country, Sex, Age_Recr, Agexit_Stom, Agexit_Uadt, Esoph_Mal_Adeno, Esoph_Mal_Scc, HN_Mal_SiteCode, 
                      HN_Mal_vitd_SiteCode, HN_Mal_INHANCE_2007_SiteCode, HN_Mal_INHANCE_2020_SiteCode, Larynx_Mal_INHANCE_2007_SiteCode, 
                      Pharynx_Mal_vitd_SiteCode, Oral_Mal_vitd_SiteCode, Oropharynx_Mal_vitd_SiteCode, Hypopharynx_Mal_INHANCE_2007_SiteCode, 
                      Oral_cavity_Mal_INHANCE_2007_SiteCode, Oropharynx_Mal_INHANCE_2007_SiteCode, Hypopharynx_Mal_INHANCE_2007_SiteCode, 
                      Unspecified_Mal_INHANCE_2007_SiteCode, L_School, Cntr_C, Age_Recr_1y, Pa_Index, Height_C, Smoke_Stat, Alc_Re, 
                      QE_P_M_N4, QE_M_N4, QE_P_M_N3, QE_M_N3, QE_P_M_N2, QE_M_N2, QE_P_M_N1, QE_M_N1)

#categorical variables as factors
data.preimp$L_School <- as.factor(data.preimp$L_School)
data.preimp$Country <- as.factor(data.preimp$Country)
data.preimp$Sex <- as.factor(data.preimp$Sex)
data.preimp$Cntr_C <- as.factor(data.preimp$Cntr_C)
data.preimp$Pa_Index <- as.factor(data.preimp$Pa_Index)
data.preimp$Smoke_Stat <- as.factor(data.preimp$Smoke_Stat)

#perform imputation
set.seed(123)
data.imp <- mice(data.preimp, m = 5, method = 'pmm')

#check one of the imputed datasets
head(complete(data.imp, 3))

#check imputed variables
data.imp$imp

#check completeness
data.comp <- complete(data.imp, "long", include = TRUE)
table(data.comp$.imp)
data.comp$Smoke_Stat.NA <- cci(data$Smoke_Stat)
head(data.comp[, c("Smoke_Stat", "Smoke_Stat.NA")])
# ggplot(data.comp, 
#        aes(x = .imp, y = Smoke_Stat, color = Smoke_Stat.NA)) + 
#   geom_jitter(show.legend = FALSE, 
#               width = .1)



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
cox_table_func <- function(cox_object, model_number, exposure_var, outcome_var) {
          x_table <- tidy(cox_object) %>% 
            filter(., row_number()==1) %>%
            mutate(.,"estimate_10" = estimate*10, "std.error_10" = std.error*10) %>% 
            mutate(., "HR_10" = exp(estimate_10), "conf.low_10" = exp(estimate_10-1.96*std.error_10), "conf.high_10" = exp(estimate_10+1.96*std.error_10)) %>%
            add_column(., "model" = model_number, "exposure" = exposure_var, "outcome" = outcome_var) %>%
            mutate(conf.low_10=format(round(conf.low_10, 2), nsmall = 2), conf.high_10=format(round(conf.high_10, 2), nsmall = 2)) %>%
            mutate("95% CI" = paste(conf.low_10, "-", conf.high_10, sep = "")) %>%
            select(., model, exposure, outcome, HR_10, "95% CI", p.value) %>%
            mutate(HR_10=round(HR_10,2), p.value=p_format(p_round(p.value), digits=3, accuracy = 0.001)) %>%
            rename(., "Model" = model, "Exposure" = exposure, "Outcome" = outcome, "HR" = HR_10, "P-value" = p.value)
          print(x_table)
}


cox_func <- function(ageexit, outcome, exposure, exposure_var, outcome_var) {
  #write function to run cox regression model 1
  x1 <- with(data.imp, survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y)))
  pooled_x1 <- pool(x1)
  #write function to run cox regression model 2
  x2 <- with(data.imp, survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat) +
                          factor(L_School) + factor(Pa_Index) + Height_C))
  pooled_x2 <- pool(x2)
  #write function to run cox regression model 3
  x3 <- with(data.imp, survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat) +
                          factor(L_School) + factor(Pa_Index) + Height_C + Alc_Re))
  pooled_x3 <- pool(x3)
  
  cox_table1 <- cox_table_func(pooled_x1, 1, exposure_var, outcome_var)
  cox_table2 <- cox_table_func(pooled_x2, 2, exposure_var, outcome_var)
  cox_table3 <- cox_table_func(pooled_x3, 3, exposure_var, outcome_var)
  
  cox_table_all <- rbind(cox_table1, cox_table2, cox_table3)
  
  print(cox_table_all)
  
}

cox_Esoph_Mal_Adeno_table <- cox_func(data$Agexit_Stom, data$Esoph_Mal_Adeno, data$QE_P_M_N4, "UPF intake in %g/d", "Oesophageal adenocarcinoma")

#create table for oesophageal cancer cox results 
Esoph_Mal_morp_table <- cox_Esoph_Mal_Adeno_table %>% flextable(.)


cox_HN_Mal_INHANCE_2007_SiteCode_table <- cox_func(data$Agexit_Uadt, data$HN_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4, "UPF intake in %g/d", "Head and neck cancer (INHANCE_2007)")

#create table for head and neck cancer cox results 
HN_Mal_table <- cox_HN_Mal_INHANCE_2007_SiteCode_table %>% flextable(.)


cox_Larynx_Mal_INHANCE_2007_table <- cox_func(data$Agexit_Uadt, data$Larynx_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4, "UPF intake in %g/d", "Larynx cancer (INHANCE_2007)")
cox_Oral_Mal_INHANCE_2007_table <- cox_func(data$Agexit_Uadt, data$Oral_cavity_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4, "UPF intake in %g/d", "Oral cancer (INHANCE_2007)")
cox_Oropharynx_Mal_INHANCE_2007_SiteCode_table <- cox_func(data$Agexit_Uadt, data$Oropharynx_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4, "UPF intake in %g/d", "Oropharynx cancer (INHANCE_2007)")
cox_Hypopharynx_Mal_INHANCE_2007_SiteCode_table <- cox_func(data$Agexit_Uadt, data$Hypopharynx_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4, "UPF intake in %g/d", "Hypopharynx cancer (INHANCE_2007)")
cox_Unspecified_Mal_INHANCE_2007_SiteCode_table <- cox_func(data$Agexit_Uadt, data$Unspecified_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4, "UPF intake in %g/d", "Unspecified/overlapping cancer (INHANCE_2007)")

#create table for head and neck cancer cox results 
HN_subtypes_Mal_table <- rbind(cox_Larynx_Mal_INHANCE_2007_table, cox_Oral_Mal_INHANCE_2007_table, cox_Oropharynx_Mal_INHANCE_2007_SiteCode_table,
                               cox_Hypopharynx_Mal_INHANCE_2007_SiteCode_table, cox_Unspecified_Mal_INHANCE_2007_SiteCode_table) %>% flextable(.)

#-----------------------------------------------------------------------------#
#                                 Save tables                                 #
#-----------------------------------------------------------------------------#
#export tables to word
sect_properties <- prop_section(page_size = page_size(orient = "landscape"))

#save %g/d tables
Esoph_Mal_morp_table %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)

HN_Mal_table %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)

HN_subtypes_Mal_table %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)













############ export to excel file %g/d
cox_func_excel <- function(cox_object, model_num, exposure_var, outcome_var) {
  x <- tidy(cox_object) %>% 
    filter(., row_number()==1) %>%
    mutate(.,"b" = estimate*10, "se" = std.error*10) %>% 
    add_column(., "model" = model_num, "exposure" = exposure_var, "outcome" = outcome_var) %>%
    select(., model, exposure, outcome, b, "se", p.value)
  print(x)
}

excel_func <- function(ageexit, outcome, exposure, exposure_var, outcome_var) {
  #write function to run cox regression model 1
  x1 <- with(data.imp, survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y)))
  pooled_x1 <- pool(x1)
  #write function to run cox regression model 2
  x2 <- with(data.imp, survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat) +
                                         factor(L_School) + factor(Pa_Index) + Height_C))
  pooled_x2 <- pool(x2)
  #write function to run cox regression model 3
  x3 <- with(data.imp, survival::coxph(Surv(Age_Recr, ageexit, outcome) ~ exposure + strata(Sex, Cntr_C, Age_Recr_1y) + factor(Smoke_Stat) +
                                         factor(L_School) + factor(Pa_Index) + Height_C + Alc_Re))
  pooled_x3 <- pool(x3)
  
  cox_table1 <- cox_func_excel(pooled_x1, 1, exposure_var, outcome_var)
  cox_table2 <- cox_func_excel(pooled_x2, 2, exposure_var, outcome_var)
  cox_table3 <- cox_func_excel(pooled_x3, 3, exposure_var, outcome_var)
  
  cox_table_all <- rbind(cox_table1, cox_table2, cox_table3)
  
  print(cox_table_all)
  
}

cox_HN_Mal_excel <- excel_func(data$Agexit_Uadt, data$HN_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4, "Ultra-processed foods", "Head and neck cancer")
cox_Esoph_Mal_Adeno_excel <- excel_func(data$Agexit_Stom, data$Esoph_Mal_Adeno, data$QE_P_M_N4, "Ultra-processed foods", "Oesophageal adenocarcinoma")


HN_O_Mal_excel <- rbind(cox_HN_Mal_excel,
                        cox_Esoph_Mal_Adeno_excel)

write.xlsx(HN_O_Mal_excel, file = "FILE.xlsx")





cox_HN_oral_Mal_excel <- excel_func(data$Agexit_Uadt, data$Oral_cavity_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4, "Ultra-processed foods", "Oral cavity cancer")
cox_HN_oropharynx_Mal_excel <- excel_func(data$Agexit_Uadt, data$Oropharynx_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4, "Ultra-processed foods", "Oropharynx cancer")
cox_HN_hypopharynx_Mal_excel <- excel_func(data$Agexit_Uadt, data$Hypopharynx_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4, "Ultra-processed foods", "Hypopharynx cancer")
cox_HN_larynx_Mal_excel <- excel_func(data$Agexit_Uadt, data$Larynx_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4, "Ultra-processed foods", "Larynx cancer")
cox_HN_unspecified_Mal_excel <- excel_func(data$Agexit_Uadt, data$Unspecified_Mal_INHANCE_2007_SiteCode, data$QE_P_M_N4, "Ultra-processed foods", "Unspecified/overlapping cancer")

HN_subtypes_excel <- rbind(cox_HN_oral_Mal_excel,
                           cox_HN_oropharynx_Mal_excel, 
                           cox_HN_hypopharynx_Mal_excel, 
                           cox_HN_larynx_Mal_excel,
                           cox_HN_unspecified_Mal_excel)

write.xlsx(HN_subtypes_excel, file = "FILE.xlsx")









