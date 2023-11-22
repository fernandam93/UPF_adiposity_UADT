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
p_load_gh("BS1125/CMAverse")

#-----------------------------------------------------------------------------#
#                                  Read dataset                               #
#-----------------------------------------------------------------------------#

data <- readRDS(file = "FILE.RDS")

#-----------------------------------------------------------------------------#
#                               Mediation analysis                            #
#-----------------------------------------------------------------------------#

#create a smaller dataset
z <- data[, c("Sex", "Cntr_C", "Age_Recr_1y", "QE_P_M_N4", "Height_C", "Pa_Index_C", "L_School_C", 
              "Smoke_Stat_C", "Bmi_Adj", "Whr_Adj", "Hypert", "MRMed_Score_C", "QE_US_ENERGY", 
              "Alc_Re", "HN_Mal_INHANCE_2007_SiteCode", "Esoph_Mal_Adeno", "Esoph_Mal_Scc", "Py_Uadt", "Py_Stom")]

#rescale UPF intake so that HR is per 10% increase
z$QE_P_M_N4 <- z$QE_P_M_N4/10

#set categorical variables as factors
z$Cntr_C <- as.factor(z$Cntr_C)
z$L_School_C <- as.factor(z$L_School_C)
z$Smoke_Stat_C <- as.factor(z$Smoke_Stat_C)
z$Pa_Index_C <- as.factor(z$Pa_Index_C)

######################################################
####### MODELS FOR CONTINUOUS EXPOSURES ##############
######################################################

############Head and neck cancer

#whr 
med_Whr_Adj_HN_0.76 <- CMAverse::cmest(data = z[complete.cases(z[, 10]), ], exposure = "QE_P_M_N4", mediator = "Whr_Adj", 
                                           basec = c("Sex", "Cntr_C", "Age_Recr_1y", "Height_C", "Pa_Index_C", "L_School_C", 
                                                     "Smoke_Stat_C", "Alc_Re"), yreg = "coxph", full = T, casecontrol = F,
                                           mreg = list("linear"), outcome = "Py_Uadt", astar = 0, a=1, inference = "bootstrap", nboot = 1000, estimation = "imputation",
                                           EMint = T, mval = list(0.76), model = "rb", event = "HN_Mal_INHANCE_2007_SiteCode")

############Oesophageal adenocarcinoma

  #bmi
med_Bmi_Adj_adeno_22.68 <- CMAverse::cmest(data = z[complete.cases(z[, 9]), ], exposure = "QE_P_M_N4", mediator = "Bmi_Adj", 
                                     basec = c("Sex", "Cntr_C", "Age_Recr_1y", "Height_C", "Pa_Index_C", "L_School_C", 
                                               "Smoke_Stat_C", "Alc_Re"), yreg = "coxph", full = T, casecontrol = F,
                                     mreg = list("linear"), outcome = "Py_Stom", astar = 0, a=1, inference = "bootstrap", nboot = 1000, estimation = "imputation",
                                     EMint = T, mval = list(22.68), model = "rb", event = "Esoph_Mal_Adeno")

#whr
med_Whr_Adj_adeno_0.76 <- CMAverse::cmest(data = z[complete.cases(z[, 10]), ], exposure = "QE_P_M_N4", mediator = "Whr_Adj", 
                                            basec = c("Sex", "Cntr_C", "Age_Recr_1y", "Height_C", "Pa_Index_C", "L_School_C", 
                                                      "Smoke_Stat_C", "Alc_Re"), yreg = "coxph", full = T, casecontrol = F,
                                            mreg = list("linear"), outcome = "Py_Stom", astar = 0, a=1, inference = "bootstrap", nboot = 1000, estimation = "imputation",
                                            EMint = T, mval = list(0.76), model = "rb", event = "Esoph_Mal_Adeno")

#-----------------------------------------------------------------------------#
#                                 Create tables                                 #
#-----------------------------------------------------------------------------#

set_flextable_defaults(table.layout = "autofit")

#create tables to then save
mediation_table_func <- function(x, exposure_var, mediator_var, outcome_var) {
  effect_df <- tibble(Effect = factor(names(x$effect.pe), levels = names(x$effect.pe)),
                        Estimate = x$effect.pe, CI_lower = x$effect.ci.low,
                        CI_upper = x$effect.ci.high, P = x$effect.pval)  
  effect_df %>% add_column(., "Exposure" = exposure_var, "Mediator" = mediator_var, "Outcome" = outcome_var, "Mval" = as.numeric(x$ref[[3]]), .before = "Effect") %>%
                mutate(Estimate=round(Estimate,2), CI_upper=format(round(CI_upper, digits=2), nsmall=2), CI_lower=format(round(CI_lower, digits=2), nsmall = 2) , 
                      P=p_format(p_round(P), digits=3, accuracy = 0.001)) %>% 
    mutate("95% CI" = paste(CI_lower, "-", CI_upper, sep = "")) %>%
    select(., Exposure, Mediator, Outcome, Mval, Effect, Estimate, "95% CI", P) %>% flextable(.) %>% 
    footnote(., i = 1, j = 4,
             value = as_paragraph(
               c("Rcde: controlled direct effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnde: total natural direct effect hazard ratio; Rpnie: pure natural indirect effect hazard ratio; Rtnie: total natural indirect effect hazard ratio; Rte: total effect hazard ratio; ERcde: excess relative hazard due to controlled direct effect; ERintref: excess relative hazard due to reference interaction; ERintmed: excess relative hazard due to mediated interaction; ERpnie: excess relative hazard due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated")
          ),
          ref_symbols = c(""),
          part = "header")
  
}



#continuous exposure
#hn
med_int_Whr_Adj_HN_table <- mediation_table_func(med_Whr_Adj_HN_0.76, "10% increase in share of UPFs (%g/d)", "Waist-to-hip ratio", "Head and neck cancer") 

#adeno
med_int_Bmi_Adj_adeno_table <- mediation_table_func(med_Bmi_Adj_adeno_22.68, "10% increase in share of UPFs (%g/d)", "Body mass index", "Oesophageal adenocarcinoma") 
med_int_Whr_Adj_adeno_table <- mediation_table_func(med_Whr_Adj_adeno_0.76, "10% increase in share of UPFs (%g/d)", "Waist-to-hip ratio", "Oesophageal adenocarcinoma") 


#-----------------------------------------------------------------------------#
#                                 Save tables                                 #
#-----------------------------------------------------------------------------#
#export tables to word
sect_properties <- prop_section(page_size = page_size(orient = "landscape"))

#save %g/d tables
med_int_Whr_Adj_HN_table %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)
med_int_Bmi_Adj_adeno_table %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)
med_int_Whr_Adj_adeno_table %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)


############ export to excel file

mediation_func_excel <- function(med_object, exposure_var, mediator_var, outcome_var) {
  effect_df <- tibble(Effect = factor(names(med_object$effect.pe), levels = names(med_object$effect.pe)),
                      Estimate = med_object$effect.pe, SE = med_object$effect.se, CI_lower = med_object$effect.ci.low,
                      CI_upper = med_object$effect.ci.high, P = med_object$effect.pval)  
  effect_df %>% add_column(., "Exposure" = exposure_var, "Mediator" = mediator_var, "Outcome" = outcome_var, "Mval" = as.numeric(med_object$ref[[3]]), .before = "Effect") %>% 
    select(., Exposure, Mediator, Outcome, Mval, Effect, Estimate, SE, CI_lower, CI_upper, P) 
  
}


#hn
med_int_Whr_Adj_HN_excel <- mediation_func_excel(med_Whr_Adj_HN_0.76, "Ultra-processed foods", "Waist-to-hip ratio", "Head and neck cancer") 

#adeno
med_int_Bmi_Adj_adeno_excel <- mediation_func_excel(med_Bmi_Adj_adeno_22.68, "Ultra-processed foods", "Body mass index", "Oesophageal adenocarcinoma") 
med_int_Whr_Adj_adeno_excel <- mediation_func_excel(med_Whr_Adj_adeno_0.76, "Ultra-processed foods", "Waist-to-hip ratio", "Oesophageal adenocarcinoma") 

HN_O_Mal_excel <- rbind(med_int_Whr_Adj_HN_excel, med_int_Bmi_Adj_adeno_excel, 
                        med_int_Whr_Adj_adeno_excel)

write.xlsx(HN_O_Mal_excel, file = "FILE.xlsx")
