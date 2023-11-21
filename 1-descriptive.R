###############################################################################
#                UPF - HEAD & NECK AND OESOPHAGEAL CANCERS PROJECT            #
###############################################################################
#last modified: 5 October 2023

#-----------------------------------------------------------------------------#
#                                  Housekeeping                               #
#-----------------------------------------------------------------------------#

rm(list=ls()) #Remove any existing objects in R 

setwd("WORKING_DIR") 

if (!require("pacman")) install.packages("pacman")
pacman::p_load("haven", "tidyverse", "dplyr", "data.table", "expss", "table1", "dplyr", "flextable", "magrittr", "officer", "finalfit")

#-----------------------------------------------------------------------------#
#                                  Read dataset                               #
#-----------------------------------------------------------------------------#

data <- readRDS(file = "DATA.RDS")

#-----------------------------------------------------------------------------#
#                            Descriptive tables                               #
#-----------------------------------------------------------------------------#
#https://cran.r-project.org/web/packages/table1/vignettes/table1-examples.html

#factor categorical variables
data$QE_P_M_N4_quartile <- 
  factor(data$QE_P_M_N4_quartile, levels=c(1,2,3,4),
         labels=c("Q1", 
                  "Q2",
                  "Q3",
                  "Q4"))

data$QE_P_M_N3_quartile <- 
  factor(data$QE_P_M_N3_quartile, levels=c(1,2,3,4),
         labels=c("Q1", 
                  "Q2",
                  "Q3",
                  "Q4"))

data$QE_P_M_N2_quartile <- 
  factor(data$QE_P_M_N2_quartile, levels=c(1,2,3,4),
         labels=c("Q1", 
                  "Q2",
                  "Q3",
                  "Q4"))

data$QE_P_M_N1_quartile <- 
  factor(data$QE_P_M_N1_quartile, levels=c(1,2,3,4),
         labels=c("Q1", 
                  "Q2",
                  "Q3",
                  "Q4"))

data$QE_P_M_N4_quartile_SS <- 
  factor(data$QE_P_M_N4_quartile_SS, levels=c(1,2,3,4),
         labels=c("Q1", 
                  "Q2",
                  "Q3",
                  "Q4"))

data$QE_P_M_N3_quartile_SS <- 
  factor(data$QE_P_M_N3_quartile_SS, levels=c(1,2,3,4),
         labels=c("Q1", 
                  "Q2",
                  "Q3",
                  "Q4"))

data$QE_P_M_N2_quartile_SS <- 
  factor(data$QE_P_M_N2_quartile_SS, levels=c(1,2,3,4),
         labels=c("Q1", 
                  "Q2",
                  "Q3",
                  "Q4"))

data$QE_P_M_N1_quartile_SS <- 
  factor(data$QE_P_M_N1_quartile_SS, levels=c(1,2,3,4),
         labels=c("Q1", 
                  "Q2",
                  "Q3",
                  "Q4"))

data$Sex <- 
  factor(data$Sex, levels=c(1,2),
         labels=c("Male", 
                  "Female"))

data$Smoke_Stat <- 
  factor(data$Smoke_Stat, levels=c(1,2,3),
         labels=c("Never", 
                  "Former",
                  "Smoker"))

data$Smoke_Stat_C <- 
  factor(data$Smoke_Stat_C, levels=c(1,2,3),
         labels=c("Never", 
                  "Former",
                  "Smoker"))

data$L_School <- 
  factor(data$L_School, levels=c(0,1,2,3,4),
         labels=c("None", 
                  "Primary school completed",
                  "Technical/Professional school",
                  "Secondary school",
                  "Longer education"))
data$L_School_C <- 
  factor(data$L_School_C, levels=c(0,1,2,3,4),
         labels=c("None", 
                  "Primary school completed",
                  "Technical/Professional school",
                  "Secondary school",
                  "Longer education"))

data$Pa_Index <- 
  factor(data$Pa_Index, levels=c(1,2,3,4),
         labels=c("Inactive", 
                  "Moderately inactive",
                  "Moderately active",
                  "Active"))

data$Pa_Index_C <- 
  factor(data$Pa_Index_C, levels=c(1,2,3,4),
         labels=c("Inactive", 
                  "Moderately inactive",
                  "Moderately active",
                  "Active"))

data$MRMed_Score_C <- 
  factor(data$MRMed_Score_C, levels=c(1,2,3),
         labels=c("Low", 
                  "Medium",
                  "High"))

data$Country <- 
  factor(data$Country, levels=c(1,2,3,4,5,7,8,9,"B"),
         labels=c("France", 
                  "Italy",
                  "Spain",
                  "United Kingdom",
                  "The Netherlands",
                  
                  "Germany",
                  "Sweden",
                  "Denmark",
                  "Norway"))


units(data$Age_Recr) <- "years"
units(data$Height_C) <- "cm"
units(data$Bmi_Adj) <- "kg/m2"
units(data$Alc_Re) <- "g/d"
units(data$QE_US_CHOCDF) <- "g/d"
units(data$QE_US_NA) <- "mg/d"
units(data$QE_US_FAT) <- "g/d"
units(data$QE_US_FIBTG) <- "g/d"
units(data$QE_US_PROCNT) <- "g/d"
units(data$QE_US_ENERGY) <- "kcal/d"
units(data$QE_P_M_N4) <- "%g/d"
units(data$QE_P_M_N3) <- "%g/d"
units(data$QE_P_M_N2) <- "%g/d"
units(data$QE_P_M_N1) <- "%g/d"
units(data$QE_M_N4) <- "g/d"
units(data$QE_M_N3) <- "g/d"
units(data$QE_M_N2) <- "g/d"
units(data$QE_M_N1) <- "g/d"



table1::label(data$Smoke_Stat) <- "Smoking status"
table1::label(data$Smoke_Stat_C) <- "Smoking status"
table1::label(data$Age_Recr) <- "Age at recruitment"
table1::label(data$Sex) <- "Sex"
table1::label(data$Height_C) <- "Height"
table1::label(data$Bmi_Adj) <- "Measured BMI"
table1::label(data$L_School) <- "Educational level"
table1::label(data$Pa_Index) <- "Physical activity level"
table1::label(data$L_School_C) <- "Educational level"
table1::label(data$Pa_Index_C) <- "Physical activity level"
table1::label(data$Whr_Adj) <- "Measured WHR"
table1::label(data$Country) <- "Country"
table1::label(data$Alc_Re) <- "Alcohol intake"
table1::label(data$QE_US_ENERGY) <- "Energy intake"
table1::label(data$QE_US_CHOCDF) <- "Carbohydrate intake"
table1::label(data$QE_US_NA) <- "Sodium intake"
table1::label(data$QE_US_FAT) <- "Fat intake"
table1::label(data$QE_US_PROCNT) <- "Protein intake"
table1::label(data$QE_US_FIBTG) <- "Fibre intake"
table1::label(data$MRMed_Score_C) <- "Modified relative Mediterranean diet score"
table1::label(data$QE_P_M_N4_quartile) <- "Quartiles of share of ultra-processed foods in the diet"
table1::label(data$QE_P_M_N3_quartile) <- "Quartiles of share of processed foods in the diet"
table1::label(data$QE_P_M_N2_quartile) <- "Quartiles of share of processed ingredients in the diet"
table1::label(data$QE_P_M_N1_quartile) <- "Quartiles of share of unprocessed and minimally processed foods in the diet"
table1::label(data$QE_P_M_N4) <- "Relative intake of ultra-processed foods in the diet"
table1::label(data$QE_P_M_N3) <- "Relative intake of processed foods in the diet"
table1::label(data$QE_P_M_N2) <- "Relative intake of processed ingredients in the diet"
table1::label(data$QE_P_M_N1) <- "Relative intake of unprocessed and minimally processed foods in the diet"
table1::label(data$QE_M_N4) <- "Absolute intake of ultra-processed foods in the diet"
table1::label(data$QE_M_N3) <- "Absolute intake of processed foods in the diet"
table1::label(data$QE_M_N2) <- "Absolute intake of processed ingredients in the diet"
table1::label(data$QE_M_N1) <- "Absolute intake of unprocessed and minimally processed foods in the diet"

#summary exposures
sum_N4_SS <- tapply(data$QE_P_M_N4, data$Sex, 
                    function(x) format(summary(x), scientific = F))
sum_N4 <- summary(data$QE_P_M_N4)

sum_N3_SS <- tapply(data$QE_P_M_N3, data$Sex, 
                    function(x) format(summary(x), scientific = F))
sum_N3 <- summary(data$QE_P_M_N3)

sum_N2_SS <- tapply(data$QE_P_M_N2, data$Sex, 
                    function(x) format(summary(x), scientific = F))
sum_N2 <- summary(data$QE_P_M_N2)

sum_N1_SS <- tapply(data$QE_P_M_N1, data$Sex, 
                    function(x) format(summary(x), scientific = F))
sum_N1 <- summary(data$QE_P_M_N1)


# sum_N4[[2]] #Q1
# sum_N4[[3]] #Q2
# sum_N4[[5]] #Q3
# sum_N4[[6]] #Q4
# 
# sum_N4_SS$Male[[2]] #Q1
# sum_N4_SS$Male[[3]] #Q2
# sum_N4_SS$Male[[5]] #Q3
# sum_N4_SS$Male[[6]] #Q4


#create tables
#sex specific quartiles
x4_SS <- table1::table1(~Age_Recr + Sex + Country + Height_C + Bmi_Adj + Whr_Adj + 
                       L_School + Pa_Index + Smoke_Stat + Alc_Re + 
                       QE_US_ENERGY + QE_US_CHOCDF + QE_US_FAT + 
                       QE_US_FIBTG + QE_US_NA + QE_P_M_N4 + + QE_M_N4 + QE_P_M_N3 + QE_M_N3 + QE_P_M_N1 + QE_M_N1 | QE_P_M_N4_quartile_SS, 
                     data = data, footnote = paste("Baseline characteristics have been stratified by sex-specific quartiles of relative intake of ultra-processed foods (%g/d).
Abbreviations: BMI, Body mass index; Q1, quartile 1; Q2, quartile 2; Q3, quartile 3; Q4, quartile 4.
Quartile cut-offs were defined as follows: Male Q1 = ", format(round(as.numeric(sum_N4_SS$Male[[2]]), 2), nsmall = 2), " %g/d, Q2 = ", format(round(as.numeric(sum_N4_SS$Male[[3]]), 2), nsmall = 2), " %g/d, Q3 = ", format(round(as.numeric(sum_N4_SS$Male[[5]]), 2), nsmall = 2), " %g/d, Q4 = ", format(round(as.numeric(sum_N4_SS$Male[[6]]), 2), nsmall = 2), " %g/d; Female Q1 = ", format(round(as.numeric(sum_N4_SS$Female[[2]]), 2), nsmall = 2), " %g/d, Q2 = ", format(round(as.numeric(sum_N4_SS$Female[[3]]), 2), nsmall = 2), " %g/d, Q3 = ", format(round(as.numeric(sum_N4_SS$Female[[5]]), 2), nsmall = 2), " %g/d, Q4 = ", format(round(as.numeric(sum_N4_SS$Female[[6]]), 2), nsmall = 2), " %g/d.", sep = ""))


#-----------------------------------------------------------------------------#
#                                 Save tables                                 #
#-----------------------------------------------------------------------------#

sect_properties <- prop_section(page_size = page_size(orient = "landscape"))

#sex specific quartiles
t1flex(x4_SS) %>% 
  save_as_docx(path="FILE.docx", pr_section = sect_properties)

################



#-----------------------------------------------------------------------------#
#                           Histogram UPF intake                              #
#-----------------------------------------------------------------------------#

#histogram
hist(data$QE_P_M_N4, xlab = "Proportion of ultra-processed foods in the diet (%g/d)", breaks = 50, col = "skyblue2", border = "black", xlim = c(0,60), main = title(main = ""))

#save
  png("FILE.png", res=300, height=1500, width=3000)
  hist(data$QE_P_M_N4, xlab = "Proportion of ultra-processed foods in the diet (%g/d)", breaks = 50, col = "skyblue2", border = "black", xlim = c(0,60), main = title(main = ""))
  dev.off()

#calculate summary statistics for each NOVA group
  summary(data$QE_P_M_N4)
  summary(data$QE_P_M_N3)
  summary(data$QE_P_M_N1)
  