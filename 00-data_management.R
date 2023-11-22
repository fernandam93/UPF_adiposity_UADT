###############################################################################
#                UPF - HEAD & NECK AND OESOPHAGEAL CANCERS PROJECT            #
###############################################################################
#last modified: 12 Mar 2023

#-----------------------------------------------------------------------------#
#                                  Housekeeping                               #
#-----------------------------------------------------------------------------#

rm(list=ls()) #Remove any existing objects in R 

setwd("WORKING_DIR") 

if (!require("pacman")) install.packages("pacman")
pacman::p_load("haven", "tidyverse", "dplyr", "data.table", "expss", "table1", "dplyr", "flextable", "magrittr", "officer")

#-----------------------------------------------------------------------------#
#                                  Read dataset                               #
#-----------------------------------------------------------------------------#

data <- readRDS(file = "ORIGINAL_FILE_PREMANAGEMENT.rds")

#-----------------------------------------------------------------------------#
#                         Recode missing values                               #
#-----------------------------------------------------------------------------#
#deal with missing values

#educational level
table(data$L_School, useNA = "always")
data$L_School[data$L_School==5] <- NA

#smoking status
table(data$Smoke_Stat, useNA = "always")
data$Smoke_Stat[data$Smoke_Stat==4] <- NA

#physical activity index
table(data$Pa_Index, useNA = "always")
data$Pa_Index[data$Pa_Index==5] <- NA

#-----------------------------------------------------------------------------#
#                         Imputation for missing values                       #
#-----------------------------------------------------------------------------#
#education has 3.7% missing, PA 2% missing, smoking 1.9% missing

#We can impute education, physical activity and smoking with modal value 
table(data$L_School, useNA = "always")
table(data$Pa_Index, useNA = "always")
table(data$Smoke_Stat, useNA = "always")

# Create the function.
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#find modes
getmode(data$L_School) #1-Primary school completed
getmode(data$Pa_Index) #2-Moderately inactive
getmode(data$Smoke_Stat) #1-Never

#create new variables and impute missing with mode
data$L_School_C <- data$L_School
data$L_School_C[is.na(data$L_School_C)] <- getmode(data$L_School) 

data$Pa_Index_C <- data$Pa_Index
data$Pa_Index_C[is.na(data$Pa_Index_C)] <- getmode(data$Pa_Index) 

data$Smoke_Stat_C <- data$Smoke_Stat
data$Smoke_Stat_C[is.na(data$Smoke_Stat_C)] <- getmode(data$Smoke_Stat) 

#-----------------------------------------------------------------------------#
#                 Create new exposure variables                               #
#-----------------------------------------------------------------------------#

#create UPF quartiles
data$QE_P_M_N4_quartile <- ntile(data$QE_P_M_N4, 4) 
#create quartiles for other NOVA categories
data$QE_P_M_N3_quartile <- ntile(data$QE_P_M_N3, 4)
data$QE_P_M_N2_quartile <- ntile(data$QE_P_M_N2, 4)
data$QE_P_M_N1_quartile <- ntile(data$QE_P_M_N1, 4)

#create sex specific quartiles for UPFs and other NOVA groups
data <- data %>%
  group_by(Sex) %>%
  mutate(QE_P_M_N4_quartile_SS = ntile(QE_P_M_N4, 4)) %>%
  mutate(QE_P_M_N3_quartile_SS = ntile(QE_P_M_N3, 4)) %>%
  mutate(QE_P_M_N2_quartile_SS = ntile(QE_P_M_N2, 4)) %>%
  mutate(QE_P_M_N1_quartile_SS = ntile(QE_P_M_N1, 4))

#-----------------------------------------------------------------------------#
#                  Create new Oesophageal cancer variables                    #
#-----------------------------------------------------------------------------#

#CREATE SCC AND ADENO CODES BASED ON CLASSIFICATIONS
adeno_codes <- c("8140/3", "8144/3", "8480/3", "8481/3", "8490/3") 

#CREATE OESOPHAGEAL CANCER CODES
oesophagus_codes <- c("C150", "C151", "C152", "C153", "C154", "C155", "C158", "C159")


#USE CODES TO CREATE OESOPHAGEAL CANCER VARIABLES
#overall oesophageal cancer
data$Esoph_Mal_SiteCode <- ifelse(data$Sitestom %in% oesophagus_codes, 1, 0)
#scc and adenocarcinoma of oesophagus
data$Esoph_Mal_Adeno_SiteCode <- ifelse((data$Sitestom %in% oesophagus_codes & data$Morpstom %in% adeno_codes), 1, 0) #MATCHES VARAIBLE ALREADY AVAILABLE IN EPIC

table(data$Esoph_Mal_Adeno, data$Esoph_Mal_Adeno_SiteCode2) #we can use the original EPIC variable 'Esoph_Mal_Adeno'

#-----------------------------------------------------------------------------#
#                        Create new HNC variables                             #
#-----------------------------------------------------------------------------#

#CREATE HNC CODES
#larynx
larynx_codes <- c("C320", "C321", "C322", "C323", "C328", "C329")
#oral cavity
lip_codes <- c("C000", "C001", "C002", "C003", "C004", "C005", "C006", "C008", "C009")
tongue_codes <- c("C019", "C020", "C021", "C022", "C023", "C024", "C028", "C029")
gum_codes <- c("C030", "C031", "C039")
floorofmouth_codes <- c("C040", "C041", "C048", "C049")
palate_codes <- c("C050", "C051", "C052", "C058", "C059")
othermouth_codes <- c("C060", "C061", "C062", "C068", "C069")
salivarygland_codes <- c("C079", "C080", "C081", "C088", "C089")
#pharynx
tonsil_codes <- c("C090", "C091", "C098", "C099") #part of oropharynx
oropharynx_codes <- c("C100", "C101", "C102", "C103", "C104", "C108", "C109") 
nasopharynx_codes <- c("C110", "C111", "C112", "C113", "C118", "C119")
pyriformsinus_code <- c("C129") #part of hypopharynx
hypopharynx_codes <- c("C130", "C131", "C132", "C138", "C139")
#other
otherHN_codes <- c("C140", "C142", "C148")


#INHANCE coding for HNC
INHANCE_2007_oral_codes <- c("C003", "C004", "C005", "C006", "C008", "C009", "C020", "C021", "C022", "C023",
                        gum_codes, floorofmouth_codes, "C050", othermouth_codes)
INHANCE_2007_oropharynx_codes <- c("C019", "C024", "C051", "C052", tonsil_codes, oropharynx_codes)
INHANCE_2007_hypopharynx_codes <- c(pyriformsinus_code, hypopharynx_codes)
INHANCE_2007_larynx_codes <- larynx_codes
INHANCE_2007_unspecified_codes <- c("C028", "C029", "C058", "C059", otherHN_codes)

INHANCE_2007_codes <- c(INHANCE_2007_oral_codes, INHANCE_2007_oropharynx_codes, INHANCE_2007_hypopharynx_codes,
                        INHANCE_2007_unspecified_codes, INHANCE_2007_larynx_codes) 

#USE CODES TO CREATE HNC VARIABLES
#overall HN cancer
data$HN_Mal_INHANCE_2007_SiteCode <- ifelse(data$Siteuadt %in% INHANCE_2007_codes, 1, 0) 
#INHANCE subtypes
data$Oral_cavity_Mal_INHANCE_2007_SiteCode <- ifelse(data$Siteuadt %in% INHANCE_2007_oral_codes, 1, 0)
data$Oropharynx_Mal_INHANCE_2007_SiteCode <- ifelse(data$Siteuadt %in% INHANCE_2007_oropharynx_codes, 1, 0)
data$Hypopharynx_Mal_INHANCE_2007_SiteCode <- ifelse(data$Siteuadt %in% INHANCE_2007_hypopharynx_codes, 1, 0)
data$Larynx_Mal_INHANCE_2007_SiteCode <- ifelse(data$Siteuadt %in% INHANCE_2007_larynx_codes, 1, 0)
data$Unspecified_Mal_INHANCE_2007_SiteCode <- ifelse(data$Siteuadt %in% INHANCE_2007_unspecified_codes, 1, 0)


#-----------------------------------------------------------------------------#
#                        Create other variables                               #
#-----------------------------------------------------------------------------#

#create age at recruitment in 1-year categories variable 
data$Age_Recr_1y <- trunc(data$Age_Recr)

#-----------------------------------------------------------------------------#
#                             Save new dataset                                #
#-----------------------------------------------------------------------------#

#save edited dataset
saveRDS(data, file = "FILE.rds") 


