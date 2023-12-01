#######################################################################
#               FOOD PROCESSING HNC SUBTYPE HETEROGENEITY             #
#######################################################################
# R version 4.2.3 (2023-03-15)
# Last modified: 13 March 2023

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("DIRECTORY") 

if (!require("pacman")) install.packages("pacman")
pacman::p_load("meta", "metafor", "tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot", "openxlsx", "cowplot")


#---------------------------------------------------------------------#
#                            Read results                             #----
#---------------------------------------------------------------------#

results_N1 <- read.xlsx("RESULTS_FILE1.xlsx")
results_N3 <- read.xlsx("RESULTS_FILE2.xlsx")
results_N4 <- read.xlsx("RESULTS_FILE3.xlsx")

results <- rbind(results_N1, results_N3, results_N4
                 )

#---------------------------------------------------------------------#
#                      Prepare for analyses                           #----
#---------------------------------------------------------------------#

results$food_labels[results$exposure=="Unprocessed/minimally processed foods"] <- "NOVA 1"
results$food_labels[results$exposure=="Processed foods"] <- "NOVA 3"
results$food_labels[results$exposure=="Ultra-processed foods"] <- "NOVA 4"

results <- results[results$model==3,]

#---------------------------------------------------------------------#
#                     Function for meta-analysis                       #----
#---------------------------------------------------------------------#

meta_func <- function(exp_varname, out1, out2="", out3="", out4="", 
                      out5="", out6="", out7="", out8="", out9="", out10="", 
                      out11="", out12="", out13="", out14="", out15="", 
                      out16="", out17="", out18="", out19="", out20="", out21="")
{
  input <- results
  input <- input[which(input$exposure==exp_varname),]
  input <- input[which(input$outcome==out1 | input$outcome==out2 | input$outcome==out3 | 
                         input$outcome==out4 | input$outcome==out5 | input$outcome==out6 | 
                         input$outcome==out7 | input$outcome==out8 | input$outcome==out9 | 
                         input$outcome==out10 | input$outcome==out11 | input$outcome==out12 |
                         input$outcome==out13 | input$outcome==out14 | input$outcome==out15 | 
                         input$outcome==out16 | input$outcome==out17 | input$outcome==out18 | 
                         input$outcome==out19 | input$outcome==out20 | input$outcome==out21),]
  #meta-analysis
  a <- metagen(TE = b, seTE = se, data = input, 
               studlab = paste(outcome), sm = "HR",
               hakn = FALSE, byvar = c(exposure),
               method.tau="DL", comb.fixed = T, comb.random = F) 
  print(a)

  #extract values from meta output
  TE.tibble <- as_tibble(a$TE.fixed.w)
  se.tibble <- as_tibble(a$seTE.fixed.w)
  p.tibble <- as_tibble(a$pval.fixed.w)
  bylevs.tibble <- as_tibble(a$bylevs)
  #combine tibbles and change column names
  tibble <- cbind(TE.tibble, se.tibble, p.tibble, bylevs.tibble)
  colnames(tibble) <- c("b", "se", "pval", "exposure")
  #add columns for outcomes and het test
  tibble$outcomes <- paste(out1, out2, out3, out4, out5, sep = ", ")
  tibble$het_test <- a$pval.Q
  tibble
  
}

#---------------------------------------------------------------------#
#                Run function to calculate heterogeneity              #----
#---------------------------------------------------------------------#

het_test_N1 <- meta_func("Unprocessed/minimally processed foods", "Oral cavity cancer", "Hypopharynx cancer" , "Oropharynx cancer", "Larynx cancer")
het_test_N3 <- meta_func("Processed foods", "Oral cavity cancer", "Hypopharynx cancer" , "Oropharynx cancer", "Larynx cancer")
het_test_N4 <- meta_func("Ultra-processed foods", "Oral cavity cancer", "Hypopharynx cancer" , "Oropharynx cancer", "Larynx cancer")

het_test_N1_unspecified <- meta_func("Unprocessed/minimally processed foods", "Oral cavity cancer", "Hypopharynx cancer" , "Oropharynx cancer", "Larynx cancer", "Unspecified/overlapping cancer")
het_test_N3_unspecified <- meta_func("Processed foods", "Oral cavity cancer", "Hypopharynx cancer" , "Oropharynx cancer", "Larynx cancer", "Unspecified/overlapping cancer")
het_test_N4_unspecified <- meta_func("Ultra-processed foods", "Oral cavity cancer", "Hypopharynx cancer" , "Oropharynx cancer", "Larynx cancer", "Unspecified/overlapping cancer")

# 
# het_test_N1_pharynx<- meta_func("Unprocessed/minimally processed foods", "Hypopharynx cancer", "Oropharynx cancer")
# het_test_N3_pharynx<- meta_func("Processed foods", "Hypopharynx cancer", "Oropharynx cancer")
# het_test_N4_pharynx<- meta_func("Ultra-processed foods", "Hypopharynx cancer", "Oropharynx cancer")

all_combined <- rbind(het_test_N1, het_test_N3, het_test_N4, 
                      het_test_N1_unspecified, het_test_N3_unspecified, het_test_N4_unspecified 
                      #het_test_N1_pharynx, het_test_N3_pharynx, het_test_N4_pharynx
                      )

#---------------------------------------------------------------------#
#                             Save to excel                           #----
#---------------------------------------------------------------------#

write.xlsx(all_combined, file = "FILE.xlsx")



