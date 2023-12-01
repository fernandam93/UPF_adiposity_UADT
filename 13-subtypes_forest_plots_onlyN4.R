#######################################################################
#               FOOD PROCESSING HNC SUBTYPE FOREST PLOTS              #
#######################################################################
# R version 4.2.3 (2023-03-15)
# Last modified: 29 April 2023

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("DIRECTORY") 

if (!require("pacman")) install.packages("pacman")
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot", "openxlsx")

#---------------------------------------------------------------------#
#                            Read results                             #----
#---------------------------------------------------------------------#

results <- read.xlsx("RESULTS_FILE.xlsx")

#---------------------------------------------------------------------#
#                         Create forest plot                          #----
#---------------------------------------------------------------------#

results$model <- factor(results$model, levels = c("3", "2", "1"))

as.factor(results$outcome)
results$outcome <- factor(results$outcome, levels = c("Oral cavity cancer", "Oropharynx cancer", "Hypopharynx cancer", "Larynx cancer", "Unspecified/overlapping cancer"))

myforestplot <- function(df, title, xlab)
{
  x <- forestplot(
    df = df,
    estimate = b,
    se = se,
    pvalue = p.value,
    name = outcome,
    logodds = T,
    colour = model,
    title = title,
    xlab = xlab,
    #xlim= c(0.5,1.5)
  ) 
  colours_BP <- c("#E69F00", "#56B4E9", "#009E73")
  x <- x + scale_color_manual(values=colours_BP)
  x <- x + scale_x_continuous(breaks = c(0.60, 0.8, 1.00, 1.2, 1.5, 1.7, 2))
  print(x)
}

plot <- myforestplot(results, "", "Hazard ratio (95% CI) per 10% g/d higher UPF intake")

#---------------------------------------------------------------------#
#                               Save plot                             #----
#---------------------------------------------------------------------#

#save
save_func <- function(file_name, plot_name)
{
  png(file_name, res=300, height=2000, width=2000)
  print(plot_name)
  dev.off()
}
save_func('FILE.png', plot)

