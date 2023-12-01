#######################################################################
#               FOOD PROCESSING HNC AND OC FOREST PLOTS                #
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

results <- read.xlsx("RESULTS_FILE1.xlsx")

results_grams <- read.xlsx("RESULTS_FILE2.xlsx")

results_kcal_p <- read.xlsx("RESULTS_FILE3.xlsx")

results_kcal <- read.xlsx("RESULTS_FILE4.xlsx")

#---------------------------------------------------------------------#
#                         Create forest plots                         #----
#---------------------------------------------------------------------#

########################### grams% ##########################

results$model <- factor(results$model, levels = c("3", "2", "1"))

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
    xlim= c(0.8,1.5)
  ) 
  colours_BP <- c("#E69F00", "#56B4E9", "#009E73")
  x <- x + scale_color_manual(values=colours_BP)
  x <- x + scale_x_continuous(breaks = c(0.8, 1.0, 1.2, 1.4, 1.6))
  print(x)
}

plot <- myforestplot(results, "", "HR (95% CI) per 10% g/d higher UPF intake")


########################### grams ##########################

results_grams$model <- factor(results_grams$model, levels = c("3", "2", "1"))

myforestplot_grams <- function(df, title, xlab)
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
    xlim= c(0.95,1.13)
  ) 
  colours_BP <- c("#E69F00", "#56B4E9", "#009E73")
  x <- x + scale_color_manual(values=colours_BP)
  x <- x + scale_x_continuous(breaks = c(0.95, 1.00, 1.05, 1.1))
  print(x)
}

plot_grams <- myforestplot_grams(results_grams, "", "HR (95% CI) per 100 g higher UPF intake")


######################## kcal% ###########################

results_kcal_p$model <- factor(results_kcal_p$model, levels = c("3", "2", "1"))

myforestplot_kcal_p <- function(df, title, xlab)
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
    xlim= c(0.8,1.5)
  ) 
  colours_BP <- c("#E69F00", "#56B4E9", "#009E73")
  x <- x + scale_color_manual(values=colours_BP)
  x <- x + scale_x_continuous(breaks = c(0.8, 1.0, 1.2, 1.4, 1.6))
  print(x)
}

plot_kcal_p <- myforestplot_kcal_p(results_kcal_p, "", "HR (95% CI) per 10% kcal/d higher UPF intake")


########################### kcal ##########################

results_kcal$model <- factor(results_kcal$model, levels = c("3", "2", "1"))

myforestplot_kcal <- function(df, title, xlab)
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
    xlim= c(0.95,1.13)
  ) 
  colours_BP <- c("#E69F00", "#56B4E9", "#009E73")
  x <- x + scale_color_manual(values=colours_BP)
  x <- x + scale_x_continuous(breaks = c(0.95, 1.00, 1.05, 1.1))
  print(x)
}

plot_kcal <- myforestplot_kcal(results_kcal, "", "HR (95% CI) per 100 kcal higher UPF intake")


#---------------------------------------------------------------------#
#                           Combine plots                             #----
#---------------------------------------------------------------------#


all <- ggarrange(plot, plot_grams, plot_kcal_p, plot_kcal, ncol = 2, nrow=2, legend = "bottom", labels = c("A", "B", "C", "D"), common.legend = T)


#---------------------------------------------------------------------#
#                               Save plot                             #----
#---------------------------------------------------------------------#

#save
save_func2 <- function(file_name, plot_name)
{
  png(file_name, res=300, height=2000, width=4000)
  print(plot_name)
  dev.off()
}
save_func2('FILE.png', all)






