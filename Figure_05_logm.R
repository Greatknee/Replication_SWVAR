# Load necessary libraries
library(abind)
library(dplyr)
library(ggplot2)
library(psych)
library(sparsevar)
library(vars)
library(MASS) 
library(forecast)
library(patchwork)
library(glasso)
#Set the directory
setwd("C://Users//greatknee//Desktop//Mortality//Reproducibility//Replication_SWVAR")

# Source functions
source("functions/utils.R")
source("functions/datapre.R")
source("functions/model_benchmark.R")
source("functions/model_SWVAR.R")
source("functions/func_figure_05_logm.R")




# Load data
data1 <- datapre_out(group = 1) 
data2 <- datapre_out(group = 3) 

# Run analysis
g1 = func_figure_05_logm(data1)
g2 = func_figure_05_logm(data2)
g1 + g2
