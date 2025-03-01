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

#Set the directory
setwd("C://Users//greatknee//Desktop//Mortality//Reproducibility//Replication_SWVAR")

# Source functions
source("functions/utils.R")
source("functions/datapre.R")
source("functions/model_benchmark.R")
source("functions/model_SWVAR.R")
source("functions/func_figure_01_coef.R")

# Load data
data <- datapre_in(group = 1) 

# Run analysis
out = func_figure_01_coef(data)
out$g1+out$g2+out$g3+out$g4
