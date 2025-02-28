# Load necessary libraries
library(abind)
library(dplyr)
library(ggplot2)
library(psych)
library(sparsevar)
library(vars)
library(MASS) 
library(forecast)
library(vital)
library(demography)
library(tsibble)
library(dplyr)
library(CVXR)
library(patchwork)

#Set the directory
setwd("C://Users//greatknee//Desktop//Mortality//Reproducibility//Replication_SWVAR")

# Source functions
source("functions/utils.R")
source("functions/datapre.R")
source("functions/model_benchmark.R")
source("functions/model_SWVAR.R")
source("functions/func_figure_06_le_1.R")
source("functions/func_figure_06_le_2.R")


# Load data
data <- datapre_out(group = 1) 
# Figure 6.a
f1 = func_figure_07_multijump()

f1
