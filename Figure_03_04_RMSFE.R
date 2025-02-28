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
source("functions/func_figure_03_prefit.R")
source("functions/func_figure_03_draw.R")
source("functions/func_figure_04_draw.R")

# Load data
data <- datapre_out(group = 1) 

# Run analysis
fit = func_figure_03_prefit(data)

# Figure 3
f1 = func_figure_03_draw(fit)
# Figure 4
f2 = func_figure_04_draw(fit)

f1
f2
