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

# Source functions
source("functions/utils.R")
source("functions/datapre.R")
source("functions/model_benchmark.R")
source("functions/model_SWVAR.R")
source("functions/func_figure_02_corr.R")




# Load data
data1 <- datapre_in(group = 1) 
data2 <- datapre_in(group = 2) 

# Run analysis
g1 = func_figure_02_corr(data1)
g2 = func_figure_02_corr(data2)

g1+g2

#save
pdf('output/result_figure_021.pdf', width=8, height=6)
print(g1)
dev.off()
pdf('output/result_figure_022.pdf', width=8, height=6)
print(g2)
dev.off()