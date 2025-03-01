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
library(tseries)
library(feasts)

# Source functions
source("functions/utils.R")
source("functions/datapre.R")
source("functions/model_benchmark.R")
source("functions/model_SWVAR.R")
source("functions/func_figure_09_le_special.R")


start_time <- Sys.time()
# Figure 9
out = func_figure_09_le_special()

pdf('output/result_figure_091.pdf', width=15, height=6)
print(out$p1 + out$p2)
dev.off()

end_time <- Sys.time()
print(end_time - start_time)

#