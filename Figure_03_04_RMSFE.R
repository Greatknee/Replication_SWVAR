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
library(feasts)

# Source functions
source("functions/utils.R")
source("functions/datapre.R")
source("functions/model_benchmark.R")
source("functions/model_SWVAR.R")
source("functions/func_figure_03_prefit.R")
source("functions/func_figure_03_draw.R")
source("functions/func_figure_04_draw.R")

# # Sample
# # Load data
# data <- datapre_out(group = 1)
# # Run analysis
# fit = func_figure_03_prefit(data)
# # Figure 3
# f1 = func_figure_03_draw(fit)
# # Figure 4
# f2 = func_figure_04_draw(fit)
# 
# f1
# f2


#forloop
time1 <- Sys.time()
group = c(1,2,3,5)
for (i in 1:length(group)) {
  data <- datapre_out(group = group[i])
  fit = func_figure_03_prefit(data)
  f1 = func_figure_03_draw(fit)
  if (group[i] == 1) {
    f2 = func_figure_04_draw(fit)
    pdf('output/result_figure_041.pdf', width=15, height=6)
    print(f2)
    dev.off()
  }
  if (group[i] == 3) {
    f2 = func_figure_04_draw(fit)
    pdf('output/result_figure_042.pdf', width=15, height=6)
    print(f2)
    dev.off()
  }
  
  # save output
  pdf(paste('output/result_figure_03',i,'.pdf',sep = ''), width=8, height=6)
  print(f1)
  dev.off()
}

time2 <- Sys.time()
print(time2 - time1)

# Time difference of 10.66827 mins