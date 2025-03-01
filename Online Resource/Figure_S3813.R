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


# Source functions
source("functions/utils.R")
source("functions/datapre.R")
source("functions/model_benchmark.R")
source("functions/model_SWVAR.R")
source("functions/func_figure_03_prefit.R")
source("Online Resource/sfunctions/func_figure_s3.R")

time1 <- Sys.time()
# Load data
data <- datapre_out(group = 1)
# Run analysis
fit = func_figure_03_prefit(data)
# Figure 4
f1 = func_figure_s31(fit)
f1[[1]] + f1[[2]]
#save output
pdf(paste('Online Resource/output/figure_s',3,1,'.pdf',sep = ''), width=8, height=10)
print(f1[[1]])  
dev.off()
pdf(paste('Online Resource/output/figure_s',3,2,'.pdf',sep = ''), width=8, height=7.5)
print(f1[[2]])  
dev.off()

#forloop
ggroup = c(2,3)
for (i in 1:length(ggroup)) {
  data <- datapre_out(group = ggroup[i]) 
  # Run analysis
  fit = func_figure_03_prefit(data)
  out = func_figure_s32(fit)
  pdf(paste('Online Resource/output/figure_s',3+5*(ggroup[i]-1),1,'.pdf',sep = ''), width=8, height=7.5)
  print(out)  
  dev.off()
}

time2 <- Sys.time()
print(time2 - time1)

# Time difference of 8.99442 mins


