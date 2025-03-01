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
source("Online Resource/sfunctions/func_figure_s4.R")


time1 <- Sys.time()

# Load data
data <- datapre_out(group = 1)
# Figure 4
f1 = func_figure_s41(data)
f1[[1]] + f1[[2]]
#save output
pdf(paste('Online Resource/output/figure_s',4,1,'.pdf',sep = ''), width=8, height=10)
print(f1[[1]])  
dev.off()
pdf(paste('Online Resource/output/figure_s',4,2,'.pdf',sep = ''), width=8, height=7.5)
print(f1[[2]])  
dev.off()

#forloop
ggroup = c(2,3)
for (i in 1:length(ggroup)) {
  data <- datapre_out(group = ggroup[i]) 
  
  # Run analysis
  out = func_figure_s42(data)
  pdf(paste('Online Resource/output/figure_s',4+5*(ggroup[i]-1),1,'.pdf',sep = '') , width=8, height=7.5)
  print(out)  
  dev.off()
}

time2 <- Sys.time()
print(time2 - time1)

