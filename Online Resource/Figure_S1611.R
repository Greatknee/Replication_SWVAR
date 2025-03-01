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


# Source functions
source("functions/utils.R")
source("functions/datapre.R")
source("functions/model_benchmark.R")
source("functions/model_SWVAR.R")
source("Online Resource/sfunctions/func_figure_s1.R")


# Load data
ggroup = 1
data <- datapre_in(group = ggroup) 

# Run analysis
out = func_figure_s11(data)
print(out$g1+out$g3+out$g5+out$g7)#SVAR  
print(out$g2+out$g4+out$g6+out$g8)#SWVAR


#save output
for (i in 1:length(out)) {
  pdf(paste('Online Resource/output/figure_s',1+5*(ggroup-1),i,'.pdf',sep = ''), width=7, height=6)
  print(out[[i]])  
  dev.off()
}

#######################################
# Load data
ggroup = c(2,3)
for (i in 1:length(ggroup)) {
  data <- datapre_in(group = ggroup[i]) 
  
  # Run analysis
  out = func_figure_s12(data)
  print(out$g1+out$g2)#SVAR+SWVAR

  
  #save output
  for (ii in 1:length(out)) {
    pdf(paste('Online Resource/output/figure_s',1+5*(ggroup[i]-1),ii,'.pdf',sep = ''), width=7, height=6)
    print(out[[ii]])  
    dev.off()
  }
}
