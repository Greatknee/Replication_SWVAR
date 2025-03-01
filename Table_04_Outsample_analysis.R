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
#Set the directory
setwd("C://Users//greatknee//Desktop//Mortality//Reproducibility//Replication_SWVAR")
#remove all existing variables
rm()
# Source functions
source("functions/utils.R")
source("functions/datapre.R")
source("functions/model_benchmark.R")
source("functions/model_SWVAR.R")
source("functions/func_table_04_outsample.R")
# # Load data
# data <- datapre_out(group = 5)
# # Run analysis
# result <- func_table_04_outsample(data)
# result
# Formal
start_time <- Sys.time()
group = c(1,2,3,5)
result_all = matrix(0,5*length(group),8)
rownames(result_all) = rep(c('RMSFE','sigma','max','min','MAFE'),length(group))
colnames(result_all) = c('Li-Lee','FDM','STAR1','STAR2','VAR','SVAR1','SVAR2','SWVAR')
for (i in 1:length(group)) {
  data <- datapre_out(group = group[i])
  result_all[(1:5)+(i-1)*5,] <- func_table_04_outsample(data)
}
result_all
# Save output
write.csv(result_all, "output/result_table_04.csv")
end_time <- Sys.time()
print(end_time - start_time)
#Time difference of 12.5 mins


