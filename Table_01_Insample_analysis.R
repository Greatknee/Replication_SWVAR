# Load necessary libraries
library(dplyr)
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

# Source functions
source("functions/utils.R")
source("functions/datapre.R")
source("functions/model_benchmark.R")
source("functions/model_SWVAR.R")
source("functions/func_table_01_insample.R")

# # Sample
# # Load data
# data <- datapre_in(group = 5)
# 
# # Run analysis
# result <- func_table_01_insample_m(data)
# #write.csv(result, "output/result_table_01_temp.csv")
# result

# Formal
start_time <- Sys.time()
result_all = matrix(0,5*5,8)
rownames(result_all) = rep(c('RMSE','edf','logL','AIC','BIC'),5)
colnames(result_all) = c('Li-Lee','Co-FDM','STAR1','STAR2','VAR','SVAR1','SVAR2','SWVAR')
for (i in 1:5) {
  data <- datapre_in(group = i)
  result_all[(1:5)+(i-1)*5,] <- func_table_01_insample_m(data)
}
result_all
# Save output
write.csv(result_all, "output/result_table_01_compare.csv")
end_time <- Sys.time()
print(end_time - start_time)
#Time difference of 10.87298 mins


# Remarks: The warning message (NAs introduced by coercion to integer range) originates from the CVXR function while solving the objective function of the STAR model. 
# This indicates that the solution provided by CVXR is suboptimal. Such issues commonly arise when handling large-scale optimization tasks and cannot be resolved simply by increasing the number of iterations.
# See also: https://yetanothermathprogrammingconsultant.blogspot.com/2022/08/large-sparse-transportation-model-with.html



