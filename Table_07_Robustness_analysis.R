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
library(tseries)
library(feasts)


# Source functions
source("functions/utils.R")
source("functions/datapre.R")
source("functions/model_benchmark.R")
source("functions/model_SWVAR.R")
source("functions/func_table_07_robust_g1_5.R")
source("functions/func_table_07_robust_g2_1.R")

# Sample
# Load data
#group = '1'/'2','3','5',setting = '1'(5 year age group)/'2'(1 year age group,middle age), gender = 'male'/'female'
# dataf <- datapre_robust(group = 1, setting = '1', gender = 'female')
# 
# # Run analysis
# resultf <- func_table_07_robust_g1_5(dataf)
# 
# datam <- datapre_robust(group = 1, setting = '1', gender = 'male')
# resultm <- func_table_07_robust_g1_5(datam)
# 
# resultf
# resultm

######################################
#The first half part for Setting 1 (5-year-age group)
#You can change the parameter in main function "func_table_07_robust_g1_5(~,fitstar)" to decide whether carry out it, 
#the default in the reproduction package is "fitstar:FALSE" means not run it.

#Forloop
start_time1 <- Sys.time()
result_set1 = matrix(0,4*2,2*5)
gen = c('female','male')
group = c(1,2,3,5)
rownames(result_set1) = rep(c('RMSFE','MAFE'),length(group))
colnames(result_set1) = rep(c('Li-Lee','CoFDM','STAR','SVAR','SWVAR'),2)
for (i in 1:length(gen)) {
  gi = gen[i]
  for (j in 1:length(group)) {
    gj = group[j]
    data <- datapre_robust(group = gj, setting = '1', gender = gi) 
    result_set1[((1:2)+(j-1)*2),((1:5)+(i-1)*5)] <- func_table_07_robust_g1_5(data)
  }
}
result_set1
# Run analysis

# Save output
write.csv(result_set1, "output/result_table_07_set1.csv")
end_time1 <- Sys.time()
print(end_time1 - start_time1)

#Time difference of 15.19703 mins

######################################
######################################
#The second half part for Setting 2 (1-year-age group).
#Now there are 45 age groups, fitting the STAR model will cost more then 1 hours.
#You can change the parameter in main function "func_table_07_robust_g2_1(~,fitstar= TRUE)" to decide whether carry out it, 
#the default in the reproduction package is not run it.


#Forloop
start_time2 <- Sys.time()
result_set2 = matrix(0,4*2,2*5)
gen = c('female','male')
group = c(1,2,3,5)
rownames(result_set2) = rep(c('RMSFE','MAFE'),length(group))
colnames(result_set2) = rep(c('Li-Lee','CoFDM','STAR','SVAR','SWVAR'),2)
for (i in 1:length(gen)) {
  gi = gen[i]
  for (j in 1:length(group)) {
    gj = group[j]
    data <- datapre_robust(group = gj, setting = '2', gender = gi) 
    result_set2[((1:2)+(j-1)*2),((1:5)+(i-1)*5)] <- func_table_07_robust_g2_1(data,star = FALSE)
  }
}
result_set2
# Save output
write.csv(result_set2, "output/result_table_07_set2.csv")
end_time2 <- Sys.time()
print(end_time2 - start_time2)

# Remarks: The warning message (NAs introduced by coercion to integer range) originates from the CVXR function while solving the objective function of the STAR model. 
# This indicates that the solution provided by CVXR is suboptimal. Such issues commonly arise when handling large-scale optimization tasks and cannot be resolved simply by increasing the number of iterations.
# See also: https://yetanothermathprogrammingconsultant.blogspot.com/2022/08/large-sparse-transportation-model-with.html
