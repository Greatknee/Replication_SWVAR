# Load necessary libraries
library(dplyr)
library(ggplot2)

#Set the directory
setwd("C://Users//greatknee//Desktop//Mortality//Reproducibility//Replication_SWVAR")

# Source functions
source("functions/utils.R")
source("functions/datapre.R")
source("functions/model_benchmark.R")
source("functions/model_SWVAR.R")
source("functions/model_vartest.R")
source("functions/func_table_0203_test.R")



# # Sample
# # Load data
# data <- datapre_in(group = 1)
# # Run analysis
# result <- func_table_0203_test(data)
# result

# Formal
resultC = matrix(0,nrow = 5,ncol = 4)
resultJ = matrix(0,nrow = 5,ncol = 2)

for (i in 1:5) {
  data <- datapre_in(group = i)

  # Run analysis
  result <- func_table_0203_test(data)
  resultC[i,] = result$Ctest
  resultJ[i,] = result$Jtest
}

resultC
resultJ

#Save output
write.csv(resultC, "output/result_table02.csv")
write.csv(resultJ, "output/result_table03.csv")



