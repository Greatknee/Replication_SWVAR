# Load necessary libraries
library(dplyr)
library(ggplot2)
library(MASS)
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


start_time = Sys.time()
# Formal
resultC = matrix(0,nrow = 5,ncol = 4)
resultJ = matrix(0,nrow = 5,ncol = 2)
colnames(resultC) = c("Estimation", "Std","t statistic","pvalue")
colnames(resultJ) = c("Lambda statistic","pvalue")

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
write.csv(resultC, "output/result_table_02.csv")
write.csv(resultJ, "output/result_table_03.csv")

end_time <- Sys.time()
print(end_time - start_time)


