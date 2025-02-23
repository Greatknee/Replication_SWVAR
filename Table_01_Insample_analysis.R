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
source("functions/func_table_01_insample.R")

# Load data
data <- datapre_in(group = 5)

# Run analysis
result <- func_table_01_insample(data)

result

# Save output
#write.csv(result, "output/result.csv")