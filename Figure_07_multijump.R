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

#Set the directory
setwd("C://Users//greatknee//Desktop//Mortality//Reproducibility//Replication_SWVAR")

# Source functions
source("functions/utils.R")
source("functions/utils_multijump.R")
source("functions/datapre.R")
source("functions/model_benchmark.R")
source("functions/model_SWVAR.R")
source("functions/func_figure_07_multijump.R")


####################################
#Now there are 45 age groups, fitting the STAR model will cost more then 1 hours.
#You can change the parameter in main function "func_table_07_robust_g2_1(~,fitstar= TRUE)" to decide whether carry out it, 
#the default in the reproduction package is not run it.


start_time = Sys.time()
group = c(1,2,3,5)
flist = list()
for (i in 1:length(group)) {
  # Load data
  data <- datapre_out(group = group[i])
  # Figure 6.a
  f1 = func_figure_07_multijump(data)
  flist[[i]] = f1
}

flist[[i]]
pdf('output/result_figure_071.pdf', width=10, height=8)
print(f1)
dev.off()

end_time <- Sys.time()
print(end_time - start_time)