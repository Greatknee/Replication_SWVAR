# Load necessary libraries
library(abind)
library(ggplot2)
library(patchwork)

#Set the directory
setwd("C://Users//greatknee//Desktop//Mortality//Reproducibility//Replication_SWVAR")

# Source functions
source("functions/datapre.R")
source("functions/func_figure_08_lilee.R")


# Fit Li-Lee
out <- func_figure_08_lilee() 

#
glok = out$data$glok
gloi = out$data$gloi
glot = out$data$glot
coulist = out$data$coulist
# Figure 8.a
par(mfrow=c(1,2))
plot(c(0,1,seq(5,99,5)),out$lB,type = 'l',xlim=c(0,100),xlab = 'Age',ylab = 'Bx')
plot(1950:(1950+glot-1),out$lK,type = 'l',xlim=c(1950,2000),xlab = 'Year',ylab = 'Kt')

# Figure 8.b
dfb = data.frame(bx = as.vector(out$listb))
dfb$age = rep(1:glok,each = gloi)
dfb$pop = rep(coulist, glok)
g1 = ggplot(dfb,aes(x = age,y = bx,color = pop))+geom_point()+geom_line()

dfk = data.frame(kt = as.vector(out$listk))
dfk$year = rep(1:glot,each = gloi)
dfk$pop = rep(coulist, glot)
g2 = ggplot(dfk,aes(x = year,y = kt,color = pop))+geom_point()+geom_line()

g1+g2



