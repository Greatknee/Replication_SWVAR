setwd("C://Users//greatknee//Desktop//Mortality")


# #data pre
# coulist = c('AUT','CAN','DNK','ENG','FIN','FRA','ITA','JPN','NLD','NOR','ESP','SWE','SWI','USA')
# #coulist = c('AUT','ENG','FRA','ITA','NLD','SWI')
# gloi = length(coulist)
# years = 1950
# yeare = 2018
# PYear = years:yeare
# glot = length(PYear)
# glok = 21

# #data pre group2
# coulist = c('BGR','CZE','HUN','LTU','RUS')
# gloi = length(coulist)
# years = 1959
# yeare = 2014
# PYear = years:yeare
# glot = length(PYear)
# glok = 21
# #  
# #data pre group3
# coulist = c('BLR','EST','LVA','LTU','RUS','UKR')
# gloi = length(coulist)
# years = 1959
# yeare = 2013
# PYear = years:yeare
# glot = length(PYear)
# glok = 21

#data pre group4
coulist = c('BGR','CZE','HUN','BLR','EST','LVA','LTU','RUS','UKR')
gloi = length(coulist)
years = 1959
yeare = 2013
PYear = years:yeare
glot = length(PYear)
glok = 21

Mltaddress = paste('data/',coulist,'/STATS/mltper_5x1.txt',sep = '')


#用male 总表算的
TurnRawData1 <-  function(rdatad,ys=1950,ye=2019) {
  dataD = rdatad[-1,1:10]
  colnames(dataD) = rdatad[1,1:10]
  PYear = ys:ye
  NumY = length(PYear)
  NumA = glok
  datayear = as.numeric(dataD$Year)
  dataD1 = subset(dataD, dataD$Year >= ys &dataD$Year <= ye )
  
  #--------------D,E,M,m----------------------
  D = matrix(as.numeric(dataD1$dx),nrow = 24)[1:glok,]
  E = matrix(as.numeric(dataD1$lx),nrow = 24)[1:glok,]
  M = matrix(as.numeric(dataD1$mx),nrow = 24)[1:glok,]
  m = log(M)
  Output = list("death" = D, "exposure" = E,'mortality'=M,"logmortality"=m)
  return(Output)
}

#------------------B,F,G,S----------------------------
#---------------------A,Y,Pop-------------------------------
data = array(0,dim = c(glok,glot,gloi))
tempdmat = matrix(0,glok,glot)
tempemat = matrix(0,glok,glot)


for (i in 1:gloi) {
  tempdata = read.table(Mltaddress[i],header = T,fill = T)
  tempO = TurnRawData1(tempdata,years,yeare)
  data[,,i] = tempO$logmortality  
  
  tempdmat = tempdmat+tempO$death
  tempemat = tempemat+tempO$exposure
}

datar = data
data_star_low = log(tempdmat/tempemat)


# 
# ##################################################
# #For female
# Mltaddress = paste('data/',coulist,'/STATS/fltper_5x1.txt',sep = '')
# 
# #------------------B,F,G,S----------------------------
# #---------------------A,Y,Pop-------------------------------
# data = array(0,dim = c(glok,glot,gloi))
# tempdmat = matrix(0,glok,glot)
# tempemat = matrix(0,glok,glot)
# 
# 
# for (i in 1:gloi) {
#   tempdata = read.table(Mltaddress[i],header = T,fill = T)
#   tempO = TurnRawData1(tempdata,years,yeare)
#   data[,,i] = tempO$logmortality  
#   
#   tempdmat = tempdmat+tempO$death
#   tempemat = tempemat+tempO$exposure
# }
# 
# datar_f = data
# data_star_low_f = log(tempdmat/tempemat)

