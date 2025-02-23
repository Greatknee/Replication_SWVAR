

datapre_in = function(group){
  if (group == 1) {
    coulist = c('AUT','CAN','DNK','ENG','FIN','FRA','ITA','JPN','NLD','NOR','ESP','SWE','SWI','USA')
    gloi = length(coulist)
    years = 1950
    yeare = 2018
    PYear = years:yeare
    glot = length(PYear)
    glok = 21
  }else if (group == 2) {
    #data pre group2
    coulist = c('BGR','CZE','HUN','LTU','RUS')
    gloi = length(coulist)
    years = 1959
    yeare = 2014
    PYear = years:yeare
    glot = length(PYear)
    glok = 21
  }else if (group == 3) {
    #data pre group3
    coulist = c('BLR','EST','LVA','LTU','RUS','UKR')
    gloi = length(coulist)
    years = 1959
    yeare = 2013
    PYear = years:yeare
    glot = length(PYear)
    glok = 21
  }else if (group == 4) {
    #data pre group4
    coulist = c('BGR','CZE','HUN','BLR','EST','LVA','LTU','RUS','UKR')
    gloi = length(coulist)
    years = 1959
    yeare = 2013
    PYear = years:yeare
    glot = length(PYear)
    glok = 21
  }else if (group == 5) {
    coulist = c('AUT','ENG','FRA','ITA','NLD','SWI')
    gloi = length(coulist)
    years = 1950
    yeare = 2018
    PYear = years:yeare
    glot = length(PYear)
    glok = 21
  }
  
  Mltaddress = paste('data/',coulist,'/mltper_5x1.txt',sep = '')
  
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
  return(list(group = group, coulist = coulist, datar = datar,data_star_low = data_star_low,years = years,yeare = yeare,glok = glok,glot = glot, gloi = gloi))
}



datapre_out = function(group){
  library(abind)
  if (group == 1) {
    coulist = c('AUT','CAN','DNK','ENG','FIN','FRA','ITA','JPN','NLD','NOR','ESP','SWE','SWI','USA')
    gloi = length(coulist)
    ystr = 1950
    yetr = 1999
    yste = 2000
    yete = 2019
    PYear = ystr:yetr
    glot = length(PYear)
    glok = 21
  }else if (group == 2) {
    #data pre group2
    coulist = c('BGR','CZE','HUN','LTU','RUS')
    gloi = length(coulist)
    ystr = 1959
    yetr = 1999
    yste = 2000
    yete = 2013
    glot = length(PYear)
    glok = 21
  }else if (group == 3) {
    #data pre group3
    coulist = c('BLR','EST','LVA','LTU','RUS','UKR')
    gloi = length(coulist)
    ystr = 1959#1959
    yetr = 1999#1999
    yste = 2000
    yete = 2013
    PYear = ystr:yetr
    glot = length(PYear)
    glok = 21
  }else if (group == 4) {
    return("We did not conduct out-of-sample analysis to Group 4")
  }else if (group == 5) {
    coulist = c('AUT','ENG','FRA','ITA','NLD','SWI')
    gloi = length(coulist)
    ystr = 1950
    yetr = 1999
    yste = 2000
    yete = 2019
    PYear = ystr:yetr
    glot = length(PYear)
    glok = 21
  }
  
  glote = length(yste:yete)
  Mltaddress = paste('data/',coulist,'/mltper_5x1.txt',sep = '')
  
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
    D = matrix(as.numeric(dataD1$dx),nrow = 24)[1:21,]
    E = matrix(as.numeric(dataD1$lx),nrow = 24)[1:21,]
    M = matrix(as.numeric(dataD1$mx),nrow = 24)[1:21,]
    #M = D/E
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
    tempO = TurnRawData1(tempdata,ys = ystr,ye = yetr)
    data[,,i] = tempO$logmortality  
    #data[,,i] = log(tempO$death/tempO$exposure)
    tempdmat = tempdmat+tempO$death
    tempemat = tempemat+tempO$exposure
  }
  
  datatr = data
  #processing for li-lee
  dim(datatr)
  data_star_low_tr = log(tempdmat/tempemat)
  data_star_low_tr[is.na(data_star_low_tr)] = min(data_star_low_tr)
  dim(data_star_low_tr)
  
  
  #create test dataset
  datate = array(0,dim = c(glok,glote,gloi))
  
  
  for (i in 1:gloi) {
    tempdata = read.table(Mltaddress[i],header = T,fill = T)
    tempO = TurnRawData1(tempdata,ys = yste,ye = yete)
    #datate[,,i] = log(tempO$death/tempO$exposure)
    datate[,,i] = tempO$logmortality  
    
  }
  
  dataall = abind(datatr, datate, along = 2)
  return(list(group = group, coulist = coulist, datatr = datatr, datate = datate, data_star_low_tr = data_star_low_tr, years = years,yeare = yeare,glok = glok,glot = glot, gloi = gloi))
}
