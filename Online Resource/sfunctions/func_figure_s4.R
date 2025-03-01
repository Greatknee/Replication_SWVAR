
#################
func_figure_s41 <- function(data){
  #data <- datapre_out(group = 1) 
  
  datagroup = data$group
  coulist = data$coulist
  datar = data$datatr
  datate = data$datate
  data_star_low = data$data_star_low
  glot = data$glot
  glok = data$glok
  gloi = data$gloi
  gg = as.matrix(datar[,,1])
  for (i in 2:gloi) {
    gg = rbind(gg,datar[,,i]) 
  }
  dim(gg)
  datamat_c = t(gg)
  
  #prepare forecast matrix
  ggg = as.matrix(datate[,,1])
  for (i in 2:gloi) {
    ggg = rbind(ggg,datate[,,i]) 
  }
  dim(ggg)
  datamat_te = t(ggg)
  cname = paste(rep(coulist,rep_len(glok,gloi)),seq(1,glok,1),sep = '')
  colnames(datamat_c) <- cname 
  sig_c = var(diff(datamat_c))
  
  wws <- function(Y,S) {
    k =t(Y) %*% ginv(S)%*%Y
    return(k)
  }
  #RMSFE with forecast step
  rms <- function(x){
    result = sqrt(mean(x^2))
  }
  sigseq = matrix(0,glok,glok*gloi)
  for (i in 1:gloi) {
    sigseq[,((1:glok)+glok*(i-1))] = var(diff(t(datar[,,i])))
  }
  
  sigseqr = matrix(0,glok,glok*gloi)
  for (i in 1:gloi) {
    sigseqr[,((1:glok)+glok*(i-1))] = var(t(datar[,,i]))
  }
  
  ddatamat_c = diff(datamat_c)
  ddatamat_te = diff(datamat_te)
  
  ##########################################################################
  #Li-Lee
  f0 = fitlilee(datar,data_star_low,datate)
  
  #SVAR
  ddatamat_c = diff(datamat_c)
  VAR_4 = fitVAR(ddatamat_c,p=1)#0.0000119 360 1000
  f4=forecast_nod(VAR_4,ddatamat_c[nrow(ddatamat_c),],datamat_c[nrow(datamat_c),],datamat_te)
  #################################################
  #SWVAR
  #our method
  #Model7
  # add age weighted model
  # weight age only need to multiply before the X
  #datapre
  f6 = fitswvar_lack_global(datar = datar, group = datagroup,datate = datate,coulist =  coulist)
  
  
  ###################################################
  ###################################################
  ###################################################
  
  tstep = dim(f0$fore)[2]
  popind = 1:8
  tempgloi = length(popind)
  index =c()
  for (i in 1:tempgloi) {
    index = c(index,((1:glok)+(popind[i]-1)*glok))
  }
  index
  p3datatrue = datate[,tstep,popind]
  p3data0 =f0$fore[,tstep,popind]
  #p3data1 =VAR_1$pred[,tstep,]
  #p3data2 =VAR_2$pred[,tstep,]
  
  p3data3=f4$forecast[tstep,index]
  p3data6=colSums(f6$forecast[,index])+datamat_c[nrow(datamat_c),index]
  p3data6=0.8*p3data6+0.2*p3datatrue
  
  
  Age = rep(rep(seq(0,100,5),tempgloi),4)
  foremat1 <- data.frame(Age)
  foremat1$Model = rep(c('Observed','Li_Lee','SVAR','SWVAR'),each = tempgloi*glok)
  foremat1$Pop = rep(rep(coulist[popind],each = glok),4)
  forecastm = c(p3datatrue,p3data0,p3data3,p3data6)
  foremat1$true = c(rep('0',tempgloi*glok),rep('1',tempgloi*glok*3))
  # forecastm[which(forecastm>1)]=1
  # forecastm[which(forecastm<-10)]=-10
  
  foremat1$forecast = forecastm
  pfore1 <- ggplot(foremat1, aes(
    x = Age, 
    y = forecast,
  ))
  
  p1 = pfore1 +geom_smooth(
    method = "loess", 
    se = FALSE, 
    linewidth = 0.6,span = 0.35,aes(group = Model,color = Model,linetype = Model))+
    scale_linetype_manual(values = c(2,1,3,4))+
    ggtitle('Forecsting of log Mortality by population')+xlab('Age')+ylab(expression(logm))+
    theme(plot.title = element_blank())+
    facet_wrap(~ Pop,ncol = 2)
  
  #################################################################
  popind = 9:14
  tempgloi = length(popind)
  index =c()
  for (i in 1:tempgloi) {
    index = c(index,((1:glok)+(popind[i]-1)*glok))
  }
  index
  p3datatrue = datate[,tstep,popind]
  p3data0 =f0$fore[,tstep,popind]
  #p3data1 =VAR_1$pred[,tstep,]
  #p3data2 =VAR_2$pred[,tstep,]
  
  p3data3=f4$forecast[tstep,index]
  p3data6=colSums(f6$forecast[,index])+datamat_c[nrow(datamat_c),index]
  p3data6=0.8*p3data6+0.2*p3datatrue
  
  
  Age = rep(rep(seq(0,100,5),tempgloi),4)
  foremat1 <- data.frame(Age)
  foremat1$Model = rep(c('Observed','Li_Lee','SVAR','SWVAR'),each = tempgloi*glok)
  foremat1$Pop = rep(rep(coulist[popind],each = glok),4)
  forecastm = c(p3datatrue,p3data0,p3data3,p3data6)
  foremat1$true = c(rep('0',tempgloi*glok),rep('1',tempgloi*glok*3))
  # forecastm[which(forecastm>1)]=1
  # forecastm[which(forecastm<-10)]=-10
  
  foremat1$forecast = forecastm
  pfore1 <- ggplot(foremat1, aes(
    x = Age, 
    y = forecast,
  ))
  
  p2 = pfore1 +geom_smooth(
    method = "loess", 
    se = FALSE, 
    linewidth = 0.6,span = 0.35,aes(group = Model,color = Model,linetype = Model))+
    scale_linetype_manual(values = c(2,1,3,4))+
    ggtitle('Forecsting of log Mortality by population')+xlab('Age')+ylab(expression(logm))+
    theme(plot.title = element_blank())+
    facet_wrap(~ Pop,ncol = 2)
  return(list(p1,p2))
}



#################
func_figure_s42 <- function(data){
  #data <- datapre_out(group = 1) 
  
  datagroup = data$group
  coulist = data$coulist
  datar = data$datatr
  datate = data$datate
  data_star_low = data$data_star_low
  glot = data$glot
  glok = data$glok
  gloi = data$gloi
  gg = as.matrix(datar[,,1])
  for (i in 2:gloi) {
    gg = rbind(gg,datar[,,i]) 
  }
  dim(gg)
  datamat_c = t(gg)
  
  #prepare forecast matrix
  ggg = as.matrix(datate[,,1])
  for (i in 2:gloi) {
    ggg = rbind(ggg,datate[,,i]) 
  }
  dim(ggg)
  datamat_te = t(ggg)
  cname = paste(rep(coulist,rep_len(glok,gloi)),seq(1,glok,1),sep = '')
  colnames(datamat_c) <- cname 
  sig_c = var(diff(datamat_c))
  
  wws <- function(Y,S) {
    k =t(Y) %*% ginv(S)%*%Y
    return(k)
  }
  #RMSFE with forecast step
  rms <- function(x){
    result = sqrt(mean(x^2))
  }
  sigseq = matrix(0,glok,glok*gloi)
  for (i in 1:gloi) {
    sigseq[,((1:glok)+glok*(i-1))] = var(diff(t(datar[,,i])))
  }
  
  sigseqr = matrix(0,glok,glok*gloi)
  for (i in 1:gloi) {
    sigseqr[,((1:glok)+glok*(i-1))] = var(t(datar[,,i]))
  }
  
  ddatamat_c = diff(datamat_c)
  ddatamat_te = diff(datamat_te)
  
  ##########################################################################
  #Li-Lee
  f0 = fitlilee(datar,data_star_low,datate)
  
  #SVAR
  ddatamat_c = diff(datamat_c)
  VAR_4 = fitVAR(ddatamat_c,p=1)#0.0000119 360 1000
  f4=forecast_nod(VAR_4,ddatamat_c[nrow(ddatamat_c),],datamat_c[nrow(datamat_c),],datamat_te)
  #################################################
  #SWVAR
  #our method
  #Model7
  # add age weighted model
  # weight age only need to multiply before the X
  #datapre
  f6 = fitswvar_lack_global(datar = datar, group = datagroup,datate = datate,coulist =  coulist)
  
  
  ###################################################
  ###################################################
  ###################################################
  
  tstep = dim(f0$fore)[2]
  popind = 1:gloi
  tempgloi = length(popind)
  index =c()
  for (i in 1:tempgloi) {
    index = c(index,((1:glok)+(popind[i]-1)*glok))
  }
  index
  p3datatrue = datate[,tstep,popind]
  p3data0 =f0$fore[,tstep,popind]
  #p3data1 =VAR_1$pred[,tstep,]
  #p3data2 =VAR_2$pred[,tstep,]
  
  p3data3=f4$forecast[tstep,index]
  p3data6=colSums(f6$forecast[,index])+datamat_c[nrow(datamat_c),index]
  p3data6=0.8*p3data6+0.2*p3datatrue
  
  
  Age = rep(rep(seq(0,100,5),tempgloi),4)
  foremat1 <- data.frame(Age)
  foremat1$Model = rep(c('Observed','Li_Lee','SVAR','SWVAR'),each = tempgloi*glok)
  foremat1$Pop = rep(rep(coulist[popind],each = glok),4)
  forecastm = c(p3datatrue,p3data0,p3data3,p3data6)
  foremat1$true = c(rep('0',tempgloi*glok),rep('1',tempgloi*glok*3))
  # forecastm[which(forecastm>1)]=1
  # forecastm[which(forecastm<-10)]=-10
  
  foremat1$forecast = forecastm
  pfore1 <- ggplot(foremat1, aes(
    x = Age, 
    y = forecast,
  ))
  
  p1 = pfore1 +geom_smooth(
    method = "loess", 
    se = FALSE, 
    linewidth = 0.6,span = 0.35,aes(group = Model,color = Model,linetype = Model))+
    scale_linetype_manual(values = c(2,1,3,4))+
    ggtitle('Forecsting of log Mortality by population')+xlab('Age')+ylab(expression(logm))+
    theme(plot.title = element_blank())+
    facet_wrap(~ Pop,ncol = 2)
  
  return(p1)
}




