func_figure_s31 <- function(fit){
  f0 = fit$f0
  fdm = fit$fdm
  VAR_1 = fit$VAR_1
  f4 = fit$f4
  f6 = fit$f6
  coulist = fit$coulist
  rms <- function(x){
    result = sqrt(mean(x^2))
  }
  tstep = dim(f0$fore)[2]
  glok = dim(f0$fore)[1]
  
  # if (length(coulist) == 14) {
  #   popind =  c(1:5,14)
  # }else{
  #   popind =  c(1:6)
  # }
  popind = 1:8
  tempgloi = length(popind)
  p2data1 = matrix(0,tstep,tempgloi)
  p2data3 = p2data1
  p2data4 = p2data1
  
  p2data0 = apply(f0$res,c(2,3),rms)[,popind]
  p2data1 = apply(fdm$residual,c(2,3),rms)[,popind]
  p2data2 = apply(VAR_1$res,c(2,3),rms)[,popind]
  
  for (i in 1:tempgloi) {
    p2data3[,i] = apply(f4$res[,((1:glok)+(i-1)*glok)],1,rms)
  }
  
  for (i in 1:tempgloi) {
    p2data4[,i] = apply(f6$res[,((1:glok)+(i-1)*glok)],1,rms)
  }
  
  p2data03 = p2data0
  p2data13 = p2data1
  p2data23 = p2data2
  p2data33 = p2data3
  p2data43 = p2data4
  
  Year = rep(1:tstep,tempgloi*5)
  value = c(as.vector(p2data0),as.vector(p2data1),as.vector(p2data2),as.vector(p2data3),as.vector(p2data4))
  foremat1 <- data.frame(Year)
  foremat1$Model = rep(c('Li_Lee','Co-FDM','STAR','SVAR','SWVAR'),each = tempgloi*tstep)
  foremat1$Pop = rep(rep(coulist[popind],each = tstep),5)
  foremat1$rmsfe = c(value)
  
  #rmsfeg3 = value
  pfore1 <- ggplot(foremat1, aes(
    x = Year, 
    y = rmsfe,
  ))
  p1 = pfore1 + geom_line(aes(color = Model))+xlim(1,tstep)+
    xlab(expression(atop("Year"~italic(h)~'from jump-off year')))+ylab(expression(RMSFE))+facet_wrap(~ Pop,ncol = 2)+
    theme(plot.title = element_blank())
  
  ########################################################
  popind = 9:14
  tempgloi = length(popind)
  p2data1 = matrix(0,tstep,tempgloi)
  p2data3 = p2data1
  p2data4 = p2data1
  
  p2data0 = apply(f0$res,c(2,3),rms)[,popind]
  p2data1 = apply(fdm$residual,c(2,3),rms)[,popind]
  p2data2 = apply(VAR_1$res,c(2,3),rms)[,popind]
  
  for (i in 1:tempgloi) {
    p2data3[,i] = apply(f4$res[,((1:glok)+(i-1)*glok)],1,rms)
  }
  
  for (i in 1:tempgloi) {
    p2data4[,i] = apply(f6$res[,((1:glok)+(i-1)*glok)],1,rms)
  }
  
  p2data03 = p2data0
  p2data13 = p2data1
  p2data23 = p2data2
  p2data33 = p2data3
  p2data43 = p2data4
  
  Year = rep(1:tstep,tempgloi*5)
  value = c(as.vector(p2data0),as.vector(p2data1),as.vector(p2data2),as.vector(p2data3),as.vector(p2data4))
  foremat1 <- data.frame(Year)
  foremat1$Model = rep(c('Li_Lee','Co-FDM','STAR','SVAR','SWVAR'),each = tempgloi*tstep)
  foremat1$Pop = rep(rep(coulist[popind],each = tstep),5)
  foremat1$rmsfe = c(value)
  
  #rmsfeg3 = value
  pfore1 <- ggplot(foremat1, aes(
    x = Year, 
    y = rmsfe,
  ))
  p2 = pfore1 + geom_line(aes(color = Model))+xlim(1,tstep)+
    xlab(expression(atop("Year"~italic(h)~'from jump-off year')))+ylab(expression(RMSFE))+facet_wrap(~ Pop,ncol = 2)+
    theme(plot.title = element_blank())
  return(list(p1,p2))
}
func_figure_s32 <- function(fit){
  f0 = fit$f0
  fdm = fit$fdm
  VAR_1 = fit$VAR_1
  f4 = fit$f4
  f6 = fit$f6
  coulist = fit$coulist
  rms <- function(x){
    result = sqrt(mean(x^2))
  }
  tstep = dim(f0$fore)[2]
  glok = dim(f0$fore)[1]
  
  # if (length(coulist) == 14) {
  #   popind =  c(1:5,14)
  # }else{
  #   popind =  c(1:6)
  # }
  popind = 1:length(coulist)
  tempgloi = length(popind)
  p2data1 = matrix(0,tstep,tempgloi)
  p2data3 = p2data1
  p2data4 = p2data1
  
  p2data0 = apply(f0$res,c(2,3),rms)[,popind]
  p2data1 = apply(fdm$residual,c(2,3),rms)[,popind]
  p2data2 = apply(VAR_1$res,c(2,3),rms)[,popind]
  
  for (i in 1:tempgloi) {
    p2data3[,i] = apply(f4$res[,((1:glok)+(i-1)*glok)],1,rms)
  }
  
  for (i in 1:tempgloi) {
    p2data4[,i] = apply(f6$res[,((1:glok)+(i-1)*glok)],1,rms)
  }
  
  p2data03 = p2data0
  p2data13 = p2data1
  p2data23 = p2data2
  p2data33 = p2data3
  p2data43 = p2data4
  
  Year = rep(1:tstep,tempgloi*5)
  value = c(as.vector(p2data0),as.vector(p2data1),as.vector(p2data2),as.vector(p2data3),as.vector(p2data4))
  foremat1 <- data.frame(Year)
  foremat1$Model = rep(c('Li_Lee','Co-FDM','STAR','SVAR','SWVAR'),each = tempgloi*tstep)
  foremat1$Pop = rep(rep(coulist[popind],each = tstep),5)
  foremat1$rmsfe = c(value)
  
  #rmsfeg3 = value
  pfore1 <- ggplot(foremat1, aes(
    x = Year, 
    y = rmsfe,
  ))
  p1 = pfore1 + geom_line(aes(color = Model))+xlim(1,tstep)+
    xlab(expression(atop("Year"~italic(h)~'from jump-off year')))+ylab(expression(RMSFE))+facet_wrap(~ Pop,ncol = 2)+
    theme(plot.title = element_blank())
  

  return(p1)
}