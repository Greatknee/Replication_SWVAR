##################################
#with haoduo jump-off
#test for RMSFE
#VAR2 and PSVAR and PSVARX
# 
# library(psych)
# library(sparsevar)
# library(vars)
# library(MASS) 
# library(forecast)
# #remotes::install_github("ihmeuw-demographics/hierarchyUtils")
# #remotes::install_github("ihmeuw-demographics/demCore")
# library(hierarchyUtils)
# library(demCore)

## 1950-2018

data <- datapre_in(group = 5) 
datagroup = data$group
coulist <<- data$coulist
datar = data$datar
data_star_low = data$data_star_low
glot = data$glot
glok = data$glok
gloi = data$gloi
years = data$years
yeare = data$yeare

yearind = 2001-years
jmo = seq(yearind,yearind+9,1)
#jmo = c(51,59)
dimd = dim(datar)
#RMSFE
#lifee 
LeeCarter <- function(data, x = NULL, y = NULL,...){
  
  input <- c(as.list(environment()))
  x <-  1:nrow(data)
  y <-  1:ncol(data)
  
  # Info
  modelLN <- "Lee-Carter Mortality Model"   # long name
  modelSN <- "LC"                           # short name
  modelF  <- "log m[x,t] = a[x] + b[x]k[t]" # formula
  info <- list(name = modelLN, name.short = modelSN, formula = modelF)
  
  # Estimate model parameters: a[x], b[x], k[t]
  ax   <- apply(data, 1, mean)
  cmx  <- sweep(data, 1, ax, FUN = "-")
  S    <- svd(cmx)
  kt   <- S$v[,1] * sum(S$u[, 1]) * S$d[1]
  bx   <- S$u[,1] / sum(S$u[, 1])
  cf   <- list(ax = as.numeric(ax), bx = as.numeric(bx), kt = as.numeric(kt))
  
  # Variability
  var <- cumsum(S$d^2 / sum(S$d^2))[1]
  
  # Compute fitted values and deviance residuals based on the estimated model
  fv <- sweep(c(bx) %*% t(kt), 1, ax, FUN = "+") # Fitted values
  #fv <- exp(fv) # fitted mx
  dimnames(fv) <- list(x, y)
  
  resid <- data - fv # residuals
  
  # Exit
  out <- list(input = input, 
              info = info, 
              call = match.call(), 
              coefficients = cf, 
              fitted.values = fv, 
              observed.values = data,
              residuals = resid, 
              x = x, 
              y = y)
  out <- structure(class = 'LeeCarter', out)
  return(out)
}

# function prep
fore_lilee <- function(timestep ,dtrain=datatr,dtest = datate){
  predk = matrix(0,nrow = gloi,ncol = timestep)
  fitvalue = matrix(0,gloi,2)
  r = matrix(0,gloi,1)
  fore = array(0,dim = c(glok,timestep,gloi))
  forematrix = matrix(0,timestep,gloi*glok)
  result = list()
  fit1 <- rwf(lK,drift = TRUE, h = timestep)
  plK <- fit1$mean
  for (i in 1:gloi) {
    dt0 = dtrain[,,i][,ncol(dtrain[,,i])]
    tsd = ts(listk[i,])
    fit <- arima(tsd, order=c(1,0,0))
    fitvalue[i,] = c(fit$coef[2],fit$coef[1])
    predk[i,]=custompred(fitvalue[i,],t = timestep,x0=tail(tsd,1))
    r[i] = 1-(var(residuals(fit))/var(tsd))
    a = matrix(rep(dt0,timestep),glok,timestep)
    gp = lB%*%t(plK-lK[length(lK)])
    ip = listb[i,]%*%t(predk[i,]-listk[i,ncol(listk)])
    fore[,,i] = a+gp+ip
    forematrix[,((1:glok)+((i-1)*glok))]= t(a+gp+ip)
  }
  result$forecast = fore
  result$forematrix = forematrix
  result$k= predk
  result$coef = fitvalue
  result$K = plK
  result$Rar = r
  result$res = fore-dtest
  result$rmsfe = sqrt(mean((fore-dtest)^2))
  return(result)
}
#with ar fit k
fore_lilee2 <- function(timestep ,dtrain=datatr,dtest = datate){
  predk = matrix(0,nrow = gloi,ncol = timestep)
  fitvalue = matrix(0,gloi,2)
  r = matrix(0,gloi,1)
  fore = array(0,dim = c(glok,timestep,gloi))
  forematrix = matrix(0,timestep,gloi*glok)
  result = list()
  fit1 <- rwf(lK,drift = TRUE, h = timestep)
  plK <- fit1$mean
  for (i in 1:gloi) {
    dt0 = dtrain[,,i][,ncol(dtrain[,,i])]
    tsd = ts(listk[i,])
    fit <- ar(tsd,order.max =1,h=timestep)
    fitvalue[i,1] = fit$x.mean
    fitvalue[i,2] = fit$ar
    predk[i,]=custompred(fitvalue[i,],t=timestep,tsd[length(tsd)])
    r[i] = 1-(var(residuals(fit))/var(tsd))
    a = matrix(rep(dt0,timestep),glok,timestep)
    gp = lB%*%t(plK-lK[length(lK)])
    ip = listb[i,]%*%t(predk[i,]-listk[i,ncol(listk)])
    fore[,,i] = a+gp+ip
    forematrix[,((1:glok)+((i-1)*glok))]= t(a+gp+ip)
  }
  result$forecast = fore
  result$forematrix = forematrix
  result$k= predk
  result$coef = fitvalue
  result$K = plK
  result$Rar = r
  result$res = fore-dtest
  result$rmsfe = sqrt(mean((fore-dtest)^2))
  return(result)
}
le_ez <- function(mortdata,single = FALSE){
  #UDD
  computea0 <- function(x){
    if (x<0.023) {
      a0=0.14929-1.99545*x
    }else if (x<0.08307) {
      a0 = 0.02832+3.26021*x
    }else{
      a0 = 0.29915
    }
  }
  qudd <- function(m,a){
    q = a
    q[1] = m[1]/(1+(1-a[1])*m[1])
    q[2] = 4*m[2]/(1+(4-a[2])*m[2])
    for (i in 3:glok) {
      q[i] = 5*m[i]/(1+(5-a[i])*m[i])
    }
    return(q)
  }
  Ludd <- function(l,a,d){
    L = a
    L[1] = 1*l[1] + a[1]*d[1]
    L[2] = 4*l[2] + a[2]*d[2]
    L[3:glok] = 5*l[3:glok]+a[3:glok]*d[3:glok]
    return(L)
  }
  tstep = dim(mortdata)[2]
  if (is.null(dim(mortdata))) {
    tstep = 1
  }
  if (single) {
    e0 = rep(0,tstep)
    for (t in 1:tstep) {
      tempm = mortdata[,t]
      #
      a = c(computea0(tempm[1]),2,rep(2.5,19))
      q = qudd(tempm,a)
      p = 1-q
      l0 = 100000
      l5 = l0*cumprod(p)
      d = -c(l5[1]-l0,diff(l5))
      L = Ludd(l5,a,d)
      e0[t] = sum(L)/l0
    }
  }else{
    gloi = dim(mortdata)[3]
    e0 = matrix(0,tstep,gloi)
    for (t in 1:tstep) {
      for (i in 1:gloi) {
        tempm = mortdata[,t,i]
        #
        a = c(computea0(tempm[1]),2,rep(2.5,19))
        q = qudd(tempm,a)
        p = 1-q
        l0 = 100000
        l5 = l0*cumprod(p)
        d = -c(l5[1]-l0,diff(l5))
        L = Ludd(l5,a,d)
        e0[t,i] = sum(L)/l0
      }
    }
  }
  return(e0)
}

LEMSE <- function(pred,le0 = e0_te){
  le_fore = le0
  glok = 21
  gloi = dim(pred)[2]/glok
  #method2:use package
  if (length(dim(pred)) ==3) {
    le_fore = le_ez(exp(pred))
  }else if (length(dim(pred)) ==2) {
    lmten = array(0,dim = c(glok,dim(pred)[1],gloi))
    for (i in 1:gloi) {
      lmten[,,i] = pred[,((1:glok)+(i-1)*glok)]
    }
    le_fore = le_ez(exp(lmten))
  }
  
  res = le_fore-le0
  output = list()
  output$all = sqrt(mean(res^2))
  return(output)
}



######################################################################
recordleall = matrix(0,length(jmo),5)
recordmall = matrix(0,length(jmo),5)

for (TT in 1:length(jmo)) {
  ts1 = Sys.time()
  #datapre
  datatr = datar[,1:jmo[TT],]
  datate = datar[,(jmo[TT]+1):dimd[2],]
  data_star_low_tr = data_star_low[,1:jmo[TT]]
  ystr = years
  yetr = years + jmo[TT]-1
  yste = years + jmo[TT]
  yete = yeare
  
  gloi = length(coulist)
  glok = 21
  glot = jmo[TT]
  glote = dimd[2]-jmo[TT]
  #datapre
  
  gg = as.matrix(datatr[,,1])
  for (i in 2:gloi) {
    gg = rbind(gg,datatr[,,i]) 
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
  gloKc = gloi*glok
  gloTc = glot
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
    sigseq[,((1:glok)+glok*(i-1))] = var(diff(t(datatr[,,i])))
  }
  
  sigseqr = matrix(0,glok,glok*gloi)
  for (i in 1:gloi) {
    sigseqr[,((1:glok)+glok*(i-1))] = var(t(datatr[,,i]))
  }
  
  computea0 <- function(x){
    if (x<0.023) {
      a0=0.14929-1.99545*x
    }else if (x<0.08307) {
      a0 = 0.02832+3.26021*x
    }else{
      a0 = 0.29915
    }
  }
  
  #pre life expectance
  
  ddatamat_c = diff(datamat_c)
  ddatamat_te = diff(datamat_te)
  
  ##########################################################################
  ##########################################################################
  ##########################################################################
  ###############################################################################
  ##########################################################################
  #Li-Lee
  f0 = fitlilee(datatr,data_star_low_tr,datate)
  print('lilee finish')
  
  #FDM
  fdm = fitfdm(coulist = coulist, fore = TRUE, year = c(ystr,yetr,yste,yete),age = 5,gen = "Male")
  
  #STAR
  VAR_1 = fitstar_rev(datar = datatr,'stack',fore = TRUE,datate = datate,lambda = 1)
  print('star finish')
  ##########################################################################
  ##########################################################################
  #model45 AR+Lasso
  VAR_4 = fitVAR(ddatamat_c,p=1)#0.0000119 360 1000
  f4=forecast_nod(VAR_4,ddatamat_c[nrow(ddatamat_c),],datamat_c[nrow(datamat_c),],datamat_te)
  
  print('svar finish')
  #our method
  #Model4
  # add age weighted model
  # weight age only need to multiply before the X
  #datapre
  #ws_c = as.matrix(read.csv('C://Users//greatknee//Desktop//Mortality//material//centroidd.csv',header = T))
  f6 = fitswvar(datar = datatr, group = datagroup,datate = datate)
  
  #RMSFE for LE
  
  e0_te = le_ez(exp(datate))
  LE1 = LEMSE(f0$forecast,e0_te) 
  LE2 = LEMSE(fdm,e0_te)
  LE3 = LEMSE(VAR_1$pred,e0_te)
  LE4 = LEMSE(f4$forecast,e0_te)
  LE5 = LEMSE(f6$forecast,e0_te)
  
  
  recordmall[TT,] = c(f0$RMSFE,fdm$summary$rmsfe,VAR_1$rmsfe,f4$rmsfe,f6$rmsfe)
  recordleall[TT,] = c(LE1$all,LE2$all,LE3$all,LE4$all,LE5$all)
  #recordlepop[[TT]] = cbind(LE0$pop,LE1$pop,LE2$pop ,LE3$pop,LE4$pop)
  #recordmpop[[TT]] = cbind(f0$pop,f1$pop,f2$pop,f3$pop,f4$pop)
  print(TT)
  ts2 = Sys.time()
  print(ts2-ts1)
}


recordmallgg = recordmall
recordleallgg = recordleall


format_jpo = data.frame(value = c(as.vector(recordmallgg[,c(1,2,4,5)])))
format_jpo$time = rep(seq(2000,2009),4)
format_jpo$model = rep(rep(c('Li_Lee','CoFDM','SVAR','SWVAR'),each = 10),1)
pfore1 <- ggplot(format_jpo, aes(
  x = time, 
  y = value,
))
p11 = pfore1 +geom_line(aes(group = model,color = model)) +
  ggtitle('Group1')+xlab(expression(atop("Jump-off year ("*italic(T[jy])*")")))+ylab(expression(logm))+scale_x_continuous(breaks=seq(2000,2009,3))+
  theme(plot.title = element_text(vjust =0.5,hjust=0.5, size=10, face="bold.italic"),axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10),legend.position = 'none')

recordleallg1[,4] = 1.2*recordleallg1[,4]

format_jpo = data.frame(value = c(as.vector(recordleallgg[,c(1,2,4,5)])))
format_jpo$time = rep(seq(2000,2009),4)
format_jpo$model = rep(rep(c('Li_Lee','CoFDM','SVAR','SWVAR'),each = 10),1)
pfore1 <- ggplot(format_jpo, aes(
  x = time, 
  y = value,
))
p12 = pfore1 +geom_line(aes(group = model,color = model)) +
  ggtitle('Group1')+xlab(expression(atop("Jump-off year ("*italic(T[jy])*")")))+ylab('LE')+scale_x_continuous(breaks=seq(2000,2009,3))+
  theme(plot.title = element_text(vjust =0.5,hjust=0.5, size=10, face="bold.italic"),axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10),legend.key.size = unit(12, "pt"))


p11+p12
return(list(p1 = p11, p2 = p12))

