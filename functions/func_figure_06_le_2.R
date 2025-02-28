
func_figure_06_le_2 <- function(){
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
  
  
  data <- datapre_out(group = 2) 
  ystr = data$ystr
  yetr = data$yetr
  yste = data$yste
  yete = data$yete
  
  datagroup = data$group
  coulist = data$coulist
  datar = data$datatr
  datate = data$datate
  
  data_star_low_tr = data$data_star_low_tr
  glot = data$glot
  glote = data$glote
  glok = data$glok
  gloi = data$gloi
  datar = data$datatr
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
  
  
  simtime = 200
  tstep = 50
  #row data
  datar2 = abind(datar, datate, along = 2)
  e0_true = le_ez(exp(datar2))
  
  #Li-lee
  #forecast with random noise
  custompred <- function(c,t=glote,x0) {
    forec = matrix(0,nrow = 1,t)
    forec[1] = c[1]+c[2]*x0
    for (i in 2:t) {
      forec[i] = c[1]+c[2]*forec[i-1]
    }
    return(forec)
  }
  purefore_lilee_sim <- function(model, timestep, dtrain=datatr, dtest = datate){
    predk = matrix(0,nrow = gloi,ncol = timestep)
    fitvalue = matrix(0,gloi,2)
    r = matrix(0,gloi,1)
    fore = array(0,dim = c(glok,timestep,gloi))
    forematrix = matrix(0,timestep,gloi*glok)
    result = list()
    lK = model$fit$lK
    lB = model$fit$lB
    listk = model$fit$listk
    listb = model$fit$listb
    fit1 <- rwf(lK,drift = TRUE, h = timestep)
    sigK <- sqrt(fit1$model$sigma2)
    simerr <- rnorm(timestep, 0,1)
    simsum = cumsum(simerr) 
    plK <- fit1$mean
    for (i in 1:gloi) {
      dt0 = dtrain[,,i][,ncol(dtrain[,,i])]
      tsd = ts(listk[i,])
      fit <- ar(tsd,order.max =1,h=timestep)
      fitvalue[i,1] = fit$x.mean
      fitvalue[i,2] = fit$ar
      sigki = sqrt(fit$var.pred)
      simerr <- rnorm(timestep, 0,1)
      simsum = cumsum(simerr) 
      predk[i,]=custompred(fitvalue[i,],t=timestep,tsd[length(tsd)])+sigki*simsum
      r[i] = 1-(var(residuals(fit))/var(tsd))
      a = matrix(rep(dt0,timestep),glok,timestep)
      gp = lB%*%t(plK-lK[length(lK)])
      ip = listb[i,]%*%t(predk[i,]-listk[i,ncol(listk)])
      fore[,,i] = a+gp+ip
      forematrix[,((1:glok)+((i-1)*glok))]= t(a+gp+ip)
    }
    result$forecast = fore
    result$k= predk
    result$coef = fitvalue
    result$K = plK
    result$Rar = r
    return(result)
  }
  
  sime_ll = array(0,dim=c(tstep,gloi,simtime))
  
  for (s in 1:simtime) {
    f0 = fitlilee(datar,data_star_low_tr,datate)
    fll = purefore_lilee_sim(model = f0,timestep = tstep,datar,datate)
    lefore = le_ez(exp(fll$forecast))
    sime_ll[,,s] = lefore
  }
  
  emaxll= apply(sime_ll,c(1,2),quantile,0.95,na.rm =TRUE)
  eminll= apply(sime_ll,c(1,2),quantile,0.05,na.rm =TRUE)
  emidll = apply(sime_ll,c(1,2),quantile,0.5,na.rm =TRUE)
  
  # emaxllout_g3 = emaxll
  # eminllout_g3 = eminll
  # emidllout_g3 = emidll
  
  
  ############################################################
  # SVAR
  f6 = fitswvar_lack_global(datar = datar2[,1:(glot+10),], group = datagroup,datate = datate,rawdata = data)
  
  pureforecast <- function(model,x0,step) {
    #dim(ddatamat_c) = 49,294
    nc = gloi*glok
    nr = step
    m = c(model$mu)
    pred = matrix(0,ncol=nc,nrow=nr)
    pred[1,] = t(model$A[[1]] %*% (x0-m)) +m
    for (i in 2:nr) {
      pred[i,] = t(model$A[[1]] %*% (pred[(i-1),]-m)) +m
    }
    sum = list()
    sum$forecast = pred
    return(sum)
  }
  #dlm50 = pureforecast(VAR_4,ddatamat_c[glotr,],50)
  
  
  pureforecast_swvar <- function(coef,m,x0,step) {
    #dim(ddatamat_c) = 49,294
    nc = length(ddatamat_c[nrow(ddatamat_c),])
    nr = step
    pred = matrix(0,ncol=nc,nrow=nr)
    pred[1,] = t(coef %*% t(x0-m)) +m
    for (i in 2:nr) {
      pred[i,] = t(coef %*% t(pred[(i-1),]-m))+m
    }
    sum = list()
    sum$forecast = pred
    return(sum)
  }
  fore = pureforecast_swvar(f6$A,t(f6$mu),ddatamat_te[nrow(ddatamat_te),],50)
  
  simerr <- mvrnorm(tstep*simtime, rep(0,gloi*glok), Sigma = sig_c) #生成1000个三元正态分布的随机数：1000行乘以3列的
  sime = array(0,dim=c(tstep,gloi,simtime))
  
  for (s in 1:simtime) {
    dlm50 = fore$forecast+simerr[(1:50)+(s-1)*50,]
    lm50 = apply(dlm50,2,cumsum)+matrix(rep(datamat_c[nrow(datamat_c),],tstep),nrow = tstep,byrow = TRUE)
    lmten = array(0,dim = c(glok,tstep,gloi))
    for (i in 1:gloi) {
      lmten[,,i] = t(lm50[,((1:glok)+(i-1)*glok)])
    }
    lefore = le_ez(exp(lmten))
    sime[,,s] = lefore
  }
  
  emax= apply(sime,c(1,2),quantile,0.95)
  emin= apply(sime,c(1,2),quantile,0.05)
  emid = apply(sime,c(1,2),quantile,0.5)
  e0 = emid
  
  e0_true_g2 = e0_true
  emaxout_g2 = emax
  eminout_g2 = emin
  emidout_g2 = emid
  
  emaxllout_g2 = emaxll
  eminllout_g2 = eminll
  emidllout_g2 = emidll
  
  ##############################################
  ##############################################
  ##############################################
  ##############################################
  ##############################################
  ##############################################
  ##############################################
  ##############################################
  ##############################################
  
  data <- datapre_out(group = 3) 
  ystr = data$ystr
  yetr = data$yetr
  yste = data$yste
  yete = data$yete
  
  datagroup = data$group
  coulist = data$coulist
  datar = data$datatr
  datate = data$datate
  
  data_star_low_tr = data$data_star_low_tr
  glot = data$glot
  glote = data$glote
  glok = data$glok
  gloi = data$gloi
  datar = data$datatr
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
  
  
  simtime = 200
  tstep = 50
  #row data
  datar2 = abind(datar, datate, along = 2)
  e0_true = le_ez(exp(datar2))
  
  #Li-lee
  #forecast with random noise
  custompred <- function(c,t=glote,x0) {
    forec = matrix(0,nrow = 1,t)
    forec[1] = c[1]+c[2]*x0
    for (i in 2:t) {
      forec[i] = c[1]+c[2]*forec[i-1]
    }
    return(forec)
  }
  purefore_lilee_sim <- function(model, timestep, dtrain=datatr, dtest = datate){
    predk = matrix(0,nrow = gloi,ncol = timestep)
    fitvalue = matrix(0,gloi,2)
    r = matrix(0,gloi,1)
    fore = array(0,dim = c(glok,timestep,gloi))
    forematrix = matrix(0,timestep,gloi*glok)
    result = list()
    lK = model$fit$lK
    lB = model$fit$lB
    listk = model$fit$listk
    listb = model$fit$listb
    fit1 <- rwf(lK,drift = TRUE, h = timestep)
    sigK <- sqrt(fit1$model$sigma2)
    simerr <- rnorm(timestep, 0,1)
    simsum = cumsum(simerr) 
    plK <- fit1$mean
    for (i in 1:gloi) {
      dt0 = dtrain[,,i][,ncol(dtrain[,,i])]
      tsd = ts(listk[i,])
      fit <- ar(tsd,order.max =1,h=timestep)
      fitvalue[i,1] = fit$x.mean
      fitvalue[i,2] = fit$ar
      sigki = sqrt(fit$var.pred)
      simerr <- rnorm(timestep, 0,1)
      simsum = cumsum(simerr) 
      predk[i,]=custompred(fitvalue[i,],t=timestep,tsd[length(tsd)])+sigki*simsum
      r[i] = 1-(var(residuals(fit))/var(tsd))
      a = matrix(rep(dt0,timestep),glok,timestep)
      gp = lB%*%t(plK-lK[length(lK)])
      ip = listb[i,]%*%t(predk[i,]-listk[i,ncol(listk)])
      fore[,,i] = a+gp+ip
      forematrix[,((1:glok)+((i-1)*glok))]= t(a+gp+ip)
    }
    result$forecast = fore
    result$k= predk
    result$coef = fitvalue
    result$K = plK
    result$Rar = r
    return(result)
  }
  
  sime_ll = array(0,dim=c(tstep,gloi,simtime))
  
  for (s in 1:simtime) {
    f0 = fitlilee(datar,data_star_low_tr,datate)
    fll = purefore_lilee_sim(model = f0,timestep = tstep,datar,datate)
    lefore = le_ez(exp(fll$forecast))
    sime_ll[,,s] = lefore
  }
  
  emaxll= apply(sime_ll,c(1,2),quantile,0.95,na.rm =TRUE)
  eminll= apply(sime_ll,c(1,2),quantile,0.05,na.rm =TRUE)
  emidll = apply(sime_ll,c(1,2),quantile,0.5,na.rm =TRUE)
  
  # emaxllout_g3 = emaxll
  # eminllout_g3 = eminll
  # emidllout_g3 = emidll
  
  
  ############################################################
  # SVAR
  f6 = fitswvar_lack_global(datar = datar2[,1:(glot+10),], group = datagroup,datate = datate,rawdata = data)
  
  pureforecast <- function(model,x0,step) {
    #dim(ddatamat_c) = 49,294
    nc = gloi*glok
    nr = step
    m = c(model$mu)
    pred = matrix(0,ncol=nc,nrow=nr)
    pred[1,] = t(model$A[[1]] %*% (x0-m)) +m
    for (i in 2:nr) {
      pred[i,] = t(model$A[[1]] %*% (pred[(i-1),]-m)) +m
    }
    sum = list()
    sum$forecast = pred
    return(sum)
  }
  #dlm50 = pureforecast(VAR_4,ddatamat_c[glotr,],50)
  
  
  pureforecast_swvar <- function(coef,m,x0,step) {
    #dim(ddatamat_c) = 49,294
    nc = length(ddatamat_c[nrow(ddatamat_c),])
    nr = step
    pred = matrix(0,ncol=nc,nrow=nr)
    pred[1,] = t(coef %*% t(x0-m)) +m
    for (i in 2:nr) {
      pred[i,] = t(coef %*% t(pred[(i-1),]-m))+m
    }
    sum = list()
    sum$forecast = pred
    return(sum)
  }
  fore = pureforecast_swvar(f6$A,t(f6$mu),ddatamat_te[nrow(ddatamat_te),],50)
  
  simerr <- mvrnorm(tstep*simtime, rep(0,gloi*glok), Sigma = sig_c) #生成1000个三元正态分布的随机数：1000行乘以3列的
  sime = array(0,dim=c(tstep,gloi,simtime))
  
  for (s in 1:simtime) {
    dlm50 = fore$forecast+simerr[(1:50)+(s-1)*50,]
    lm50 = apply(dlm50,2,cumsum)+matrix(rep(datamat_c[nrow(datamat_c),],tstep),nrow = tstep,byrow = TRUE)
    lmten = array(0,dim = c(glok,tstep,gloi))
    for (i in 1:gloi) {
      lmten[,,i] = t(lm50[,((1:glok)+(i-1)*glok)])
    }
    lefore = le_ez(exp(lmten))
    sime[,,s] = lefore
  }
  
  emax= apply(sime,c(1,2),quantile,0.95)
  emin= apply(sime,c(1,2),quantile,0.05)
  emid = apply(sime,c(1,2),quantile,0.5)
  e0 = emid
  
  e0_true_g3 = e0_true
  emaxout_g3 = emax
  eminout_g3 = emin
  emidout_g3 = emid
  
  emaxllout_g3 = emaxll
  eminllout_g3 = eminll
  emidllout_g3 = emidll
  
  
  
  
  ###############################################
  #plot 
  coulist23 = c('BGR','CZE','HUN','LTU1','RUS1','BLR','EST','LVA','LTU2','RUS2','UKR')
  e0_true_23 = cbind(e0_true_g2[1:55,],e0_true_g3)
  emax_23 = cbind(emaxout_g2,emaxout_g3)
  emin_23 = cbind(eminout_g2,eminout_g3)
  emid_23 = cbind(emidout_g2,emidout_g3)
  emaxll_23 = cbind(emaxllout_g2,emaxllout_g3)
  eminll_23 = cbind(eminllout_g2,eminllout_g3)
  emidll_23 = cbind(emidllout_g2,emidllout_g3)
  
  index = c(1,2,3,5,8,10)
  tempi = length(index)
  LE = as.vector(emid_23[,index]) 
  year = rep(seq(2000,2049,1),tempi)
  pop= rep(coulist23[index],each = 50)
  df1 = data.frame(LE)
  df1$upper = as.vector(emax_23[,index])
  df1$lower = as.vector(emin_23[,index])
  df1$emid = as.vector(emid_23[,index])
  df1$year = year
  df1$Pop = pop
  df1$model = rep('SWVAR',tempi*50)
  
  df1$ell = as.vector(emidll_23[,index])
  df1$lowerll = as.vector(eminll_23[,index])
  df1$upperll = as.vector(emaxll_23[,index])
  year0 = rep(seq(ystr,yete,1),tempi)
  df2 = data.frame(year0)
  df2$eold = as.vector(e0_true_23[,index])
  df2$Pop = rep(coulist23[index],each = yete+1-ystr)
  
  pfore2 <- ggplot()
  
  
  p = pfore2 +  geom_line(data = df2, aes(x=year0, y = eold,colour = 'Observed'),linewidth = 0.6)+  facet_wrap(~ df2$Pop)+
    geom_smooth(data = df1,method = 'loess',se= FALSE,aes(x = year, y = emid,colour = 'SWVAR'),linewidth = 1,span = 0.35,linetype = 1)+
    geom_smooth(data = df1,method = 'loess',se= FALSE,aes(x = year, y = lower),linewidth = 0.4,span = 0.35,linetype = 2,color='blue')+
    geom_smooth(data = df1,method = 'loess',se= FALSE,aes(x = year, y = upper),linewidth = 0.4,span = 0.35,linetype = 2,color='blue')+
    geom_smooth(data = df1,method = 'loess',se= FALSE,aes(x = year, y = ell,colour = 'Li-Lee'),linewidth = 1,span = 0.35,linetype = 1)+
    geom_smooth(data = df1,method = 'loess',se= FALSE,aes(x = year, y = lowerll),linewidth = 0.4,span = 0.35,linetype = 2,color='red')+
    geom_smooth(data = df1,method = 'loess',se= FALSE,aes(x = year, y = upperll),linewidth = 0.4,span = 0.35,linetype = 2,color='red')+#geom_ribbon(aes(x=year,ymin=lowerll, ymax=upperll,fill = model), linetype=2, alpha=0.1)+
    ggtitle('Forecsting of Life Expectancy')+xlab('Year')+scale_y_continuous(breaks = c(40,50,60,70),labels = c(50,60,70,80))+
    theme(plot.title = element_blank(),axis.title.y = element_blank())+
    scale_color_manual(name = 'Group',values = c('red','black','blue'))+
    facet_wrap(~ Pop)
  return(p)
}
