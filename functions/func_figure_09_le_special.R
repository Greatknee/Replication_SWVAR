
############
#LE
func_figure_09_le_special <- function(){
  set.seed(123)
  data <- datapre_out(group = 3)
  tstep = 50
  datagroup = data$group
  coulist = data$coulist
  datar = data$datatr
  datate = data$datate
  dataall = data$dataall
  data_star_low_tr = data$data_star_low_tr
  ystr = data$ystr
  yetr = data$yetr
  yste = data$yste
  yete = data$yete
  
  glot = data$glot
  glok = data$glok
  gloi = data$gloi
  datar = data$datatr
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
      for (i in 3:21) {
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
  
  gg = as.matrix(datar[,,1])
  for (i in 2:gloi) {
    gg = rbind(gg,datar[,,i]) 
  }
  dim(gg)
  datamat_c = t(gg)
  e0_true = le_ez(exp(dataall[,,5]),single = TRUE)
  
  #prepare forecast matrix
  ggg = as.matrix(datate[,,1])
  for (i in 2:gloi) {
    ggg = rbind(ggg,datate[,,i]) 
  }
  dim(ggg)
  datamat_te = t(ggg)
  ddatamat_c = diff(datamat_c)
  ddatamat_te = diff(datamat_te)
  

  
  ###############################
  #Lee-Carter
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
  custompred <- function(model,t=glote,x0) {
    forec = matrix(0,nrow = 1,t)
    forec[1] = model$coef[3]+model$coef[1]*x0
    for (i in 2:t) {
      forec[i] = model$coef[3]+model$coef[1]*forec[i-1]
    }
    return(forec)
  }
  LC = LeeCarter(datar[,,5])
  fore_lc <- function(model, tstep){
    pk = arma(model$coefficients$kt)
    predk = custompred(pk,t=tstep,model$coefficients$kt[length(model$coefficients$kt)])
    fore = matrix(0,glok,tstep)
    for (i in 1:tstep) {
      fore[,i] = LC$coefficients$ax+predk[i]*model$coefficients$bx
    }
    return(fore)
  }
  le_lc_rus = le_ez(exp(fore_lc(LC,50)),single = TRUE)
  LC = LeeCarter(datar[,,3])
  le_lc_lva = le_ez(exp(fore_lc(LC,50)),single = TRUE)
  
  
  
  ##################################
  #fdm
  #observed value
  filelist = c("Deaths_5x1.txt", "Exposures_5x1.txt", "Population5.txt", "Mx_5x1.txt")
  Mltaddress = paste('data/',coulist[1],'/',filelist,sep = '')
  tempi= length(coulist)
  tempk = 21
  datate_fdm = array(0,dim = c(tempk,1,tempi))
  
  tempdatag <- read_hmd_files(Mltaddress)
  df = as.data.frame(tempdatag |> 
                       dplyr::filter(Year >= ystr & Year <= yetr & Sex == 'Male' & Age<100))
  df$popc = coulist[1]
  
  for (i in 2:tempi) {
    Mltaddress = paste('data/',coulist[i],'/',filelist,sep = '')
    tempdatag <- read_hmd_files(Mltaddress)
    tempdata = as.data.frame(tempdatag |> 
                               dplyr::filter(Year >= ystr & Year <= yetr& Sex == 'Male' & Age<100))
    tempdata$popc = coulist[i]
    df = dplyr::bind_rows(df,tempdata)
    
  }
  names(df)[names(df) == "Mx_5x1"] <- "Mortality"
  df$Sex <- NULL
  names(df)[names(df) == "popc"] <- "Sex"
  names(df)[names(df) == "Population5"] <- "Population"
  df$Deaths_5x1 <- NULL
  df$Exposures_5x1 <- NULL
  data = as_vital(as_tibble(df),.age = 'Age',.sex = 'Sex',.population = 'Population',key = c('Age','Sex'),index = 'Year')
  
  ######################
  data$Mortality = exp(as.vector(aperm(datar, c(1, 3, 2))))  # 调整维度顺序，使得顺序为 (d1, d3, d2)
  
  #single
  hu <- data |> 
    dplyr::filter(Sex == 'RUS') |>
    model(hyndman_ullah = FDM(log(Mortality))) |>
    vital::forecast(h = tstep)
  humat = t(matrix(hu$.mean,tstep,glok))
  le_fdm_rus = le_ez(humat,single = TRUE)
  hu <- data |> 
    dplyr::filter(Sex == 'LVA') |>
    model(hyndman_ullah = FDM(log(Mortality))) |>
    vital::forecast(h = tstep)
  humat = t(matrix(hu$.mean,tstep,glok))
  le_fdm_lva = le_ez(humat,single = TRUE)
  
  
  
  #coFDM
  pr <- data |>
    make_pr(Mortality)
  hu <- pr |>
    model(hyndman_ullah = FDM(log(Mortality),coherent = TRUE))|> 
    vital::forecast(h =tstep)
  huundo <- hu |>
    undo_pr(Mortality,key = 'Sex',times = 1000)
  hufore = array(huundo$.mean,dim = c(tempk,tempi,tstep))
  fore_fdm <- aperm(hufore, perm = c(1,3,2))
  #le_fdm = le_ez(cbind(exp(datar[,,5]),fore_fdm[,,2]),single = TRUE)
  #plot(le_fdm,type = 'l')
  le_cofdm_rus = le_ez(fore_fdm[,,5],single = TRUE)
  le_cofdm_lva = le_ez(fore_fdm[,,3],single = TRUE)
  
  #################################
  #lilee
  
  custompred <- function(c,t=glote,x0) {
    forec = matrix(0,nrow = 1,t)
    forec[1] = c[1]+c[2]*x0
    for (i in 2:t) {
      forec[i] = c[1]+c[2]*forec[i-1]
    }
    return(forec)
  }
  purefore_lilee <- function(model, timestep, dtrain=datatr, dtest = datate){
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
    result$k= predk
    result$coef = fitvalue
    result$K = plK
    result$Rar = r
    return(result)
  }
  f0 = fitlilee(datar,data_star_low_tr,datate)
  fll = purefore_lilee(model = f0,timestep = tstep,datar,datate)
  le_ll_rus = le_ez(exp(fll$forecast[,,5]),single = TRUE)
  le_ll_lva = le_ez(exp(fll$forecast[,,3]),single = TRUE)
  #plot(le_ll[,1],type = 'l')
  
  ##########################################################
  #STAR
  purefore_star <- function(model,tstep,x0){
    tempk = dim(model$fore)[1]
    tempi = dim(model$fore)[3]
    coef = model$A
    mm = model$mu
    predmat = array(0,c(tempk*tempi,tstep))
    predmat[,1] = coef%*% x0 +mm
    
    for (t in 1:(tstep-1)){
      predmat[,(t+1)] = coef%*% predmat[,t] + mm 
    }
    pred = array(predmat,dim = c(tempk,tempi,tstep))
    pred = aperm(pred,c(1,3,2))
    return(pred)
  }
  VAR_1 = fitstar(datar = datar,'stack',fore = TRUE,datate = datate,lambda = 1,kappa = 0)
  lmstar = purefore_star(model = VAR_1,tstep,x0 = matrix(datar[,glot,],ncol = 1))
  # lmstar[,,5] = 0.1*lmstar[,,5]+0.9*lmten[,,5]
  # lmstar[,,3] = 0.05*lmstar[,,3]+0.95*lmten[,,3]
  le_star_rus = le_ez(exp(lmstar[,,5]),single = TRUE)
  le_star_lva = le_ez(exp(lmstar[,,3]),single = TRUE)
  
  ##########################################################

  ##########################################################
  #VAR
  pureforecast_s <- function(model,x0,step) {
    #dim(ddatamat_c) = 49,294
    nc = length(ddatamat_c[nrow(ddatamat_c),])
    nr = step
    m = c(model$mu)
    pred = matrix(0,ncol=nc,nrow=nr)
    pred[1,] = t(model$A[[1]] %*% (x0-m)) +m
    for (i in 2:nr) {
      pred[i,] = t(model$A[[1]] %*% (pred[(i-1),]))+0.985*m
    }
    sum = list()
    sum$forecast = pred
    return(sum)
  }
  VAR_4 = fitVAR(datar,require_diff = TRUE,p=1)#0.0000119 360 1000
  dlm50_s = pureforecast_s(VAR_4,ddatamat_c[nrow(ddatamat_c),],50)
  lm50_s = apply(dlm50_s$forecast,2,cumsum)+matrix(rep(datamat_c[nrow(datamat_c),],tstep),nrow = tstep,byrow = TRUE)
  lmten_s = array(0,dim = c(glok,tstep,gloi))
  for (i in 1:gloi) {
    lmten_s[,,i] = t(lm50_s[,((1:glok)+(i-1)*glok)])
  }
  #le_swvar = le_ez(exp(abind(datar[,,4:5],lmten[,,4:5],along = 2)))
  le_svar_rus = le_ez(exp(lmten_s[,,5]),single = TRUE)
  le_svar_lva = le_ez(exp(lmten_s[,,3]),single = TRUE)
  
  ###########################################################
  #SWVAR
  # bestp = f6$best
  # write.csv(bestp,file ='C:/Users/greatknee/Desktop/Mortality/material_VAR/swvar_param/group2')
  # bestpp = read.csv('C:/Users/greatknee/Desktop/Mortality/material_VAR/swvar_param/group2')
  # bestppp = as.matrix(bestpp[,-1])
  f6 = fitswvar_lack_global(datar = dataall, group =datagroup,datate = datate,coulist = coulist)
  pureforecast_swvar <- function(coef,m,x0,step) {
    #dim(ddatamat_c) = 49,294
    nc = length(ddatamat_c[nrow(ddatamat_c),])
    nr = step
    pred = matrix(0,ncol=nc,nrow=nr)
    pred[1,] = t(coef %*% t(x0-m)) +m
    for (i in 2:nr) {
      pred[i,] = t(coef %*% t(pred[(i-1),]-m))+0.84*m
    }
    sum = list()
    sum$forecast = pred
    return(sum)
  }
  
  dlm50 = pureforecast_swvar(f6$A,VAR_4$mu,ddatamat_c[nrow(ddatamat_c),],50)
  lm50 = apply(dlm50$forecast,2,cumsum)+matrix(rep(datamat_c[nrow(datamat_c),],tstep),nrow = tstep,byrow = TRUE)
  lmten = array(0,dim = c(glok,tstep,gloi))
  for (i in 1:gloi) {
    lmten[,,i] = t(lm50[,((1:glok)+(i-1)*glok)])
  }
  #le_swvar = le_ez(exp(abind(datar[,,4:5],lmten[,,4:5],along = 2)))
  le_swvar_rus = le_ez(exp(lmten[,,5]),single = TRUE)
  le_swvar_lva = le_ez(exp(lmten[,,3]),single = TRUE)
  
  
  ##############################################################
  ##############################################################
  ##############################################################
  ##############################################################
  #plot
  LE = as.vector(c(le_lc_rus,le_ll_rus,le_fdm_rus,le_cofdm_rus,le_star_rus,le_svar_rus,le_swvar_rus)) 
  #pop= rep(coulist,each = 50)
  df1 = data.frame(LE)
  modelist = c('Lee-Carter','Li-Lee','Hyndman-Ullah','Co-FDM','STAR','SVAR','SWVAR')
  df1$Model = rep(modelist,each = 50)
  year = rep(seq(2000,2049,1),length(modelist))
  df1$year = year
  
  year0 = seq(ystr,yete,1)
  df2 = data.frame(year0)
  e0_true = le_ez(exp(dataall[,,5]),single = TRUE)
  df2$eold = as.vector(e0_true)
  
  pfore1 <- ggplot()
  
  p1 = pfore1 +  geom_line(data = df2, aes(x=year0, y = eold),linewidth = 1, colour = 'black')+
    geom_line(data = df1, aes(x=year, y = LE, colour = Model),linewidth = 0.6)+ 
    ggtitle('RUS')+xlab('Year')+ylab('Life expectancy')+theme(legend.position="none")#theme(plot.title = element_blank())
  
  ######################################################################
  LE = as.vector(c(le_lc_lva,le_ll_lva,le_fdm_lva,le_cofdm_lva,le_star_lva,le_svar_lva,le_swvar_lva)) 
  df1 = data.frame(LE)
  df1$Model = rep(modelist,each = 50)
  year = rep(seq(2000,2049,1),length(modelist))
  df1$year = year
  
  year0 = seq(ystr,yete,1)
  df2 = data.frame(year0)
  e0_true = le_ez(exp(dataall[,,3]),single = TRUE)
  df2$eold = as.vector(e0_true)
  
  pfore1 <- ggplot()
  
  p2 = pfore1 +  geom_line(data = df2, aes(x=year0, y = eold),linewidth = 1, colour = 'black')+
    geom_line(data = df1, aes(x=year, y = LE, colour = Model),linewidth = 0.6)+ 
    ggtitle('LVA')+xlab('Year')+ylab('Life expectancy')
  return(list(p1 = p1, p2 =p2))
}
  
