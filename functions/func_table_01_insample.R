#prediction
#如果SWVAR报错，加载utilsVAR
#14 low mortality pop
library(psych)
library(sparsevar)
library(vars)
library(MASS) 
library(forecast)



func_table_01_insample <- function(data){
  
  #datapre
  datagroup = data$group
  coulist = data$coulist
  datar = data$datar
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
    sigseq[,((1:glok)+glok*(i-1))] = var(diff(t(datar[,,i])))
  }
  
  sigseqr = matrix(0,glok,glok*gloi)
  for (i in 1:gloi) {
    sigseqr[,((1:glok)+glok*(i-1))] = var(t(datar[,,i]))
  }
  
  #############################################################################
  #Model1,2:li-lee(5 age group/ 1 age group)
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
  
  #1age
  #--------------------li_lee---------------------------
  Global = LeeCarter(data_star_low)#K*T
  lA = Global$coefficients$ax#[age]
  lB = Global$coefficients$bx
  lK = Global$coefficients$kt
  #----------------LCperpopulation---------------------------
  
  #S1,get all lm_{x,i,t}
  dims = dim(datar)
  lista = array(0,dim = c(gloi,glok))
  listb = array(0,dim = c(gloi,glok))
  listk = array(0,dim = c(gloi,glot))
  BK = lB%*%t(lK)
  
  for (i in 1:dims[3]) {
    dataraw = datar[,,i]
    
    ax  <- apply(dataraw, 1, mean)
    A.cmx <- sweep(dataraw, 1, ax, FUN = "-")
    cmx <- A.cmx - BK
    S   <- svd(cmx)
    kt  <- S$v[, 1] * sum(S$u[, 1]) * S$d[1]
    bx  <- S$u[, 1] / sum(S$u[, 1])
    lista[i,] = ax 
    listb[i,] = bx#[pop,age]
    listk[i,] = kt#[pop,year]
    
    #data = sweep(dataraw,1,Global$coefficients$ax, FUN = "-")-BK
    # model = LeeCarter(data)
    # lista[i,] = model$coefficients$ax 
    # listb[i,] = model$coefficients$bx#[pop,age]
    # listk[i,] = model$coefficients$kt#[pop,year]
  }
  #lista = sweep(lista,2,Global$coefficients$ax, FUN = "+")#[pop,age]
  out = list(lista,listb,listk)
  
  #in sample
  predy = array(0,dim = c(glok,glot,gloi))
  for (i in 1:gloi) {
    db = datar[,,i]
    temppred = matrix(0,glok,glot)
    for (t in 1:glot) {
      for (k in 1:glok) {
        temppred[k,t] = lista[i,k]+lB[k]*lK[t] +listb[i,k]*listk[i,t] 
      }
    }
    predy[,,i] = temppred
  }
  
  #number of parameters
  orin = length(lista)*gloi+length(listb)+length(listk)+length(lB)+length(lK)
  orin
  
  dfy = datar-predy
  
  diffpredy = aperm(apply(datar, c(1, 3), diff),c(2,1,3)) - aperm(apply(predy, c(1, 3), diff),c(2,1,3))
  
  
  SSE = sum((dfy^2))
  MSE = mean(dfy^2)
  MSE^(0.5)
  
  l=0
  diffpredy = aperm(apply(datar, c(1, 3), diff),c(2,1,3)) - aperm(apply(predy, c(1, 3), diff),c(2,1,3))
  for (i in 1:gloi) {
    sig = sigseq[,((1:glok)+glok*(i-1))]
    templ = sapply(as.data.frame(diffpredy[,,i]),wws,S=sig)
    l = l  -1/2*sum(templ)
  }
  
  #loopresdidual = matrix(0,glot,gloi*glok)
  # for (i in 1:gloi) {
  #   loopresdidual[,(1:glok)+(i-1)*glok] = t(dfy[,,i])
  # }
  # templ = sapply(as.data.frame(t(loopresdidual[-1,])),wws,S=sig_c)
  # sum(templ)
  # l =-1/2*sum(templ)
  l
  AIC = 2*l-2*orin
  BIC = 2*l-log(glot)*orin
  
  AIC
  BIC
  cu1 = list()
  cu1$RMSE = MSE^(0.5)
  cu1$edf = orin
  cu1$l = l
  cu1$AIC = AIC
  cu1$BIC = BIC
  ##########################################################################
  ##########################################################################
  ##########################################################################
  ##########################################################################
  ##########################################################################
  ##########################################################################
  #co FDM
  fdm = fitfdm(coulist = coulist, year = c(years,yeare),order = 6)
  dfy = fdm$resfit
  
  l=0
  diffpredy = aperm(apply(datar, c(1, 3), diff),c(2,1,3)) - aperm(apply(fdm$pred, c(1, 3), diff),c(2,1,3))
  for (i in 1:gloi) {
    sig = var(t(log(fdm$datar[,,i])))
    templ = sapply(as.data.frame(diffpredy[,,i]),wws,S=sig)
    l = l  -1/2*sum(templ)
  }
  fdm$sl = l
  fdm$sedf =  fdm$edf 
  fdm$AIC =  2*l - fdm$edf*2
  fdm$BIC = 2*l - fdm$edf*log(glot)
  
  
  #STAR sep-model2 stack-model3
  
  VAR_1 = fitstar_revv(datar = datar,'stack')
  VAR_1$summary
  # x <- 1:(glok*gloi)
  # y <- (gloi*glok):1
  # data <- expand.grid(X=x, Y=y)
  # AA = t(VAR_1$A)
  # sum(AA<= -0.1)
  # data$Z <- c(t(VAR_1$A))
  # 
  # # Heatmap 
  # ggplot(data, aes(X, Y, fill= Z)) +   geom_tile() + 
  #   scale_fill_gradient(low="white", high="blue")
  
  VAR_2 = fitstar_revv(datar = datar,'sep')
  VAR_2$summary
  
  
  ##########################################################################
  ##########################################################################
  ##########################################################################
  ##########################################################################
  ##########################################################################
  ##########################################################################
  ##########################################################################
  
  #model3  VAR
  ddatamat_c = diff(datamat_c)
  VAR_3 = fitVAR(ddatamat_c,p=1,lambdas_list=c(0))
  #model7 SVAR
  VAR_4 = fitVAR(ddatamat_c,p=1)#0.0000119 360 1000
  
  #customsummary(VAR_1,K = gloKc,T = gloTc,s = sig_c)
  #cu3 = customsummary2(VAR_3,K = glok,T = glot,gloi,s = sigseq)
  cu3p = customsummary_BICinrate(VAR_3,K = glok,T = glot,gloi,rawdata = datamat_c,s = sigseqr)
  cu3p$sl
  cu4p = customsummary_BICinrate(VAR_4,K = glok,T = glot,gloi,rawdata = datamat_c,s = sigseqr)
  cu4p$BIC
  
  
  
  ###################################################
  #Model6 seperate SVAR
  l =matrix(0,gloi,1)
  edf = matrix(0,gloi,1)
  residual4pop = matrix(0,glot-2,gloi*glok)
  cu5 = list()
  for (i in 1:gloi) {
    data = t(datar[,,i])
    ddata = diff(data)
    VAR6i = fitVAR(ddata,p=1)
    s = sigseqr[,((1:glok)+(i-1)*glok)]
    residual4pop[,(1:glok)+(i-1)*glok] = VAR6i$residuals[-1,]
    tempdfy = VAR6i$pred+data[1:(glot-1),] - data[2:glot,]
    templ = sapply(as.data.frame(t(tempdfy)),wws,S=s)
    l[i] =  - 1/2*sum(templ)
    edf[i] =  sum(!VAR6i$coef ==0)
  }
  sl=sum(l)
  sedf = sum(edf)
  AIC = 2*sl-2*sedf
  BIC = 2*sl - log(glot)*sedf
  AIC
  BIC
  cu5$RMSE = sqrt(mean(residual4pop^2))
  cu5$edf = sedf
  cu5$l = sl
  cu5$AIC = AIC
  cu5$BIC = BIC
  
  
  
  #our method
  #Model7
  # add age weighted model
  # weight age only need to multiply before the X
  #datapre
  if (datagroup == 1 | datagroup == 5) {
    ws_c = as.matrix(read.csv('data//centroidd.csv',header = T))
  }else if(datagroup == 2 | datagroup == 3 | datagroup == 4) {
    ws_c = as.matrix(read.csv('data//centroidd2.csv',header = T))
  }
  ws_c = ws_c[match(coulist,colnames(ws_c)),match(coulist,colnames(ws_c))]
  listl = matrix(0,gloi,1)
  templist = 1:gloi
  besta = matrix(0,gloi,1)
  bests = matrix(0,gloi,1)
  predmse = matrix(0,gloi,1)
  cu7 = list()
  #grid search
  for (i in 1:gloi) {
    tsdata = t(datar[,,i])
    dtsdata = diff(tsdata)
    tempind = templist[templist!=i]
    spamat1 = as.matrix(datar[,,tempind[1]])
    for (ni in 2:length(tempind)) {
      spamat1 = rbind(spamat1,datar[,,tempind[ni]])
    }
    spamat = t(spamat1)
    dspamat = diff(spamat)
    colnames(dspamat) =  paste('aug',seq(1,dim(spamat)[2]),sep = '')
    
    sig = sigseqr[,((1:glok)+(i-1)*glok)]
    aseq = 1.5#seq(1,50,5)
    sseq = 1#seq(1,6,0.5)
    #penalized term 
    tempws = ws_c[i,][ws_c[i,]!=0]
    gridsc = matrix(0,length(aseq),length(sseq))
    gridsAIC = matrix(0,length(aseq),length(sseq))
    gridsl =matrix(0,length(aseq),length(sseq))
    for (ii in 1:length(aseq)) { 
      for (jj in 1:length(sseq)) {
        wa = matrix(0,glok,glok)
        for (i3 in 1:glok) {
          for (j3 in 1:glok) {
            wa[i3,j3] = exp(abs(i3-j3)/aseq[ii])
          }
        }
        scaledws = c(min(tempws),sseq[jj]*tempws)/sum(tempws)
        result_list = list()
        for (k in 1:gloi) {
          result_list[[k]] <-  scaledws[k] *wa
        }
        ws = do.call(cbind, result_list)
        ws = as.vector(t(ws))
        
        VAR = fitVAR_weighted(dtsdata,dspamat,weight = ws,p=1)
        edf = sum(!VAR$coef ==0)
        dfy = VAR$pred + tsdata[1:(glot-1),] - tsdata[2:glot,]
        templ = sapply(as.data.frame(t(dfy)),wws,S=sig)
        l = -1/2*sum(templ)
        AIC = 2*l - edf
        gridsc[ii,jj] = edf
        gridsAIC[ii,jj] = AIC
        gridsl[ii,jj] = l
      }
    }
    max = which(gridsAIC == max(gridsAIC), arr.ind = TRUE)  
    besta[i] = max[1,][1]
    bests[i] = max[1,][2]
  }
  
  
  
  besta
  bests
  bestagroup2 = aseq[besta]
  bestsgroup2 = sseq[bests]
  loopresidual = matrix(0,glot-2,gloi*glok)
  loopcoef = matrix(0,gloi*glok,gloi*glok)
  loopedf = matrix(0,gloi,1)
  looppred = matrix(0,glot-1,gloi*glok)
  #fitting
  for (i in 1:gloi) {
    tsdata = t(datar[,,i])
    dtsdata = diff(tsdata)
    tempind = templist[templist!=i]
    spamat1 = as.matrix(datar[,,tempind[1]])
    for (ni in 2:length(tempind)) {
      spamat1 = rbind(spamat1,datar[,,tempind[ni]])
    }
    spamat = t(spamat1)
    dspamat = diff(spamat)
    colnames(dspamat) =  paste('aug',seq(1,dim(spamat)[2]),sep = '')
    tempws = ws_c[i,][ws_c[i,]!=0]
    
    sig = sigseq[,((1:glok)+(i-1)*glok)]
    wa = matrix(0,glok,glok)
    for (i3 in 1:glok) {
      for (j3 in 1:glok) {
        wa[i3,j3] = exp(abs(i3-j3)/bestagroup2[i])
      }
    }
    
    scaledws = c(min(tempws),bestsgroup2[i]*tempws)/sum(tempws)
    result_list = list()
    for (k in 1:gloi) {
      result_list[[k]] <-  scaledws[k] *wa
    }
    ws = do.call(cbind, result_list)
    ws = as.vector(t(ws))
    VAR_8 = fitVAR_weighted(dtsdata,dspamat,weight = ws,p=1)#加了alpha = 0.9
    VAR_8 = fitVAR_weighted(dtsdata,dspamat,weight = ws,p=1,lambdas_list= 1.2*VAR_8$lambda)#加了alpha = 0.9
    looppred[,(1:glok)+(i-1)*glok] = VAR_8$pred
    loopresidual[,(1:glok)+(i-1)*glok] = VAR_8$residuals[-1,]
    loopedf[i] = sum(!VAR_8$coef ==0)
    loopcoef[(1:glok)+(i-1)*glok,] = matrix(VAR_8$coef[-1],nrow = glok,ncol = glok*gloi,byrow = TRUE)
    
  }
  
  #prepare for coefmat
  coefmat = loopcoef
  for (i in 2:(gloi)) {
    
    coefmat[((1:glok)+(i-1)*glok),(1:((i-1)*glok))] = loopcoef[((1:glok)+(i-1)*glok),((glok+1):(i*glok))]
    coefmat[((1:glok)+(i-1)*glok),(((i-1)*glok+1):(i*glok))] = loopcoef[((1:glok)+(i-1)*glok),(1:glok)]
  }
  coefmat
  
  
  
  #in sample
  besta
  bests
  sqrt(mean(loopresidual^2))
  
  
  #separate loss
  dfy = loopresidual
  
  mse = mean(dfy^2)
  sqrt(mse)
  edf = sum(loopedf)
  edf
  
  l4 = matrix(0,gloi,1)
  for (i in 1:gloi) {
    s = sigseqr[,((1:glok)+(i-1)*glok)]
    tempdfy = looppred[,((1:glok)+(i-1)*glok)] + t(datar[,1:(glot-1),i]) -  t(datar[,2:glot,i])
    templ = sapply(as.data.frame(t(tempdfy)),wws,S=s)
    l4[i] =  - 1/2*sum(templ)
  }
  l4
  sum(l4)
  sedf = sum(loopedf)
  AIC = 2*sum(l4) - sedf*2
  BIC = 2*sum(l4) - sedf*log(glot)
  AIC
  BIC
  cu7$RMSE = sqrt(mse)
  cu7$edf = sedf
  cu7$l = sum(l4)
  cu7$AIC = AIC
  cu7$BIC = BIC
  
  
  conclu = matrix(0,5,8)
  conclu[,1] = c(cu1$RMSE,cu1$edf,cu1$l,cu1$AIC,cu1$BIC)
  conclu[,2] = c(fdm$rmse,fdm$edf,fdm$sl,fdm$AIC,fdm$BIC)
  conclu[,3] = c(VAR_1$summary$rmse,VAR_1$summary$edf,VAR_1$summary$sl,VAR_1$summary$AIC,VAR_1$summary$BIC)
  conclu[,4] = c(VAR_2$summary$rmse,VAR_2$summary$edf,VAR_2$summary$sl,VAR_2$summary$AIC,VAR_2$summary$BIC)
  conclu[,5] = c(cu3p$rmse,cu3p$edf,cu3p$sl,cu3p$AIC,cu3p$BIC)
  conclu[,6] = c(cu4p$rmse,cu4p$edf,cu4p$sl,cu4p$AIC,cu4p$BIC)
  conclu[,7] = c(cu5$RMSE,cu5$edf,cu5$l,cu5$AIC,cu5$BIC)
  conclu[,8] = c(cu7$RMSE,cu7$edf,cu7$l,cu7$AIC,cu7$BIC)
  
  rownames(conclu) = c('rmse','edf','logL','AIC','BIC')
  colnames(conclu) = c('Li-Lee','FDM','STAR1','STAR2','VAR','SVAR1','SVAR2','SWVAR')
  conclu[1,] = round(conclu[1,],4)
  conclu[2,] = round(conclu[2,],0)
  conclu[2:5,] = round(conclu[2:5,],0)
  return(conclu)
}
