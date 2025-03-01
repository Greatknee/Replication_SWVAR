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

