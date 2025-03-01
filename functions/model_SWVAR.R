fitVAR_weighted <- function(data,dataplus,weight, p = 1, penalty = "ENET", method = "cv",age_weight = FALSE,spatial_weight = FALSE, ...) {
  opt <- list(...)
  
  # convert data to matrix
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  if (!is.matrix(dataplus)) {
    dataplus <- as.matrix(dataplus)
  }
  opt$dataplus = dataplus
  opt$dp = NULL
  colnames(data) = as.character(1:dim(data)[2])
  cnames <- cbind(t(colnames(data)) ,t(colnames(dataplus)))
  if (!is.null(weight)){
    pf = weight
    opt$pf = pf
    opt$spatial_weight = TRUE
  }
  if (method == "cv") {
    
    # use CV to find lambda
    opt$method <- "cv"
    out <- cvVAR_aug_w(data, p, penalty, opt)
  } else {
    # error: unknown method
    stop("Unknown method. Possible values are \"cv\" or \"timeSlice\"")
  }
  
  # Add the names of the variables to the matrices
  if (!is.null(cnames)) {
    for (k in 1:length(out$fitmat$A)) {
      colnames(out$fitmat$A[[k]]) <- cnames[1:dim(data)[2]]
      rownames(out$fitmat$A[[k]]) <- cnames[1:dim(data)[2]]
    }
    for (k in 1:length(out$fitmat$B)) {
      colnames(out$fitmat$B[[k]]) <- cnames[(1:dim(data)[2])+(k*dim(data)[2])]
      rownames(out$fitmat$B[[k]]) <- cnames[1:dim(data)[2]]
    }
    
  }
  return(out)
}

cvVAR_aug_w <- function(data, p, penalty = "ENET", opt = NULL) {
  nc <- ncol(data)
  nr <- nrow(data)
  
  picasso <- ifelse(!is.null(opt$picasso), opt$picasso, FALSE)
  threshold <- ifelse(!is.null(opt$threshold), opt$threshold, FALSE)
  
  threshold_type <- ifelse(!is.null(opt$threshold_type),
                           opt$threshold_type, "soft"
  )
  
  return_fit <- ifelse(!is.null(opt$return_fit), opt$return_fit, FALSE)
  
  
  # transform the dataset
  tr_dt <- transformData(data, p, opt)
  
  if (penalty == "ENET") {
    
    # fit the ENET model
    t <- Sys.time()
    fit <- cvVAR_ENET_w(tr_dt$X, tr_dt$y, nvar = nc, opt)
    elapsed <- Sys.time() - t
    
    # extract what is needed
    lambda <- ifelse(is.null(opt$lambda), "lambda.min", opt$lambda)
    
    # extract the coefficients and reshape the matrix
    Avector <- stats::coef(fit, s = lambda)
    A <- matrix(Avector[2:length(Avector)],
                nrow = nc,
                byrow = TRUE
    )
    
    
    mse <- min(fit$cvm)
  }else {
    
    # Unknown penalty error
    stop("Unkown penalty. Available penalties are: ENET, SCAD, MCP.")
  }
  
  # Get back the list of VAR matrices (of length p)
  fitmat <- splitMatrix(A, p, aug=1)
  
  # Now that we have the matrices compute the residuals
  res <- computeResiduals_aug(tr_dt$series, tr_dt$spatialserials,fitmat$A,fitmat$B,m = tr_dt$mu,mp = tr_dt$muplus)
  
  # To extract the sd of mse
  if (penalty == "ENET") {
    ix <- which(fit$cvm == min(fit$cvm))
    mse_sd <- fit$cvsd[ix]
  } else {
    ix <- which(fit$cve == min(fit$cve))
    mse_sd <- fit$cvse[ix]
  }
  
  # Create the output
  output <- list()
  output$fitmat = fitmat
  output$mu <- tr_dt$mu
  output$mup <- tr_dt$muplus
  output$A <- A
  
  # Do you want the fit?
  if (return_fit == TRUE) {
    output$fit <- fit
  }
  
  # Return the "best" lambda
  output$lambda <- fit$lambda.min
  output$pred <- tr_dt$series+res
  output$mse <- mse
  output$mse_sd <- mse_sd
  output$time <- elapsed
  output$series <- tr_dt$series
  output$residuals <- res
  
  # Variance/Covariance estimation
  output$sigma <- estimateCovariance(res)
  
  output$penalty <- penalty
  output$method <- "cv"
  
  #addition
  index = which(fit$lambda %in% fit$lambda.min)
  tLL <- fit$glmnet.fit$nulldev *  fit$glmnet.fit$dev.ratio[index]
  k <- fit$glmnet.fit$df[index]
  nn <- fit$glmnet.fit$nobs
  AIC <- tLL-2*k 
  BIC<-tLL-log(nn)*k
  output$logL = tLL
  output$AIC = AIC
  output$BIC = BIC
  output$coef = Avector
  output$Deviance = (1-fit$glmnet.fit$dev.ratio[index])*fit$glmnet.fit$nulldev
  output$pf = opt$pf
  
  attr(output, "class") <- "var"
  attr(output, "type") <- "fit"
  return(output)
}

cvVAR_ENET_w <- function(X, y, nvar, opt) {
  a <- ifelse(is.null(opt$alpha), 1, opt$alpha)
  nl <- ifelse(is.null(opt$nlambda), 100, opt$nlambda)
  tm <- ifelse(is.null(opt$type.measure), "mse", opt$type.measure)
  nf <- ifelse(is.null(opt$nfolds), 10, opt$nfolds)
  parall <- ifelse(is.null(opt$parallel), FALSE, opt$parallel)
  ncores <- ifelse(is.null(opt$ncores), 1, opt$ncores)
  
  
  
  # Vector of lambdas to work on
  if (!is.null(opt$lambdas_list)) {
    lambdas_list <- opt$lambdas_list
  } else {
    lambdas_list <- c(0)
  }
  
  # Assign ids to the CV-folds (useful for replication of results)
  if (is.null(opt$folds_ids)) {
    folds_ids <- numeric(0)
  } else {
    nr <- nrow(X)
    folds_ids <- rep(sort(rep(seq(nf), length.out = nr / nvar)), nvar)
  }
  
  
  if (length(folds_ids) == 0) {
    if (length(lambdas_list) < 2) {
      cvfit <- glmnet::cv.glmnet(X, y,
                                 alpha = a, nlambda = nl,
                                 type.measure = tm, nfolds = nf,
                                 parallel = FALSE, standardize = FALSE , penalty.factor = opt$pf,intercept=FALSE
      )
    } else {
      cvfit <- glmnet::cv.glmnet(X, y,
                                 alpha = a, lambda = lambdas_list,
                                 type.measure = tm, nfolds = nf,
                                 parallel = FALSE, standardize = FALSE, penalty.factor = opt$pf,intercept=FALSE
      )
    }
  } else {
    if (length(lambdas_list) < 2) {
      cvfit <- glmnet::cv.glmnet(X, y,
                                 alpha = a, nlambda = nl,
                                 type.measure = tm, foldid = folds_ids,
                                 parallel = FALSE, standardize = FALSE, penalty.factor = opt$pf,intercept=FALSE
      )
    } else {
      cvfit <- glmnet::cv.glmnet(X, y,
                                 alpha = a, lambda = lambdas_list,
                                 type.measure = tm, foldid = folds_ids,
                                 parallel = FALSE, standardize = FALSE, penalty.factor = opt$pf,intercept=FALSE
      )
    }
  }
  
  return(cvfit)
}

fitswvar <- function(datar, datate, group, best = NULL,sg = FASLE){
  glok = dim(datar)[1]
  glot = dim(datar)[2]
  glote = dim(datate)[2]
  gloi = dim(datar)[3]
  if (group == 1 | group == 5) {
    ws_c = as.matrix(read.csv('data//centroidd.csv',header = T))
  }else if(group == 2 | group == 3) {
    ws_c = as.matrix(read.csv('data//centroidd2.csv',header = T))
  }
  ws_c = ws_c[match(coulist,colnames(ws_c)),match(coulist,colnames(ws_c))]
  sigseq = matrix(0,glok,glok*gloi)
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
  ddatamat_c = diff(datamat_c)
  ddatamat_te = diff(datamat_te)
  for (i in 1:gloi) {
    sigseq[,((1:glok)+glok*(i-1))] = var(diff(t(datar[,,i])))
  }
  wws <- function(Y,S) {
    k =t(Y) %*% ginv(S)%*%Y
    return(k)
  }
  templist = 1:gloi
  besta = matrix(0,gloi,1)
  bests = matrix(0,gloi,1)
  loopresidual = matrix(0,glot-2,gloi*glok)
  loopcoef = matrix(0,gloi*glok,gloi*glok)
  loopedf = matrix(0,gloi,1)
  predmse = matrix(0,gloi,1)
  loopm = matrix(0,gloi*glok,1)
  loopforer = matrix(0,glote-1,gloi*glok)
  #grid search
  if (is.null(best)) {
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
      
      sig = sigseq[,((1:glok)+(i-1)*glok)]
      aseq = 1.5#10c(2,5,10,20,50,100)
      sseq = 1#c(1,1.5,2,3,5)
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
          wa[upper.tri(wa)] = max(wa)
          scaledws = c(min(log(tempws)),sseq[jj]*log(tempws))/sum(log(tempws))
          result_list = list()
          for (k in 1:gloi) {
            result_list[[k]] <-  scaledws[k] *wa
          }
          ws = do.call(cbind, result_list)
          ws = as.vector(t(ws))
          if(sg){
            VAR_6 = fitVAR_weighted_sgl(dtsdata,dspamat,weight = ws,p=1)#加了alpha = 0.9
          }else{
            VAR_6 = fitVAR_weighted(dtsdata,dspamat,weight = ws,p=1)#加了alpha = 0.9
          }
          edf = sum(!VAR$coef ==0)
          dfy = VAR$residuals[-1,]
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
    besta4fgroup_x = aseq[besta]
    bests4fgroup_x = sseq[bests]
  }else{
    besta4fgroup_x = best[,2]
    bests4fgroup_x = best[,3]
  }
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
        wa[i3,j3] = exp(abs(i3-j3)/besta4fgroup_x[i])
      }
    }
    result_list = list()
    scaledws = c(min(log(tempws)),bests4fgroup_x[i]*log(tempws))/sum(log(tempws))
    for (k in 1:gloi) {
      result_list[[k]] <-  scaledws[k] *wa
    }
    ws = do.call(cbind, result_list)
    ws = as.vector(t(ws))
    if(sg){
      VAR_6 = fitVAR_weighted_sgl(dtsdata,dspamat,weight = ws,p=1)#加了alpha = 0.9
    }else{
      VAR_6 = fitVAR_weighted(dtsdata,dspamat,weight = ws,p=1)#加了alpha = 0.9
    }
    loopresidual[,(1:glok)+(i-1)*glok] = VAR_6$residuals[-1,]
    loopedf[i] = sum(!VAR_6$coef ==0)
    loopm[(1:glok)+(i-1)*glok] = VAR_6$mu
    loopcoef[(1:glok)+(i-1)*glok,] = matrix(VAR_6$coef[-1],nrow = glok,ncol = glok*gloi,byrow = TRUE)
  }
  
  #prepare for coefmat
  coefmat = loopcoef
  for (i in 2:(gloi)) {
    
    coefmat[((1:glok)+(i-1)*glok),(1:((i-1)*glok))] = loopcoef[((1:glok)+(i-1)*glok),((glok+1):(i*glok))]
    coefmat[((1:glok)+(i-1)*glok),(((i-1)*glok+1):(i*glok))] = loopcoef[((1:glok)+(i-1)*glok),(1:glok)]
  }
  
  f6=forecast_our2_nod(coefmat,mean = loopm,ddatamat_c[nrow(ddatamat_c),],datamat_te[1,],datamat_te)
  output = f6
  output$A = coefmat
  output$mu = loopm
  output$best = rbind(besta4fgroup_x, bests4fgroup_x)
  return(output)
}

fitswvar_lack_global <- function(datar, datate, group, best = NULL, coulist,gd = FALSE){
  set.seed(123)
  glok = dim(datar)[1]
  glot = dim(datar)[2]
  glote = dim(datate)[2]
  gloi = dim(datar)[3]
  coulist = coulist
  if (group == 1 | group == 5) {
    ws_c = as.matrix(read.csv('data//centroidd.csv',header = T))
  }else if(group == 2 | group == 3) {
    ws_c = as.matrix(read.csv('data//centroidd2.csv',header = T))
  }
  ws_c = ws_c[match(coulist,colnames(ws_c)),match(coulist,colnames(ws_c))]
  sigseq = matrix(0,glok,glok*gloi)
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
  ddatamat_c = diff(datamat_c)
  ddatamat_te = diff(datamat_te)
  for (i in 1:gloi) {
    sigseq[,((1:glok)+glok*(i-1))] = var(diff(t(datar[,,i])))
  }
  wws <- function(Y,S) {
    k =t(Y) %*% ginv(S)%*%Y
    return(k)
  }
  templist = 1:gloi
  besta = matrix(0,gloi,1)
  bests = matrix(0,gloi,1)
  loopresidual = matrix(0,glot-2,gloi*glok)
  loopcoef = matrix(0,gloi*glok,gloi*glok)
  loopedf = matrix(0,gloi,1)
  predmse = matrix(0,gloi,1)
  loopm = matrix(0,gloi*glok,1)
  loopforer = matrix(0,glote-1,gloi*glok)
  #grid search
  if (is.null(best)) {
    if (isTRUE(gd)) {
      aseq = c(1.5,2,5,10,20,50,100)
      sseq = c(1,1.5,2,3,5)
    }else{
      aseq = 1.5
      sseq = 1
    }
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
      
      sig = sigseq[,((1:glok)+(i-1)*glok)]
      
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
          wa[upper.tri(wa)] = max(wa)
          scaledws = c(min(log(tempws)),sseq[jj]*log(tempws))/sum(log(tempws))
          result_list = list()
          for (k in 1:gloi) {
            result_list[[k]] <-  scaledws[k] *wa
          }
          ws = do.call(cbind, result_list)
          ws = as.vector(t(ws))
          VAR = fitVAR_weighted(dtsdata,dspamat,weight = ws,p=1)
          edf = sum(!VAR$coef ==0)
          dfy = VAR$residuals[-1,]
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
    besta4fgroup_x = aseq[besta]
    bests4fgroup_x = sseq[bests]
  }else{
    besta4fgroup_x = best[,2]
    bests4fgroup_x = best[,3]
  }
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
        wa[i3,j3] = exp(abs(i3-j3)/besta4fgroup_x[i])
      }
    }
    result_list = list()
    scaledws = c(min(log(tempws)),bests4fgroup_x[i]*log(tempws))/sum(log(tempws))
    for (k in 1:gloi) {
      result_list[[k]] <-  scaledws[k] *wa
    }
    ws = do.call(cbind, result_list)
    ws = as.vector(t(ws))
    VAR_6 = fitVAR_weighted(dtsdata,dspamat,weight = ws,p=1)#加了alpha = 0.9
    loopresidual[,(1:glok)+(i-1)*glok] = VAR_6$residuals[-1,]
    loopedf[i] = sum(!VAR_6$coef ==0)
    loopm[(1:glok)+(i-1)*glok] = VAR_6$mu
    loopcoef[(1:glok)+(i-1)*glok,] = matrix(VAR_6$coef[-1],nrow = glok,ncol = glok*gloi,byrow = TRUE)
  }
  
  #prepare for coefmat
  coefmat = loopcoef
  for (i in 2:(gloi)) {
    
    coefmat[((1:glok)+(i-1)*glok),(1:((i-1)*glok))] = loopcoef[((1:glok)+(i-1)*glok),((glok+1):(i*glok))]
    coefmat[((1:glok)+(i-1)*glok),(((i-1)*glok+1):(i*glok))] = loopcoef[((1:glok)+(i-1)*glok),(1:glok)]
  }
  
  f6=forecast_our2_nod(coefmat,mean = loopm,ddatamat_c[nrow(ddatamat_c),],datamat_te[1,],datamat_te,tempk = glok)
  output = f6
  output$A = coefmat
  output$mu = loopm
  output$best = rbind(besta4fgroup_x, bests4fgroup_x)
  return(output)
}

##################################################
##################################################
#sparsegl
fitVAR_weighted_sgl <- function(data,dataplus,group,weightspatial,weightage, p = 1, penalty = "ENET", method = "cv",age_weight = FALSE,spatial_weight = FALSE, ...) {
  opt <- list(...)
  
  # convert data to matrix
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  if (!is.matrix(dataplus)) {
    dataplus <- as.matrix(dataplus)
  }
  opt$dataplus = dataplus
  opt$dp = NULL
  colnames(data) = as.character(1:dim(data)[2])
  cnames <- cbind(t(colnames(data)) ,t(colnames(dataplus)))
  #prepare for sparsegl
  tempX = ncol(data)
  tempI = 1+(ncol(dataplus)/tempX)
  opt$group = group
  opt$ws = weightspatial
  opt$wa = weightage
  opt$spatial_weight = TRUE
  if (method == "cv") {
    
    # use CV to find lambda
    opt$method <- "cv"
    out <- cvVAR_aug_sgl(data, p, penalty, opt)
  } else {
    # error: unknown method
    stop("Unknown method. Possible values are \"cv\" or \"timeSlice\"")
  }
  
  # Add the names of the variables to the matrices
  if (!is.null(cnames)) {
    for (k in 1:length(out$fitmat$A)) {
      colnames(out$fitmat$A[[k]]) <- cnames[1:dim(data)[2]]
      rownames(out$fitmat$A[[k]]) <- cnames[1:dim(data)[2]]
    }
    for (k in 1:length(out$fitmat$B)) {
      colnames(out$fitmat$B[[k]]) <- cnames[(1:dim(data)[2])+(k*dim(data)[2])]
      rownames(out$fitmat$B[[k]]) <- cnames[1:dim(data)[2]]
    }
    
  }
  return(out)
}

cvVAR_aug_sgl <- function(data, p, penalty = "ENET", opt = NULL) {
  nc <- ncol(data)
  nr <- nrow(data)
  
  picasso <- ifelse(!is.null(opt$picasso), opt$picasso, FALSE)
  threshold <- ifelse(!is.null(opt$threshold), opt$threshold, FALSE)
  
  threshold_type <- ifelse(!is.null(opt$threshold_type),
                           opt$threshold_type, "soft"
  )
  
  return_fit <- ifelse(!is.null(opt$return_fit), opt$return_fit, FALSE)
  
  
  # transform the dataset
  tr_dt <- transformData(data, p, opt)
  
  if (penalty == "ENET") {
    
    # fit the ENET model
    t <- Sys.time()
    fit <- cvVAR_ENET_sgl(tr_dt$X, tr_dt$y, nvar = nc, opt)
    elapsed <- Sys.time() - t
    
    # extract what is needed
    lambda <- ifelse(is.null(opt$lambda), "lambda.min", opt$lambda)
    
    # extract the coefficients and reshape the matrix
    Avector <- stats::coef(fit, s = lambda)
    A <- matrix(Avector[2:length(Avector)],
                nrow = nc,
                byrow = TRUE
    )
    
    
    mse <- min(fit$cvm)
  }else {
    
    # Unknown penalty error
    stop("Unkown penalty. Available penalties are: ENET, SCAD, MCP.")
  }
  
  # Get back the list of VAR matrices (of length p)
  fitmat <- splitMatrix(A, p, aug=1)
  
  # Now that we have the matrices compute the residuals
  res <- computeResiduals_aug(tr_dt$series, tr_dt$spatialserials,fitmat$A,fitmat$B,m = tr_dt$mu,mp = tr_dt$muplus)
  
  # To extract the sd of mse
  if (penalty == "ENET") {
    ix <- which(fit$cvm == min(fit$cvm))
    mse_sd <- fit$cvsd[ix]
  } else {
    ix <- which(fit$cve == min(fit$cve))
    mse_sd <- fit$cvse[ix]
  }
  
  # Create the output
  output <- list()
  output$fitmat = fitmat
  output$mu <- tr_dt$mu
  output$mup <- tr_dt$muplus
  output$A <- A
  
  # Do you want the fit?
  if (return_fit == TRUE) {
    output$fit <- fit
  }
  
  # Return the "best" lambda
  output$lambda <- fit$lambda.min
  
  output$mse <- mse
  output$mse_sd <- mse_sd
  output$time <- elapsed
  output$series <- tr_dt$series
  output$residuals <- res
  
  # Variance/Covariance estimation
  output$sigma <- estimateCovariance(res)
  
  output$penalty <- penalty
  output$method <- "cv"
  
  #addition
  index = which(fit$lambda %in% fit$lambda.min)
  tLL <- fit$sparsegl.fit$nulldev *  fit$sparsegl.fit$dev.ratio[index]
  k <- fit$sparsegl.fit$df[index]
  nn <- fit$sparsegl.fit$nobs
  AIC <- tLL-2*k 
  BIC<-tLL-log(nn)*k
  output$logL = tLL
  output$AIC = AIC
  output$BIC = BIC
  output$coef = Avector
  output$Deviance = (1-fit$sparsegl.fit$dev.ratio[index])*fit$sparsegl.fit$nulldev
  output$pf = opt$pf
  
  attr(output, "class") <- "var"
  attr(output, "type") <- "fit"
  return(output)
}

cvVAR_ENET_sgl <- function(X, y, nvar, opt) {
  a <- ifelse(is.null(opt$alpha), 1, opt$alpha)
  nl <- ifelse(is.null(opt$nlambda), 100, opt$nlambda)
  tm <- ifelse(is.null(opt$type.measure), "mse", opt$type.measure)
  nf <- ifelse(is.null(opt$nfolds), 10, opt$nfolds)
  parall <- ifelse(is.null(opt$parallel), FALSE, opt$parallel)
  ncores <- ifelse(is.null(opt$ncores), 1, opt$ncores)
  
  
  
  # Vector of lambdas to work on
  if (!is.null(opt$lambdas_list)) {
    lambdas_list <- opt$lambdas_list
  } else {
    lambdas_list <- c(0)
  }
  
  # Assign ids to the CV-folds (useful for replication of results)
  if (is.null(opt$folds_ids)) {
    folds_ids <- numeric(0)
  } else {
    nr <- nrow(X)
    folds_ids <- rep(sort(rep(seq(nf), length.out = nr / nvar)), nvar)
  }
  
  
  
  cvfit <- sparsegl::cv.sparsegl(X, y, group = opt$group, nlambda = nl,
                                 nfolds = nf, standardize = FALSE,
                                 pf_group = opt$ws, pf_sparse = opt$wa, intercept=FALSE)
  
  return(cvfit)
}

