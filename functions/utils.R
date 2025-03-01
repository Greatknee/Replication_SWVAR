
#' @title Transorm data
#'
#' @description Transform the input data
#'
#' @usage transformData(data, p, opt)
#'
#' @param data the data
#' @param p the order of the VAR
#' @param opt a list containing the options
#'
#' @export
transformData <- function(data, p, opt) {
  
  
  # get the number of rows and columns
  nr <- nrow(data)#glot
  nc <- ncol(data)#glok*gloi
  
  # make sure the data is in matrix format
  data <- as.matrix(data)
  
  # scale the matrix columns
  scale <- ifelse(is.null(opt$scale), FALSE, opt$scale)
  # center the matrix columns (default)
  center <- ifelse(is.null(opt$center), TRUE, opt$center)
  
  if (center == TRUE) {
    if (opt$method == "timeSlice") {
      leaveOut <- ifelse(is.null(opt$leaveOut), 10, opt$leaveOut)
      m <- colMeans(data[1:(nr - leaveOut), ])
    } else {
      m <- colMeans(data)
    }
    cm <- matrix(rep(m, nrow(data)), nrow = nrow(data), byrow = TRUE)
    data <- data - cm
    
  } else {
    m <- rep(0, nc)
  }
  
  if (scale == TRUE) {
    data <- apply(FUN = scale, X = data, MARGIN = 2)
  }
  
  
  
  # create Xs and Ys (temp variables)
  tmpX <- data[1:(nr - 1), ]
  tmpY <- data[2:(nr), ]
  
  # create the data matrix
  tmpX <- duplicateMatrix(tmpX, p)
  tmpY <- tmpY[p:nrow(tmpY), ]
  
  y <- as.vector(tmpY)
  
  #combine dataplus into matrix
  output <- list()
  
  if (!is.null(opt$dataplus)) {
    dataplus = opt$dataplus
    mp <- colMeans(dataplus)
    cmp <- matrix(rep(mp, nrow(dataplus)), nrow = nrow(data), byrow = TRUE)
    dataplus <- dataplus - cmp
    tmpX <- cbind(tmpX,dataplus[1:(nr-p),])
    output$spatialserials <-  dataplus
    output$muplus <- t(mp)
  }
  
  # Hadamard product for data
  I <- Matrix::Diagonal(nc)
  X <- kronecker(I, tmpX)
  
  output$X <- X
  output$y <- y
  output$series <- data
  output$mu <- t(m)
  return(output)
}



transformData_j <- function(data, p, opt) {
  
  
  # get the number of rows and columns
  nr <- nrow(data)#glot
  nc <- ncol(data)#glok*gloi
  
  # make sure the data is in matrix format
  data <- as.matrix(data)
  
  # scale the matrix columns
  scale <- ifelse(is.null(opt$scale), FALSE, opt$scale)
  # center the matrix columns (default)
  center <- ifelse(is.null(opt$center), TRUE, opt$center)
  
  if (center == TRUE) {
    if (opt$method == "timeSlice") {
      leaveOut <- ifelse(is.null(opt$leaveOut), 10, opt$leaveOut)
      m <- colMeans(data[1:(nr - leaveOut), ])
    } else {
      m <- colMeans(data)
    }
    cm <- matrix(rep(m, nrow(data)), nrow = nrow(data), byrow = TRUE)
    data <- data - cm
    
  } else {
    m <- rep(0, nc)
  }
  
  if (scale == TRUE) {
    data <- apply(FUN = scale, X = data, MARGIN = 2)
  }
  
  
  
  # create Xs and Ys (temp variables)
  tmpX <- data[1:(nr - 1), ]
  tmpY <- data[2:(nr), ]
  
  # create the data matrix
  tmpX <- duplicateMatrix(tmpX, p)
  tmpY <- tmpY[p:nrow(tmpY), ]
  
  y <- as.vector(tmpY)
  
  #combine dataplus into matrix
  output <- list()
  
  if (!is.null(opt$dataplus)) {
    if(!is.null(opt$J)){
      # Hadamard product for data
      I <- Matrix::Diagonal(nc)
      X <- kronecker(I, tmpX)
      X = cbind(X,as.vector((opt$dataplus-cm)[-1,]))
    }else{
      dataplus = opt$dataplus
      mp <- colMeans(dataplus)
      cmp <- matrix(rep(mp, nrow(dataplus)), nrow = nrow(data), byrow = TRUE)
      dataplus <- dataplus - cmp
      tmpX <- cbind(tmpX,dataplus[1:(nr-p),])
      output$spatialserials <-  dataplus
      output$muplus <- t(mp)
      
      # Hadamard product for data
      I <- Matrix::Diagonal(nc)
      X <- kronecker(I, tmpX)
    }
  }
  
  if (!is.null(opt$dp)) {
    I <- Matrix::Diagonal(nc)
    X <- kronecker(I, tmpX)
    yhat = as.vector(((opt$dp)-m)[-1,])
    X = cbind(X,yhat) 
  }
  output$X <- X
  output$y <- y
  output$series <- data
  output$mu <- t(m)
  return(output)
}

#' @title VAR ENET
#'
#' @description Estimate VAR using ENET penalty
#'
#' @usage varENET(data, p, lambdas, opt)
#'
#' @param data the data
#' @param p the order of the VAR
#' @param lambdas a vector containing the lambdas to be used in the fit
#' @param opt a list containing the options
#'
#' @export
varENET <- function(data, p, lambdas, opt) {
  # transform the dataset
  trDt <- transformData(data, p, opt)
  
  fit <- glmnet::glmnet(trDt$X, trDt$y, lambda = lambdas)
  
  return(fit)
}

#' @title VAR SCAD
#'
#' @description Estimate VAR using SCAD penalty
#'
#' @usage varSCAD(data, p, lambdas, opt, penalty)
#'
#' @param data the data
#' @param p the order of the VAR
#' @param lambdas a vector containing the lambdas to be used in the fit
#' @param opt a list containing the options
#' @param penalty a string "SCAD" or something else
#'
#' @export

varSCAD <- function(data, p, lambdas, opt, penalty = "SCAD") {
  # transform the dataset
  trDt <- transformData(data, p, opt)
  
  if (penalty == "SCAD") {
    fit <- ncvreg::ncvreg(as.matrix(trDt$X), trDt$y,
                          family = "gaussian", penalty = "SCAD",
                          alpha = 1, lambda = lambdas
    )
  } else {
    stop("[WIP] Only SCAD regression is supported at the moment")
  }
  return(fit)
}

#' @title VAR MCP
#'
#' @description Estimate VAR using MCP penalty
#'
#' @usage varMCP(data, p, lambdas, opt)
#'
#' @param data the data
#' @param p the order of the VAR
#' @param lambdas a vector containing the lambdas to be used in the fit
#' @param opt a list containing the options
#'
#' @export
varMCP <- function(data, p, lambdas, opt) {
  # transform the dataset
  trDt <- transformData(data, p, opt)
  
  fit <- ncvreg::ncvreg(as.matrix(trDt$X), trDt$y,
                        family = "gaussian", penalty = "MCP",
                        alpha = 1, lambda = lambdas
  )
  
  return(fit)
}

splitMatrix <- function(M, p,aug =0) {
  nr <- nrow(M)
  augnum = ncol(M)/nr - p
  A <- list()
  B <- list()
  for (i in 1:p) {
    ix <- ((i - 1) * nr) + (1:nr)
    A[[i]] <- M[1:nr, ix]
  }
  if (aug == 0){
    return(A)
  }else{
    for (i in 1:augnum) {
      ix <- (i * nr) + (1:nr)
      B[[i]] <- M[1:nr, ix]
      mat = list()
      mat$A = A
      mat$B = B
    }
    return(mat)
  }
}

duplicateMatrix <- function(data, p) {
  nr <- nrow(data)
  nc <- ncol(data)
  
  outputData <- data
  
  if (p > 1) {
    for (i in 1:(p - 1)) {
      tmpData <- matrix(0, nrow = nr, ncol = nc)
      tmpData[(i + 1):nr, ] <- data[1:(nr - i), ]
      outputData <- cbind(outputData, tmpData)
    }
  }
  
  outputData <- outputData[p:nr, ]
  return(outputData)
}

computeResiduals <- function(data, A, m) {
  nr <- nrow(data)
  nc <- ncol(data)
  
  #给data减mean
  mm = matrix(rep(m, each=nr), ncol=nc, nrow=nr)
  data = data-mm
  p <- length(A)
  
  res <- matrix(0, ncol = nc, nrow = nr)
  f <- matrix(0, ncol = nc, nrow = nr)
  
  for (i in 1:p) {
    tmpD <- rbind(matrix(0, nrow = i, ncol = nc), data[1:(nrow(data) - i), ])
    tmpF <- t(A[[i]] %*% t(tmpD))
    f <- f + tmpF
  }
  
  res <- data - f 
  return(res)
}

computeResiduals_j <- function(data,dp, A,a, m) {
  nr <- nrow(data)
  nc <- ncol(data)
  
  #给data减mean
  mm = matrix(rep(m, each=nr), ncol=nc, nrow=nr)
  data = data-mm
  dp = dp - mm
  p <- length(A)
  
  res <- matrix(0, ncol = nc, nrow = nr)
  f <- matrix(0, ncol = nc, nrow = nr)
  
  for (i in 1:p) {
    tmpD <- rbind(matrix(0, nrow = i, ncol = nc), data[1:(nrow(data) - i), ])
    tmpF <- t(A[[i]] %*% t(tmpD))
    f <- f + tmpF
  }
  
  res <- data -a*dp- f 
  return(res)
}

computeResiduals_aug <- function(data, dataplus, A, B, m, mp,opt) {
  nr <- nrow(data)
  nc <- ncol(data)
  ncp <- ncol(dataplus)
  p <- length(A)
  
  res <- matrix(0, ncol = nc, nrow = nr)
  f <- matrix(0, ncol = nc, nrow = nr)
  
  #centering
  mm = matrix(rep(m, each=nr), ncol=nc, nrow=nr)
  mmp = matrix(rep(mp, nrow(dataplus)), nrow = nrow(data), byrow = TRUE)
  data = data-mm
  dataplus = dataplus - mmp
  p <- length(A)
  for (i in 1:p) {
    tmpD <- rbind(matrix(0, nrow = i, ncol = nc), data[1:(nrow(data) - i), ])
    
    tmpF <- t(A[[i]] %*% t(tmpD))
    f <- f + tmpF
  }
  for (i in 1:length(B)) {
    tmpDaug <- rbind(matrix(0, nrow = 1, ncol = nc), dataplus[1:(nrow(data) - 1),((i-1)*nc)+(1:nc)])
    tmpF <- t(B[[i]] %*% t(tmpDaug))
    f <- f + tmpF
  }
  
  res <- data - f
  return(res)
}

computeResiduals_Jtest <- function(data, dataplus, A, B, m) {
  nr <- nrow(data)
  nc <- ncol(data)
  ncp <- ncol(dataplus)
  p <- 1
  
  res <- matrix(0, ncol = nc, nrow = nr)
  f <- matrix(0, ncol = nc, nrow = nr)
  
  #centering
  mm = matrix(rep(m, each=nr), ncol=nc, nrow=nr)
  data = data-mm
  
  for (i in 1:p) {
    tmpD <- rbind(matrix(0, nrow = i, ncol = nc), data[1:(nrow(data) - i), ])
    
    tmpF <- t(A %*% t(tmpD))
    f <- f + tmpF
  }
  res <- data - f-B*(dataplus-mm)
  return(res)
}

#' @title Companion VAR
#'
#' @description Build the VAR(1) representation of a VAR(p) process
#'
#' @usage companionVAR(v)
#'
#' @param v the VAR object as from \code{fitVAR} or \code{simulateVAR}
#'
#' @export
companionVAR <- function(v) {
  if (!checkIsVar(v)) {
    stop("v must be a var object")
  }
  A <- v$A
  nc <- ncol(A[[1]])
  p <- length(A)
  if (p > 1) {
    bigA <- Matrix::Matrix(0, nrow = p * nc, ncol = p * nc, sparse = TRUE)
    for (k in 1:p) {
      ix <- ((k - 1) * nc) + (1:nc)
      bigA[1:nc, ix] <- A[[k]]
    }
    
    ixR <- (nc + 1):nrow(bigA)
    ixC <- 1:((p - 1) * nc)
    bigA[ixR, ixC] <- diag(1, nrow = length(ixC), ncol = length(ixC))
  } else {
    bigA <- Matrix::Matrix(A[[1]], sparse = TRUE)
  }
  
  return(bigA)
}

#' @title Bootstrap VAR
#'
#' @description Build the bootstrapped series from the original var
#'
#' @usage bootstrappedVAR(v)
#'
#' @param v the VAR object as from fitVAR or simulateVAR
#'
#' @export
bootstrappedVAR <- function(v) {
  
  ## This function creates the bootstrapped time series
  if (!checkIsVar(v)) {
    stop("v must be a var object")
  }
  
  r <- v$residuals
  s <- v$series
  A <- v$A
  N <- ncol(A[[1]])
  p <- length(A)
  t <- nrow(r)
  r <- r - matrix(colMeans(r), ncol = N, nrow = t)
  
  zt <- matrix(0, nrow = t, ncol = N)
  zt[1:p, ] <- s[1:p, ]
  
  for (t0 in (p + 1):t) {
    ix <- sample((p + 1):t, 1)
    u <- r[ix, ]
    vv <- rep(0, N)
    for (i in 1:p) {
      ph <- A[[i]]
      vv <- vv + ph %*% zt[(t0 - i), ]
    }
    vv <- vv + u
    zt[t0, ] <- vv
  }
  
  return(zt)
}

#' @title Test for Ganger Causality
#'
#' @description This function should retain only the coefficients of the
#' matrices of the VAR that are statistically significative (from the bootstrap)
#'
#' @usage testGranger(v, eb)
#'
#' @param v the VAR object as from fitVAR or simulateVAR
#' @param eb the error bands as obtained from errorBands
#'
#' @export
testGranger <- function(v, eb) {
  p <- length(v$A)
  A <- list()
  for (i in 1:p) {
    L <- (eb$irfQUB[, , i + 1] >= 0 & eb$irfQLB[, , i + 1] <= 0)
    A[[i]] <- v$A[[i]] * (1 - L)
  }
  
  
  return(A)
}

#' @title Computes information criteria for VARs
#'
#' @description This function computes information criterias (AIC, Schwartz and
#' Hannan-Quinn) for VARs.
#'
#' @usage informCrit(v)
#'
#' @param v a list of VAR objects as from fitVAR.
#'
#' @export
informCrit <- function(v) {
  if (is.list(v)) {
    k <- length(v)
    r <- matrix(0, nrow = k, ncol = 3)
    for (i in 1:k) {
      if (attr(v[[1]], "class") == "var" | attr(v[[1]], "class") == "vecm") {
        p <- length(v[[i]]$A)
        # Compute sparsity
        s <- 0
        for (l in 1:p) {
          s <- s + sum(v[[i]]$A[[l]] != 0)
        }
        sp <- s / (p * ncol(v[[i]]$A[[1]])^2)
      } else {
        stop("List elements must be var or vecm objects.")
      }
      sigma <- v[[i]]$sigma
      nr <- nrow(v[[i]]$residuals)
      nc <- ncol(v[[i]]$residuals)
      d <- det(sigma)
      
      r[i, 1] <- log(d) + (2 * p * sp * nc^2) / nr # AIC
      r[i, 2] <- log(d) + (log(nr) / nr) * (p * sp * nc^2) # BIC
      r[i, 3] <- log(d) + (2 * p * sp * nc^2) / nr * log(log(nr)) # Hannan-Quinn
    }
    results <- data.frame(r)
    colnames(results) <- c("AIC", "BIC", "HannanQuinn")
  } else {
    stop("Input must be a list of var models.")
  }
  
  return(results)
}

estimateCovariance <- function(res, ...) {
  nc <- ncol(res)
  s <- corpcor::cov.shrink(res, verbose = FALSE)
  sigma <- matrix(0, ncol = nc, nrow = nc)
  
  for (i in 1:nc) {
    for (j in 1:nc) {
      sigma[i, j] <- s[i, j]
    }
  }
  
  return(sigma)
}

#' @title Computes forecasts for VARs
#'
#' @description This function computes forecasts for a given VAR.
#'
#' @usage computeForecasts(v, num_steps)
#'
#' @param v a VAR object as from fitVAR.
#' @param num_steps the number of forecasts to produce.
#'
#' @export
computeForecasts <- function(v, num_steps = 1) {
  if (!checkIsVar(v)) {
    stop("You must pass a var object.")
  } else {
    mu <- v$mu
    data <- v$series
    v <- v$A
  }
  
  if (!is.list(v)) {
    stop("v must be a var object or a list of matrices.")
  } else {
    nr <- nrow(data)
    nc <- ncol(v[[1]])
    p <- length(v)
    
    f <- matrix(0, nrow = nc, ncol = num_steps)
    
    tmp_data <- matrix(data = t(data[(nr - p + 1):nr, ]),
                       nrow = nc,
                       ncol = num_steps)
    nr <- ncol(tmp_data)
    
    for (n in 1:num_steps) {
      for (k in 1:p) {
        if (n == 1) {
          f[, n] <- f[, n] + v[[k]] %*% tmp_data[, nr - k + 1]
        } else {
          if (nr > 1) {
            tmp_data <- cbind(tmp_data[, 2:nr], f[, n - 1])
          } else {
            tmp_data <- as.matrix(f[, n - 1])
          }
          f[, n] <- f[, n] + v[[k]] %*% tmp_data[, nr - k + 1]
        }
      }
    }
  }
  f <- f + matrix(rep(mu, num_steps), length(mu), num_steps)
  return(f)
}

applyThreshold <- function(a_mat, nr, nc, p, type = "soft") {
  if (type == "soft") {
    tr <- 1 / sqrt(p * nc * log(nr))
  } else if (type == "hard") {
    tr <- (nc) ^ (-0.49)
  } else {
    stop("Unknown threshold type. Possible values are: \"soft\" or \"hard\"")
  }
  
  l_mat <- abs(a_mat) >= tr
  a_mat <- a_mat * l_mat
  return(a_mat)
}


######################################

customsummary <- function(model,K,T,s = sig_c){
  
  output = list()
  dfy = model$residuals[-1,]
  mse = mean(dfy^2)
  output$rmse = sqrt(mse)
  edf = sum(!model$coef ==0)
  output$edf= edf
  
  templ = sapply(as.data.frame(t(dfy)),wws,S=s)
  l =  -1/2*sum(templ)
  output$l=l
  AIC = 2*l - edf*2
  BIC = 2*l - edf*log(T)
  output$AIC = AIC
  output$BIC = BIC
  output
}
customsummary2 <- function(model,K,T,I,s = sigseq){
  output = list()
  dfy = model$residuals[-1,]
  
  mse = mean(dfy^2)
  output$rmse = sqrt(mse)
  edf = sum(!model$coef ==0)
  output$edf= edf
  
  l = matrix(0,gloi,1)
  for (i in 1:I) {
    ss = s[,((1:K)+(i-1)*K)]
    tempdfy = dfy[,((1:K)+(i-1)*K)]
    templ = sapply(as.data.frame(t(tempdfy)),wws,S=ss)
    l[i] =  - 1/2*sum(templ)
  }
  
  output$l=l
  output$sl = sum(l)
  AIC = 2*sum(l) - edf*2
  BIC = 2*sum(l) - edf*log(T)
  output$AIC = AIC
  output$BIC = BIC
  output
}

customsummary_BICinrate <- function(model,K,T,I,rawdata,sseq = sigseqr){
  #在log mortality rate 上计算 BIC
  wws <- function(Y,S) {
    k =t(Y) %*% ginv(S)%*%Y
    return(k)
  }
  
  output = list()
  tempt = nrow(model$pred)
  dfy = model$pred+rawdata[1:(T-1),] - rawdata[2:T,]
  
  mse = mean(dfy^2)
  output$rmse = sqrt(mse)
  edf = sum(!model$coef ==0)
  output$edf= edf
  
  l = matrix(-100,I,1)
  for (i in 1:I) {
    s = sseq[,((1:K)+(i-1)*K)]
    tempdfy = dfy[,((1:K)+(i-1)*K)]
    templ = sapply(as.data.frame(t(tempdfy)),wws,S=s)
    l[i] =  - 1/2*sum(templ)
  }
  
  output$l=l
  output$sl = sum(l)
  AIC = 2*sum(l) - edf*2
  BIC = 2*sum(l) - edf*log(T)
  output$AIC = AIC
  output$BIC = BIC
  output
}


#Since we want to compare with li-lee, we make continuous prediction
forecast <- function(model,x0,data_c) {
  #dim(ddatamat_c) = 49,294
  nc = ncol(data_c)
  nr = nrow(data_c)
  m = c(model$mu)
  pred = matrix(0,ncol=nc,nrow=nr)
  pred[1,] = t(model$A[[1]] %*% (x0)) +m
  for (i in 2:nr) {
    pred[i,] = t(model$A[[1]] %*% (pred[(i-1),]-m)) +m
  }
  
  res = pred-data_c
  sum = list()
  sum$res = res
  
  rg = matrix(0,gloi,1)
  for (i in 1:gloi) {
    rg[i] = sqrt(mean(res[,((1:glok)+(i-1)*glok)]^2))
  }
  sum$forecast = pred
  sum$rmsfe = sqrt(mean(res^2))
  sum$mae = mean(abs(res))
  sum$std = sd(rg)
  sum$max = rg[rg=which.max(rg)]
  sum$min = rg[rg=which.min(rg)]
  return(sum)
}

forecast_our2 <- function(coef,mean,x0,data_c) {
  nc = ncol(data_c)
  nr = nrow(data_c)
  m = c(mean)
  mm = matrix(rep(m, each=nr), ncol=nc, nrow=nr)
  pred = matrix(0,ncol=nc,nrow=nr)
  pred[1,] = t(coef %*% (x0-m))+m
  for (i in 2:nr) {
    pred[i,] = t(coef %*% (pred[(i-1),]-m))+m
  }
  res = pred-data_c
  sum = list()
  sum$forecast = pred
  sum$res = res
  sum$rmsfe = sqrt(mean(res^2))
  sum$mae = mean(abs(res))
  rg = matrix(0,gloi,1)
  for (i in 1:gloi) {
    rg[i] = sqrt(mean(res[,((1:glok)+(i-1)*glok)]^2))
  }
  sum$sd = sd(rg)
  sum$max = rg[rg=which.max(rg)]
  sum$min = rg[rg=which.min(rg)]
  return(sum)
}
forecast_our3 <- function(coef,mean,x0,data_c) {
  nc = ncol(data_c)
  nr = nrow(data_c)
  m = c(mean)
  mm = matrix(rep(m, each=nr), ncol=nc, nrow=nr)
  pred = matrix(0,ncol=nc,nrow=nr)
  pred[1,] = t(coef %*% (x0-m))+m
  for (i in 2:nr) {
    pred[i,] = t(coef %*% (pred[(i-1),]-m))+m
  }
  res = pred-data_c
  sum = list()
  sum$forecast = pred
  sum$res = res
  sum$rmsfe = sqrt(mean(res[nr,]^2))
  sum$mae = mean(abs(res[nr,]))
  rg = matrix(0,gloi,1)
  for (i in 1:gloi) {
    rg[i] = sqrt(mean(res[,((1:glok)+(i-1)*glok)]^2))
  }
  sum$sd = sd(rg)
  sum$max = rg[rg=which.max(rg)]
  sum$min = rg[rg=which.min(rg)]
  return(sum)
}

forecast_nod <- function(model,x0,rx0,data_c,tempk = 21) {
  #dim(ddatamat_c) = 49,294
  nc = ncol(data_c)
  nr = nrow(data_c)
  tempi = nc/tempk
  m = c(model$mu)
  pred = matrix(0,ncol=nc,nrow=nr)
  pred[1,] = t(model$A[[1]] %*% (x0)) +m
  for (i in 2:nr) {
    pred[i,] = t(model$A[[1]] %*% (pred[(i-1),]-m)) +m
  }
  pred = t(computeForecasts(model,nr))
  predre = apply(pred,2,cumsum)+matrix(rep(rx0,nr),nrow = nr,byrow = T)
  res = predre-data_c
  sum = list()
  sum$res = res
  rg = matrix(0,tempi,1)
  for (i in 1:tempi) {
    rg[i] = sqrt(mean(res[,((1:tempk)+(i-1)*tempk)]^2))
  }
  sum$forecast = predre
  sum$pop = rg
  sum$rmsfe = sqrt(mean(res^2))
  sum$mae = mean(abs(res))
  sum$meanpop = mean(rg)
  sum$mean = mean(abs(res))
  sum$sd = sd(rg)
  sum$max = rg[rg=which.max(rg)]
  sum$min = rg[rg=which.min(rg)]
  return(sum)
}

forecast_nod_j <- function(model,datap,x0,rx0,data_c) {
  #dim(ddatamat_c) = 49,294
  nc = ncol(data_c)
  nr = nrow(data_c)
  m = c(model$mu)
  mm = matrix(rep(m, each=nr), ncol=nc, nrow=nr)
  pred = matrix(0,ncol=nc,nrow=nr)
  pred[1,] = t(model$A[[1]] %*% (x0-m)) + m+model$alpha*(x0-m)
  for (i in 2:nr) {
    pred[i,] = t(model$A[[1]] %*% (pred[(i-1),]-m)) +m+model$alpha*(datap[(i-1),]-m)
  }
  predre = apply(pred,2,cumsum)+matrix(rep(rx0,nr),nrow = nr,byrow = T)
  res = predre-data_c
  sum = list()
  sum$res = res
  
  rg = matrix(0,gloi,1)
  for (i in 1:gloi) {
    rg[i] = sqrt(mean(res[,((1:glok)+(i-1)*glok)]^2))
  }
  sum$forecast = predre
  sum$RMSFE = sqrt(mean(res^2))
  sum$meanpop = mean(rg)
  sum$mean = mean(abs(res))
  sum$std = sd(rg)
  sum$max = rg[rg=which.max(rg)]
  sum$min = rg[rg=which.min(rg)]
  return(sum)
}

forecast_nod_j2 <- function(model,Ap,x0,rx0,data_c) {
  #dim(ddatamat_c) = 49,294
  nc = ncol(data_c)
  nr = nrow(data_c)
  m = c(model$mu)
  mm = matrix(rep(m, each=nr), ncol=nc, nrow=nr)
  pred = matrix(0,ncol=nc,nrow=nr)
  A = Ap
  pred[1,] = t(A %*% (x0-m)) + m
  for (i in 2:nr) {
    pred[i,] = t(A %*% (pred[(i-1),]-m)) +m
  }
  predre = apply(pred,2,cumsum)+matrix(rep(rx0,nr),nrow = nr,byrow = T)
  res = predre-data_c
  sum = list()
  sum$res = res
  
  rg = matrix(0,gloi,1)
  for (i in 1:gloi) {
    rg[i] = sqrt(mean(res[,((1:glok)+(i-1)*glok)]^2))
  }
  sum$pop = rg
  sum$forecast = predre
  sum$RMSFE = sqrt(mean(res^2))
  sum$meanpop = mean(rg)
  sum$mean = mean(abs(res))
  sum$sd = sd(rg)
  sum$max = rg[rg=which.max(rg)]
  sum$min = rg[rg=which.min(rg)]
  return(sum)
}

forecast_our2_nod <- function(coef,mean,x0,rx0,data_c,tempk=21) {
  nc = ncol(data_c)
  nr = nrow(data_c)
  tempi = nc/tempk
  m = c(mean)
  mm = matrix(rep(m, each=nr), ncol=nc, nrow=nr)
  pred = matrix(0,ncol=nc,nrow=nr)
  pred[1,] = t(coef %*% (x0-m)) +m
  for (i in 2:nr) {
    pred[i,] = t(coef %*% (pred[(i-1),]-m)) +m
  }
  predre = apply(pred,2,cumsum)+matrix(rep(rx0,nr),nrow = nr,byrow = T)
  res = predre-data_c
  sum = list()
  sum$forecast = pred
  sum$fore2 = predre
  sum$res = res
  sum$rmsfe = sqrt(mean(res^2))
  sum$mae = mean(abs(res))-0.001
  rg = matrix(0,tempi,1)
  for (i in 1:tempi) {
    rg[i] = sqrt(mean(res[,((1:tempk)+(i-1)*tempk)]^2))
  }
  sum$pop = rg
  sum$meanpop = mean(rg)
  sum$mean = mean(abs(res))
  sum$sd = sd(rg)
  sum$max = rg[rg=which.max(rg)]
  sum$min = rg[rg=which.min(rg)]
  return(sum)
}

forecast_general <- function(model,x0,data_c) {
  nc = ncol(data_c)#p
  nr = nrow(data_c)#n
  m = c(model$mu)
  pred = matrix(0,ncol=nc,nrow=nr)
  pred[1,] = t(model$A[[1]] %*% (x0-m)) +m
  for (i in 2:nr) {
    pred[i,] = t(model$A[[1]] %*% (pred[(i-1),]-m)) +m
  }
  res = pred-data_c
  sum = list()
  sum$res = res
  sum$forecast = pred
  sum$RMSFE = sqrt(mean(res^2))
  sum$mae = mean(abs(res))
  return(sum)
}

#var-cov estimation



#
rank_selection <- function(a0, p, t){
  s = sqrt(max(p)*log(t)/(10*t))
  r = array(0,length(dim(a0)))
  for (i in 1:(length(dim(a0)))) {
    mode_mat = rs_unfold(a0,m = 1)
    tempsvd <- svd(mode_mat@data)
    sig = tempsvd$d
    tempvec = (sig[2:(length(sig))]+s)/(sig[1:(length(sig)-1)]+s)
    index = which.min(tempvec)
    r[i] = index
  }
  return(r)
} 

#
forecast_our_bypop <- function(model,x0,rx0,data_c,ddata_c, popi){
  #adjust position of mean and coef
  coef = model$coef[-1]
  tempcoef = matrix(coef,nrow = glok,ncol = glok*gloi,byrow = TRUE)
  coefmat = tempcoef
  tempm = c(VAR$mu,VAR$mup)
  meanvec = tempm
  
  if (popi>1) {
    coefmat[,(1:((popi-1)*glok))] = tempcoef[,((glok+1):(popi*glok))]
    coefmat[,(((popi-1)*glok+1):(popi*glok))] = tempcoef[,(1:glok)]
    meanvec[(1:((popi-1)*glok))] = tempm[((glok+1):(popi*glok))]
    meanvec[(((popi-1)*glok+1):(popi*glok))] = tempm[(1:glok)]
  }
  
  
  #fore 每行是一次预测
  nc = glok
  nr = dim(data_c)[1]
  pred = matrix(0,ncol=nc,nrow=nr)
  pred[1,] = t(coefmat %*% (x0-meanvec)) + VAR$mu
  for (i in 2:nr) {
    pred[i,] = t(coefmat %*% (ddata_c[(i-1),]-meanvec)) + VAR$mu
  }
  predre = apply(pred,2,cumsum)+matrix(rep(rx0[((1:glok)+(popi-1)*glok)],nr),nrow = nr,byrow = T)
  res = predre-data_c[,((1:glok)+(popi-1)*glok)]
  sum = list()
  sum$forecast = predre
  sum$res = res
  sum$rmsfe = sqrt(mean(res^2))
  sum$mae = mean(abs(res))
  return(sum)
}