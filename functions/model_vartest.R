fitVAR_Jtest2 <- function(data,datap, p = 1, penalty = "ENET", method = "cv", ...) {
  opt <- list(...)
  
  # convert data to matrix
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  opt$dp = datap
  cnames <- colnames(data)
  
  if (method == "cv") {
    
    # use CV to find lambda
    opt$method <- "cv"
    out <- cvVAR_j(data, p, penalty, opt)
  } 
  else if (method == "timeSlice") {
    
    # use timeslice to find lambda
    opt$method <- "timeSlice"
    out <- timeSliceVAR(data, p, penalty, opt)
  } 
  else {
    
    # error: unknown method
    stop("Unknown method. Possible values are \"cv\" or \"timeSlice\"")
  }
  
  # Add the names of the variables to the matrices
  if (!is.null(cnames)) {
    for (k in 1:length(out$A)) {
      colnames(out$A[[k]]) <- cnames
      rownames(out$A[[k]]) <- cnames
    }
  }
  
  return(out)
}

cvVAR_j <- function(data, p, penalty = "ENET", opt = NULL) {
  nc <- ncol(data)
  nr <- nrow(data)
  
  picasso <- ifelse(!is.null(opt$picasso), opt$picasso, FALSE)
  threshold <- ifelse(!is.null(opt$threshold), opt$threshold, FALSE)
  
  threshold_type <- ifelse(!is.null(opt$threshold_type),
                           opt$threshold_type, "soft"
  )
  
  return_fit <- ifelse(!is.null(opt$return_fit), opt$return_fit, TRUE)
  
  if (picasso) {
    stop("picasso available only with timeSlice method.")
  }
  # transform the dataset
  tr_dt <- transformData_j(data, p, opt)
  
  if (penalty == "ENET") {
    
    # fit the ENET model
    t <- Sys.time()
    fit <- cvVAR_ENET(tr_dt$X, tr_dt$y, nvar = nc, opt)
    elapsed <- Sys.time() - t
    
    # extract what is needed
    lambda <- ifelse(is.null(opt$lambda), "lambda.min", opt$lambda)
    
    # extract the coefficients and reshape the matrix
    Avector <- stats::coef(fit, s = lambda)
    A <- matrix(Avector[2:(length(Avector)-1)],
                nrow = nc, ncol = nc * p,
                byrow = TRUE
    )
    alphaa = Avector[length(Avector)]
    
    mse <- min(fit$cvm)
  } else if (penalty == "SCAD") {
    
    # convert from sparse matrix to std matrix (SCAD does not work with sparse
    # matrices)
    tr_dt$X <- as.matrix(tr_dt$X)
    
    # fit the SCAD model
    t <- Sys.time()
    fit <- cvVAR_SCAD(tr_dt$X, tr_dt$y, opt)
    elapsed <- Sys.time() - t
    
    # extract the coefficients and reshape the matrix(每行是对单个age group的估计的coef)
    Avector <- stats::coef(fit, s = "lambda.min")
    A <- matrix(Avector[2:length(Avector)],
                nrow = nc, ncol = nc * p,
                byrow = TRUE
    )
    mse <- min(fit$cve)
  }
  else if (penalty == "MCP") {
    
    # convert from sparse matrix to std matrix (MCP does not work with sparse
    # matrices)
    tr_dt$X <- as.matrix(tr_dt$X)
    
    # fit the MCP model
    t <- Sys.time()
    fit <- cvVAR_SCAD(tr_dt$X, tr_dt$y, opt)
    elapsed <- Sys.time() - t
    
    # extract the coefficients and reshape the matrix
    Avector <- stats::coef(fit, s = "lambda.min")
    A <- matrix(Avector[2:length(Avector)],
                nrow = nc, ncol = nc * p,
                byrow = TRUE
    )
    mse <- min(fit$cve)
  } else {
    
    # Unknown penalty error
    stop("Unkown penalty. Available penalties are: ENET, SCAD, MCP.")
  }
  
  # If threshold = TRUE then set to zero all the entries that are smaller than
  # the threshold
  if (threshold == TRUE) {
    A <- applyThreshold(A, nr, nc, p, type = threshold_type)
  }
  
  # Get back the list of VAR matrices (of length p)
  A <- splitMatrix(A, p)
  
  # Now that we have the matrices compute the residuals
  res <- computeResiduals_j(tr_dt$series,opt$dp, A,alphaa, tr_dt$mu)
  
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
  output$mu <- tr_dt$mu
  output$A <- A
  output$alpha = alphaa
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
  #output$sigma <- estimateCovariance(res)
  
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
  output$datap = opt$dp
  output$X = tr_dt$X
  output$Y = tr_dt$y
  attr(output, "class") <- "var"
  attr(output, "type") <- "fit"
  return(output)
}

cvVAR_ENET <- function(X, y, nvar, opt) {
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
  
  if (parall == TRUE) {
    if (ncores < 1) {
      stop("The number of cores must be > 1")
    } else {
      cl <- doParallel::registerDoParallel(cores = ncores)
      
      if (length(folds_ids) == 0) {
        if (length(lambdas_list) < 2) {
          cvfit <- glmnet::cv.glmnet(X, y,
                                     alpha = a, nlambda = nl,
                                     type.measure = tm, nfolds = nf,
                                     parallel = TRUE, standardize = FALSE,intercept=FALSE
          )
        } else {
          cvfit <- glmnet::cv.glmnet(X, y,
                                     alpha = a, lambda = lambdas_list,
                                     type.measure = tm, nfolds = nf,
                                     parallel = TRUE, standardize = FALSE,intercept=FALSE
          )
        }
      } else {
        if (length(lambdas_list) < 2) {
          cvfit <- glmnet::cv.glmnet(X, y,
                                     alpha = a, nlambda = nl,
                                     type.measure = tm, foldid = folds_ids,
                                     parallel = TRUE, standardize = FALSE,intercept=FALSE
          )
        } else {
          cvfit <- glmnet::cv.glmnet(X, y,
                                     alpha = a, lambda = lambdas_list,
                                     type.measure = tm, foldid = folds_ids,
                                     parallel = TRUE, standardize = FALSE,intercept=FALSE
          )
        }
      }
    }
  } else {
    if (length(folds_ids) == 0) {
      if (length(lambdas_list) < 2) {
        cvfit <- glmnet::cv.glmnet(X, y,
                                   alpha = a, nlambda = nl,
                                   type.measure = tm, nfolds = nf,
                                   parallel = FALSE, standardize = FALSE,intercept=FALSE
        )
      } else {
        cvfit <- glmnet::cv.glmnet(X, y,
                                   alpha = a, lambda = lambdas_list,
                                   type.measure = tm, nfolds = nf,
                                   parallel = FALSE, standardize = FALSE,intercept=FALSE
        )
      }
    } else {
      if (length(lambdas_list) < 2) {
        cvfit <- glmnet::cv.glmnet(X, y,
                                   alpha = a, nlambda = nl,
                                   type.measure = tm, foldid = folds_ids,
                                   parallel = FALSE, standardize = FALSE,intercept=FALSE
        )
      } else {
        cvfit <- glmnet::cv.glmnet(X, y,
                                   alpha = a, lambda = lambdas_list,
                                   type.measure = tm, foldid = folds_ids,
                                   parallel = FALSE, standardize = FALSE,intercept=FALSE
        )
      }
    }
  }
  
  return(cvfit)
}

cvVAR_SCAD <- function(X, y, opt) {
  e <- ifelse(is.null(opt$eps), 0.01, opt$eps)
  nf <- ifelse(is.null(opt$nfolds), 10, opt$nfolds)
  parall <- ifelse(is.null(opt$parallel), FALSE, opt$parallel)
  ncores <- ifelse(is.null(opt$ncores), 1, opt$ncores)
  picasso <- ifelse(is.null(opt$picasso), FALSE, TRUE)
  
  if (!picasso) {
    if (parall == TRUE) {
      if (ncores < 1) {
        stop("The number of cores must be > 1")
      } else {
        cl <- parallel::makeCluster(ncores)
        cvfit <- ncvreg::cv.ncvreg(X, y,
                                   nfolds = nf, penalty = "SCAD",
                                   eps = e, cluster = cl
        )
        parallel::stopCluster(cl)
      }
    } else {
      cvfit <- ncvreg::cv.ncvreg(X, y, nfolds = nf, penalty = "SCAD", eps = e)
    }
  } else {
    cvfit <- picasso::picasso(X, y, method = "scad")
  }
  
  return(cvfit)
}

cvVAR_MCP <- function(X, y, opt) {
  e <- ifelse(is.null(opt$eps), 0.01, opt$eps)
  nf <- ifelse(is.null(opt$nfolds), 10, opt$nfolds)
  parall <- ifelse(is.null(opt$parallel), FALSE, opt$parallel)
  ncores <- ifelse(is.null(opt$ncores), 1, opt$ncores)
  
  if (parall == TRUE) {
    if (ncores < 1) {
      stop("The number of cores must be > 1")
    } else {
      cl <- parallel::makeCluster(ncores)
      cvfit <- ncvreg::cv.ncvreg(X, y,
                                 nfolds = nf, penalty = "MCP",
                                 eps = e, cluster = cl
      )
      parallel::stopCluster(cl)
    }
  } else {
    cvfit <- ncvreg::cv.ncvreg(X, y, nfolds = nf, penalty = "MCP", eps = e)
  }
  
  return(cvfit)
}