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

