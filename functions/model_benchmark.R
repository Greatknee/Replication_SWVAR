#all benchmark models
fitfdm <- function(coulist, fore = FALSE, year,age = 5,gen = 'Male',order = 1){
  
  tempk = 21
  tempi = length(coulist)
  if (gen == 'Male') {
    genum = 2
  }else if (gen =='Female') {
    genum = 1
  }
  if (fore==FALSE) {
    filelist = c("Mx_5x1.txt", "Exposures_5x1.txt")
    ystr = year[1]
    yetr = year[2]
    tempt = yetr-ystr+1
    datar = array(0,dim=c(tempk,tempt,tempi))
    pred = array(0,dim=c(tempk,tempt,tempi))
    for (i in 1:tempi) {
      Mltaddress = paste('data/',coulist[i],'/',filelist,sep = '')
      tempdemdata <- read.demogdata(Mltaddress[1], Mltaddress[2], type="mortality", label=coulist[i])
      tempdemdata$age = c(0,1,seq(5,110,5))
      tempdemdata_mod = extract.years(extract.ages(tempdemdata,0:99,FALSE),ystr:yetr)
      tempfdm = demography::fdm(tempdemdata_mod,series = names(tempdemdata_mod$rate)[2],order = order)
      datar[,,i] = tempdemdata_mod$rate$male
      pred[,,i] = tempfdm$fitted$y
    }
    output = list()
    resfit = log(datar) - pred
    edf = order*tempt*(tempi+1) + (order + order*tempi)*tempk + tempi
    
    # l = matrix(0,tempi,1)
    # for (i in 1:tempi) {
    #   s = var(t(log(datar[,,i])))
    #   tempdfy = resfit[,,i]
    #   templ = sapply(as.data.frame(tempdfy),wws,S=s)
    #   l[i] =  - 0.5*sum(templ)
    # }
    #output$l=
    output$rmse = sqrt(mean(resfit^2))
    output$edf= edf
    #output$S = s
    output$datar = datar
    output$pred = pred
    output$resfit = resfit
    #output$sl = sum(l)
    #AIC = 2*sum(l) - edf*2
    #BIC = 2*sum(l) - edf*log(glot)
    #output$AIC = AIC
    #output$BIC = BIC
    return(output)
  }else if (fore) {
    ###outout rmse
    filelist = c("Mx_5x1.txt", "Exposures_5x1.txt")
    ystr = year[1]
    yetr = year[2]
    tempt = yetr-ystr+1
    datar = array(0,dim=c(tempk,tempt,tempi))
    pred = array(0,dim=c(tempk,tempt,tempi))
    for (i in 1:tempi) {
      Mltaddress = paste('data/',coulist[i],'/',filelist,sep = '')
      tempdemdata <- read.demogdata(Mltaddress[1], Mltaddress[2], type="mortality", label=coulist[i])
      tempdemdata$age = c(0,1,seq(5,110,5))
      tempdemdata_mod = extract.years(extract.ages(tempdemdata,0:99,FALSE),ystr:yetr)
      tempfdm = demography::fdm(tempdemdata_mod,series = names(tempdemdata_mod$rate)[genum],order = 1)
      datar[,,i] = tempdemdata_mod$rate$male
      pred[,,i] = tempfdm$fitted$y
    }
    res = log(datar) - pred
    rmse = sqrt(mean(res^2))
    mae = mean(abs(res))
    ############
    #forecast
    ystr = year[1]
    yetr = year[2]
    tempt = yetr-ystr+1
    yste = year[3]
    yete = year[4]
    tempte = yete-yste+1
    if (age == 1) {
      tempk = 49
      filelist = c("Deaths_1x1.txt", "Exposures_1x1.txt", "Population.txt", "Mx_1x1.txt")
      Mltaddress = paste('data/',coulist[1],'/',filelist,sep = '')
      
      datate = array(0,dim = c(tempk,tempte,tempi))
      
      tempdatag <- read_hmd_files(Mltaddress)
      df = as.data.frame(tempdatag |> 
                           dplyr::filter(Year >= ystr & Year <= yetr& Sex == gen &Age>=41& Age<=89))
      dfte1 = as.data.frame(tempdatag |> 
                              dplyr::filter(Year >= yste & Year <= yete& Sex == gen &Age>=41& Age<=89))
      datate[,,1] = matrix(dfte1$Mortality,tempk,tempte)
      df$popc = coulist[1]
      
      for (i in 2:tempi) {
        Mltaddress = paste('data/',coulist[i],'/',filelist,sep = '')
        tempdatag <- read_hmd_files(Mltaddress)
        tempdata = as.data.frame(tempdatag |> 
                                   dplyr::filter(Year >= ystr & Year <= yetr& Sex == gen  &Age>=41& Age<=89))
        dfte = as.data.frame(tempdatag |> 
                               dplyr::filter(Year >= yste & Year <= yete& Sex == gen &Age>=41& Age<=89))
        datate[,,i] = matrix(dfte$Mortality,tempk,tempte)
        tempdata$popc = coulist[i]
        df = dplyr::bind_rows(df,tempdata)
        
      }
      df$Sex <- NULL
      names(df)[names(df) == "popc"] <- "Sex"
      
      data = as_vital(as_tibble(df),.age = 'Age',.sex = 'Sex',.population = 'Population',key = c('Age','Sex'),index = 'Year')
      
      pr <- data |>
        make_pr(Mortality)
      hu <- pr |>
        model(hyndman_ullah = FDM(log(Mortality),coherent = TRUE,order = 10))|> 
        vital::forecast(h =tempte)
      huundo <- hu |>
        undo_pr(Mortality,key = 'Sex',times = 100)
      
      hufore = array(huundo$.mean,dim = c(tempk,tempi,tempte))
      fore <- aperm(hufore, perm = c(1,3,2))
      output = list()
      output$rmse = rmse
      output$mae = mae
      output$pred = pred
      output$fore = fore
      output$residual = log(fore)-log(datate)
      sum = list()
      sum$rmsfe = sqrt(mean((log(fore)-log(datate))^2))
      pop = sqrt(apply((log(fore)-log(datate))^2, 3, mean))
      sum$pop = pop
      sum$mafe = mean(abs(log(fore)-log(datate)))
      output$summary = sum
    }
    else{
      tempk = 21
      filelist = c("Deaths_5x1.txt", "Exposures_5x1.txt", "Population5.txt", "Mx_5x1.txt")
      Mltaddress = paste('data/',coulist[1],'/',filelist,sep = '')
      
      datate = array(0,dim = c(tempk,tempte,tempi))
      
      tempdatag <- read_hmd_files(Mltaddress)
      df = as.data.frame(tempdatag |> 
                           dplyr::filter(Year >= ystr & Year <= yetr& Sex == gen & Age<100))
      dfte1 = as.data.frame(tempdatag |> 
                              dplyr::filter(Year >= yste & Year <= yete& Sex == gen & Age<100))
      datate[,,1] = matrix(dfte1$Mx_5x1,tempk,tempte)
      df$popc = coulist[1]
      
      for (i in 2:tempi) {
        Mltaddress = paste('data/',coulist[i],'/',filelist,sep = '')
        tempdatag <- read_hmd_files(Mltaddress)
        tempdata = as.data.frame(tempdatag |> 
                                   dplyr::filter(Year >= ystr & Year <= yetr& Sex == gen & Age<100))
        dfte = as.data.frame(tempdatag |> 
                               dplyr::filter(Year >= yste & Year <= yete& Sex == gen & Age<100))
        datate[,,i] = matrix(dfte$Mx_5x1,tempk,tempte)
        tempdata$popc = coulist[i]
        df = dplyr::bind_rows(df,tempdata)
        
      }
      names(df)[names(df) == "Mx_5x1"] <- "Mortality"
      df$Sex <- NULL
      names(df)[names(df) == "popc"] <- "Sex"
      
      data = as_vital(as_tibble(df),.age = 'Age',.sex = 'Sex',.population = 'Population',key = c('Age','Sex'),index = 'Year')
      
      pr <- data |>
        make_pr(Mortality)
      hu <- pr |>
        model(hyndman_ullah = FDM(log(Mortality),coherent = TRUE))|> 
        vital::forecast(h =tempte)
      huundo <- hu |>
        undo_pr(Mortality,key = 'Sex',times = 200)
      
      hufore = array(huundo$.mean,dim = c(tempk,tempi,tempte))
      fore <- aperm(hufore, perm = c(1,3,2))
      
      output = list()
      output$data = abind(datar,datate,along = 2)
      output$rmse = rmse
      output$mae = mae
      output$pred = pred
      output$fore = fore
      output$residual = log(fore)-log(datate)
      sum = list()
      sum$rmsfe = sqrt(mean((log(fore)-log(datate))^2))
      pop = sqrt(apply((log(fore)-log(datate))^2, c(2,3), mean))
      sum$pop = pop
      sum$mafe = mean(abs(log(fore)-log(datate)))
      output$summary = sum
    }
    # sd(pop)
    # max(pop)
    # min(pop)
    
    return(output)
  }
}


fitfdm_bisex <- function(coulist, fore = FALSE, year,age = 5){
  tempk = 21
  tempi = length(coulist)
  
  if (fore==FALSE) {
    filelist = c("Mx_5x1.txt", "Exposures_5x1.txt")
    ystr = year[1]
    yetr = year[2]
    tempt = yetr-ystr+1
    pred = datar
    datar = array(0,dim=c(tempk,tempt,tempi))
    pred = array(0,dim=c(tempk,tempt,tempi))
    for (i in 1:tempi) {
      Mltaddress = paste('data/',coulist[i],'/',filelist,sep = '')
      tempdemdata <- read.demogdata(Mltaddress[1], Mltaddress[2], type="mortality", label=coulist[i])
      tempdemdata$age = c(0,1,seq(5,110,5))
      tempdemdata_mod = extract.years(extract.ages(tempdemdata,0:99,FALSE),ystr:yetr)
      tempfdm = demography::fdm(tempdemdata_mod,series = names(tempdemdata_mod$rate)[2],order = 2)
      datar[,,i] = tempdemdata_mod$rate$male
      pred[,,i] = tempfdm$fitted$y
    }
    output = list()
    res = log(datar) - pred
    edf = tempk*(tempi+1) + (2 + 2*tempi)*glot
    
    l = matrix(0,tempi,1)
    for (i in 1:tempi) {
      s = var(t(log(datar[,,i])))
      tempdfy = res[,,i]
      templ = sapply(as.data.frame(tempdfy),wws,S=s)
      l[i] =  - 1/2*sum(templ)
    }
    output$l=l
    output$rmse = sqrt(mean(res^2))
    output$edf= edf
    output$sl = sum(l)
    AIC = 2*sum(l) - edf*2
    BIC = 2*sum(l) - edf*log(glot)
    output$AIC = AIC
    output$BIC = BIC
    return(output)
  }else if (fore) {
    ystr = year[1]
    yetr = year[2]
    tempt = yetr-ystr+1
    yste = year[3]
    yete = year[4]
    tempte = yete-yste+1
    if (age == 1) {
      tempk = 49
      filelist = c("Deaths_1x1.txt", "Exposures_1x1.txt", "Population.txt", "Mx_1x1.txt")
      Mltaddress = paste('data/',coulist[1],'/',filelist,sep = '')
      
      datate = array(0,dim = c(tempk,tempte,tempi))
      
      tempdatag <- read_hmd_files(Mltaddress)
      df = as.data.frame(tempdatag |> 
                           dplyr::filter(Year >= ystr & Year <= yetr& Sex == gen &Age>=41& Age<=89))
      dfte1 = as.data.frame(tempdatag |> 
                              dplyr::filter(Year >= yste & Year <= yete& Sex == gen &Age>=41& Age<=89))
      datate[,,1] = matrix(dfte1$Mortality,tempk,tempte)
      df$popc = coulist[1]
      
      for (i in 2:tempi) {
        Mltaddress = paste('data/',coulist[i],'/',filelist,sep = '')
        tempdatag <- read_hmd_files(Mltaddress)
        tempdata = as.data.frame(tempdatag |> 
                                   dplyr::filter(Year >= ystr & Year <= yetr& Sex == gen  &Age>=41& Age<=89))
        dfte = as.data.frame(tempdatag |> 
                               dplyr::filter(Year >= yste & Year <= yete& Sex == gen &Age>=41& Age<=89))
        datate[,,i] = matrix(dfte$Mortality,tempk,tempte)
        tempdata$popc = coulist[i]
        df = dplyr::bind_rows(df,tempdata)
        
      }
      df$Sex <- NULL
      names(df)[names(df) == "popc"] <- "Sex"
      
      data = as_vital(as_tibble(df),.age = 'Age',.sex = 'Sex',.population = 'Population',key = c('Age','Sex'),index = 'Year')
      
      pr <- data |>
        make_pr(Mortality)
      hu <- pr |>
        model(hyndman_ullah = FDM(log(Mortality),coherent = TRUE,order = 1))|> 
        vital::forecast(h =tempte)
      huundo <- hu |>
        undo_pr(Mortality,key = 'Sex',times = 100)
      
      hufore = array(huundo$.mean,dim = c(tempk,tempi,tempte))
      fore <- aperm(hufore, perm = c(1,3,2))
      output = list()
      output$fore = fore
      output$residual = log(fore)-log(datate)
      sum = list()
      sum$rmsfe = sqrt(mean((log(fore)-log(datate))^2))
      pop = sqrt(apply((log(fore)-log(datate))^2, c(2,3), mean))
      sum$pop = pop
      sum$mae = mean(abs(log(fore)-log(datate)))
      output$summary = sum
    }
    else{
      tempk = 21
      filelist = c("Deaths_5x1.txt", "Exposures_5x1.txt", "Population5.txt", "Mx_5x1.txt")
      Mltaddress = paste('data/',coulist[1],'/',filelist,sep = '')
      datate = array(0,dim = c(tempk,tempte,2*tempi))
      tempdatag <- read_hmd_files(Mltaddress)
      df1 = as.data.frame(tempdatag |> 
                            dplyr::filter(Year >= ystr & Year <= yetr& Sex == 'Male' & Age<100))
      df2 = as.data.frame(tempdatag |> 
                            dplyr::filter(Year >= ystr & Year <= yetr& Sex == 'Female' & Age<100))
      dfte1 = as.data.frame(tempdatag |>
                              dplyr::filter(Year >= yste & Year <= yete& Sex == 'Male'& Age<100))
      dfte2 = as.data.frame(tempdatag |>
                              dplyr::filter(Year >= yste & Year <= yete& Sex == 'Female'& Age<100))
      datate[,,1] = matrix(dfte2$Mx_5x1,tempk,tempte) #和后面对应 F在前
      datate[,,2] = matrix(dfte1$Mx_5x1,tempk,tempte)
      df1$popc = paste(coulist[1],'M')
      df2$popc = paste(coulist[1],'F')
      df = dplyr::bind_rows(df1,df2)
      
      for (i in 2:tempi) {
        Mltaddress = paste('data/',coulist[i],'/',filelist,sep = '')
        tempdatag <- read_hmd_files(Mltaddress)
        tempdata1 = as.data.frame(tempdatag |> 
                                    dplyr::filter(Year >= ystr & Year <= yetr& Sex == 'Male' & Age<100))
        tempdata2 = as.data.frame(tempdatag |> 
                                    dplyr::filter(Year >= ystr & Year <= yetr& Sex == 'Female' & Age<100))
        dfte1 = as.data.frame(tempdatag |>
                                dplyr::filter(Year >= yste & Year <= yete& Sex == 'Male'& Age<100))
        dfte2 = as.data.frame(tempdatag |>
                                dplyr::filter(Year >= yste & Year <= yete& Sex == 'Female'& Age<100))
        datate[,,(2*i-1)] = matrix(dfte2$Mx_5x1,tempk,tempte)
        datate[,,(2*i)] = matrix(dfte1$Mx_5x1,tempk,tempte)
        
        tempdata1$popc = paste(coulist[i],'M')
        tempdata2$popc = paste(coulist[i],'F')
        df = dplyr::bind_rows(df,tempdata1)
        df = dplyr::bind_rows(df,tempdata2)
      }
      names(df)[names(df) == "Mx_5x1"] <- "Mortality"
      df$Sex <- NULL
      names(df)[names(df) == "popc"] <- "Sex"
      
      data = as_vital(as_tibble(df),.age = 'Age',.sex = 'Sex',.population = 'Population',key = c('Age','Sex'),index = 'Year')
      
      pr <- data |>
        make_pr(Mortality)
      hu <- pr |>
        model(hyndman_ullah = FDM(log(Mortality),coherent = TRUE,order = 1))|> 
        vital::forecast(h =tempte)
      huundo <- hu |>
        undo_pr(Mortality,key = 'Sex',times = 200)
      
      hufore = array(huundo$.mean,dim = c(tempk,2*tempi,tempte))
      fore <- aperm(hufore, perm = c(1,3,2))
      output = list()
      output$fore = fore
      output$residual = log(fore)-log(datate)
      sum = list()
      sum$rmsfe = sqrt(mean((log(fore)-log(datate))^2))
      pop = sqrt(apply((log(fore)-log(datate))^2, c(2,3), mean))
      sum$pop = pop
      sum$mae = mean(abs(log(fore)-log(datate)))
      output$summary = sum
    }
    # sd(pop)
    # max(pop)
    # min(pop)
    return(output)
  }
}

####################STAR

fitstar = function(datar, method, fore= FALSE, datate = NULL, lambda = 1 , kappa = 0){
  result = list()
  tempk = dim(datar)[1]
  tempt = dim(datar)[2]
  tempi = dim(datar)[3]
  wws <- function(Y,S) {
    k =t(Y) %*% ginv(S)%*%Y
    return(k)
  }
  if (fore == FALSE) {
    pred = array(0,c(tempk,tempt-1,tempi))
    if (method == 'sep') {
      edf = tempi*(tempk+tempk-1+tempk-2)
      #We do not grid search
      lambda_m = lambda
      lambda_a = lambda
      lambda_b = lambda
      cvxl = seq(0,tempi)
      for (II in 1:tempi) {
        X <- datar[,,II]
        Y <- X
        # 定义P
        m <- Variable(tempk)
        A <- Variable(tempk,tempk)
        
        #加常数项
        B = cbind(m,A)
        Y_add = rbind(1,Y)
        # 定义优化目标：最小化残差平方和
        objective <- Minimize(sum_squares(Y[,2:tempt]  - B %*% Y_add[,1:(tempt-1)])+lambda_m*sum_squares(diff(m))
                              +lambda_a*sum_squares(diff(diag(A[2:tempk,1:(tempk-1)])))
                              +lambda_b*sum_squares(diff(diag(A[3:tempk,1:(tempk-2)]))))
        
        #
        zero_constraint = c(A[upper.tri(A,diag = FALSE)]==0)
        row_sum_constraint <- sapply(1:tempk, function(i) sum(A[i, ]) == 1)
        AA = A[4:tempk,1:(tempk-3)]
        row_sum_constraint2 <- c(AA[lower.tri(AA,diag = TRUE)]==0)
        
        nonnegative_constraint1 <- sapply(2:tempk, function(i) A[i,(i-1)] >= 0)
        nonnegative_constraint2 <- sapply(3:tempk, function(i) A[i,(i-2)] >= 0)
        constraints <- c(zero_constraint,row_sum_constraint2,row_sum_constraint,nonnegative_constraint1,nonnegative_constraint2)
        
        # 定义问题
        problem <- Problem(objective,constraints = constraints)
        
        # 求解优化问题
        a = solve(problem)
        coef <- a$getValue(A)
        coef
        cvxl[II] = a$value
        m <- a$getValue(m)
        for (t in 1:(tempt-1)) {
          pred[,t,II] = coef%*% X[,t]+m
        }
      }
      resf = datar[,2:tempt,] - pred
    }
    else if (method == 'stack') {
      edf = tempk*(tempi^2)-tempi
      lamalist = rep(lambda,tempi)
      lammlist=rep(lambda,tempi)
      lamrholist=matrix(lambda,tempi,tempi)
      
      # 创建一个示例的3D数组 X
      X <- datar
      
      # 合并 d1 和 d3 维度并重排顺序
      X_new <- aperm(X, c(1, 3, 2))  # 调整维度顺序，使得顺序为 (d1, d3, d2)
      
      # 将其转换为矩阵，合并 d1 和 d3 维度
      Y <- matrix(X_new, nrow = tempk*tempi, ncol = tempt)
      
      # 定义P
      m <- Variable(tempk*tempi)
      BBmat = Variable(tempk*tempi,tempk*tempi)
      #加常数项
      A = cbind(m,BBmat)
      Y_add = rbind(1,Y)
      #penalized term
      normm = 0
      norma = 0
      normrho = 0
      for (ni in 1:tempi) {
        mk = m[(1:tempk)+(ni-1)*tempk]
        diffm = diff(mk)
        normm = normm + lamalist[ni] * sum_squares(diffm[2:length(diffm)])
      }
      for (ni in 1:tempi) {
        ak = diag(BBmat[((2:tempk)+(ni-1)*tempk),((1:(tempk-1))+(ni-1)*tempk)])
        diffa = diff(ak)
        norma = norma + lamalist[ni] * sum_squares(diffm[2:length(diffa)])
      }
      for (ni in 1:tempi) {
        for (nj in 1:tempi) {
          if (ni!=nj) {
            BBtemp = BBmat[((1:tempk)+(ni-1)*tempk),((1:tempk)+(nj-1)*tempk)]
            rhok = diag(BBtemp)
            diffrho = diff(rhok)
            normrho = normrho+normrho
          }
        }
      }
      # 定义优化目标：最小化残差平方和
      objective <- Minimize(sum_squares(Y[,2:tempt]  - A %*% Y_add[,1:(tempt-1)])+normm+norma+normrho)
      
      #
      constraints <- c()
      for (i in 1:tempi) {
        for (j in 1:tempi) {
          if (i==j) {
            tempBB = BBmat[((1:tempk)+(i-1)*tempk),((1:tempk)+(j-1)*tempk)]
            tempBBlow = tempBB[2:tempk,1:(tempk-1)]
            #row_sum_constraint <- sapply(2:tempk, function(ii) tempBB[ii,ii-1]+tempBB[ii,ii] == 1)
            constraints = c(constraints,
                            diag(tempBB)>=0, tempBB[upper.tri(tempBB)]==0,
                            tempBBlow[lower.tri(tempBBlow)]==0)
          }else{
            tempBB = BBmat[((1:tempk)+(i-1)*tempk),((1:tempk)+(j-1)*tempk)]
            zero_constraint <- c(tempBB[upper.tri(tempBB)]==0,tempBB[lower.tri(tempBB)]==0,diag(tempBB)>=0)
            constraints = c(constraints,zero_constraint)
          }
        }
      }
      row_sum_constraint_p <- sapply(1:(tempi*tempk), function(i) sum(BBmat[i, ]) == 1)
      constraints = c(constraints,row_sum_constraint_p,BBmat<=1)
      
      # 定义问题
      problem <- Problem(objective,constraints = constraints)
      
      
      # 求解优化问题
      a = solve(problem)
      coef = a$getValue(BBmat)
      cvxl = a$value
      mm = a$getValue(m)
      predmat = array(0,c(tempk*tempi,(tempt-1)))
      for (t in 1:(tempt-1)){
        predmat[,t] = coef%*% Y[,t]+mm
      }
      pred = array(predmat,dim = c(tempk,tempi,(tempt-1)))
      pred = aperm(pred,c(1,3,2))
    }
    res = datar[,2:tempt,] - pred
    result$pred = pred
    result$residual = res
    result$status = a$status
    result$A = coef
    #in-sample-summary
    in_summary = list()
    l = matrix(0,tempi,1)
    
    for (i in 1:tempi) {
      tempres = res[,,i]
      s = var(t(datar[,,i]))
      templ = sapply(as.data.frame(tempres),wws,S=s)
      l[i] =  - 1/2*sum(templ)
    }
    in_summary$l = l
    in_summary$rmse = sqrt(mean(res^2))
    in_summary$edf = edf
    in_summary$sl = sum(l)
    in_summary$AIC = 2*sum(l) - 2*edf
    in_summary$BIC = 2*sum(l) - log(tempt)*edf
    result$summary = in_summary
  }
  else if (fore == TRUE& !is.null(datate)) {
    tempte = dim(datate)[2]
    fore = array(0,c(tempk,tempte,tempi))
    coef_A = matrix(0,tempk*tempi,tempk*tempi)
    if (method == 'sep') {
      edf = tempi*(tempk+tempk-1+tempk-2)
      #We do not grid search
      lambda_m = lambda
      lambda_a = lambda
      lambda_b = lambda
      cvxl = seq(0,tempi)
      mu = rep(0,tempi*tempk)
      for (II in 1:tempi) {
        X <- datar[,,II]
        Y <- X
        # 定义P
        m <- Variable(tempk)
        A <- Variable(tempk,tempk)
        
        #加常数项
        B = cbind(m,A)
        Y_add = rbind(1,Y)
        # 定义优化目标：最小化残差平方和
        objective <- Minimize(sum_squares(Y[,2:tempt]  - B %*% Y_add[,1:(tempt-1)])+lambda_m*sum_squares(diff(m))
                              +lambda_a*sum_squares(diff(diag(A[2:tempk,1:(tempk-1)])))
                              +lambda_b*sum_squares(diff(diag(A[3:tempk,1:(tempk-2)]))))
        
        #
        zero_constraint = c(A[upper.tri(A,diag = FALSE)]==0)
        row_sum_constraint <- sapply(1:tempk, function(i) sum(A[i, ]) == 1)
        AA = A[4:tempk,1:(tempk-3)]
        row_sum_constraint2 <- c(AA[lower.tri(AA,diag = TRUE)]==0)
        
        nonnegative_constraint1 <- sapply(2:tempk, function(i) A[i,(i-1)] >= 0)
        nonnegative_constraint2 <- sapply(3:tempk, function(i) A[i,(i-2)] >= 0)
        constraints <- c(zero_constraint,row_sum_constraint2,row_sum_constraint,nonnegative_constraint1,nonnegative_constraint2)
        
        # 定义问题
        problem <- Problem(objective,constraints = constraints)
        
        # 求解优化问题
        a = solve(problem)
        coef <- a$getValue(A)
        mm = a$getValue(m)
        coef[coef<= kappa] = 0
        coef_A[((1:tempk)+(II-1)*tempk),((1:tempk)+(II-1)*tempk)] = coef
        cvxl[II] = a$value
        fore[,1,II] = coef%*% X[,tempt]+mm
        for (t in 1:(tempte-1)) {
          fore[,(t+1),II] = coef%*% fore[,t,II]+mm
        }
        mu[(1:tempk)+(II-1)*tempi] = mm
      }
      result$A = coef_A
      result$mu = mu
    }
    else if (method == 'stack') {
      edf = tempi*(tempk^2)+tempk*tempi
      lamalist = rep(lambda,tempi)
      lammlist=rep(lambda,tempi)
      lamrholist=matrix(lambda,tempi,tempi)
      
      #centra
      # 创建一个示例的3D数组 X
      X <- datar
      
      # 合并 d1 和 d3 维度并重排顺序
      X_new <- aperm(X, c(1, 3, 2))  # 调整维度顺序，使得顺序为 (d1, d3, d2)
      
      # 将其转换为矩阵，合并 d1 和 d3 维度
      Y <- matrix(X_new, nrow = tempk*tempi, ncol = tempt)
      # 定义P
      BBmat = Variable(tempk*tempi,tempk*tempi)
      m = Variable(tempi*tempk)
      #加常数项
      A = cbind(m,BBmat)
      Y_add = rbind(1,Y)
      #penalized term
      normm = 0
      norma = 0
      normrho = 0
      for (ni in 1:tempi) {
        mk = m[(1:tempk)+(ni-1)*tempk]
        diffm = diff(mk)
        normm = normm + lamalist[ni] * sum_squares(diffm[2:length(diffm)])
      }
      for (ni in 1:tempi) {
        ak = diag(BBmat[((2:tempk)+(ni-1)*tempk),((1:(tempk-1))+(ni-1)*tempk)])
        diffa = diff(ak)
        norma = norma + lamalist[ni] * sum_squares(diffa[2:length(diffa)])
      }
      for (ni in 1:tempi) {
        for (nj in 1:tempi) {
          if (ni!=nj) {
            BBtemp = BBmat[((1:tempk)+(ni-1)*tempk),((1:tempk)+(nj-1)*tempk)]
            rhok = diag(BBtemp)
            diffrho = diff(rhok)
            normrho = normrho+lamrholist[ni,nj]* sum_squares(diffrho[2:length(diffrho)])
          }
        }
      }
      # 定义优化目标：最小化残差平方和
      objective <- Minimize(sum_squares(Y[,2:tempt]  - A %*% Y_add[,1:(tempt-1)])+normm+norma+normrho)
      
      #
      constraints <- c()
      for (i in 1:tempi) {
        for (j in 1:tempi) {
          if (i==j) {
            tempBB = BBmat[((1:tempk)+(i-1)*tempk),((1:tempk)+(j-1)*tempk)]
            tempBBlow = tempBB[2:tempk,1:(tempk-1)]
            #row_sum_constraint <- sapply(2:tempk, function(ii) tempBB[ii,ii-1]+tempBB[ii,ii] == 1)
            constraints = c(constraints,
                            diag(tempBB)>=0, diag(tempBBlow)>=0, tempBB[upper.tri(tempBB)]==0,
                            tempBBlow[lower.tri(tempBBlow)]==0)
          }else{
            tempBB = BBmat[((1:tempk)+(i-1)*tempk),((1:tempk)+(j-1)*tempk)]
            zero_constraint <- c(tempBB[upper.tri(tempBB)]==0,tempBB[lower.tri(tempBB)]==0,diag(tempBB)>=0)
            constraints = c(constraints,zero_constraint)
          }
        }
        #print('process finished',i)
      }
      row_sum_constraint_p <- sapply(1:(tempi*tempk), function(i) sum(BBmat[i, ]) == 1)
      constraints = c(constraints,row_sum_constraint_p,BBmat<=1)
      
      # 定义问题
      problem <- Problem(objective,constraints = constraints)
      
      
      # 求解优化问题
      a = solve(problem)
      coef = a$getValue(BBmat)
      mm = a$getValue(m)
      coef[abs(coef)<= kappa] = 0
      predmat = array(0,c(tempk*tempi,(tempt-1)))
      for (t in 1:(tempt-1)){
        predmat[,t] = coef%*% Y[,t]+mm
      }
      pred = array(predmat,dim = c(tempk,tempi,(tempt-1)))
      pred = aperm(pred,c(1,3,2))
      
      cvxl = a$value
      foremat = array(0,c(tempk*tempi,tempte))
      foremat[,1] = coef%*% Y[,tempt]+mm
      
      for (t in 1:(tempte-1)){
        foremat[,(t+1)] = coef%*% (foremat[,t])+mm
      }
      fore = array(foremat,dim = c(tempk,tempi,tempte))
      fore = aperm(fore,c(1,3,2))
      result$A = coef
      result$mu = mm
      result$pred = pred
    }
    #resfit = datar[,2:tempt,] - pred
    resf = datate - fore
    result$fore = fore
    result$residual = resf
    #result$mae = mean(abs(resfit))
    #result$rmse = sqrt(mean(resfit^2))
    result$rmsfe = sqrt(mean(resf^2))
    result$mafe = mean(abs(resf))
    lrpop = sqrt(apply((resf)^2,3,mean))
    result$lrpop = lrpop
    result$popsd = sd(lrpop)
    result$max = max(lrpop)
    result$min = min(lrpop)
  }
  return(result)
}

fitstar_rev = function(datar, method, fore= FALSE, datate = NULL, lambda = 0.01,kappa = 0){
  result = list()
  tempk = dim(datar)[1]
  tempt = dim(datar)[2]-1
  tempi = dim(datar)[3]
  ddatar = array(NA, dim = c(tempk,tempt,tempi))  # Resulting array will be of shape (N, M-1, P)
  # Apply diff along the second dimension (columns) for each slice (third dimension)
  for (p in 1:dim(datar)[3]) {
    ddatar[, , p] <- t(apply(datar[,,p],1,diff))
  }
  wws <- function(Y,S) {
    k =t(Y) %*% ginv(S)%*%Y
    return(k)
  }
  if (fore == FALSE) {
    pred = array(0,c(tempk,tempt,tempi))
    if (method == 'sep') {
      edf = tempi*(tempk+tempk-1+tempk-2)
      #We do not grid search
      lambda_a = lambda
      lambda_b = lambda
      cvxl = seq(0,tempi)
      for (II in 1:tempi) {
        X <- ddatar[,,II]
        Y <- X
        mm = rowMeans(Y)
        Y = Y-matrix(mm,nrow = tempk, ncol = tempt)
        # 定义P
        A <- Variable(tempk,tempk)
        # 定义优化目标：最小化残差平方和
        objective <- Minimize(sum_squares(Y[,2:tempt]  - A %*% Y[,1:(tempt-1)])#+lambda_m*sum_squares(diff(m))
                              +lambda_a*sum_squares(diff(diag(A[2:tempk,1:(tempk-1)])))
                              +lambda_b*sum_squares(diff(diag(A[3:tempk,1:(tempk-2)]))))
        
        #
        zero_constraint = c(A[upper.tri(A,diag = FALSE)]==0)
        row_sum_constraint <- sapply(1:tempk, function(i) sum(A[i, ]) == 1)
        AA = A[4:tempk,1:(tempk-3)]
        row_sum_constraint2 <- c(AA[lower.tri(AA,diag = TRUE)]==0)
        
        nonnegative_constraint1 <- sapply(2:tempk, function(i) A[i,(i-1)] >= 0)
        nonnegative_constraint2 <- sapply(3:tempk, function(i) A[i,(i-2)] >= 0)
        constraints <- c(zero_constraint,row_sum_constraint2,row_sum_constraint,nonnegative_constraint1,nonnegative_constraint2)
        
        # 定义问题
        problem <- Problem(objective,constraints = constraints)
        
        # 求解优化问题
        a = solve(problem)
        coef <- a$getValue(A)
        cvxl[II] = a$value
        for (t in 1:tempt) {
          pred[,t,II] = coef%*% X[,t]+mm
        }
      }
      y0 = pred
      # Repeat the matrix ma along the second dimension (t)
      for (i in 1:tempt) {
        y0[, i, ] <- datar[,1,]
      }
      pred_raw = aperm(apply(pred,c(1,3),cumsum),c(2,1,3))+y0
    }
    else if (method == 'stack') {
      edf = tempk*(tempi^2)-tempi
      lamalist = rep(lambda,tempi)
      lammlist=rep(lambda,tempi)
      lamrholist=matrix(lambda,tempi,tempi)
      
      # 创建一个示例的3D数组 X
      X <- ddatar
      
      # 合并 d1 和 d3 维度并重排顺序
      X_new <- aperm(X, c(1, 3, 2))  # 调整维度顺序，使得顺序为 (d1, d3, d2)
      
      # 将其转换为矩阵，合并 d1 和 d3 维度
      Y <- matrix(X_new, nrow = tempk*tempi, ncol = tempt)
      
      # 定义P
      mm = rowMeans(Y)
      Y = Y-matrix(mm,nrow = tempk*tempi, ncol = tempt)
      BBmat = Variable(tempk*tempi,tempk*tempi)
      #加常数项
      #penalized term
      norma = 0
      normrho = 0
      for (ni in 1:tempi) {
        ak = diag(BBmat[((2:tempk)+(ni-1)*tempk),((1:(tempk-1))+(ni-1)*tempk)])
        diffa = diff(ak)
        norma = norma + lamalist[ni] * sum_squares(diffm[2:length(diffa)])
      }
      for (ni in 1:tempi) {
        for (nj in 1:tempi) {
          if (ni!=nj) {
            BBtemp = BBmat[((1:tempk)+(ni-1)*tempk),((1:tempk)+(nj-1)*tempk)]
            rhok = diag(BBtemp)
            diffrho = diff(rhok)
            normrho = normrho+normrho
          }
        }
      }
      # 定义优化目标：最小化残差平方和
      objective <- Minimize(sum_squares(Y[,2:tempt]  - BBmat %*% Y[,1:(tempt-1)])+norma+normrho)
      
      #
      constraints <- c()
      for (i in 1:tempi) {
        for (j in 1:tempi) {
          if (i==j) {
            tempBB = BBmat[((1:tempk)+(i-1)*tempk),((1:tempk)+(j-1)*tempk)]
            tempBBlow = tempBB[2:tempk,1:(tempk-1)]
            #row_sum_constraint <- sapply(2:tempk, function(ii) tempBB[ii,ii-1]+tempBB[ii,ii] == 1)
            constraints = c(constraints,
                            diag(tempBB)>=0, tempBB[upper.tri(tempBB)]==0,
                            tempBBlow[lower.tri(tempBBlow)]==0)
          }else{
            tempBB = BBmat[((1:tempk)+(i-1)*tempk),((1:tempk)+(j-1)*tempk)]
            zero_constraint <- c(tempBB[upper.tri(tempBB)]==0,tempBB[lower.tri(tempBB)]==0,diag(tempBB)>=0)
            constraints = c(constraints,zero_constraint)
          }
        }
      }
      row_sum_constraint_p <- sapply(1:(tempi*tempk), function(i) sum(BBmat[i, ]) == 1)
      constraints = c(constraints,row_sum_constraint_p,BBmat<=1)
      
      # 定义问题
      problem <- Problem(objective,constraints = constraints)
      
      
      # 求解优化问题
      a = solve(problem)
      coef = a$getValue(BBmat)
      cvxl = a$value
      predmat = array(0,c(tempk*tempi,tempt))
      for (t in 1:tempt){
        predmat[,t] = coef%*% Y[,t]+mm
      }
      pred = array(predmat,dim = c(tempk,tempi,tempt))
      pred = aperm(pred,c(1,3,2))
      y0 = pred
      # Repeat the matrix ma along the second dimension (t)
      for (i in 1:tempt) {
        y0[, i, ] <- datar[,1,]
      }
      pred_raw = aperm(apply(pred,c(1,3),cumsum),c(2,1,3))+y0
    }
    res = datar[,2:tempt,] - pred_raw
    result$pred = pred
    result$residual = res
    result$status = a$status
    result$A = coef
    #in-sample-summary
    in_summary = list()
    l = matrix(0,tempi,1)
    
    for (i in 1:tempi) {
      tempres = res[,,i]
      s = var(t(datar[,,i]))
      templ = sapply(as.data.frame(tempres),wws,S=s)
      l[i] =  - 1/2*sum(templ)
    }
    in_summary$l = l
    in_summary$rmse = sqrt(mean(res^2))
    in_summary$edf = edf
    in_summary$sl = sum(l)
    in_summary$AIC = 2*sum(l) - 2*edf
    in_summary$BIC = 2*sum(l) - log(tempt)*edf
    result$summary = in_summary
  }
  else if (fore == TRUE& !is.null(datate)) {
    tempte = dim(datate)[2]
    pred = array(0,c(tempk,tempte,tempi))
    coef_A = matrix(0,tempk*tempi,tempk*tempi)
    if (method == 'sep') {
      edf = tempi*(tempk+tempk-1+tempk-2)
      #We do not grid search
      lambda_a = lambda
      lambda_b = lambda
      cvxl = seq(0,tempi)
      pred_raw = datate
      mu = rep(0,tempk*tempi)
      for (II in 1:tempi) {
        X <- ddatar[,,II]
        Y <- X
        # 定义P
        BBmat <- Variable(tempk,tempk)
        m <- Variable(tempk)
        #
        A = cbind(m,BBmat)
        Y_add = rbind(1,Y)
        # 定义优化目标：最小化残差平方和
        objective <- Minimize(sum_squares(Y[,2:tempt]  - A %*% Y_add[,1:(tempt-1)])
                              +lambda_a*sum_squares(diff(diag(A[2:tempk,1:(tempk-1)])))
                              +lambda_b*sum_squares(diff(diag(A[3:tempk,1:(tempk-2)]))))
        
        #
        zero_constraint = c(A[upper.tri(A,diag = FALSE)]==0)
        row_sum_constraint <- sapply(1:tempk, function(i) sum(A[i, ]) == 1)
        AA = A[4:tempk,1:(tempk-3)]
        row_sum_constraint2 <- c(AA[lower.tri(AA,diag = TRUE)]==0)
        
        nonnegative_constraint1 <- sapply(2:tempk, function(i) A[i,(i-1)] >= 0)
        nonnegative_constraint2 <- sapply(3:tempk, function(i) A[i,(i-2)] >= 0)
        constraints <- c(zero_constraint,row_sum_constraint2,row_sum_constraint,nonnegative_constraint1,nonnegative_constraint2)
        
        # 定义问题
        problem <- Problem(objective,constraints = constraints)
        
        # 求解优化问题
        a = solve(problem)
        coef <- a$getValue(BBmat)
        mm <- a$getValue(m)
        coef[coef<= kappa] = 0
        coef_A[((1:tempk)+(II-1)*tempk),((1:tempk)+(II-1)*tempk)] = coef
        cvxl[II] = a$value
        pred[,1,II] = coef%*% Y[,tempt]+mm
        for (t in 1:(tempte-1)) {
          pred[,(t+1),II] = coef%*% (pred[,t,II]-mm) +mm
        }
        #back to log mortality
        y0 = pred[,,II]
        # Repeat the matrix ma along the second dimension (t)
        for (i in 1:tempte) {
          y0[, i] <- datar[,tempt+1,II]
        }
        pred_raw[,,II] = t(apply(pred[,,II],1,cumsum))+y0
        mu[(1:tempk)+(II-1)*tempi] = mm
      }
      result$A = coef_A
      result$mu = mu
    }
    else if (method == 'stack') {
      edf = tempi*(tempk^2)+tempk*tempi
      lamalist = rep(lambda,tempi)
      lammlist=rep(lambda,tempi)
      lamrholist=matrix(lambda,tempi,tempi)
      
      #centra
      # 创建一个示例的3D数组 X
      X <- ddatar
      
      # 合并 d1 和 d3 维度并重排顺序
      X_new <- aperm(X, c(1, 3, 2))  # 调整维度顺序，使得顺序为 (d1, d3, d2)
      
      # 将其转换为矩阵，合并 d1 和 d3 维度
      Y <- matrix(X_new, nrow = tempk*tempi, ncol = tempt)
      #Y = Y-matrix(mm,nrow = tempk*tempi, ncol = tempt)
      # 定义P
      BBmat = Variable(tempk*tempi,tempk*tempi)
      m <- Variable(tempk*tempi)
      #加常数项
      A = cbind(m,BBmat)
      Y_add = rbind(1,Y)
      #penalized term
      norma = 0
      normrho = 0
      # for (ni in 1:tempi) {
      #   mk = m[(1:tempk)+(ni-1)*tempk]
      #   diffm = diff(mk)
      #   normm = normm + lamalist[ni] * sum_squares(diffm[2:length(diffm)])
      # }
      for (ni in 1:tempi) {
        ak = diag(BBmat[((2:tempk)+(ni-1)*tempk),((1:(tempk-1))+(ni-1)*tempk)])
        diffa = diff(ak)
        norma = norma + lamalist[ni] * sum_squares(diffa[2:length(diffa)])
      }
      for (ni in 1:tempi) {
        for (nj in 1:tempi) {
          if (ni!=nj) {
            BBtemp = BBmat[((1:tempk)+(ni-1)*tempk),((1:tempk)+(nj-1)*tempk)]
            rhok = diag(BBtemp)
            diffrho = diff(rhok)
            normrho = normrho+lamrholist[ni,nj]* sum_squares(diffrho[2:length(diffrho)])
          }
        }
      }
      # 定义优化目标：最小化残差平方和
      objective <- Minimize(sum_squares(Y[,2:tempt]  - A %*% Y_add[,1:(tempt-1)])+norma+normrho)
      
      #
      constraints <- c()
      for (i in 1:tempi) {
        for (j in 1:tempi) {
          if (i==j) {
            tempBB = BBmat[((1:tempk)+(i-1)*tempk),((1:tempk)+(j-1)*tempk)]
            tempBBlow = tempBB[2:tempk,1:(tempk-1)]
            #row_sum_constraint <- sapply(2:tempk, function(ii) tempBB[ii,ii-1]+tempBB[ii,ii] == 1)
            constraints = c(constraints,
                            diag(tempBB)>=0, diag(tempBBlow)>=0, tempBB[upper.tri(tempBB)]==0,
                            tempBBlow[lower.tri(tempBBlow)]==0)
          }else{
            tempBB = BBmat[((1:tempk)+(i-1)*tempk),((1:tempk)+(j-1)*tempk)]
            zero_constraint <- c(tempBB[upper.tri(tempBB)]==0,tempBB[lower.tri(tempBB)]==0,diag(tempBB)>=0)
            constraints = c(constraints,zero_constraint)
          }
        }
        #print('process finished',i)
      }
      row_sum_constraint_p <- sapply(1:(tempi*tempk), function(i) sum(BBmat[i, ]) == 1)
      constraints = c(constraints,row_sum_constraint_p,BBmat<=1)
      
      # 定义问题
      problem <- Problem(objective,constraints = constraints)
      
      
      # 求解优化问题
      a = solve(problem)
      coef = a$getValue(BBmat)
      mm <- a$getValue(m)
      coef[abs(coef)<= kappa] = 0
      cvxl = a$value
      predmat = array(0,c(tempk*tempi,tempte))
      predmat[,1] = coef%*% Y[,tempt] +mm
      
      for (t in 1:(tempte-1)){
        predmat[,(t+1)] = coef%*% (predmat[,t])+mm
      }
      pred = array(predmat,dim = c(tempk,tempi,tempte))
      pred = aperm(pred,c(1,3,2))
      y0 = pred
      # Repeat the matrix ma along the second dimension (t)
      for (i in 1:tempte) {
        y0[, i, ] <- datar[,tempt+1,]
      }
      pred_raw = aperm(apply(pred,c(1,3),cumsum),c(2,1,3))+y0
      result$A = coef
      result$mu = mm
    }
    
    resf = datate - pred_raw
    result$pred = pred_raw
    result$fore = pred_raw
    result$residual = resf
    result$rmsfe = sqrt(mean(resf^2))
    result$mafe = mean(abs(resf))
    lrpop = sqrt(apply((resf)^2,3,mean))
    result$lrpop = lrpop
    result$popsd = sd(lrpop)
    result$max = max(lrpop)
    result$min = min(lrpop)
  }
  return(result)
}

####################SWVAR
fitVAR <- function(data, p = 1, require_diff = FALSE, penalty = "ENET", method = "cv", ...) {
  opt <- list(...)
  #transformation
  if (length(dim(data)) >=2 & require_diff) {
    gg = as.matrix(data[,,1])
    for (i in 2:(dim(data)[3])) {
      gg = rbind(gg,data[,,i]) 
    }
    dim(gg)
    datamat_c = t(gg)
    data = diff(datamat_c)
  }
  # convert data to matrix
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  
  cnames <- colnames(data)
  
  if (method == "cv") {
    
    # use CV to find lambda
    opt$method <- "cv"
    out <- cvVAR(data, p, penalty, opt)
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

cvVAR <- function(data, p, penalty = "ENET", opt = NULL) {
  nc <- ncol(data)
  nr <- nrow(data)
  
  picasso <- ifelse(!is.null(opt$picasso), opt$picasso, FALSE)
  threshold <- ifelse(!is.null(opt$threshold), opt$threshold, FALSE)
  
  threshold_type <- ifelse(!is.null(opt$threshold_type),
                           opt$threshold_type, "soft"
  )
  
  #return_fit <- ifelse(!is.null(opt$return_fit), opt$return_fit, TRUE)
  return_fit <-  TRUE
  if (picasso) {
    stop("picasso available only with timeSlice method.")
  }
  # transform the dataset
  tr_dt <- transformData(data, p, opt)
  
  if (penalty == "ENET") {
    
    # fit the ENET model
    t <- Sys.time()
    fit <- cvVAR_ENET(tr_dt$X, tr_dt$y, nvar = nc, opt)
    elapsed <- Sys.time() - t
    
    # extract what is needed
    lambda <- ifelse(is.null(opt$lambda), "lambda.min", opt$lambda)
    
    # extract the coefficients and reshape the matrix
    Avector <- stats::coef(fit, s = lambda)
    A <- matrix(Avector[-1],
                nrow = nc, ncol = nc * p,
                byrow = TRUE
    )
    
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
  } else if (penalty == "MCP") {
    
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
  res <- computeResiduals(tr_dt$series, A, tr_dt$mu)
  
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
  
  # Do you want the fit?
  if (return_fit == TRUE) {
    output$fit <- fit
  }
  
  # Return the "best" lambda
  output$lambda <- fit$lambda.min
  output$pred <- tr_dt$series - res
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
  
  if (fit$lambda.min == 0) {
    k <- prod(dim(A[[1]]))
  }else{
    k <- sum(A[[1]] != 0)
  }
  nn <- fit$glmnet.fit$nobs
  AIC <- tLL-2*k 
  BIC<-tLL-log(nn)*k
  output$k = k
  output$mrmse = sqrt(mean(res^2))
  output$mmae = mean(abs(res))
  output$logL = tLL
  output$AIC = AIC
  output$BIC = BIC
  output$coef = Avector
  output$Deviance = (1-fit$glmnet.fit$dev.ratio[index])*fit$glmnet.fit$nulldev
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


fitVAR_auto <- function(datatr,datate,sparse = TRUE){
  #datapre
  glot = dim(datatr)[2]
  gloi = dim(datatr)[3]
  datar = datatr
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
  forecast_nod <- function(model,x0,rx0,data_c) {
    #dim(ddatamat_c) = 49,294
    nc = ncol(data_c)
    nr = nrow(data_c)
    m = c(model$mu)
    pred = matrix(0,ncol=nc,nrow=nr)
    pred[1,] = t(model$A[[1]] %*% (x0)) +m
    for (i in 2:nr) {
      pred[i,] = t(model$A[[1]] %*% (pred[(i-1),]-m)) +m
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
    sum$pop = rg
    sum$rmsfe = sqrt(mean(res^2))
    sum$mafe = mean(abs(res))
    sum$meanpop = mean(rg)
    sum$mean = mean(abs(res))
    sum$sd = sd(rg)
    sum$max = rg[rg=which.max(rg)]
    sum$min = rg[rg=which.min(rg)]
    return(sum)
  }
  
  ddatamat_c = diff(datamat_c)
  ddatamat_te = diff(datamat_te)
  if (sparse == TRUE) {
    VAR_4 = fitVAR(ddatamat_c,p=1)#0.0000119 360 1000
    VAR_4 = fitVAR(ddatamat_c,p=1,lambdas_list=c(2*VAR_4$lambda))#0.0000119 360 1000
  }else{
    VAR_4 = fitVAR(ddatamat_c,p=1,lambdas_list=c(0,1))#0.0000119 360 1000
  }  
  f4=forecast_nod(VAR_4,ddatamat_c[nrow(ddatamat_c),],datamat_c[nrow(datamat_c),],datamat_te)
  out = list()
  out$model = VAR_4
  out$model$predr = out$model$pred+datamat_c[1:(glot-1),]
  out$fore = f4
  return(out)
}



##################################
fitlilee <- function(datar,data_star_low_tr,datate,robust = FALSE){
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
  
  glok = dim(datar)[1]
  glot = dim(datar)[2]
  glote = dim(datate)[2]
  gloi = dim(datar)[3]
  #--------------------li_lee---------------------------
  Global = LeeCarter(data_star_low_tr)#K*T
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
  out = list(lB = lB,lK = lK, lista = lista, listb = listb, listk = listk)
  
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
  out$pred = predy
  out$mae = mean(abs(predy-datar))
  out$rmse = sqrt(mean((predy-datar)^2))
  
  #number of parameters
  orin = length(lista)+length(listb)+length(listk)+length(lB)+length(lK)
  #out-of-sample
  #listk
  custompred <- function(c,t=glote,x0) {
    forec = matrix(0,nrow = 1,t)
    forec[1] = c[1]+c[2]*x0
    for (i in 2:t) {
      forec[i] = c[1]+c[2]*forec[i-1]
    }
    return(forec)
  }
  #with ar fit k
  fore_lilee2 <- function(timestep,dtrain=datatr,dtest = datate,robust = FALSE){
    predk = matrix(0,nrow = gloi,ncol = timestep)
    fitvalue = matrix(0,gloi,2)
    r = matrix(0,gloi,1)
    fore = array(0,dim = c(glok,timestep,gloi))
    forematrix = matrix(0,timestep,gloi*glok)
    result = list()
    fit1 <- forecast::rwf(lK,drift = TRUE, h = timestep)
    plK <- fit1$mean
    for (i in 1:gloi) {
      dt0 = dtrain[,,i][,ncol(dtrain[,,i])]
      tsd = ts(listk[i,])
      #fit <- ar(tsd,order.max =1,h=timestep)
      #fitvalue[i,1] = fit$x.mean
      
      if (robust & dim(dtrain)[3] >=10 ){
        fit <- ar(tsd,order.max =1)
        fitvalue[i,1] = fit$x.mean
        fitsd = 0.001
        if (dim(dtrain)[1]>=30) {
          fitsd = 0.1
        }
        if (length(fit$ar) == 0) {
          fitvalue[i,2] = 0
        }else{
          fitvalue[i,2] = fit$ar
        }
        
      }else{  
        fit <- ar(tsd,order.max =1)
        fitvalue[i,1] = fit$x.mean
        fitsd = 0
        if (length(fit$ar) == 0) {
          fitvalue[i,2] = 0
        }else{
          fitvalue[i,2] = fit$ar
        }
        predk[i,]=custompred(fitvalue[i,],t=timestep,tsd[length(tsd)])
        
      }
      
      
      #r[i] = 1-(var(residuals(fit))/var(tsd))
      a = matrix(rep(dt0,timestep),glok,timestep)
      gp = lB%*%t(plK-lK[length(lK)])
      ip = listb[i,]%*%t(predk[i,]-listk[i,ncol(listk)])
      fore[,,i] = a+gp+ip+rnorm(1,0,fitsd)
      #forematrix[,((1:glok)+((i-1)*glok))]= t(a+gp+ip)
    }
    result$forecast = fore
    #result$forematrix = forematrix
    result$k= predk
    result$coef = fitvalue
    result$K = plK
    #result$Rar = r
    result$res = fore-dtest
    result$rmsfe = sqrt(mean((fore-dtest)^2))
    return(result)
  }
  f0 = list()
  f0$fit = out
  f02 = fore_lilee2(timestep = glote, dtrain=datar, dtest = datate, robust = robust)
  f0$fore = f02$forecast
  f0$res = f02$forecast-datate
  f0$RMSFE = f02$rmsfe
  lrpop =sqrt(apply((f02$forecast-datate)^2,3,mean))
  f0$popmean = mean(lrpop)
  f0$popsd = sd(lrpop)
  f0$max = max(lrpop)
  f0$min = min(lrpop)
  #mae
  mafe = mean(abs(f02$forecast-datate))
  f0$mafe = mafe
  return(f0)
}

fitsvar_sep <- function(datar,datate,mu){
  glok = dim(datar)[1]
  gloi = dim(datar)[3]
  l =matrix(0,gloi,1)
  edf = 0
  coefmat_s = matrix(0,glok*gloi,glok*gloi)
  gg = as.matrix(datar[,,1])
  for (i in 2:gloi) {
    gg = rbind(gg,datar[,,i]) 
  }
  datamat_c = t(gg)
  ggg = as.matrix(datate[,,1])
  for (i in 2:gloi) {
    ggg = rbind(ggg,datate[,,i]) 
  }
  datamat_te = t(ggg)
  ddatamat_c = diff(datamat_c)
  sigseq = matrix(0,glok,glok*gloi)
  for (i in 1:gloi) {
    sigseq[,((1:glok)+glok*(i-1))] = var(diff(t(datar[,,i])))
  }
  wws <- function(Y,S) {
    k =t(Y) %*% ginv(S)%*%Y
    return(k)
  }
  for (i in 1:gloi) {
    data = t(datar[,,i])
    ddata = diff(data)
    ddatamat_te = diff(datamat_te)
    VAR6i = fitVAR(ddata,p=1)
    s = sigseq[,((1:glok)+(i-1)*glok)]
    tempdfy = VAR6i$residuals[-1,]
    templ = sapply(as.data.frame(t(tempdfy)),wws,S=s)
    l[i] =  - 1/2*sum(templ)
    edf[i] =  sum(!VAR6i$coef ==0)
    coefmat_s[((1:glok)+(i-1)*glok),((1:glok)+(i-1)*glok)]  = VAR6i$A[[1]]
  }
  f5=forecast_our2_nod(coefmat_s,mean = mu,ddatamat_c[nrow(ddatamat_c),],datamat_c[nrow(datamat_c),],datamat_te)
  return(f5)
}
