
func_figure_08_lilee <- function(){
  data <- datapre_ll(group = 1) 
  ystr = data$ystr
  yetr = data$yetr
  yste = data$yste
  yete = data$yete
  
  datagroup = data$group
  coulist = data$coulist
  datar = data$datatr
  
  data_star_low_tr = data$data_star_low_tr
  glot = data$glot
  glote = data$glote
  glok = data$glok
  gloi = data$gloi
  
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
  
  #--------------------li_lee---------------------------
  Global = LeeCarter(data_star_low_tr)
  lA = Global$coefficients$ax#[age]
  lB = Global$coefficients$bx
  lK = Global$coefficients$kt
  #----------------LCperpopulation---------------------------
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
    #listb[13,] = c(-0.7,-0.65,0.6,0.26,0.45,-0.45,-1.26,-1.31,-0.42,0.1,0.35,0.55,0.62,0.67,0.7,0.64,0.5,0.35,0.2,0.1,0.03)
    # model = LeeCarter(data)
    # lista[i,] = model$coefficients$ax 
    # listb[i,] = model$coefficients$bx#[pop,age]
    # listk[i,] = model$coefficients$kt#[pop,year]
  }
  #lista = sweep(lista,2,Global$coefficients$ax, FUN = "+")#[pop,age]
  listb[13,] = -listb[13,]/400
  out = list(lB = lB,lK = lK, lista = lista, listb = listb, listk = listk,data = data)
}
  