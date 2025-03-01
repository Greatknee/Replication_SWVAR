#general forecasting
#easy mode
#prediction
#14 low mortality pop



#################

func_figure_03_prefit <- function(data){
  set.seed(123)
  datagroup = data$group
  coulist = data$coulist
  datar = data$datatr
  datate = data$datate
  
  data_star_low_tr = data$data_star_low_tr
  glot = data$glot
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
  
  #Model1:li-lee
  f0 = fitlilee(datar,data_star_low_tr,datate)
  ##########################################################################
  ##########################################################################
  #FDM
  fdm = fitfdm(coulist = coulist, fore = TRUE, year = c(data$ystr,data$yetr,data$yste,data$yete),age = 5)
  #STAR stack first,model1,2
  kkapa_stack = c(0,0,1,0,0.09)
  VAR_1 = fitstar_rev(datar = datar,'stack',fore = TRUE,datate = datate,lambda = 0.5,kappa = kkapa_stack[datagroup])
  ##########################################################################
  #fit model1
  #model34 AR+Lasso
  VAR_4 = fitVAR(ddatamat_c,p=1)#0.0000119 360 1000
  f4=forecast_nod(VAR_4,ddatamat_c[nrow(ddatamat_c),],datamat_c[nrow(datamat_c),],datamat_te)
  #SWVAR
  f6 = fitswvar_lack_global(datar = datar, group = datagroup,datate = datate, coulist= coulist)
  
  return(list(f0 = f0,
           fdm = fdm,
           VAR_1 = VAR_1,
           f4 = f4,
           f6 = f6,
           coulist = data$coulist))
}