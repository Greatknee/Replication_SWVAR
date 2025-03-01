#################
func_table_04_outsample <- function(data){
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
  VAR_1$rmsfe
  VAR_1$mafe
  VAR_2 = fitstar(datar = datar,'sep',fore = TRUE,datate = datate,kappa = 0.015)
  VAR_2$rmsfe
  ##########################################################################
  #fit model1
  #model34 AR+Lasso
  VAR_3 = fitVAR(ddatamat_c,p=1,lambdas_list=c(0,0.5,1))
  VAR_4 = fitVAR(ddatamat_c,p=1)#0.0000119 360 1000
  f3=forecast_nod(VAR_3,ddatamat_c[nrow(ddatamat_c),],datamat_c[nrow(datamat_c),],datamat_te)
  f4=forecast_nod(VAR_4,ddatamat_c[nrow(ddatamat_c),],datamat_c[nrow(datamat_c),],datamat_te)
  #################################################
  #VAR5 fit every pop
  f5 = fitsvar_sep(datar,datate,mu = VAR_3$mu)
  #################################################
  #VAR5 fit every pop
  #our method
  #Model6
  # add age weighted model
  # weight age only need to multiply before the X
  #datapre
  #ws_c = as.matrix(read.csv('C://Users//greatknee//Desktop//Mortality//material_VAR//centroidd.csv',header = T))
  #ws_c = as.matrix(read.csv('C://Users//greatknee//Desktop//Mortality//material_VAR//centroidd2.csv',header = T))
  #read.csv(file ='C:/Users/greatknee/Desktop/Mortality/material_VAR/swvar_param/group5_1950_male.csv')
  f6 = fitswvar_lack_global(datar = datar, group = datagroup,datate = datate,coulist = coulist)
  
  #write.csv(f6$best, file ='C:/Users/greatknee/Desktop/Mortality/material_VAR/swvar_param/group1_1950_male.csv')
  # f62 = fitswvar_ori(datar = datar, spaw = 1)
  # f63 = fitswvar_reweight(datar = datar,weight = VAR_1$A, spaw = 1)
  ####################################################
  
  conclu = matrix(0,5,8)
  conclu[1,] = c(f0$RMSFE,fdm$summary$rmsfe,VAR_1$rmsfe,VAR_2$rmsfe,f3$rmsfe,f4$rmsfe,f5$rmsfe,f6$rmsfe)
  conclu[2,] = c(f0$popsd,sd(fdm$summary$pop),VAR_1$popsd,VAR_2$popsd,f3$sd,f4$sd,f5$sd,f6$sd)
  conclu[3,] = c(f0$max,max(fdm$summary$pop),VAR_1$max,VAR_2$max,f3$max,f4$max,f5$max,f6$max)
  conclu[4,] = c(f0$min,min(fdm$summary$pop),VAR_1$min,VAR_2$min,f3$min,f4$min,f5$min,f6$min)
  conclu[5,] = c(f0$mafe,fdm$summary$mafe,VAR_1$mafe,VAR_2$mafe,f3$mae,f4$mae,f5$mae,f6$mae)
  
  
  rownames(conclu) = c('RMSFE','sigma','max','min','MAFE')
  colnames(conclu) = c('Li-Lee','FDM','STAR1','STAR2','VAR','SVAR1','SVAR2','SWVAR')
  
  conclu = round(conclu,4)
  
  #write.csv(conclu, "output/result.csv")
  return(conclu)
}