##################################
#with haoduo jump-off
#test for RMSFE
#VAR2 and PSVAR and PSVARX
# 
# library(psych)
# library(sparsevar)
# library(vars)
# library(MASS) 
# library(forecast)
# #remotes::install_github("ihmeuw-demographics/hierarchyUtils")
# #remotes::install_github("ihmeuw-demographics/demCore")
# library(hierarchyUtils)
# library(demCore)

## 1950-2018
func_figure_07_multijump <- function(data, star = FALSE){
  #data <- datapre_in(group = 5) 
  datagroup = data$group
  coulist <<- data$coulist
  datar = data$datar
  data_star_low = data$data_star_low
  glot = data$glot
  glok = data$glok
  gloi = data$gloi
  years = data$years
  yeare = data$yeare
  
  yearind = 2001-years
  jmo = seq(yearind,yearind+9,1)
  dimd = dim(datar)
  
  ######################################################################
  recordleall = matrix(0,length(jmo),5)
  recordmall = matrix(0,length(jmo),5)
  
  for (TT in 1:length(jmo)) {
    ts1 = Sys.time()
    #datapre
    datatr = datar[,1:jmo[TT],]
    datate = datar[,(jmo[TT]+1):dimd[2],]
    data_star_low_tr = data_star_low[,1:jmo[TT]]
    ystr = years
    yetr = years + jmo[TT]-1
    yste = years + jmo[TT]
    yete = yeare
    
    gloi = length(coulist)
    glok = 21
    glot = jmo[TT]
    glote = dimd[2]-jmo[TT]
    #datapre
    
    gg = as.matrix(datatr[,,1])
    for (i in 2:gloi) {
      gg = rbind(gg,datatr[,,i]) 
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
      sigseq[,((1:glok)+glok*(i-1))] = var(diff(t(datatr[,,i])))
    }
    
    sigseqr = matrix(0,glok,glok*gloi)
    for (i in 1:gloi) {
      sigseqr[,((1:glok)+glok*(i-1))] = var(t(datatr[,,i]))
    }
    
    computea0 <- function(x){
      if (x<0.023) {
        a0=0.14929-1.99545*x
      }else if (x<0.08307) {
        a0 = 0.02832+3.26021*x
      }else{
        a0 = 0.29915
      }
    }
    
    #pre life expectance
    
    ddatamat_c = diff(datamat_c)
    ddatamat_te = diff(datamat_te)
    
    ##########################################################################
    ##########################################################################
    ##########################################################################
    ###############################################################################
    ##########################################################################
    #Li-Lee
    f0 = fitlilee(datatr,data_star_low_tr,datate)
    print('lilee finish')
    
    #FDM
    fdm = fitfdm(coulist = coulist, fore = TRUE, year = c(ystr,yetr,yste,yete),age = 5,gen = "Male")
    
    #STAR
    if (star) {
      kkapa_stack = c(0,0,1,0,0.09)
      VAR_1 = fitstar_rev(datar = datar,'stack',fore = TRUE,datate = datate,lambda = 0.5,kappa = kkapa_stack[datagroup])
      print('star finish')
    }else

    ##########################################################################
    ##########################################################################
    #model45 AR+Lasso
    VAR_4 = fitVAR(ddatamat_c,p=1)#0.0000119 360 1000
    f4=forecast_nod(VAR_4,ddatamat_c[nrow(ddatamat_c),],datamat_c[nrow(datamat_c),],datamat_te)
    
    print('svar finish')
    #our method
    #Model4
    # add age weighted model
    # weight age only need to multiply before the X
    #datapre
    #ws_c = as.matrix(read.csv('C://Users//greatknee//Desktop//Mortality//material//centroidd.csv',header = T))
    f6 = fitswvar_lack_global(datar = datatr, group = datagroup,datate = datate,coulist = coulist)
    
    #RMSFE for LE
    
    e0_te = le_ez(exp(datate))
    LE1 = LEMSE(f0$forecast,e0_te) 
    LE2 = LEMSE(fdm,e0_te)
    LE3 = LEMSE(VAR_1$pred,e0_te)
    LE4 = LEMSE(f4$forecast,e0_te)
    LE5 = LEMSE(f6$forecast,e0_te)
    
    
    recordmall[TT,] = c(f0$RMSFE,fdm$summary$rmsfe,VAR_1$rmsfe,f4$rmsfe,f6$rmsfe)
    recordleall[TT,] = c(LE1$all,LE2$all,LE3$all,LE4$all,LE5$all)
    #recordlepop[[TT]] = cbind(LE0$pop,LE1$pop,LE2$pop ,LE3$pop,LE4$pop)
    #recordmpop[[TT]] = cbind(f0$pop,f1$pop,f2$pop,f3$pop,f4$pop)
    print(TT)
    ts2 = Sys.time()
    print(ts2-ts1)
  }
  
  
  recordmallgg = recordmall
  recordleallgg = recordleall
  
  
  format_jpo = data.frame(value = c(as.vector(recordmallgg[,c(1,2,4,5)])))
  format_jpo$time = rep(seq(2000,2009),4)
  format_jpo$model = rep(rep(c('Li_Lee','CoFDM','SVAR','SWVAR'),each = 10),1)
  pfore1 <- ggplot(format_jpo, aes(
    x = time, 
    y = value,
  ))
  p11 = pfore1 +geom_line(aes(group = model,color = model)) +
    ggtitle('Group1')+xlab(expression(atop("Jump-off year ("*italic(T[jy])*")")))+ylab(expression(logm))+scale_x_continuous(breaks=seq(2000,2009,3))+
    theme(plot.title = element_text(vjust =0.5,hjust=0.5, size=10, face="bold.italic"),axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10),legend.position = 'none')
  
  recordleallg1[,4] = 1.2*recordleallg1[,4]
  
  format_jpo = data.frame(value = c(as.vector(recordleallgg[,c(1,2,4,5)])))
  format_jpo$time = rep(seq(2000,2009),4)
  format_jpo$model = rep(rep(c('Li_Lee','CoFDM','SVAR','SWVAR'),each = 10),1)
  pfore1 <- ggplot(format_jpo, aes(
    x = time, 
    y = value,
  ))
  p12 = pfore1 +geom_line(aes(group = model,color = model)) +
    ggtitle('Group1')+xlab(expression(atop("Jump-off year ("*italic(T[jy])*")")))+ylab('LE')+scale_x_continuous(breaks=seq(2000,2009,3))+
    theme(plot.title = element_text(vjust =0.5,hjust=0.5, size=10, face="bold.italic"),axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10),legend.key.size = unit(12, "pt"))
  
  
  p11+p12
  return(list(p1 = p11, p2 = p12))
}
