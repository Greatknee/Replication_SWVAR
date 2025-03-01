# fit the model
#general forecasting
#easy mode
#prediction
#14 low mortality pop




#################
func_figure_01_coef <- function(data){
  datagroup = data$group
  coulist = data$coulist
  datar = data$datar
  data_star_low = data$data_star_low
  glot = data$glot
  glok = data$glok
  gloi = data$gloi
  gg = as.matrix(datar[,,1])
  for (i in 2:gloi) {
    gg = rbind(gg,datar[,,i]) 
  }
  dim(gg)
  datamat_c = t(gg)
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
    sigseq[,((1:glok)+glok*(i-1))] = var(diff(t(datar[,,i])))
  }
  
  sigseqr = matrix(0,glok,glok*gloi)
  for (i in 1:gloi) {
    sigseqr[,((1:glok)+glok*(i-1))] = var(t(datar[,,i]))
  }
  
  ##########################################################################
  #SVAR
  ddatamat_c = diff(datamat_c)
  VAR_4 = fitVAR(ddatamat_c,p=1)#0.0000119 360 1000
  #################################################
  #SWVAR
  #our method
  #Model7
  # add age weighted model
  # weight age only need to multiply before the X
  #datapre
  if (datagroup == 1 | datagroup == 5) {
    ws_c = as.matrix(read.csv('data//centroidd.csv',header = T))
  }else if(datagroup == 2 | datagroup == 3 | datagroup == 4) {
    ws_c = as.matrix(read.csv('data//centroidd2.csv',header = T))
  }
  ws_c = ws_c[match(coulist,colnames(ws_c)),match(coulist,colnames(ws_c))]
  listl = matrix(0,gloi,1)
  templist = 1:gloi
  besta = matrix(0,gloi,1)
  bests = matrix(0,gloi,1)
  predmse = matrix(0,gloi,1)
  cu7 = list()
  #grid search
  aseq = 1.5#seq(1,50,5)
  sseq = 1#seq(1,6,0.5)
  
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
    
    sig = sigseqr[,((1:glok)+(i-1)*glok)]
  
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
        scaledws = c(min(tempws),sseq[jj]*tempws)/sum(tempws)
        result_list = list()
        for (k in 1:gloi) {
          result_list[[k]] <-  scaledws[k] *wa
        }
        ws = do.call(cbind, result_list)
        ws = as.vector(t(ws))
        
        VAR = fitVAR_weighted(dtsdata,dspamat,weight = ws,p=1)
        edf = sum(!VAR$coef ==0)
        dfy = VAR$pred + tsdata[1:(glot-1),] - tsdata[2:glot,]
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
  
  
  
  besta
  bests
  bestagroup2 = aseq[besta]
  bestsgroup2 = sseq[bests]
  loopresidual = matrix(0,glot-2,gloi*glok)
  loopcoef = matrix(0,gloi*glok,gloi*glok)
  loopedf = matrix(0,gloi,1)
  looppred = matrix(0,glot-1,gloi*glok)
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
        wa[i3,j3] = exp(abs(i3-j3)/bestagroup2[i])
      }
    }
    
    scaledws = c(min(tempws),bestsgroup2[i]*tempws)/sum(tempws)
    result_list = list()
    for (k in 1:gloi) {
      result_list[[k]] <-  scaledws[k] *wa
    }
    ws = do.call(cbind, result_list)
    ws = as.vector(t(ws))
    VAR_8 = fitVAR_weighted(dtsdata,dspamat,weight = ws,p=1)#加了alpha = 0.9
    looppred[,(1:glok)+(i-1)*glok] = VAR_8$pred
    loopresidual[,(1:glok)+(i-1)*glok] = VAR_8$residuals[-1,]
    loopedf[i] = sum(!VAR_8$coef ==0)
    loopcoef[(1:glok)+(i-1)*glok,] = matrix(VAR_8$coef[-1],nrow = glok,ncol = glok*gloi,byrow = TRUE)
    
  }
  
  #prepare for coefmat
  coefmat = loopcoef
  for (i in 2:(gloi)) {
    
    coefmat[((1:glok)+(i-1)*glok),(1:((i-1)*glok))] = loopcoef[((1:glok)+(i-1)*glok),((glok+1):(i*glok))]
    coefmat[((1:glok)+(i-1)*glok),(((i-1)*glok+1):(i*glok))] = loopcoef[((1:glok)+(i-1)*glok),(1:glok)]
  }
  coefmat
  
  
  
  
  
  
  
  #四个部分
  tempcou = c('AUT','NLD','NOR','SWE')
  index = which(coulist %in% tempcou)
  indm = rep(1:21,length(tempcou))+rep((index-1)*21,each = 21)
  x <- VAR_4$A[[1]][indm,indm]
  x <-  x /(max(x) - min(x))
  
  fk = 21
  fi = length(tempcou)
  xrev = apply(x, 2, rev)
  df1 = data.frame(value = as.vector(x))
  aglist = seq(0,100,5)#as.character(seq(1,fk,1))
  aglistinv = seq(100,0,-5)#as.character(seq(fk,1,-1))
  df1$x = rep(rep(aglist,each = fk*fi),fi)
  df1$y = rep(rep(aglistinv,fi),fk*fi)
  df1$pop_r = rep(rep(tempcou,each = fk),fk*fi)
  df1$pop_e = rep(tempcou,each = fk^2*fi)
  
  #longData<-longData[longData$value!=0,]
  
  g1=ggplot(df1, aes(x = x, y = y)) + 
    geom_tile(aes(fill=value)) + facet_grid(cols = vars(pop_e),rows = vars(pop_r))+
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)+
    labs(x="Age Group(Explanatory Variable)", y="Response Variable", title="SVAR") +
    theme(panel.grid=element_blank())+ scale_x_continuous(breaks=seq(0,100,20))+scale_y_continuous(breaks = seq(0,100,20),labels = as.character(seq(100,0,-20)))+
    theme(axis.text.x=element_text(size=4, angle=0, vjust=0.3),
          axis.text.y=element_text(size=4),
          plot.title=element_blank())
  
  x <- coefmat[indm,indm]
  x <-  x /(max(x) - min(x))
  
  
  xrev = apply(x, 2, rev)
  df2 = data.frame(value = as.vector(x))
  aglist = seq(0,100,5)#as.character(seq(1,fk,1))
  aglistinv = seq(100,0,-5)#as.character(seq(fk,1,-1))
  df2$x = rep(rep(aglist,each = fk*fi),fi)
  df2$y = rep(rep(aglistinv,fi),fk*fi)
  df2$pop_r = rep(rep(tempcou,each = fk),fk*fi)
  df2$pop_e = rep(tempcou,each = fk^2*fi)
  
  #longData<-longData[longData$value!=0,]
  
  g2=ggplot(df2, aes(x = x, y = y)) + 
    geom_tile(aes(fill=value)) + facet_grid(cols = vars(pop_e),rows = vars(pop_r))+
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)+
    labs(x="Age Group(Explanatory Variable)", y="Response Variable", title="SWVAR") +
    theme(panel.grid=element_blank()) +scale_x_continuous(breaks=seq(0,100,20))+scale_y_continuous(breaks = seq(0,100,20),labels = as.character(seq(100,0,-20)))+
    theme(axis.text.x=element_text(size=4, angle=0, vjust=0.3),
          axis.text.y=element_text(size=4),
          plot.title=element_blank())
  #g1+g2
  
  
  
  
  #四个部分
  tempcou = c('CAN','ENG','FRA','USA')
  index = which(coulist %in% tempcou)
  indm = rep(1:21,length(tempcou))+rep((index-1)*21,each = 21)
  x <- VAR_4$A[[1]][indm,indm]
  x <-  x /(max(x) - min(x))
  
  fk = 21
  fi = length(tempcou)
  xrev = apply(x, 2, rev)
  df1 = data.frame(value = as.vector(x))
  aglist = seq(0,100,5)#as.character(seq(1,fk,1))
  aglistinv = seq(100,0,-5)#as.character(seq(fk,1,-1))
  df1$x = rep(rep(aglist,each = fk*fi),fi)
  df1$y = rep(rep(aglistinv,fi),fk*fi)
  df1$pop_r = rep(rep(tempcou,each = fk),fk*fi)
  df1$pop_e = rep(tempcou,each = fk^2*fi)
  
  #longData<-longData[longData$value!=0,]
  
  g3=ggplot(df1, aes(x = x, y = y)) + 
    geom_tile(aes(fill=value)) + facet_grid(cols = vars(pop_e),rows = vars(pop_r))+
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)+
    labs(x="Age Group(Explanatory Variable)", y="Response Variable", title="SVAR") +
    theme(panel.grid=element_blank())+ scale_x_continuous(breaks=seq(0,100,20))+scale_y_continuous(breaks = seq(0,100,20),labels = as.character(seq(100,0,-20)))+
    theme(axis.text.x=element_text(size=4, angle=0, vjust=0.3),
          axis.text.y=element_text(size=4),
          plot.title=element_blank())
  
  
  x <- coefmat[indm,indm]
  x <-  x /(max(x) - min(x))
  xrev = apply(x, 2, rev)
  df2 = data.frame(value = as.vector(x))
  aglist = seq(0,100,5)#as.character(seq(1,fk,1))
  aglistinv = seq(100,0,-5)#as.character(seq(fk,1,-1))
  df2$x = rep(rep(aglist,each = fk*fi),fi)
  df2$y = rep(rep(aglistinv,fi),fk*fi)
  df2$pop_r = rep(rep(tempcou,each = fk),fk*fi)
  df2$pop_e = rep(tempcou,each = fk^2*fi)
  
  #longData<-longData[longData$value!=0,]
  
  g4=ggplot(df2, aes(x = x, y = y)) + 
    geom_tile(aes(fill=value)) + facet_grid(cols = vars(pop_e),rows = vars(pop_r))+
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)+
    labs(x="Age Group(Explanatory Variable)", y="Response Variable", title="SWVAR") +
    theme(panel.grid=element_blank()) +scale_x_continuous(breaks=seq(0,100,20))+scale_y_continuous(breaks = seq(0,100,20),labels = as.character(seq(100,0,-20)))+
    theme(axis.text.x=element_text(size=4, angle=0, vjust=0.3),
          axis.text.y=element_text(size=4),
          plot.title=element_blank())
  
  #save output
  pdf('output/result_figure_011.pdf', width=12, height=6)
  print(g1+g2)  
  dev.off()
  pdf('output/result_figure_012.pdf', width=12, height=6)
  print(g3+g4)  
  dev.off()
  return(list(g1 = g1, g2=g2, g3 = g3, g4 = g4))
}