data <- datapre(group = 2)

#datapre
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
#model3  VAR
ddatamat_c = diff(datamat_c)
VAR_4 = fitVAR(ddatamat_c,p=1)#0.0000119 360 1000

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
  aseq = 1.5#seq(1,50,5)
  sseq = 1#seq(1,6,0.5)
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
  VAR_8 = fitVAR_weighted(dtsdata,dspamat,weight = ws,p=1,lambdas_list= 1.2*VAR_8$lambda)#加了alpha = 0.9
  looppred[,(1:glok)+(i-1)*glok] = VAR_8$pred
  loopresidual[,(1:glok)+(i-1)*glok] = VAR_8$residuals[-1,]
  loopedf[i] = sum(!VAR_8$coef ==0)
  loopcoef[(1:glok)+(i-1)*glok,] = matrix(VAR_8$coef[-1],nrow = glok,ncol = glok*gloi,byrow = TRUE)
  
}




#1.JTest
#datapre
ghat = ddatamat_c[-1,]-loopresidual

VARJ = fitVAR_Jtest2(ddatamat_c[-1,],datap = ghat,p=1,lambdas_list=c(0,1))
sumj = customsummary2(VARJ,glok,glot,gloi)
sumj
VAR_12 = fitVAR(ddatamat_c[-1,],p=1,lambdas_list=c(0,1))
sumj0 = customsummary2(VAR_12,glok,glot,gloi)
sumj0
chi = 2*(sumj$sl-sumj0$sl)
pvalue = 1-pchisq(2*abs(sumj$sl-sumj0$sl),1)
chi
pvalue


#3.C-test
a1 = as.vector(ddatamat_c[-1,]-VAR_4$residuals[-1,])
a2 = as.vector(ddatamat_c[-1,]-loopresidual)
y = as.vector(VAR_4$residuals[-1,])
data = data.frame(y)
data$x1 = a2 - a1


fit = lm(formula = y~0+x1,data = data)
#summary(fit)
summary(fit)$coefficients
chi
pvalue
#return(list(Ctest = summary(fit)[4], Jtest = c(chi,pvalue)))
