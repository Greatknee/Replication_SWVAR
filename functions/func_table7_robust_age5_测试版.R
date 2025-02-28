#general forecasting
set.seed(123)
#datapre
#data <- datapre_robust(group = 2,setting = '1',gender = 'female')
datagroup = data$group
coulist = data$coulist
if (data$gender == 'female') {
  gen_fdm = 'Female'
}else{
  gen_fdm = 'Male'
}
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
fdm = fitfdm(coulist = coulist, fore = TRUE, year = c(data$ystr,data$yetr,data$yste,data$yete),age = 5,gen = gen_fdm)
fdm$summary$mafe
fdm$summary$rmsfe
#STAR stack first,model1,2
VAR_1 = fitstar(datar = datar,'stack',fore = TRUE,datate = datate,lambda = 0.1,kappa = 0.002)
VAR_1$rmsfe
##########################################################################
#fit model1
#model34 AR+Lasso
VAR_4 = fitVAR(ddatamat_c,p=1)#0.0000119 360 1000
f4=forecast_nod(VAR_4,ddatamat_c[nrow(ddatamat_c),],datamat_c[nrow(datamat_c),],datamat_te)
#################################################
#VAR5 fit every pop
#f5 = fitsvar_sep(datar,datate)
#our method
#Model6
# add age weighted model
# weight age only need to multiply before the X
f6 = fitswvar_lack_global(datar = datar,group = datagroup,datate = datate,coulist = coulist)

# f62 = fitswvar_ori(datar = datar, spaw = 1)
# f63 = fitswvar_reweight(datar = datar,weight = VAR_1$A, spaw = 1)
####################################################


conclu = matrix(0,2,5)
conclu[1,] = c(f0$RMSFE,fdm$summary$rmsfe,VAR_1$rmsfe,f4$rmsfe,f6$rmsfe)
conclu[2,] = c(f0$mafe,fdm$summary$mafe,VAR_1$mafe,f4$mae,f6$mae)

rownames(conclu) = c('RMSFE','MAFE')
colnames(conclu) = c('Li-Lee','FDM','STAR','SVAR','SWVAR')
conclu