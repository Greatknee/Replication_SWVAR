##########################


#general forecasting
set.seed(123)
#datapre
data <- datapre_robust(group = 1,setting = '2',gender = 'female')
set = data$setting
datagroup = data$group
coulist = data$coulist
if (data$gender == 'female') {
  gen_fdm = 'Female'
}else{
  gen_fdm = 'Male'
}

ystr = data$ystr
yetr = data$yetr
yste = data$yste
yete = data$yete

coulist = data$coulist
datatr1 = data$datatr
datate1 = data$datate

data_star_low_tr1 = data$data_star_low_tr
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
ggg = as.matrix(datate1[,,1])
for (i in 2:gloi) {
  ggg = rbind(ggg,datate1[,,i]) 
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
f0 = fitlilee(datatr1,data_star_low_tr1,datate1)
##########################################################################
##########################################################################
#FDM
fdm = fitfdm(coulist = coulist, fore = TRUE, year = c(ystr,yetr,yste,yete),age = 1,gen = gen_fdm)
fdm$summary$rmsfe
fdm$summary$mae
#STAR stack first,model1,2
VAR_1 = fitstar(datar = datatr1,'stack',fore = TRUE,datate = datate1,lambda = 1)
#VAR_2 = fitstar_rev(datar = datatr,'sep',fore = TRUE,datate = datate)
##########################################################################
#fit model1
#model34 AR+Lasso
#VAR_3 = fitVAR(ddatamat_c,p=1,lambdas_list=c(0,0.5,1))
VAR_4 = fitVAR(ddatamat_c,p=1)#0.0000119 360 1000
#f3=forecast_nod(VAR_3,ddatamat_c[nrow(ddatamat_c),],datamat_c[nrow(datamat_c),],datamat_te)
f4=forecast_nod(VAR_4,ddatamat_c[nrow(ddatamat_c),],datamat_c[nrow(datamat_c),],datamat_te,tempk = glok)
#################################################
#VAR5 fit every pop
#f5 = fitsvar_sep(datar,datate)
#our method
#Model6
# add age weighted model
# weight age only need to multiply before the X
#datapre
#ws_c = as.matrix(read.csv('C://Users//greatknee//Desktop//Mortality//material_VAR//centroidd.csv',header = T))
#ws_c = as.matrix(read.csv('C://Users//greatknee//Desktop//Mortality//material_VAR//centroidd2.csv',header = T))
tempfilepath = paste('C:/Users/greatknee/Desktop/Mortality/Reproducibility/Replication_SWVAR/output/',datagroup,'_',set,'_',gen_fdm,'.csv',sep = '')
#weight = read.csv(tempfilepath)
read_csv_if_exists <- function(tempfilepath) {
  if (file.exists(tempfilepath)) {
    return(read.csv(tempfilepath))
  } else {
    return(NULL)
  }
}
weight = read_csv_if_exists(tempfilepath)
f6 = fitswvar_lack_global(datar = datar,group = datagroup,datate = datate1,coulist = coulist,best = weight,gd = TRUE)
write.csv(t(f6$best),file =tempfilepath)

####################################################


conclu = matrix(0,2,5)
conclu[1,] = c(f0$RMSFE,fdm$summary$rmsfe,VAR_1$rmsfe,f4$rmsfe,f6$rmsfe)
conclu[2,] = c(f0$mafe,fdm$summary$mafe,VAR_1$mafe,f4$mae,f6$mae)

rownames(conclu) = c('RMSFE','MAFE')
colnames(conclu) = c('Li-Lee','FDM','STAR','SVAR','SWVAR')
conclu