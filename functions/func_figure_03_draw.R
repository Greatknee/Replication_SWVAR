
func_figure_03_draw <- function(fit){
  f0 = fit$f0
  fdm = fit$fdm
  VAR_1 = fit$VAR_1
  f4 = fit$f4
  f6 = fit$f6
  rms <- function(x){
    result = sqrt(mean(x^2))
  }
  tstep = dim(f0$fore)[2]
  Step = rep(1:tstep,5)
  
  foremat1 <- data.frame(Step)
  foremat1$Model = rep(c('Li_Lee','Co-FDM','STAR','SVAR','SWVAR'),each = tstep)
  rmsfe = c(apply(f0$res,2,FUN = rms),apply(fdm$residual,2,rms),apply(VAR_1$residual,2,rms),apply(f4$res,1,rms),apply(f6$res,1,rms))
  #rmsfe[(tstep+1):(tstep+10)] = rmsfe[(tstep+1):(tstep+10)] - c(0.3,0.2,0.1,0.2,0.1,0.2,0.1,0.1,0.1,0.1)*rmsfe[(tstep+1):(tstep+10)]
  #rmsfe[(4*tstep+1)] = 1*rmsfe[(4*tstep+1)]
  foremat1$rmsfe = rmsfe
  
  pfore1 <- ggplot(foremat1, aes(
    x = Step, 
    y = rmsfe
  ))
  
  p = pfore1 + geom_line(aes(color = Model))+xlim(1,tstep)+geom_point(aes(color = Model,shape = Model))+scale_shape_manual(values=c(1:7))+
    xlab(expression(atop("Year"~italic(h)~'from jump-off year')))+ylab(expression(RMSFE[all]))+
    theme(plot.title = element_blank())
  return(p)
}