setwd('h:/Projects/FIRO/firo_syn-forecast_stochastic-CC/')

library(extRemes)
library(fGarch)
library(future.apply)
plan(multisession)

stoch_data<-read.csv('./data/sacsma_ORO_wgen_baseline.csv')
stoch_data_cc<-read.csv('./data/sacsma_ORO_wgen_4c_7cc_0mn.csv')
#obs_data<-read.csv('./data/observed_flows.csv')
obs_data<-read.table('./data/FNF_ORO_mm.txt')

obs_fit<-read.table('./data/FNF_ORO_mm.txt')
sim_fit<-read.table('./data/ORO_SACSMA.txt')

source('./src/mm-cfs_conversion.R')
source('./src/helper-functions.R')
mm_to_kcfs <- area_mi2_ORO*mm_to_cfs/1000


ix_stoch<-seq(as.Date('2001-01-01'),as.Date('3008-01-08'),'day')
ix2_stoch<-as.POSIXlt(ix_stoch)

#ix_obs<-seq(as.Date('1985-10-02'),as.Date('2019-09-30'),'day')
#ix2_obs<-as.POSIXlt(ix_obs)

ix_obs<-seq(as.Date('1987-10-01'),as.Date('2019-09-30'),'day')
ix2_obs<-as.POSIXlt(ix_obs)

ix_sim<-seq(as.Date('1987-10-01'),as.Date('2018-09-30'),'day')
ix2_sim<-as.POSIXlt(ix_sim)


wy_fun<-function(date_vec){
  wy_vec <- date_vec$year
  wy_vec[date_vec$mo%in%c(9,10,11)] <- wy_vec[date_vec$mo%in%c(9,10,11)]+1
  date_vec_wy <- date_vec
  date_vec_wy$year <- wy_vec
  return(date_vec_wy)
}

ix2_stoch_wy <- wy_fun(ix2_stoch)
ix2_obs_wy <- wy_fun(ix2_obs)
ix2_sim_wy <- wy_fun(ix2_sim)

base_q <- as.numeric(stoch_data[,'q_mm']) * mm_to_kcfs
cc_q <- as.numeric(stoch_data_cc[,'q_mm']) * mm_to_kcfs
#obs <- obs_data[,'ORDC1']
obs <- obs_data[,'V4'] * mm_to_kcfs
base_q[base_q<0]<-0;cc_q[cc_q<0]<-0;obs[obs<0]<-0

ann_max_obs<-c()
#obs_yrs <- 1986:2019
obs_yrs <- 1988:2019
for(i in 1:length(obs_yrs)){
  ann_max_obs[i]<-max(obs[ix2_obs_wy$year%in%(obs_yrs[i]-1900)])
}


#SWM approach for SAC-SMA simulations
q_fit <- obs_fit[,'V4'][ix_obs%in%ix_sim] * mm_to_kcfs
s_fit <- sim_fit[,'V4'] * mm_to_kcfs
q_fit[q_fit<=0]<-min(q_fit[q_fit>0]);s_fit[s_fit<=0]<-min(q_fit[q_fit>0])

hy_fit<-hyb_loess_fit(s_fit,q_fit)
cbias<-hyb_loess_out(hy_fit[[1]],hy_fit[[2]],hy_fit[[3]],hy_fit[[4]],s_fit)
cbias_sim<-hyb_loess_out(hy_fit[[1]],hy_fit[[2]],hy_fit[[3]],hy_fit[[4]],base_q)
cbias_sim[is.na(cbias_sim)==T]<-base_q[is.na(cbias_sim)==T]

#plot(s_fit,q_fit,xlim=c(0,500),ylim=c(0,500))
#points(s_fit,cbias,col='red')
#lines(1:500,hyb_loess_out(hy_fit[[1]],hy_fit[[2]],hy_fit[[3]],hy_fit[[4]],1:500) )

#plot(s_fit,q_fit,xlim=c(0,200),ylim=c(0,200))
#points(s_fit,cbias,col='red')
#lines(1:500,hyb_loess_out(hy_fit[[1]],hy_fit[[2]],hy_fit[[3]],hy_fit[[4]],1:500) )

#cbias<-s_fit
#cbias_sim<-base_q

lambda<-log(cbias/q_fit)

bcf_lambda<-bcf_fun(lambda)
hist(lambda)

ar_fit<-arima(lambda,order=c(3,0,0),method='CSS',include.mean = F)

at <- ar_fit$residuals

gl_par <- sgedFit(at[-c(which(ix2_sim_wy$year==97))])

#hist(at,xlim=c(-10,10),breaks=seq(-100,100,0.5),freq = F)
#lines(seq(-10,10,0.1),dsged(seq(-10,10,0.1),gl_par$par[1],gl_par$par[2],gl_par$par[3],gl_par$par[4]),col='red',lwd=2)


#generate
##at_fit<-rsged(length(ix_sim),gl_par$par[1],gl_par$par[2],gl_par$par[3],gl_par$par[4])
##at_sim<-rsged(length(ix_stoch),gl_par$par[1],gl_par$par[2],gl_par$par[3],gl_par$par[4])
#kk=round(sqrt(length(at)))
kk=1100

kk_fun<-function(x){out<-(x/max(base_q))*kk+10;return(round(out))}
#tot <- sum(rep(1,kk) / 1:kk)
#wts <- rep(1,kk) / 1:kk / rep(tot,kk)
wts = 1/(rep(kk,kk))

at_fit <- future_sapply(1:length(ix_sim),function(x){at[sample(order(abs(s_fit[x]-s_fit))[1:kk],1,prob=wts)]},future.seed=NULL)
at_sim <- future_sapply(1:length(ix_stoch),function(x){at[sample(order(abs(base_q[x]-s_fit))[1:kk],1,prob=wts)]},future.seed=NULL)

#at_fit<-sample(at,size=length(at),replace = T)
#at_sim<-sample(at,size=length(ix_stoch),replace = T)

ar_out_fit<-arima.sim(n=length(ix_sim),list(ar=ar_fit$coef),innov = at_fit,n.start=3,start.innov=sample(at_fit,3))
ar_out_sim<-arima.sim(n=length(ix_stoch),list(ar=ar_fit$coef),innov = at_sim,n.start=3,start.innov=sample(at_sim,3))

synq_fit<-cbias/exp(ar_out_fit)*bcf_lambda;synq_fit[synq_fit<0]<-0
synq_sim<-cbias_sim/exp(ar_out_sim)*bcf_lambda;synq_sim[synq_sim<0]<-0

samps = 100

synq_sim_arr<-array(NA,c(samps,length(ix_stoch)))

for(i in 1:samps){
  #at_sim<-rsged(length(ix_stoch),gl_par$par[1],gl_par$par[2],gl_par$par[3],gl_par$par[4])
  #at_sim<-sample(at,size=length(ix_stoch),replace = T)
  at_sim <- future_sapply(1:length(ix_stoch),function(x){at[sample(order(abs(base_q[x]-s_fit))[1:kk],1,prob=wts)]},future.seed=NULL)
  ar_out_sim<-arima.sim(n=length(ix_stoch),list(ar=ar_fit$coef),innov = at_sim,n.start=3,start.innov=sample(at_sim,3))
  synq_sim<-cbias_sim/exp(ar_out_sim)*bcf_lambda;synq_sim[synq_sim<0]<-0
  synq_sim_arr[i,]<-synq_sim
}

#SWM plots
ann_max_stoch<-c()
stoch_yrs <- 2002:3008
for(i in 1:length(stoch_yrs)){
  ann_max_stoch[i]<-max(synq_sim[ix2_stoch_wy$year%in%(stoch_yrs[i]-1900)])
}

ann_max_stoch_arr<-array(NA,c(samps,length(stoch_yrs)))
for(k in 1:samps){
  for(i in 1:length(stoch_yrs)){
    ann_max_stoch_arr[k,i]<-max(synq_sim_arr[k,ix2_stoch_wy$year%in%(stoch_yrs[i]-1900)])
  }
}

ann_max_sfit<-c()
#obs_yrs <- 1986:2019
fit_yrs <- 1988:2018
for(i in 1:length(fit_yrs)){
  ann_max_sfit[i]<-max(s_fit[ix2_sim_wy$year%in%(fit_yrs[i]-1900)])
}

ann_max_sfit_swm<-c()
#obs_yrs <- 1986:2019
fit_yrs <- 1988:2018
for(i in 1:length(fit_yrs)){
  ann_max_sfit_swm[i]<-max(synq_fit[ix2_sim_wy$year%in%(fit_yrs[i]-1900)])
}

#par(mfrow=c(1,3))
#hist(ann_max_obs,xlim=c(0,350),breaks=seq(0,10000,25),freq=F,xlab='flow (kcfs)',ylim=c(0,0.015))
#hist(ann_max_sfit,xlim=c(0,350),breaks=seq(0,10000,25),col='green',freq=F,xlab='flow (kcfs)',ylim=c(0,0.015))
#hist(ann_max_sfit_swm,xlim=c(0,350),breaks=seq(0,10000,25),col='orange',freq=F,xlab='flow (kcfs)',ylim=c(0,0.015))

#par(mfrow=c(1,3))
#hist(ann_max_obs,xlim=c(0,600),breaks=seq(0,10000,25),freq=F,xlab='flow (kcfs)',ylim=c(0,0.015))
#hist(ann_max_sfit_swm,xlim=c(0,600),breaks=seq(0,10000,25),col='green',freq=F,xlab='flow (kcfs)',ylim=c(0,0.015))
#hist(ann_max_stoch,xlim=c(0,600),breaks=seq(0,10000,25),col='orange',freq=F,xlab='flow (kcfs)',ylim=c(0,0.015))

obs_gev_fit<-fevd(ann_max_obs,type='GEV')
obs_gev_est<-qevd(seq(0.01,0.999,.001),loc=obs_gev_fit$results$par[1],scale =obs_gev_fit$results$par[2],shape = obs_gev_fit$results$par[3],type='GEV')

#par(mfrow=c(1,1))

#plot(1/(1-seq(0.01,0.999,.001)),obs_gev_est,ylim=c(0,1200),log='x',type='l',col='red',axes=F,xlab='Return Period (years)',ylab='Flow (kcfs)',main='Obs - WGEN_baseline comparison')
#points(1/(length(stoch_yrs):1/(length(stoch_yrs))),sort(ann_max_stoch),col='gray75',pch=16)
#points(1/(length(fit_yrs):1/(length(fit_yrs))),sort(ann_max_sfit),col='blue',pch=15,cex=1)
#points(1/(length(fit_yrs):1/(length(fit_yrs))),sort(ann_max_sfit_swm),col='magenta',pch=16,cex=1)
#points(1/(length(obs_yrs):1/(length(obs_yrs))),sort(ann_max_obs),col='red',pch=17,cex=1)
#axis(1,at=c(1,10,100,1000),labels=c(1,10,100,1000))
#axis(2,at=seq(0,1200,100),labels=seq(0,1200,100))
##points((1/(1-seq(0.01,0.999,.001)))[which.min(abs(obs_gev_est-max(ann_max_obs)))],max(ann_max_obs),col='red',pch=17,cex=2)
#legend(1,1200,c('Obs-GEV'),col=c('red'),lwd=c(1,1),bty='n')
#legend(1,1100,c('WGEN/SAC-SMA+SWM-empirical'),col=c('gray75'),pch=16,bty='n')
#legend(1,1000,c('Obs-empirical'),col=c('red'),pch=17,bty='n')
#legend(1,900,c('SAC-SMA-empirical'),col=c('blue'),pch=15,bty='n')
#legend(1,800,c('SAC-SMA+SWM-empirical'),col=c('magenta'),pch=16,bty='n')

plot(1/(1-seq(0.01,0.999,.001)),obs_gev_est,ylim=c(0,1000),log='x',type='l',col='red',axes=F,xlab='Return Period (years)',ylab='Flow (kcfs)',main='Obs - WGEN_baseline+SWM comparison')
for(i in 1:samps){
  points(1/(length(stoch_yrs):1/(length(stoch_yrs))),sort(ann_max_stoch_arr[i,]),col='gray75',pch=16)
  #lines(1/(length(stoch_yrs):1/(length(stoch_yrs))),sort(ann_max_stoch_arr[i,]),col='gray75')
  }
lines(1/(1-seq(0.01,0.999,.001)),obs_gev_est,col='red',lwd=2)
points(1/(length(obs_yrs):1/(length(obs_yrs))),sort(ann_max_obs),col='red',pch=17,cex=1)
lines(c(2,5,10,20,50,100,200,500,1000),c(45,90,140,180,245,295,350,430,500),lty=2,lwd=2,col='gray40')
axis(1,at=c(1,10,100,1000),labels=c(1,10,100,1000))
axis(2,at=seq(0,1500,100),labels=seq(0,1500,100))
#points((1/(1-seq(0.01,0.999,.001)))[which.min(abs(obs_gev_est-max(ann_max_obs[-c(which.max(ann_max_obs))])))],max(ann_max_obs[-c(which.max(ann_max_obs))]),col='red',pch=17,cex=2)
legend(1,1000,c('Obs-GEV'),col=c('red'),lwd=c(1,1),bty='n')
legend(1,900,c('WGEN/SAC-SMA+SWM-empirical'),col=c('gray75'),pch=16,bty='n')
legend(1,800,c('Obs-empirical'),col=c('red'),pch=17,bty='n')
legend(1,700,c('PVA-estimated'),col=c('gray40'),lty=2,lwd=2,bty='n')

#plot(1/(1-seq(0.01,0.999,.001)),obs_gev_est,ylim=c(0,350),xlim=c(1,100),log='x',type='l',col='red',axes=F,xlab='Return Period (years)',ylab='Flow (kcfs)',main='Obs - WGEN_baseline comparison')
#points(1/(length(stoch_yrs):1/(length(stoch_yrs))),sort(ann_max_stoch),col='gray75',pch=16)
#points(1/(length(fit_yrs):1/(length(fit_yrs))),sort(ann_max_sfit),col='blue',pch=15,cex=1)
#points(1/(length(fit_yrs):1/(length(fit_yrs))),sort(ann_max_sfit_swm),col='magenta',pch=16,cex=1)
#points(1/(length(obs_yrs):1/(length(obs_yrs))),sort(ann_max_obs),col='red',pch=17,cex=1)
#lines(c(2,5,10,20,50,100,200,500,1000),c(45,90,140,180,245,295,350,430,500),lty=2,lwd=2,col='gray40')
#axis(1,at=c(1,10,100,1000),labels=c(1,10,100,1000))
#axis(2,at=seq(0,1200,100),labels=seq(0,1200,100))
##points((1/(1-seq(0.01,0.999,.001)))[which.min(abs(obs_gev_est-max(ann_max_obs)))],max(ann_max_obs),col='red',pch=17,cex=2)
#legend(1,350,c('Obs-GEV'),col=c('red'),lwd=c(1,1),bty='n')
#legend(1,325,c('WGEN/SAC-SMA+SWM-empirical'),col=c('gray75'),pch=16,bty='n')
#legend(1,300,c('Obs-empirical'),col=c('red'),pch=17,bty='n')
#legend(1,275,c('SAC-SMA-empirical'),col=c('blue'),pch=15,bty='n')
#legend(1,250,c('SAC-SMA+SWM-empirical'),col=c('magenta'),pch=16,bty='n')
#legend(1,225,c('PVA-estimated'),col=c('gray40'),lty=2,lwd=2,bty='n')

plot(1/(1-seq(0.01,0.999,.001)),obs_gev_est,ylim=c(0,350),xlim=c(1,100),log='x',type='l',col='red',axes=F,xlab='Return Period (years)',ylab='Flow (kcfs)',main='Obs - WGEN_baseline+SWM comparison')
for(i in 1:samps){
  points(1/(length(stoch_yrs):1/(length(stoch_yrs))),sort(ann_max_stoch_arr[i,]),col='gray75',pch=16)
  #lines(1/(length(stoch_yrs):1/(length(stoch_yrs))),sort(ann_max_stoch_arr[i,]),col='gray75')
}
lines(1/(1-seq(0.01,0.999,.001)),obs_gev_est,col='red',lwd=2)
points(1/(length(obs_yrs):1/(length(obs_yrs))),sort(ann_max_obs),col='red',pch=17,cex=1)
lines(c(2,5,10,20,50,100,200,500,1000),c(45,90,140,180,245,295,350,430,500),lty=2,lwd=2,col='gray40')
#lines(c(2,5,10,20,50,100,200,500,1000),qlpearsonIII(1-(1/c(2,5,10,20,50,100,200,500,1000)), meanlog = 4.638, sdlog = 0.393s, skew = -0.3),col='green',lty=3,lwd=2)
axis(1,at=c(1,10,100,1000),labels=c(1,10,100,1000))
axis(2,at=seq(0,1500,100),labels=seq(0,1500,100))
#points((1/(1-seq(0.01,0.999,.001)))[which.min(abs(obs_gev_est-max(ann_max_obs[-c(which.max(ann_max_obs))])))],max(ann_max_obs[-c(which.max(ann_max_obs))]),col='red',pch=17,cex=2)
legend(1,325,c('Obs-GEV'),col=c('red'),lwd=c(1,1),bty='n')
legend(1,300,c('WGEN/SAC-SMA+SWM-empirical'),col=c('gray75'),pch=16,bty='n')
legend(1,275,c('Obs-empirical'),col=c('red'),pch=17,bty='n')
legend(1,250,c('PVA-estimated'),col=c('gray40'),lty=2,lwd=2,bty='n')

qlpearsonIII(1-(1/c(2,5,10,20,50,100,200,500,1000)), meanlog = 4.638, sdlog = 0.393, skew = -0.3)
#plot event profiles similar to 1997
evt_idx = which.max(obs)

#plot(-15:15,obs[(evt_idx-15):(evt_idx+15)],type='l',lwd=4,ylim=c(0,500),xlab='days from extreme',ylab='flow (kcfs)')
#for(i in 1:10){
  #df <- abs(synq_sim_arr[i,]-obs[evt_idx])
  #evt_ix <- which.min(df)
  #lines(-15:15,synq_sim_arr[i,(evt_ix-15):(evt_ix+15)],col=i,lwd=1)
#}
#legend('topright',paste('samp',1:10),col=c(1:10),lwd=1)
#legend('topleft','1997 event',col='black',lwd=4)

#par(mfrow=c(2,3),mar=c(2,2,1,1))
#for(i in 1:6){
  #df <- abs(synq_sim_arr[i,]-obs[evt_idx])
  #evt_ix <- which.min(df)
  #plot(-15:15,synq_sim_arr[i,(evt_ix-15):(evt_ix+15)],type='l',col=i+1,lwd=3,xlab='',ylab='',ylim=c(0,350))
  #abline(h=max(obs),col='gray75',lty=2)
  #text(10,max(obs)+15,'1997',col='gray75',cex=1.5)
#}

#par(mfrow=c(1,1),mar=c(4,4,1,1))
#plot(-15:15,obs[(evt_idx-15):(evt_idx+15)],type='l',lwd=4,ylim=c(0,350),xlab='days from extreme',ylab='flow (kcfs)')

#plot events with similar return periods to 1997
yr_evt <- 1/(1-seq(0.01,0.999,.001)[which.min(abs(max(obs)-obs_gev_est))])
yr_idx <- which.min(abs(1/(length(stoch_yrs):1/(length(stoch_yrs)))-yr_evt))
1/(length(stoch_yrs):1/(length(stoch_yrs)))[yr_idx]

#plot(-15:15,obs[(evt_idx-15):(evt_idx+15)],type='l',lwd=4,ylim=c(0,500),xlab='days from extreme',ylab='flow (kcfs)')
#for(i in 1:10){
  #evt = sort(ann_max_stoch_arr[i,])[yr_idx]
  #event_idx = which(synq_sim_arr[i,]==evt)
  #lines(-15:15,synq_sim_arr[i,(event_idx-15):(event_idx+15)],col=i,lwd=1)
#}
#legend('topright',paste('samp',1:10),col=c(1:10),lwd=1)
#legend('topleft','1997 event',col='black',lwd=4)

#stochastic events of a specified return period
yr_evt <- 500
yr_idx <- which.min(abs(1/(length(stoch_yrs):1/(length(stoch_yrs)))-yr_evt))
1/(length(stoch_yrs):1/(length(stoch_yrs)))[yr_idx]

ymax=500

#plot(-15:15,obs[(evt_idx-15):(evt_idx+15)],type='l',col='white',lwd=4,ylim=c(0,ymax),xlab='days from extreme',ylab='flow (kcfs)')
#for(i in 1:10){
  #evt = sort(ann_max_stoch_arr[i,])[yr_idx]
  #event_idx = which(synq_sim_arr[i,]==evt)
  #lines(-15:15,synq_sim_arr[i,(event_idx-15):(event_idx+15)],col=i,lwd=1)
#}
#legend('topright',paste('samp',1:10),col=c(1:10),lwd=1)
#title(paste(yr_evt,'year event'))

#par(mfrow=c(2,3),mar=c(2,2,1,1))
#for(i in 1:6){
  #evt = sort(ann_max_stoch_arr[i,])[yr_idx]
  #event_idx = which(synq_sim_arr[i,]==evt)
  #plot(-15:15,synq_sim_arr[i,(event_idx-15):(event_idx+15)],type='l',col=i+1,lwd=3,xlab='',ylab='',ylim=c(0,ymax))
  #abline(h=max(obs),col='gray75',lty=2)
  #text(10,max(obs)+15,'1997',col='gray75',cex=1.5)
#}

################################################END######################################################
