setwd('h:/Projects/FIRO/firo_syn-forecast_stochastic-CC/')

library(extRemes)
library(fGarch)

stoch_data<-read.csv('./data/sacsma_ORO_wgen_baseline.csv')
stoch_data_cc<-read.csv('./data/sacsma_ORO_wgen_4c_7cc_0mn.csv')
#obs_data<-read.csv('./data/observed_flows.csv')
obs_data<-read.table('./data/FNF_ORO_mm.txt')

obs_fit<-read.table('./data/FNF_ORO_mm.txt')
sim_fit<-read.table('./data/ORO_SACSMA.txt')
obs_inp<-read.table('./data/sacsma_ORO_statvars.txt')

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

ix_inp<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2_inp<-as.POSIXlt(ix_inp)

#precip to streamflow comparison
precip <- obs_inp[,'V4']
sma_q <- sim_fit[,'V4'][ix_sim%in%ix_inp]
obs_q <- obs_data[,'V4'][ix_obs%in%ix_inp]

precip3 <- c()
for(i in 1:length(precip)){
  precip3[i]<-sum(precip[max(0,i-2):i])
}

ids <- order(obs_q,decreasing = T)[1:100]

trk <- rep(1,100)
for(i in 1:length(ids)){
  if(trk[i]==1){
    df <- abs(ids-ids[i])
    trk[df<10 & df>0]<-0}
  else{print(i)}
}

keep_ids<-ids[trk==1][1:12]
keep_dates<-ix_inp[keep_ids]

par(mfrow=c(4,3))
for(i in 1:12){
  plot(-10:5,precip[(keep_ids[i]-10):(keep_ids[i]+5)],type='l',col='blue',xlab='days from event',ylab='mm',lwd=2)
  lines(-10:5,sma_q[(keep_ids[i]-10):(keep_ids[i]+5)],col='green',lwd=2)
  lines(-10:5,obs_q[(keep_ids[i]-10):(keep_ids[i]+5)],col='black',lwd=2)
}

par(mfrow=c(1,2))
plot(-10:5,precip[(keep_ids[1]-10):(keep_ids[1]+5)],type='l',col='blue',xlab='days from event',ylab='mm',lwd=2)
lines(-10:5,sma_q[(keep_ids[1]-10):(keep_ids[1]+5)],col='green',lwd=2)
lines(-10:5,obs_q[(keep_ids[1]-10):(keep_ids[1]+5)],col='black',lwd=2)
legend('topleft',c('Precip','Obs','Sim'),lwd=c(2,2,2),col=c('blue','black','green'))

plot(-10:5,precip3[(keep_ids[1]-10):(keep_ids[1]+5)],type='l',col='blue',xlab='days from event',ylab='mm',lwd=2)
lines(-10:5,sma_q[(keep_ids[1]-10):(keep_ids[1]+5)],col='green',lwd=2)
lines(-10:5,obs_q[(keep_ids[1]-10):(keep_ids[1]+5)],col='black',lwd=2)
legend('topleft',c('Precip-3d sum','Obs','Sim'),lwd=c(2,2,2),col=c('blue','black','green'))

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

ann_max_stoch<-c()
stoch_yrs <- 2002:3008
for(i in 1:length(stoch_yrs)){
  ann_max_stoch[i]<-max(base_q[ix2_stoch_wy$year%in%(stoch_yrs[i]-1900)])
}

ann_max_stoch_cc<-c()
for(i in 1:length(stoch_yrs)){
  ann_max_stoch_cc[i]<-max(cc_q[ix2_stoch_wy$year%in%(stoch_yrs[i]-1900)])
}

ann_max_obs<-c()
#obs_yrs <- 1986:2019
obs_yrs <- 1988:2019
for(i in 1:length(obs_yrs)){
  ann_max_obs[i]<-max(obs[ix2_obs_wy$year%in%(obs_yrs[i]-1900)])
}

par(mfrow=c(1,3))
hist(ann_max_obs,xlim=c(0,350),breaks=seq(0,350,25),freq=F,ylim=c(0,0.015),xlab='flow (kcfs)')
hist(ann_max_stoch,xlim=c(0,350),breaks=seq(0,350,25),col='green',freq=F,ylim=c(0,0.015),xlab='flow (kcfs)')
hist(ann_max_stoch_cc,xlim=c(0,350),breaks=seq(0,350,25),col='orange',freq=F,ylim=c(0,0.015),xlab='flow (kcfs)')

obs_gev_fit<-fevd(ann_max_obs,type='GEV')
obs_gev_est<-qevd(seq(0.01,0.999,.001),loc=obs_gev_fit$results$par[1],scale =obs_gev_fit$results$par[2],shape = obs_gev_fit$results$par[3],type='GEV')

par(mfrow=c(1,2))

plot(1/(1-seq(0.01,0.999,.001)),obs_gev_est,ylim=c(0,1200),log='x',type='l',col='red',axes=F,xlab='Return Period (years)',ylab='Flow (kcfs)',main='Obs - WGEN_baseline comparison')
points(1/(length(stoch_yrs):1/(length(stoch_yrs))),sort(ann_max_stoch),col='gray75')
lines(1/(1-seq(0.01,0.999,.001)),obs_gev_est,col='red',lwd=2)
points(1/(length(obs_yrs):1/(length(obs_yrs))),sort(ann_max_obs),col='red',pch=17,cex=1)
axis(1,at=c(1,10,100,1000),labels=c(1,10,100,1000))
axis(2,at=seq(0,1200,100),labels=seq(0,1200,100))
#points((1/(1-seq(0.01,0.999,.001)))[which.min(abs(obs_gev_est-max(ann_max_obs)))],max(ann_max_obs),col='red',pch=17,cex=2)
legend(1,1200,c('Obs-GEV'),col=c('red'),lwd=2,bty='n')
legend(1,1100,c('WGEN-empirical'),col=c('gray75'),pch=1,bty='n')
legend(1,1000,c('Obs-empirical'),col=c('red'),pch=17,bty='n')

plot(1/(1-seq(0.01,0.999,.001)),obs_gev_est,ylim=c(0,1200),log='x',type='l',col='red',axes=F,xlab='Return Period (years)',ylab='Flow (kcfs)',main='Obs - WGEN_4c_7%cc comparison')
points(1/(length(stoch_yrs):1/(length(stoch_yrs))),sort(ann_max_stoch_cc),col='gray75')
lines(1/(1-seq(0.01,0.999,.001)),obs_gev_est,col='red',lwd=2)
points(1/(length(obs_yrs):1/(length(obs_yrs))),sort(ann_max_obs),col='red',pch=17,cex=1)
axis(1,at=c(1,10,100,1000),labels=c(1,10,100,1000))
axis(2,at=seq(0,1200,100),labels=seq(0,1200,100))
#points((1/(1-seq(0.01,0.999,.001)))[which.min(abs(obs_gev_est-max(ann_max_obs[-c(which.max(ann_max_obs))])))],max(ann_max_obs[-c(which.max(ann_max_obs))]),col='red',pch=17,cex=2)
legend(1,1200,c('Obs-GEV'),col=c('red'),lwd=2,bty='n')
legend(1,1100,c('WGEN-empirical'),col=c('gray75'),pch=1,bty='n')
legend(1,1000,c('Obs-empirical'),col=c('red'),pch=17,bty='n')



#SWM approach for SAC-SMA simulations
q_fit <- obs_fit[,'V4'][ix_obs%in%ix_sim] * mm_to_kcfs
s_fit <- sim_fit[,'V4'] * mm_to_kcfs
q_fit[q_fit<0]<-0;s_fit[s_fit<0]<-0

hy_fit<-hyb_loess_fit(s_fit,q_fit)
cbias<-hyb_loess_out(hy_fit[[1]],hy_fit[[2]],hy_fit[[3]],hy_fit[[4]],s_fit)
cbias_sim<-hyb_loess_out(hy_fit[[1]],hy_fit[[2]],hy_fit[[3]],hy_fit[[4]],base_q)
cbias_sim[is.na(cbias_sim)==T]<-base_q[is.na(cbias_sim)==T]
err<-q_fit-cbias

#err <- s_fit - q_fit

err_abs <- abs(err)

sds<-mean(err_abs)

ls_fit<-loess(err_abs~cbias,span=0.75,degree = 2,family='gaussian',control=loess.control(surface='direct'))

norm_vec<-predict(ls_fit,cbias)
norm_vec[norm_vec<0]<-sds

norm_resids<-err / norm_vec

ar_fit<-arima(norm_resids,order=c(3,0,0),method='CSS',include.mean = F)

at <- ar_fit$residuals

gl_par <- sgedFit(at)

hist(at,xlim=c(-10,10),breaks=seq(-100,100,0.5),freq = F)
lines(seq(-10,10,0.1),dsged(seq(-10,10,0.1),gl_par$par[1],gl_par$par[2],gl_par$par[3],gl_par$par[4]),col='red',lwd=2)


#generate
norm_val<-predict(ls_fit,cbias)
norm_val_sim<-predict(ls_fit,cbias_sim)
norm_val[norm_val<=0]<-sds;norm_val_sim[norm_val_sim<=0]<-sds

at_fit<-rsged(length(ix_sim),gl_par$par[1],gl_par$par[2],gl_par$par[3],gl_par$par[4])
at_sim<-rsged(length(ix_stoch),gl_par$par[1],gl_par$par[2],gl_par$par[3],gl_par$par[4])

ar_out_fit<-arima.sim(n=length(ix_sim),list(ar=ar_fit$coef),innov = at_fit,n.start=3,start.innov=sample(at_fit,3))
ar_out_sim<-arima.sim(n=length(ix_stoch),list(ar=ar_fit$coef),innov = at_sim,n.start=3,start.innov=sample(at_sim,3))

syn_resids_fit<-ar_out_fit * norm_val
syn_resids_sim<-ar_out_sim * norm_val_sim

synq_fit<-cbias-syn_resids_fit;synq_fit[synq_fit<0]<-0
synq_sim<-cbias_sim-syn_resids_sim;synq_sim[synq_sim<0]<-0

samps = 10

synq_sim_arr<-array(NA,c(samps,length(ix_stoch)))

for(i in 1:samps){
  at_sim<-rsged(length(ix_stoch),gl_par$par[1],gl_par$par[2],gl_par$par[3],gl_par$par[4])
  ar_out_sim<-arima.sim(n=length(ix_stoch),list(ar=ar_fit$coef),innov = at_sim,n.start=3,start.innov=sample(at_sim,3))
  syn_resids_sim<-ar_out_sim * norm_val_sim
  synq_sim<-cbias_sim-syn_resids_sim;synq_sim[synq_sim<0]<-0
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

par(mfrow=c(1,3))
hist(ann_max_obs,xlim=c(0,350),breaks=seq(0,350,25),freq=F,xlab='flow (kcfs)',ylim=c(0,0.015))
hist(ann_max_sfit,xlim=c(0,350),breaks=seq(0,350,25),col='green',freq=F,xlab='flow (kcfs)',ylim=c(0,0.015))
hist(ann_max_sfit_swm,xlim=c(0,350),breaks=seq(0,350,25),col='orange',freq=F,xlab='flow (kcfs)',ylim=c(0,0.015))

par(mfrow=c(1,3))
hist(ann_max_obs,xlim=c(0,600),breaks=seq(0,1000,25),freq=F,xlab='flow (kcfs)',ylim=c(0,0.015))
hist(ann_max_sfit_swm,xlim=c(0,600),breaks=seq(0,1000,25),col='green',freq=F,xlab='flow (kcfs)',ylim=c(0,0.015))
hist(ann_max_stoch,xlim=c(0,600),breaks=seq(0,1000,25),col='orange',freq=F,xlab='flow (kcfs)',ylim=c(0,0.015))

obs_gev_fit<-fevd(ann_max_obs,type='GEV')
obs_gev_est<-qevd(seq(0.01,0.999,.001),loc=obs_gev_fit$results$par[1],scale =obs_gev_fit$results$par[2],shape = obs_gev_fit$results$par[3],type='GEV')

par(mfrow=c(1,1))

plot(1/(1-seq(0.01,0.999,.001)),obs_gev_est,ylim=c(0,1200),log='x',type='l',col='red',axes=F,xlab='Return Period (years)',ylab='Flow (kcfs)',main='Obs - WGEN_baseline comparison')
points(1/(length(stoch_yrs):1/(length(stoch_yrs))),sort(ann_max_stoch),col='gray75',pch=16)
points(1/(length(fit_yrs):1/(length(fit_yrs))),sort(ann_max_sfit),col='blue',pch=15,cex=1)
points(1/(length(fit_yrs):1/(length(fit_yrs))),sort(ann_max_sfit_swm),col='magenta',pch=16,cex=1)
points(1/(length(obs_yrs):1/(length(obs_yrs))),sort(ann_max_obs),col='red',pch=17,cex=1)
axis(1,at=c(1,10,100,1000),labels=c(1,10,100,1000))
axis(2,at=seq(0,1200,100),labels=seq(0,1200,100))
#points((1/(1-seq(0.01,0.999,.001)))[which.min(abs(obs_gev_est-max(ann_max_obs)))],max(ann_max_obs),col='red',pch=17,cex=2)
legend(1,1200,c('Obs-GEV'),col=c('red'),lwd=c(1,1),bty='n')
legend(1,1100,c('WGEN/SAC-SMA+SWM-empirical'),col=c('gray75'),pch=16,bty='n')
legend(1,1000,c('Obs-empirical'),col=c('red'),pch=17,bty='n')
legend(1,900,c('SAC-SMA-empirical'),col=c('blue'),pch=15,bty='n')
legend(1,800,c('SAC-SMA+SWM-empirical'),col=c('magenta'),pch=16,bty='n')

plot(1/(1-seq(0.01,0.999,.001)),obs_gev_est,ylim=c(0,1200),log='x',type='l',col='red',axes=F,xlab='Return Period (years)',ylab='Flow (kcfs)',main='Obs - WGEN_baseline+SWM comparison')
for(i in 1:samps){
  points(1/(length(stoch_yrs):1/(length(stoch_yrs))),sort(ann_max_stoch_arr[i,]),col='gray75',pch=16)
  #lines(1/(length(stoch_yrs):1/(length(stoch_yrs))),sort(ann_max_stoch_arr[i,]),col='gray75')
  }
lines(1/(1-seq(0.01,0.999,.001)),obs_gev_est,col='red',lwd=2)
points(1/(length(obs_yrs):1/(length(obs_yrs))),sort(ann_max_obs),col='red',pch=17,cex=1)
axis(1,at=c(1,10,100,1000),labels=c(1,10,100,1000))
axis(2,at=seq(0,1200,100),labels=seq(0,1200,100))
#points((1/(1-seq(0.01,0.999,.001)))[which.min(abs(obs_gev_est-max(ann_max_obs[-c(which.max(ann_max_obs))])))],max(ann_max_obs[-c(which.max(ann_max_obs))]),col='red',pch=17,cex=2)
legend(1,1200,c('Obs-GEV'),col=c('red'),lwd=c(1,1),bty='n')
legend(1,1100,c('WGEN/SAC-SMA+SWM-empirical'),col=c('gray75'),pch=16,bty='n')
legend(1,1000,c('Obs-empirical'),col=c('red'),pch=17,bty='n')

################################################END######################################################
