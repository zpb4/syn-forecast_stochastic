setwd('h:/Projects/FIRO/firo_syn-forecast_stochastic-CC/')

library(extRemes)

stoch_data<-read.table('./data/ORO_EF-Gupta_baseline.txt')
obs_data<-read.csv('./data/observed_flows.csv')
source('./src/mm-cfs_conversion.R')
mm_to_kcfs <- area_mi2_ORO*mm_to_cfs/1000


ix_stoch<-seq(as.Date('1399-11-01'),as.Date('2017-04-30'),'day')
ix2_stoch<-as.POSIXlt(ix_stoch)

ix_obs<-seq(as.Date('1985-10-02'),as.Date('2019-09-30'),'day')
ix2_obs<-as.POSIXlt(ix_obs)

wy_fun<-function(date_vec){
  wy_vec <- date_vec$year
  wy_vec[date_vec$mo%in%c(9,10,11)] <- wy_vec[date_vec$mo%in%c(9,10,11)]+1
  date_vec_wy <- date_vec
  date_vec_wy$year <- wy_vec
  return(date_vec_wy)
}

ix2_stoch_wy <- wy_fun(ix2_stoch)
ix2_obs_wy <- wy_fun(ix2_obs)

ens_1 <- as.numeric(stoch_data[,'V10'][-c(1)]) * mm_to_kcfs
obs <- obs_data[,'ORDC1']

ens_max<-c()
idx <- 7:56
for(i in 1:50){
  ens_max[i]<-max(as.numeric(stoch_data[,paste('V',idx[i],sep='')][-c(1)]) * mm_to_kcfs)
}


ens_1 <- as.numeric(stoch_data[,paste('V',which.max(ens_max)+6,sep='')][-c(1)]) * mm_to_kcfs
mean(ens_1)
max(ens_1)

ann_max_stoch<-c()
stoch_yrs <- 1401:2016
for(i in 1:length(stoch_yrs)){
  ann_max_stoch[i]<-max(ens_1[ix2_stoch_wy$year%in%(stoch_yrs[i]-1900)])
}

ann_max_obs<-c()
obs_yrs <- 1986:2019
for(i in 1:length(obs_yrs)){
  ann_max_obs[i]<-max(obs[ix2_obs_wy$year%in%(obs_yrs[i]-1900)])
}

hist(ann_max_obs,xlim=c(0,350),breaks=seq(0,350,25))
hist(ann_max_stoch,xlim=c(0,350),breaks=seq(0,350,25))


dsgn_prob<-1 - (1/dsgn)

gev_fun<-function(x){
  fit<-fevd(x,type='GEV')
  dsn<-qevd(dsgn_prob,loc=fit$results$par[1],scale=fit$results$par[2],shape=fit$results$par[3],type='GEV')
  return(dsn)}

obs_gev_fit<-fevd(ann_max_obs,type='GEV')
obs_gev_est<-qevd(seq(0.01,0.999,.001),loc=obs_gev_fit$results$par[1],scale =obs_gev_fit$results$par[2],shape = obs_gev_fit$results$par[3],type='GEV')

plot(1/(1-seq(0.01,0.999,.001)),obs_gev_est,ylim=c(0,1200),log='x',type='l',col='red',axes=F,xlab='Return Period (years)',ylab='Flow (kcfs)',main='Obs - WGEN comparison: 1997 included')
points(1/(length(stoch_yrs):1/(length(stoch_yrs))),sort(ann_max_stoch))
points(1/(length(obs_yrs):1/(length(obs_yrs))),sort(ann_max_obs),col='red',pch=17,cex=1)
axis(1,at=c(1,10,100,1000),labels=c(1,10,100,1000))
axis(2,at=seq(0,1200,100),labels=seq(0,1200,100))
#points((1/(1-seq(0.01,0.999,.001)))[which.min(abs(obs_gev_est-max(ann_max_obs)))],max(ann_max_obs),col='red',pch=17,cex=2)
legend(1,1200,c('Obs-GEV'),col=c('red'),lwd=c(1,1),bty='n')
legend(1,1100,c('WGEN-empirical'),col=c('black'),pch=1,bty='n')
legend(1,1000,c('largest event - 1997'),col=c('red'),pch=17,bty='n')

obs_gev_fit<-fevd(ann_max_obs[-c(which.max(ann_max_obs))],type='GEV')
obs_gev_est<-qevd(seq(0.01,0.999,.001),loc=obs_gev_fit$results$par[1],scale =obs_gev_fit$results$par[2],shape = obs_gev_fit$results$par[3],type='GEV')

plot(1/(1-seq(0.01,0.999,.001)),obs_gev_est,ylim=c(0,1200),log='x',type='l',col='red',axes=F,xlab='Return Period (years)',ylab='Flow (kcfs)',main='Obs - WGEN comparison: 1997 removed')
points(1/(length(stoch_yrs):1/(length(stoch_yrs))),sort(ann_max_stoch))
points(1/((length(obs_yrs)-1):1/((length(obs_yrs)-1))),sort(ann_max_obs[-which.max(ann_max_obs)]),col='red',pch=17,cex=1)
axis(1,at=c(1,10,100,1000),labels=c(1,10,100,1000))
axis(2,at=seq(0,1200,100),labels=seq(0,1200,100))
#points((1/(1-seq(0.01,0.999,.001)))[which.min(abs(obs_gev_est-max(ann_max_obs[-c(which.max(ann_max_obs))])))],max(ann_max_obs[-c(which.max(ann_max_obs))]),col='red',pch=17,cex=2)
legend(1,1200,c('Obs-GEV'),col=c('red'),lwd=c(1,1),bty='n')
legend(1,1100,c('WGEN-empirical'),col=c('black'),pch=1,bty='n')
legend(1,1000,c('largest event - 1986'),col=c('red'),pch=17,bty='n')

################################################END######################################################
