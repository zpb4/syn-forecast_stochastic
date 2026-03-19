#setwd('h:/Projects/FIRO/firo_syn-forecast_stochastic-CC/')

#library(extRemes)
#library(fGarch)
#library(future.apply)
#plan(multicore,workers=75)

kk=106
use_bcf=F
seed=1

args = commandArgs(trailingOnly=TRUE)
print(paste('task #',args[1]))
idx = as.numeric(args[1])

library(doParallel)
parallel::detectCores()
n.cores <- parallel::detectCores()
my.cluster<-parallel::makeCluster(n.cores,type = 'FORK')
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

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
obs <- obs_data[,'V4'] * mm_to_kcfs
base_q[base_q<0]<-0;cc_q[cc_q<0]<-0;obs[obs<0]<-0

#SWM approach for SAC-SMA simulations
q_fit <- obs_fit[,'V4'][ix_obs%in%ix_sim] * mm_to_kcfs
s_fit <- sim_fit[,'V4'] * mm_to_kcfs
q_fit[q_fit<=0]<-min(q_fit[q_fit>0]);s_fit[s_fit<=0]<-min(q_fit[q_fit>0])

hy_fit<-hyb_loess_fit(s_fit,q_fit)
cbias<-hyb_loess_out(hy_fit[[1]],hy_fit[[2]],hy_fit[[3]],hy_fit[[4]],s_fit)
cbias_sim<-hyb_loess_out(hy_fit[[1]],hy_fit[[2]],hy_fit[[3]],hy_fit[[4]],base_q)
cbias_sim[is.na(cbias_sim)==T]<-base_q[is.na(cbias_sim)==T]
cbias_sim_cc<-hyb_loess_out(hy_fit[[1]],hy_fit[[2]],hy_fit[[3]],hy_fit[[4]],cc_q)
cbias_sim_cc[is.na(cbias_sim_cc)==T]<-cc_q[is.na(cbias_sim_cc)==T]

lambda<-log(cbias/q_fit)

if(use_bcf==T){bcf_lambda<-bcf_fun(lambda)}
if(use_bcf==F){bcf_lambda<-1}
#hist(lambda)

ar_fit<-arima(lambda,order=c(3,0,0),method='CSS',include.mean = F)

at <- ar_fit$residuals

#generate

#tot <- sum(rep(1,kk) / 1:kk)
#wts <- rep(1,kk) / 1:kk / rep(tot,kk)
wts = 1/(rep(kk,kk)) #equal weights

#SWM generation
print(paste('swm start',Sys.time()))

samps = 100

stoch_yrs <- unique(ix2_stoch_wy$year)+1900
fit_yrs <- unique(ix2_sim_wy$year)+1900

if(idx==1){
  print(paste('swm nocc start',Sys.time()))
  synq_sim_arr <- foreach(m = 1:samps,.combine='rbind',.inorder=F) %dopar% {
    at_sim <- sapply(1:length(ix_stoch),function(x){at[sample(order(abs(base_q[x]-s_fit))[1:kk],1,prob=wts)]})
    ar_out_sim<-arima.sim(n=length(ix_stoch),list(ar=ar_fit$coef),innov = at_sim,n.start=3,start.innov=sample(at_sim,3))
    synq_sim<-cbias_sim/exp(ar_out_sim)*bcf_lambda;synq_sim[synq_sim<0]<-0
    return(synq_sim)}

  saveRDS(synq_sim_arr,paste('./out/synq_sim_arr_use-bcf=',use_bcf,'_k=',kk,'.rds',sep=''))
  
  ann_max_stoch_arr <- foreach(m = 1:samps,.combine='rbind',.inorder=F) %dopar% {
    out <- c()
    for(i in 1:length(stoch_yrs)){
      out[i]<-max(synq_sim_arr[m,ix2_stoch_wy$year%in%(stoch_yrs[i]-1900)])
    }
    return(out)}
  
  saveRDS(ann_max_stoch_arr,paste('./out/synq_annmax_use-bcf=',use_bcf,'_k=',kk,'.rds',sep=''))

  print(paste('swm nocc end',Sys.time()))
  
  print(paste('swm sim start',Sys.time()))
  synq_sim_arr <- foreach(m = 1:samps,.combine='rbind',.inorder=F) %dopar% {
    at_sim <- sapply(1:length(ix_sim),function(x){at[sample(order(abs(s_fit[x]-s_fit))[1:kk],1,prob=wts)]})
    ar_out_sim<-arima.sim(n=length(ix_sim),list(ar=ar_fit$coef),innov = at_sim,n.start=3,start.innov=sample(at_sim,3))
    synq_sim<-cbias/exp(ar_out_sim)*bcf_lambda;synq_sim[synq_sim<0]<-0
    return(synq_sim)}
  
  saveRDS(synq_sim_arr,paste('./out/synq_arr_fit_use-bcf=',use_bcf,'_k=',kk,'.rds',sep=''))
  
  ann_max_stoch_arr <- foreach(m = 1:samps,.combine='rbind',.inorder=F) %dopar% {
    out <- c()
    for(i in 1:length(fit_yrs)){
      out[i]<-max(synq_sim_arr[m,ix2_sim_wy$year%in%(fit_yrs[i]-1900)])
    }
    return(out)}
  
  saveRDS(ann_max_stoch_arr,paste('./out/synq_annmax_fit_use-bcf=',use_bcf,'_k=',kk,'.rds',sep=''))
  
  print(paste('swm no-cc end',Sys.time()))}


if(idx==2){
  print(paste('swm cc start',Sys.time()))
  synq_sim_arr_cc <- foreach(m = 1:samps,.combine='rbind',.inorder=F) %dopar% {
    at_sim <- sapply(1:length(ix_stoch),function(x){at[sample(order(abs(base_q[x]-s_fit))[1:kk],1,prob=wts)]})
    ar_out_sim<-arima.sim(n=length(ix_stoch),list(ar=ar_fit$coef),innov = at_sim,n.start=3,start.innov=sample(at_sim,3))
    synq_sim<-cbias_sim_cc/exp(ar_out_sim)*bcf_lambda;synq_sim[synq_sim<0]<-0
    return(synq_sim)}

  saveRDS(synq_sim_arr_cc,paste('./out/synq_sim_arr_cc_use-bcf=',use_bcf,'_k=',kk,'.rds',sep=''))
  
  ann_max_stoch_arr <- foreach(m = 1:samps,.combine='rbind',.inorder=F) %dopar% {
    out <- c()
    for(i in 1:length(stoch_yrs)){
      out[i]<-max(synq_sim_arr_cc[m,ix2_stoch_wy$year%in%(stoch_yrs[i]-1900)])
    }
    return(out)}
  
  saveRDS(ann_max_stoch_arr,paste('./out/synq_annmax_cc_use-bcf=',use_bcf,'_k=',kk,'.rds',sep=''))
  
  print(paste('swm cc end',Sys.time()))}

rm(list=ls());gc()

print(paste('swm end',Sys.time()))



####################################################END####################################################
