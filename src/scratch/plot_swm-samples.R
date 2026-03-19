
#SWM samples plotting script

rm(list=ls());gc()
setwd('z:/syn-forecast_stochastic/')

library(extRemes)
library(scales)
library(tidyverse)
library(data.table)
library(gridExtra)
library(ks)
library(viridis)
library(latex2exp)
library(RColorBrewer)
clrs<-palette.colors()
show_col(clrs)

use_bcf=F
kk=106
yr_evt_cutoff=500

stoch_data<-read.csv('./data/swm/sacsma_ORO_wgen_baseline.csv')
stoch_data_cc<-read.csv('./data/swm/sacsma_ORO_wgen_4c_7cc_0mn.csv')
#obs_data<-read.csv('./data/swm/observed_flows.csv')
obs_data<-read.table('./data/swm/FNF_ORO_mm.txt')

obs_fit<-read.table('./data/swm/FNF_ORO_mm.txt')
sim_fit<-read.table('./data/swm/ORO_SACSMA.txt')

source('./src/mm-cfs_conversion.R')
source('./src/helper-functions.R')
mm_to_kcfs <- area_mi2_ORO*mm_to_cfs/1000


ix_stoch<-seq(as.Date('2001-01-01'),as.Date('3008-01-08'),'day')
ix2_stoch<-as.POSIXlt(ix_stoch)

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

ann_max_stoch<-c()
stoch_yrs <- unique(ix2_stoch_wy$year)+1900
for(i in 1:length(stoch_yrs)){
  ann_max_stoch[i]<-max(base_q[ix2_stoch_wy$year%in%(stoch_yrs[i]-1900)])
}

ann_max_stoch_cc<-c()
for(i in 1:length(stoch_yrs)){
  ann_max_stoch_cc[i]<-max(cc_q[ix2_stoch_wy$year%in%(stoch_yrs[i]-1900)])
}

synq_sim_arr<-readRDS(paste('./out/swm/synq_sim_arr_use-bcf=',use_bcf,'_k=',kk,'.rds',sep=''))
synq_sim_arr_cc<-readRDS(paste('./out/swm/synq_sim_arr_cc_use-bcf=',use_bcf,'_k=',kk,'.rds',sep=''))
synq_sim_arr_fit<-readRDS(paste('./out/swm/synq_arr_fit_use-bcf=',use_bcf,'_k=',kk,'.rds',sep=''))

ann_max_stoch_arr<-readRDS(paste('./out/swm/synq_annmax_use-bcf=',use_bcf,'_k=',kk,'.rds',sep=''))
ann_max_stoch_arr_cc<-readRDS(paste('./out/swm/synq_annmax_cc_use-bcf=',use_bcf,'_k=',kk,'.rds',sep=''))
ann_max_stoch_arr_fit<-readRDS(paste('./out/swm/synq_annmax_fit_use-bcf=',use_bcf,'_k=',kk,'.rds',sep=''))

samps <- dim(synq_sim_arr)[1]

stoch_yrs <- unique(ix2_stoch_wy$year)+1900
fit_yrs <- unique(ix2_sim_wy$year)+1900

q_fit <- obs_fit[,'V4'][ix_obs%in%ix_sim] * mm_to_kcfs
s_fit <- sim_fit[,'V4'] * mm_to_kcfs
q_fit[q_fit<=0]<-min(q_fit[q_fit>0]);s_fit[s_fit<=0]<-min(q_fit[q_fit>0])

#SWM plots
ann_max_sfit<-c()
for(i in 1:length(fit_yrs)){
  ann_max_sfit[i]<-max(s_fit[ix2_sim_wy$year%in%(fit_yrs[i]-1900)])
}

obs_gev_fit<-fevd(ann_max_obs,type='GEV')
obs_gev_est<-qevd(seq(0.01,0.999,.001),loc=obs_gev_fit$results$par[1],scale =obs_gev_fit$results$par[2],shape = obs_gev_fit$results$par[3],type='GEV')

#///////////////////////////////////////////////////////////////////////////////////
#Plots

#setup
pcntile<-.25
y_thresh<-sort(obs)[round(length(obs)*pcntile)]
x_swm = qnorm(1:(365*length(fit_yrs))/(365*length(fit_yrs)+1))

evt_idx = order(obs,decreasing=T)[1]
evt = obs[evt_idx]
evt_scale = obs[evt_idx]*1.3
est_rtns<-cbind(c(2,5,10,20,50,100,200,500,1000),c(45,90,140,180,245,295,350,430,500))
lin_intp <- function(y,rtns){
  dif_vec <- y - rtns[,2]
  lwr_id <- which.min(dif_vec[dif_vec>=0])
  slp <- (log(rtns[lwr_id+1,1])-log(rtns[lwr_id,1]))/(rtns[lwr_id+1,2]-rtns[lwr_id,2])
  est <- log(rtns[lwr_id,1])+dif_vec[lwr_id]*slp
  return(exp(est))
}

lin_intpx <- function(x,rtns){
  dif_vec <- log(x) - log(rtns[,1])
  lwr_id <- which.min(dif_vec[dif_vec>=0])
  slp <- (rtns[lwr_id+1,2]-rtns[lwr_id,2]) / (log(rtns[lwr_id+1,1])-log(rtns[lwr_id,1]))
  est <- rtns[lwr_id,2]+dif_vec[lwr_id]*slp
  return(est)
}

est_yr_evt<-lin_intp(obs[evt_idx],est_rtns)
est_yr_evt_scale<-lin_intp(obs[evt_idx]*1.3,est_rtns)

yr_evt <- est_yr_evt
yr_idx <- which.min(abs(1/(length(stoch_yrs):1/(length(stoch_yrs)))-yr_evt))

yr_evt_scale <- est_yr_evt_scale
yr_idx_scale <- which.min(abs(1/(length(stoch_yrs):1/(length(stoch_yrs)))-yr_evt_scale))

yr_idx_cutoff <- lin_intpx(yr_evt_cutoff,est_rtns)

swm_samps <- array(NA,c(samps*10,365*length(fit_yrs)))
for(i in 1:10){
  for(k in 1:samps){
    sidx<-sample(stoch_yrs,length(fit_yrs),replace = T)-1900
    swm_samps[(i-1)*samps+k,]<-sort(synq_sim_arr[k,ix2_stoch_wy$year%in%sidx][1:(365*length(fit_yrs))],decreasing = T,na.last=T)}
}

bounds <- apply(swm_samps,2,function(x){upr<-sort(x)[1000];lwr<-sort(x)[1];return(c(upr,lwr))})
bounds[is.na(bounds)==T]<-0
bounds[1,bounds[1,]==0]<-min(bounds[1,bounds[1,]>0])
bounds[2,bounds[2,]==0]<-min(bounds[2,bounds[2,]>0])

swm_samps_cc <- array(NA,c(samps*10,365*length(fit_yrs)))
for(i in 1:10){
  for(k in 1:samps){
    sidx<-sample(stoch_yrs,length(fit_yrs),replace = T)-1900
    swm_samps_cc[(i-1)*samps+k,]<-sort(synq_sim_arr_cc[k,ix2_stoch_wy$year%in%sidx][1:(365*length(fit_yrs))],decreasing = T,na.last=T)}
}

bounds_cc <- apply(swm_samps_cc,2,function(x){upr<-sort(x)[1000];lwr<-sort(x)[1];return(c(upr,lwr))})
bounds_cc[is.na(bounds_cc)==T]<-0
bounds_cc[1,bounds_cc[1,]==0]<-min(bounds_cc[1,bounds_cc[1,]>0])
bounds_cc[2,bounds_cc[2,]==0]<-min(bounds_cc[2,bounds_cc[2,]>0])

#/////////////////////////////////////////////////////////////////////////////////////////////////
#Plot 3: SWM samples
#plot 1c: SWM samples of the same return period as 1997 event
disp_samps = 5
window = 15

swm_samps <- array(NA,c(length(-window:window),disp_samps))
swm_samps_cc <- array(NA,c(length(-window:window),disp_samps))
for(i in 1:disp_samps){
  evt = sort(ann_max_stoch_arr[i,])[yr_idx]
  event_idx = which(synq_sim_arr[i,]==evt)
  swm_samps[,i] <- synq_sim_arr[i,(event_idx-15):(event_idx+15)]
  swm_samps_cc[,i] <- synq_sim_arr_cc[i,(event_idx-15):(event_idx+15)]
}

swm_samps_df <- data.table(x=-15:15,obs=obs[(evt_idx-15):(evt_idx+15)],samps=swm_samps)
colnames(swm_samps_df)<-c('x','1997',paste('samp',1:disp_samps,sep=''))
df <- melt(swm_samps_df,id='x')

#cls=c('black',viridis(5,end=0.9))
cls = palette.colors()[c(1:4,8,9)]

samps_ggp<-ggplot(df)+theme_minimal()+
  geom_line(mapping=aes(x=x,y=value,color=variable,size=variable,alpha=variable))+
  scale_color_manual(values=cls)+
  scale_size_manual(values=c(1.5,rep(.75,disp_samps)))+
  scale_alpha_manual(values=c(1,rep(0.5,disp_samps)))+
  #annotate('text',x=1,y=850,label='b)',size=6)+
  labs(x='days from event (days)',y='flow (kcfs)')+
  coord_cartesian(xlim=c(-window,window),ylim=c(0,450),expand=F)+
  theme(legend.position = 'inside',
        legend.position.inside = c(.8,.7),
        legend.box = 'horizontal',
        legend.key.spacing.y = unit(0.01,'cm'),
        legend.key.spacing.x = unit(0.01,'cm'),
        legend.title=element_blank(),
        axis.text=element_text(size=9),
        axis.title=element_text(size=12),
        legend.text = element_text(size=9),
        panel.grid.major.y = element_line(size=1))

samps_ggp

ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/swm_plots/samps-97.png',sep=''),samps_ggp,dpi=320,width=4,height=3,unit='in')


#plot 1d: SWM_4C samples of the same return period as 1997 event
show_y=F
#swm_samps <- array(NA,c(length(-window:window),disp_samps))
#for(i in 1:disp_samps){
#evt = sort(ann_max_stoch_arr_cc[i,])[yr_idx]
#event_idx = which(synq_sim_arr_cc[i,]==evt)
#swm_samps[,i] <- synq_sim_arr_cc[i,(event_idx-15):(event_idx+15)]
#}

swm_samps_df <- data.table(x=-15:15,obs=obs[(evt_idx-15):(evt_idx+15)],samps=swm_samps_cc)
colnames(swm_samps_df)<-c('x','1997',paste('samp',1:disp_samps,sep=''))
df <- melt(swm_samps_df,id='x')

#cls=c('black',inferno(5,begin=0.2,end=0.8))
cls = palette.colors()[c(1:4,8,9)]

samps_ggp_cc<-ggplot(df)+theme_minimal()+
  geom_line(mapping=aes(x=x,y=value,color=variable,size=variable,alpha=variable))+
  scale_color_manual(values=cls)+
  scale_size_manual(values=c(1.5,rep(.75,disp_samps)))+
  scale_alpha_manual(values=c(1,rep(0.5,disp_samps)))+
  #annotate('text',x=1,y=850,label='b)',size=6)+
  labs(x='days from event (days)',y='flow (kcfs)')+
  coord_cartesian(xlim=c(-window,window),ylim=c(0,450),expand=F)+
  theme(legend.position = 'inside',
        legend.position.inside = c(.8,.7),
        legend.box = 'horizontal',
        legend.key.spacing.y = unit(0.01,'cm'),
        legend.key.spacing.x = unit(0.01,'cm'),
        legend.title=element_blank(),
        axis.text=element_text(size=9),
        axis.title=element_text(size=12),
        legend.text = element_text(size=9),
        panel.grid.major.y = element_line(size=1))+
  {if(show_y==F)
    theme(axis.text.y=element_blank(),
          axis.title.y=element_blank())}

samps_ggp_cc

ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/swm_plots/samps-97_cc.png',sep=''),samps_ggp_cc,dpi=320,width=4,height=3,unit='in')


cc_gplot<-marrangeGrob(list(samps_ggp,samps_ggp_cc),nrow=1,ncol=2,top='')
ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/manuscript/SI/samps-97_comb.png',sep=''),cc_gplot,dpi=320,width=8,height=3.5,unit='in')

################################################END######################################################
