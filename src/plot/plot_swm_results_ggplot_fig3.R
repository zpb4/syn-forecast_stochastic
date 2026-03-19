#Figure 3 plotting script

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


#plot 1a: flow duration curve analysis
ffc_dt<-data.table(x=qnorm(1:length(ix_sim)/(length(ix_sim)+1))[1:length(x_swm)],obs=sort(q_fit,decreasing = T)[1:length(x_swm)],sim=sort(s_fit,decreasing = T)[1:length(x_swm)],
                   x_swm=x_swm,ymax=bounds[1,],ymin=bounds[2,],ymax_cc=bounds_cc[1,],ymin_cc=bounds_cc[2,],samp=swm_samps_cc[sample(samps,1),],swm_med=apply(swm_samps,2,median),swm_med_cc=apply(swm_samps_cc,2,median))

#log-normal plot of FDC
ffc_ggp<-ggplot(data=ffc_dt)+theme_minimal()+
  scale_y_continuous(trans='log10',
                     breaks=log10_major_break(),
                     labels=trans_format('log10', math_format(10^.x)),
                     minor_breaks= log10_minor_break(),
                     limits=c(0.1,500))+
  scale_x_continuous(breaks=qnorm(c(0.001,0.02,0.15,0.5,0.85,0.98,0.999)),
                     labels = c(0.001,0.02,0.15,0.5,0.85,0.98,0.999),
                     limits = c(-3.1,3.1))+
  #geom_ribbon(mapping=aes(x=x_swm,ymax=ymax_cc,ymin=ymin_cc,fill='swm4C_99'),alpha=0.15)+
  geom_ribbon(mapping=aes(x=x_swm,ymax=ymax,ymin=ymin,fill='swm_99'),alpha=0.15)+
  geom_line(mapping=aes(x=x,y=obs,color='obs'),linewidth=.75)+
  geom_line(mapping=aes(x=x,y=sim,color='sim'),linewidth=.75)+
  #geom_line(mapping=aes(x=x,y=samp,color='swm4C'),linewidth=.5,alpha=0.5)+
  annotate('text',x=-3.1,y=0.74*500,label='a)',size=6,hjust=-.5)+
  labs(x='exceedance probability',y='flow (kcfs)')+
  coord_cartesian(xlim=c(-3.1,1.25),ylim=c(1,500),expand=F)+
  #scale_color_manual(breaks=c('obs','sim','swm4C'),values=c(clrs[[1]],clrs[[3]],clrs[[2]]))+
  #scale_fill_manual(breaks=c('swm_99','swm4C_99'),values=c(clrs[[9]],clrs[[2]]))+
  scale_color_manual(breaks=c('obs','sim'),labels=c(TeX('$obs_{cal}$'),TeX('$sim_{cal}$')),values=c(clrs[[1]],clrs[[6]]))+
  scale_fill_manual(breaks=c('swm_99'),labels=TeX('$s.hyd^{99}_{hist}$'),values=c(clrs[[4]]))+
  theme(legend.position = 'inside',
        legend.position.inside = c(.8,.8),
        legend.key.spacing.y = unit(0.01,'cm'),
        legend.key.spacing.x = unit(0.01,'cm'),
        legend.title=element_blank(),
        axis.text=element_text(size=9),
        axis.title=element_text(size=12),
        legend.text = element_text(size=9),
        panel.grid.major.y = element_line(size=1))+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank())+
  guides(shape = guide_legend(order=2),col=guide_legend(order=1))

ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/swm_plots/ffc-annmax-samps.png',sep=''),ffc_ggp,dpi=320,width=3.5,height=3,unit='in')
ffc_ggp


#plot 1b: annual maximum log plot
swm_annmax_arr<-array(NA,c(samps,length(stoch_yrs)))
for(i in 1:samps){
  swm_annmax_arr[i,]<-sort(ann_max_stoch_arr[i,])
}

bounds_swm <- apply(swm_annmax_arr,2,function(x){upr<-sort(x)[99];lwr<-sort(x)[2];return(c(upr,lwr))})

swm_annmax_arr_cc<-array(NA,c(samps,length(stoch_yrs)))
for(i in 1:samps){
  swm_annmax_arr_cc[i,]<-sort(ann_max_stoch_arr_cc[i,])
}

bounds_swm_cc <- apply(swm_annmax_arr_cc,2,function(x){upr<-sort(x)[99];lwr<-sort(x)[2];return(c(upr,lwr))})

x_swm = 1/(length(stoch_yrs):1/(length(stoch_yrs)))
y_stoch = sort(ann_max_stoch)
y_stoch_cc = sort(ann_max_stoch_cc)
#swm_samp = swm_annmax_arr[sample(1:samps,1),]
#swm_cc_samp = swm_annmax_arr_cc[sample(1:samps,1),]
swm_annmax_med = apply(swm_annmax_arr,2,median)
swm_cc_annmax_med = apply(swm_annmax_arr_cc,2,median)
x_gev = rep(NA,length(stoch_yrs));x_gev[1:length(seq(0.01,0.999,.001))]<-1/(1-seq(0.01,0.999,.001))
p_obs = (length(obs_yrs):1) / (length(obs_yrs)+1)
x_obs = rep(NA,length(stoch_yrs));x_obs[1:length(obs_yrs)]<-1/(length(obs_yrs):1/(length(obs_yrs)))
gev_obs = rep(NA,length(stoch_yrs));gev_obs[1:length(seq(0.01,0.999,.001))]<-obs_gev_est
emp_obs = rep(NA,length(stoch_yrs));emp_obs[1:length(obs_yrs)]<-sort(ann_max_obs)
x_act = rep(NA,length(stoch_yrs));x_act[1:length(c(2,5,10,20,50,100,200,500,1000))]<-c(2,5,10,20,50,100,200,500,1000)
y_act = rep(NA,length(stoch_yrs));y_act[1:length(c(2,5,10,20,50,100,200,500,1000))]<-c(45,90,140,180,245,295,350,430,500)

annmax_dt<-data.table(x_swm=x_swm,x_gev=x_gev,x_obs=x_obs,gev_obs=gev_obs,emp_obs=emp_obs,x_act=x_act,y_act=y_act,swm_samp=swm_annmax_med,swm_cc_samp=swm_cc_annmax_med,y_stoch=y_stoch,ystoch_cc=y_stoch_cc,
                      ymax=bounds_swm[1,],ymin=bounds_swm[2,],ymax_cc=bounds_swm_cc[1,],ymin_cc=bounds_swm_cc[2,])

df<-data.frame(x1=1,x2=yr_evt,y1=evt,y2=0)
df_scale<-data.frame(x1=1,x2=yr_evt_scale,y1=evt_scale,y2=0)

annmax_ggp<-ggplot(data=annmax_dt)+theme_minimal()+
  scale_x_continuous(trans='log10',
                     breaks=log10_major_break(),
                     labels=trans_format('log10', math_format(10^.x)),
                     minor_breaks= log10_minor_break(),
                     limits=c(1,1010))+ylim(0,1000)+
  geom_point(mapping=aes(x=x_swm,y=y_stoch,color='sim'),size=.5,shape=16)+
  #geom_point(mapping=aes(x=x_swm,y=y_stoch_cc,color='sim_cc'),size=.5,shape=16)+
  #geom_ribbon(mapping=aes(x=x_swm,ymax=ymax_cc,ymin=ymin_cc,fill='swm4C_99'),alpha=0.15)+
  geom_ribbon(mapping=aes(x=x_swm,ymax=ymax,ymin=ymin,fill='swm_99'),alpha=0.15)+
  geom_line(mapping=aes(x=x_swm,y=swm_samp,color='swm-samp'),linewidth=.75)+
  #geom_line(mapping=aes(x=x_swm,y=swm_cc_samp,color='swm4C-samp'),linewidth=.75)+
  geom_line(mapping=aes(x=x_act,y=y_act,color='pva'),linewidth=.5,linetype=2)+
  geom_point(mapping=aes(x=x_obs,y=emp_obs,color='cal_obs'),size=1.5,shape=17)+
  annotate('text',x=1,y=950,label='b)',size=6,hjust=-.5)+
  labs(x='return period (years)',y='flow (kcfs)')+
  geom_segment(aes(x=x1,y=y1,xend=x2,yend=y1),data=df,color=clrs[[1]],linetype=3,linewidth=0.25)+
  geom_segment(aes(x=x2,y=y1,xend=x2,yend=y2),data=df,color=clrs[[1]],linetype=3,linewidth=0.25)+
  annotate('text',x=x_obs[length(ann_max_obs)],y=evt,label='1997',hjust=0.5,vjust=-0.5,size=4,color=clrs[[1]])+
  annotate('text',x=yr_evt,y=0,label=paste(round(yr_evt),'yr'),size=4,hjust=0.5,vjust=-0.25)+
  geom_segment(aes(x=x1,y=y1,xend=x2,yend=y1),data=df_scale,color=clrs[[1]],linetype=3,linewidth=0.25)+
  geom_segment(aes(x=x2,y=y1,xend=x2,yend=y2),data=df_scale,color=clrs[[1]],linetype=3,linewidth=0.25)+
  annotate('text',x=x_obs[length(ann_max_obs)],y=evt_scale,label='1997 (130%)',hjust=0.5,vjust=-0.5,size=4,color=clrs[[1]])+
  annotate('text',x=yr_evt_scale,y=0,label=paste(round(yr_evt_scale),'yr'),size=4,hjust=0.5,vjust=-0.25)+
  #scale_color_manual(breaks=c('cal_obs','pva','swm-samp','swm4C-samp'),values=c(clrs[[7]],clrs[[1]],clrs[[9]],clrs[[2]]))+
  #scale_fill_manual(breaks=c('swm_99','swm4C_99'),values=c(clrs[[9]],clrs[[2]]))+
  #scale_color_manual(breaks=c('cal_obs','sim','sim_cc','pva','swm-samp','swm4C-samp'),values=c(clrs[[7]],clrs[[3]],clrs[[8]],clrs[[1]],clrs[[9]],clrs[[2]]))+
  scale_fill_manual(breaks=c('swm_99'),labels=TeX('$s.hyd^{99}_{hist}$'),values=c(clrs[[4]]))+
  scale_color_manual(breaks=c('cal_obs','sim','pva','swm-samp'),labels=c(TeX('$obs_{cal}$'),TeX('$sim^{SWG}_{hist}$'),'USACE',TeX('$s.hyd^{median}_{hist}$')),values=c(clrs[[1]],clrs[[6]],clrs[[8]],clrs[[4]]))+
  coord_cartesian(xlim=c(1,1050),ylim=c(0,1000),expand=F)+
  theme(legend.position = 'inside',
        legend.position.inside = c(.4,.8),
        legend.box = 'horizontal',
        legend.key.spacing.y = unit(0.01,'cm'),
        legend.key.spacing.x = unit(0.01,'cm'),
        legend.title=element_blank(),
        axis.text=element_text(size=9),
        axis.title=element_text(size=12),
        legend.text = element_text(size=9),
        panel.grid.major.y = element_line(size=1))+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  guides(shape = guide_legend(order=2),col=guide_legend(order=1))

annmax_ggp
ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/swm_plots/log-annmax.png',sep=''),annmax_ggp,dpi=320,width=4,height=3,unit='in')

nocc_gplot<-marrangeGrob(list(ffc_ggp,annmax_ggp),nrow=1,ncol=2,top='')
ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/swm_plots/fdc-annmax_noCC.png',sep=''),nocc_gplot,dpi=320,width=8,height=3.5,unit='in')

#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#Plot 2 - Compare baseline to 4C scenario 
#plot 2a: flow duration curve analysis

#log-normal plot of FDC
ffc_ggp_cc<-ggplot(data=ffc_dt)+theme_minimal()+
  scale_y_continuous(trans='log10',
                     breaks=log10_major_break(),
                     labels=trans_format('log10', math_format(10^.x)),
                     minor_breaks= log10_minor_break(),
                     limits=c(0.1,500))+
  scale_x_continuous(breaks=qnorm(c(0.001,0.02,0.15,0.5,0.85,0.98,0.999)),
                     labels = c(0.001,0.02,0.15,0.5,0.85,0.98,0.999),
                     limits = c(-3.1,3.1))+
  geom_ribbon(mapping=aes(x=x_swm,ymax=ymax_cc,ymin=ymin_cc,fill='swm4C_99'),alpha=0.25)+
  geom_ribbon(mapping=aes(x=x_swm,ymax=ymax,ymin=ymin,fill='swm_99'),alpha=0.15)+
  #geom_line(mapping=aes(x=x,y=obs,color='obs'),linewidth=.75)+
  #geom_line(mapping=aes(x=x,y=sim,color='sim'),linewidth=.5,alpha=0.5)+
  geom_line(mapping=aes(x=x,y=swm_med,color='swm'),linewidth=.75,alpha=0.5)+
  geom_line(mapping=aes(x=x,y=swm_med_cc,color='swm_4C'),linewidth=.75,alpha=0.5)+
  annotate('text',x=-3.1,y=0.75*500,label='c)',size=6,hjust=-.5)+
  labs(x='exceedance probability',y='flow (kcfs)')+
  coord_cartesian(xlim=c(-3.1,1.25),ylim=c(1,500),expand=F)+
  #scale_color_manual(breaks=c('obs','sim','swm4C'),values=c(clrs[[1]],clrs[[3]],clrs[[2]]))+
  #scale_fill_manual(breaks=c('swm_99','swm4C_99'),values=c(clrs[[9]],clrs[[2]]))+
  scale_color_manual(breaks=c('swm','swm_4C'),labels=c(TeX('$s.hyd_{hist}^{median}$'),TeX('$s.hyd_{4C}^{median}$')),values=c(clrs[[4]],clrs[[2]]))+
  scale_fill_manual(breaks=c('swm_99','swm4C_99'),labels=c(TeX('$s.hyd^{99}_{hist}$'),TeX('$s.hyd_{4C}^{99}$')),values=c(clrs[[4]],clrs[[2]]))+
  theme(legend.position = 'inside',
        legend.position.inside = c(.8,.7),
        legend.key.spacing.y = unit(0.01,'cm'),
        legend.key.spacing.x = unit(0.01,'cm'),
        legend.title=element_blank(),
        axis.text=element_text(size=9),
        axis.title=element_text(size=12),
        legend.text = element_text(size=9),
        panel.grid.major.y = element_line(size=1))+
  guides(shape = guide_legend(order=2),col=guide_legend(order=1))

ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/swm_plots/ffc-annmax-samps_CC.png',sep=''),ffc_ggp_cc,dpi=320,width=3.5,height=3,unit='in')
ffc_ggp_cc

#plot 2b: annual maximum log plot
stoch_yr_evts <- tail(1/(length(stoch_yrs):1/(length(stoch_yrs))),100)
stoch_yr_evt_mags <- tail(swm_annmax_med,100)
rtns <- cbind(stoch_yr_evts,stoch_yr_evt_mags)

y1<-lin_intpx(yr_evt_cutoff,rtns)
df<-data.frame(x1=1,x2=yr_evt_cutoff,y1=y1,y2=0)

stoch_yr_evts <- tail(1/(length(stoch_yrs):1/(length(stoch_yrs))),100)
stoch_yr_evt_magscc <- tail(swm_cc_annmax_med,100)
rtns <- cbind(stoch_yr_evts,stoch_yr_evt_magscc)

y1<-lin_intpx(yr_evt_cutoff,rtns)
df_cc<-data.frame(x1=1,x2=yr_evt_cutoff,y1=y1,y2=0)

swm_cutoff <- df$y1
swmcc_cutoff <- df_cc$y1

swm_maxs <- apply(synq_sim_arr,1,max)
swmcc_maxs <- apply(synq_sim_arr_cc,1,max)

swm_samp_vec <- (swm_cutoff-swm_maxs);swm_samp_vec[swm_samp_vec<0]<-NA
swm_maxs_closesamps <- order(swm_samp_vec)

swmcc_samp_vec <- (swmcc_cutoff-swmcc_maxs);swmcc_samp_vec[swmcc_samp_vec<0]<-NA
swmcc_maxs_closesamps <- order(swmcc_samp_vec)

saveRDS(swm_maxs_closesamps,paste('z:/Synthetic-Forecast-v2-FIRO-DISES/data/YRS/',yr_evt_cutoff,'yr-samples_swm.rds',sep=''))
saveRDS(swmcc_maxs_closesamps,paste('z:/Synthetic-Forecast-v2-FIRO-DISES/data/YRS/',yr_evt_cutoff,'yr-samples_swmcc.rds',sep=''))

write.csv(data.frame(dist=1:length(swm_maxs_closesamps),idx=swm_maxs_closesamps),paste('z:/Synthetic-Forecast-v2-FIRO-DISES/data/YRS/',yr_evt_cutoff,'yr-samples_swm.csv',sep=''),row.names = F)
write.csv(data.frame(dist=1:length(swmcc_maxs_closesamps),idx=swmcc_maxs_closesamps),paste('z:/Synthetic-Forecast-v2-FIRO-DISES/data/YRS/',yr_evt_cutoff,'yr-samples_swmcc.csv',sep=''),row.names = F)

annmax_ggp_cc<-ggplot(data=annmax_dt)+theme_minimal()+
  scale_x_continuous(trans='log10',
                     breaks=log10_major_break(),
                     labels=trans_format('log10', math_format(10^.x)),
                     minor_breaks= log10_minor_break(),
                     limits=c(1,1010))+ylim(0,1000)+
  #geom_point(mapping=aes(x=x_swm,y=y_stoch,color='sim'),size=.5,shape=16)+
  #geom_point(mapping=aes(x=x_swm,y=y_stoch_cc,color='sim_cc'),size=.5,shape=16)+
  geom_ribbon(mapping=aes(x=x_swm,ymax=ymax_cc,ymin=ymin_cc,fill='swm4C_99'),alpha=0.25)+
  geom_ribbon(mapping=aes(x=x_swm,ymax=ymax,ymin=ymin,fill='swm_99'),alpha=0.15)+
  geom_line(mapping=aes(x=x_swm,y=swm_samp,color='swm-samp'),linewidth=.75)+
  geom_line(mapping=aes(x=x_swm,y=swm_cc_samp,color='swm4C-samp'),linewidth=.75)+
  geom_line(mapping=aes(x=x_act,y=y_act,color='pva'),linewidth=.5,linetype=2)+
  geom_segment(aes(x=x1,y=y1,xend=x2,yend=y1),data=df,color=clrs[[1]],linetype=3,linewidth=0.25)+
  geom_segment(aes(x=x2,y=y1,xend=x2,yend=y2),data=df,color=clrs[[1]],linetype=3,linewidth=0.25)+
  annotate('text',x=1.5,y=df$y1,label=paste(round(df$y1,1)),hjust=0,vjust=1.25,size=4,color=clrs[[1]])+
  geom_segment(aes(x=x1,y=y1,xend=x2,yend=y1),data=df_cc,color=clrs[[2]],linetype=3,linewidth=0.25)+
  geom_segment(aes(x=x2,y=y1,xend=x2,yend=y2),data=df_cc,color=clrs[[2]],linetype=3,linewidth=0.25)+
  annotate('text',x=1.5,y=df_cc$y1,label=paste(round(df_cc$y1,1)),hjust=0,vjust=1.25,size=4,color=clrs[[2]])+
  annotate('text',x=yr_evt_cutoff,y=0,label=paste(yr_evt_cutoff,'yr'),size=4,hjust=0.5,vjust=-0.25)+
  #geom_point(mapping=aes(x=x_obs,y=emp_obs,color='cal_obs'),size=1.5,shape=17)+
  annotate('text',x=1,y=950,label='d)',size=6,hjust=-.5)+
  labs(x='return period (years)',y='flow (kcfs)')+
  #geom_segment(aes(x=x1,y=y1,xend=x2,yend=y1),data=df,color=clrs[[1]],linetype=3,linewidth=0.25)+
  #geom_segment(aes(x=x2,y=y1,xend=x2,yend=y2),data=df,color=clrs[[1]],linetype=3,linewidth=0.25)+
  #annotate('text',x=x_obs[length(ann_max_obs)],y=evt,label='1997',hjust=1,vjust=-0.25,size=4,color=clrs[[7]])+
  #annotate('text',x=yr_evt,y=0,label=paste(round(yr_evt),'year'),size=4,hjust=0,vjust=-0.25)+
  #scale_color_manual(breaks=c('cal_obs','pva','swm-samp','swm4C-samp'),values=c(clrs[[7]],clrs[[1]],clrs[[9]],clrs[[2]]))+
  #scale_fill_manual(breaks=c('swm_99','swm4C_99'),values=c(clrs[[9]],clrs[[2]]))+
  #scale_color_manual(breaks=c('cal_obs','sim','sim_cc','pva','swm-samp','swm4C-samp'),values=c(clrs[[7]],clrs[[3]],clrs[[8]],clrs[[1]],clrs[[9]],clrs[[2]]))+
  scale_fill_manual(breaks=c('swm_99','swm4C_99'),labels=c(TeX('$s.hyd_{hist}^{99}$'),TeX('$s.hyd_{4C}^{99}$')),values=c(clrs[[4]],clrs[[2]]))+
  scale_color_manual(breaks=c('pva','swm-samp','swm4C-samp'),labels=c('USACE',TeX('$s.hyd_{hist}^{median}$'),TeX('$s.hyd_{4C}^{median}$')),values=c(clrs[[8]],clrs[[4]],clrs[[2]]))+
  coord_cartesian(xlim=c(1,1050),ylim=c(0,1000),expand=F)+
  theme(legend.position = 'inside',
        legend.position.inside = c(.4,.8),
        legend.box = 'horizontal',
        legend.key.spacing.y = unit(0.01,'cm'),
        legend.key.spacing.x = unit(0.01,'cm'),
        legend.title=element_blank(),
        axis.text=element_text(size=9),
        axis.title=element_text(size=12),
        legend.text = element_text(size=9),
        panel.grid.major.y = element_line(size=1))+
  theme(axis.title.y=element_blank())+
  guides(shape = guide_legend(order=2),col=guide_legend(order=1))

annmax_ggp_cc
ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/swm_plots/log-annmax_CC.png',sep=''),annmax_ggp_cc,dpi=320,width=4,height=3,unit='in')

cc_gplot<-marrangeGrob(list(ffc_ggp_cc,annmax_ggp_cc),nrow=1,ncol=2,top='')
ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/swm_plots/fdc-annmax_CC.png',sep=''),cc_gplot,dpi=320,width=8,height=3.5,unit='in')

#final manuscript plot
comb_gplot<-marrangeGrob(list(ffc_ggp,ffc_ggp_cc,annmax_ggp,annmax_ggp_cc),nrow=2,ncol=2,top='')
ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/manuscript/fdc-annmax_Fig3.png',sep=''),comb_gplot,dpi=320,width=8,height=6,unit='in')



################################################END########################################################