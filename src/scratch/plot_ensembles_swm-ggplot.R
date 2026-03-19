

rm(list=ls());gc()
setwd('z:/syn-forecast_stochastic//')
library(fields)
library(scales)
library(zoo)
library(scales)
library(tidyverse)
library(data.table)
library(gridExtra)
library(ks)
library(viridis)
library(RColorBrewer)
#///////////////////////////////////////////////////////////////////////////////////////////////////////////
#Data setup
syn_vers = 2
loc = 'YRS'
keysite_name = 'ORDC1'
cal_val_setup = 'cal'
pcnt_opt = 0.99
obj_pwr = 0
opt_strat = 'ecrps-dts'
same_swm = F
samp_no = 2

skill_mod = 0
skill_dcy = 0.1
skill_tail = 0.5

#swm data
use_bcf=F
k_swm=106
yr_cutoff=500
set.seed(1)

#pick event
evt_order=1
nn_order=99  #plot nn_order nearest neighbor in swm and swm_cc; 
#set to 0 to plot evt_order for swm (swm_cc inherits swm), 99 to set both swm and swm_cc to evt_order

ix_swm<-seq(as.Date('2001-01-01'),as.Date('3008-01-08'),'day')
ix2_swm<-as.POSIXlt(ix_swm)

#plot setup
disp_pcnt <- 0.99
n_samp <- 10
save_plots <- F  # T to save plots to png files, F to display in R
clrs<-palette.colors()

path = paste('z:/Synthetic-Forecast-v',syn_vers,'-FIRO-DISES/',sep='')

load(paste(path,'out/',loc,'/data_prep_rdata.RData',sep=''))
rm(hefs_forward_cumul,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid,hefs_forward_frac,hefs_inp);gc()

hefs_forward <- readRDS(paste(path,'out/',loc,'/hefs_forward_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'.rds',sep=''))
syn_hefs_forward <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',keysite_name,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_swm.rds',sep=''))
syn_hefs_forward_cc <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',keysite_name,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_swm-cc.rds',sep=''))

obs_forward_all_leads_swm <- readRDS(paste('out/',loc,'/obs-fwd_use-bcf=',use_bcf,'_k=',k_swm,'_same-swm=',same_swm,'_cutoff=',yr_cutoff,'_swm.rds',sep=''))
obs_forward_all_leads_swmcc <- readRDS(paste('out/',loc,'/obs-fwd_use-bcf=',use_bcf,'_k=',k_swm,'_same-swm=',same_swm,'_cutoff=',yr_cutoff,'_swm-cc.rds',sep=''))

ixx_gen_swm <- readRDS(paste('out/',loc,'/ixx_gen_swm.rds',sep=''))
ixx_gen <- ixx_obs_forward

shefs_fwd_swm <- syn_hefs_forward[samp_no,1,,,]
shefs_fwd_swmcc <- syn_hefs_forward_cc[samp_no,1,,,]

obs_swm <- obs_forward_all_leads_swm[samp_no,1,,1]
obs_swm <- obs_forward_all_leads_swmcc[samp_no,1,,1]

obs_fwd_swm <- obs_forward_all_leads_swm[samp_no,1,,]
obs_fwd_swmcc <- obs_forward_all_leads_swmcc[samp_no,1,,]

rm(syn_hefs_forward,syn_hefs_forward_cc,obs_forward_all_leads_swm,obs_forward_all_leads_swmcc)
gc()

source('./src/forecast_verification_functions.R')

#///////////////////////////////////////////////////////////////////////////////////////////////////////////
#Ensemble plots
keysite = which(site_names==keysite_name)

hefs_fwd <- hefs_forward[1,keysite,,,]
obs_forward_all_leads_hind <- obs_forward_all_leads_hind[keysite,,]
obs_hefs <- obs[ixx_obs%in%ixx_hefs,keysite]


#################HEFS ensemble plots################################
hind_evt_idx <- order(obs_hefs,decreasing=T)[evt_order]
hind_evt_date <- ixx_hefs[hind_evt_idx]

leads = c(1,3,5,7)
hefs_plts = vector('list',length(leads))
ylm_scale = 2
show_x = T

for(i in 1:length(leads)){
  hefs_plt_dt <- data.table(x=0:15,hefs=rbind(rep(obs_hefs[hind_evt_idx-leads[i]],dim(hefs_fwd)[1]),t(hefs_fwd[,hind_evt_idx-leads[i],])),
                          obs=c(obs_hefs[hind_evt_idx-leads[i]],obs_forward_all_leads_hind[hind_evt_idx-leads[i],]))
  
  ylm = ylm_scale * max(hefs_plt_dt$obs)
  #convert to long form
  df <- melt(hefs_plt_dt,id='x')
  
  #assign long form df to 'obs' and 'hefs' groups
  df$col <- factor(ifelse(df$variable=="obs",1,2),levels=1:2,labels=c("obs","hefs"))

  hefs_ens_plt<-ggplot(df)+theme_minimal()+
    #geom_smooth(mapping=aes(x=x,y=value,group=variable,color=grp,size=grp,alpha=grp),se=F,span=.15)+
    geom_line(mapping=aes(x=x,y=value,group=variable,color=col,linewidth=col,alpha=col))+
    scale_color_manual(values=c('obs'=clrs[[1]],'hefs'=clrs[[9]]))+
    scale_linewidth_manual(values=c('obs'=1.5,'hefs'=.75))+
    scale_alpha_manual(values=c('obs'=1,'hefs'=.35))+
    geom_vline(xintercept = leads[i],linetype='dotted',linewidth=0.25)+
    annotate('text',x=leads[i],y=0.95*ylm,label=hind_evt_date,size=4,hjust=-0.15)+
    #annotate('text',x=1,y=850,label='b)',size=6)+
    labs(x='forecast lead (days)',y='flow (kcfs)')+
    coord_cartesian(xlim=c(0,15),ylim=c(0,ylm),expand=F)+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(size=1))+
    {if(i==1)
      theme(legend.position = 'inside',
            legend.position.inside = c(.8,.7),
            legend.box = 'horizontal',
            legend.key.spacing.y = unit(0.01,'cm'),
            legend.key.spacing.x = unit(0.01,'cm'),
            legend.title=element_blank(),
            legend.text = element_text(size=12))}+
    {if(i!=1)
      theme(legend.position = 'none',
            axis.text.y=element_blank(),
            axis.title.y=element_blank())}+
    {if(show_x==F)
      theme(legend.position = 'none',
            axis.text.x=element_blank(),
            axis.title.x=element_blank())}

  hefs_plts[[i]] <- hefs_ens_plt}

mat<-matrix(1:length(leads),ncol=length(leads))
hefs_gplot<-marrangeGrob(hefs_plts,nrow=1,ncol=length(leads),top='')
#hefs_gplot

ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ind_plots/hefs-ens-plot_evt=',evt_order,'_lds=',str_flatten(leads,collapse='-'),'_skillmod=',skill_mod,'.png',sep=''),hefs_gplot,dpi=320,width=9,height=3,unit='in')


#################syn-HEFS + SWM_nocc ensemble plots################################
#pick event
if(nn_order!=0){
  diff_vec <- abs(obs_hefs[hind_evt_idx] - obs_swm)
  swm_evt_idx <- order(diff_vec,decreasing=F)[nn_order]}
if(nn_order==0|nn_order==99){
  swm_evt_idx <- order(obs_swm,decreasing=T)[evt_order]}
swm_evt_date <- ixx_gen[swm_evt_idx]

shefs_plts = vector('list',length(leads))
show_x = F

for(i in 1:length(leads)){
  shefs_plt_dt <- data.table(x=0:15,hefs=rbind(rep(obs_swm[swm_evt_idx-leads[i]],dim(hefs_fwd)[1]),t(shefs_fwd_swm[,swm_evt_idx-leads[i],])),
                            obs=c(obs_fwd_swm[swm_evt_idx-leads[i],]))
  
  ylm = ylm_scale * max(hefs_plt_dt$obs)
  #convert to long form
  df <- melt(shefs_plt_dt,id='x')
  
  #assign long form df to 'obs' and 'hefs' groups
  df$col <- factor(ifelse(df$variable=="obs",1,2),levels=1:2,labels=c("obs","shefs"))
  
  hefs_ens_plt<-ggplot(df)+theme_minimal()+
    #geom_smooth(mapping=aes(x=x,y=value,group=variable,color=grp,size=grp,alpha=grp),se=F,span=.15)+
    geom_line(mapping=aes(x=x,y=value,group=variable,color=col,linewidth=col,alpha=col))+
    scale_color_manual(values=c('obs'=clrs[[1]],'shefs'=clrs[[4]]))+
    scale_linewidth_manual(values=c('obs'=1.5,'shefs'=.75))+
    scale_alpha_manual(values=c('obs'=1,'shefs'=.25))+
    geom_vline(xintercept = leads[i],linetype='dotted',linewidth=0.25)+
    annotate('text',x=leads[i],y=0.95*ylm,label=swm_evt_date,size=4,hjust=-0.15)+
    #annotate('text',x=1,y=850,label='b)',size=6)+
    labs(x='forecast lead (days)',y='flow (kcfs)')+
    coord_cartesian(xlim=c(0,15),ylim=c(0,ylm),expand=F)+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(size=1))+
    {if(i==1)
      theme(legend.position = 'inside',
            legend.position.inside = c(.8,.7),
            legend.box = 'horizontal',
            legend.key.spacing.y = unit(0.01,'cm'),
            legend.key.spacing.x = unit(0.01,'cm'),
            legend.title=element_blank(),
            legend.text = element_text(size=12))}+
    {if(i!=1)
      theme(legend.position = 'none',
            axis.text.y=element_blank(),
            axis.title.y=element_blank())}
  {if(show_x==F)
    theme(legend.position = 'none',
          axis.text.x=element_blank(),
          axis.title.x=element_blank())}
  
  shefs_plts[[i]] <- hefs_ens_plt}

shefs_gplot<-marrangeGrob(shefs_plts,nrow=1,ncol=length(leads),top='')
#shefs_gplot

ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ind_plots/shefs_swm_ens-plot_evt=',evt_order,'_nn=',nn_order,'_lds=',str_flatten(leads,collapse='-'),'_samp=',samp_no,'_skillmod=',skill_mod,'.png',sep=''),shefs_gplot,dpi=320,width=9,height=3,unit='in')


#################syn-HEFS + SWM_cc ensemble plots################################
#pick event
if(nn_order==99){
  swm_evt_idx <- order(obs_swm_cc,decreasing=T)[evt_order]}
swm_evt_date <- ixx_gen[swm_evt_idx]

shefscc_plts = vector('list',length(leads))
show_x = F

for(i in 1:length(leads)){
  shefs_plt_dt <- data.table(x=0:15,hefs=rbind(rep(obs_swm_cc[swm_evt_idx-leads[i]],dim(hefs_fwd)[1]),t(shefs_fwd_swmcc[,swm_evt_idx-leads[i],])),
                             obs=c(obs_fwd_swmcc[swm_evt_idx-leads[i],]))
  
  ylm = ylm_scale * max(hefs_plt_dt$obs)
  #convert to long form
  df <- melt(shefs_plt_dt,id='x')
  
  #assign long form df to 'obs' and 'hefs' groups
  df$col <- factor(ifelse(df$variable=="obs",1,2),levels=1:2,labels=c("obs","shefs"))
  
  hefs_ens_plt<-ggplot(df)+theme_minimal()+
    #geom_smooth(mapping=aes(x=x,y=value,group=variable,color=grp,size=grp,alpha=grp),se=F,span=.15)+
    geom_line(mapping=aes(x=x,y=value,group=variable,color=col,linewidth=col,alpha=col))+
    scale_color_manual(values=c('obs'=clrs[[1]],'shefs'=clrs[[2]]))+
    scale_linewidth_manual(values=c('obs'=1.5,'shefs'=.75))+
    scale_alpha_manual(values=c('obs'=1,'shefs'=.25))+
    geom_vline(xintercept = leads[i],linetype='dotted',linewidth=0.25)+
    annotate('text',x=leads[i],y=0.95*ylm,label=swm_evt_date,size=4,hjust=-0.15)+
    #annotate('text',x=1,y=850,label='b)',size=6)+
    labs(x='forecast lead (days)',y='flow (kcfs)')+
    coord_cartesian(xlim=c(0,15),ylim=c(0,ylm),expand=F)+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(size=1))+
    {if(i==1)
      theme(legend.position = 'inside',
            legend.position.inside = c(.8,.7),
            legend.box = 'horizontal',
            legend.key.spacing.y = unit(0.01,'cm'),
            legend.key.spacing.x = unit(0.01,'cm'),
            legend.title=element_blank(),
            legend.text = element_text(size=12))}+
    {if(i!=1)
      theme(legend.position = 'none',
            axis.text.y=element_blank(),
            axis.title.y=element_blank())}
  {if(show_x==F)
    theme(legend.position = 'none',
          axis.text.x=element_blank(),
          axis.title.x=element_blank())}
  
  shefscc_plts[[i]] <- hefs_ens_plt}

shefscc_gplot<-marrangeGrob(shefscc_plts,nrow=1,ncol=length(leads),top='')
#shefscc_gplot

ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ind_plots/shefs_swm-cc_ens-plot_evt=',evt_order,'_nn=',nn_order,'_lds=',str_flatten(leads,collapse='-'),'_samp=',samp_no,'_skillmod=',skill_mod,'.png',sep=''),shefscc_gplot,dpi=320,width=9,height=3,unit='in')

layout_mat <- matrix(1:(3*length(leads)),nrow=3,ncol=length(leads),byrow=T)
comb_plot<-marrangeGrob(c(hefs_plts,shefs_plts,shefscc_plts),nrow=3,ncol=length(leads),top='',layout_matrix = layout_mat)
ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/forecast_plots/combined_ens-plot_evt=',evt_order,'_nn=',nn_order,'_lds=',str_flatten(leads,collapse='-'),'_samp=',samp_no,'_skillmod=',skill_mod,'.png',sep=''),comb_plot,dpi=320,width=9,height=9,unit='in')


########################################################################

