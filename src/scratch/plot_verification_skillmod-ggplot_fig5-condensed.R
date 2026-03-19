#Figure 5 plotting script

rm(list=ls());gc()
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
library(latex2exp)

setwd('z:/syn-forecast_stochastic//')
#///////////////////////////////////////////////////////////////////////////////////////////////////////////
#Data setup
syn_vers = 2
syn_vers = 2
loc = 'YRS'
site = 'ORDC1'
disp_site = 'ORDC1'
cal_val_setup = 'cal'
pcnt_opt = 0.99
obj_pwr = 0
opt_strat = 'ecrps-dts'
same_swm = F
samp_no = 2


#skillmod data
skill_mods = c(-0.4,0,0.4,0.8)
skill_dcy = 0.1
skill_tail = 0.5
lead = 5
lds = c(1,3,5,7)

#swm data
use_bcf=F
k_swm=106
yr_cutoff=500

#pick event
evt_order=1
nn_order=99  #plot nn_order nearest neighbor in swm and swm_cc; 
#set to 0 to plot evt_order for swm (swm_cc inherits swm), 99 to set both swm and swm_cc to evt_order

#plot setup
disp_pcnt <- 0.99
n_samp <- 10
save_plots <- F  # T to save plots to png files, F to display in R
clrs<-palette.colors()

path = paste('z:/Synthetic-Forecast-v',syn_vers,'-FIRO-DISES/',sep='')

load(paste('out/',loc,'/data_prep_rdata.RData',sep=''))
rm(hefs_forward_cumul,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid,hefs_forward_frac,hefs_inp);gc()

hefs_forward <- readRDS(paste(path,'out/',loc,'/hefs_forward_mod=',skill_mods[1],'_dcy=',skill_dcy,'_tail=',skill_tail,'.rds',sep=''))

obs_forward_all_leads_swm <- readRDS(paste('out/',loc,'/obs-fwd_use-bcf=',use_bcf,'_k=',k_swm,'_same-swm=',same_swm,'_cutoff=',yr_cutoff,'_swm.rds',sep=''))
obs_forward_all_leads_swmcc <- readRDS(paste('out/',loc,'/obs-fwd_use-bcf=',use_bcf,'_k=',k_swm,'_same-swm=',same_swm,'_cutoff=',yr_cutoff,'_swm-cc.rds',sep=''))

ixx_gen_swm <- readRDS(paste('out/',loc,'/ixx_gen_swm.rds',sep=''))
ixx_gen <- ixx_obs_forward

obs_swm <- obs_forward_all_leads_swm[samp_no,1,,1]
obs_swm_cc <- obs_forward_all_leads_swmcc[samp_no,1,,1]

obs_fwd_swm <- obs_forward_all_leads_swm[samp_no,1,,]
obs_fwd_swmcc <- obs_forward_all_leads_swmcc[samp_no,1,,]

rm(obs_forward_all_leads_swm,obs_forward_all_leads_swmcc)
gc()

hefs_ens_sset <- readRDS(paste('out/',loc,'/hefs_ens-sset_mod=',str_flatten(skill_mods,','),'_dcy=',skill_dcy,'_tail=',skill_tail,'_smp=',samp_no,'.rds',sep=''))
shefs_ens_sset <- readRDS(paste('out/',loc,'/shefs_ens-sset_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',str_flatten(skill_mods,','),'_ld=',lead,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_smp=',samp_no,'_swm.rds',sep=''))
shefscc_ens_sset <- readRDS(paste('out/',loc,'/hefs_ens-sset_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',str_flatten(skill_mods,','),'_ld=',lead,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_smp=',samp_no,'_swm-cc.rds',sep=''))

hefs_ecrps_sset <- readRDS(paste('./out/',loc,'/',site,'_mod=',str_flatten(skill_mods),'_lds=',str_flatten(lds),'_dcy=',skill_dcy,'_tail=',skill_tail,'_hefs-ecrps-sset.rds',sep=''))
shefs_ecrps_sset <- readRDS(paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',str_flatten(skill_mods),'_lds=',str_flatten(lds),'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefs-ecrps-sset.rds',sep=''))
shefscc_ecrps_sset <- readRDS(paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',str_flatten(skill_mods),'_lds=',str_flatten(lds),'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefscc-ecrps-sset.rds',sep=''))

hefs_ecrps_ss_sset <- readRDS(paste('./out/',loc,'/',site,'_mod=',str_flatten(skill_mods),'_lds=',str_flatten(lds),'_dcy=',skill_dcy,'_tail=',skill_tail,'_hefs-ecrps-ss-sset.rds',sep=''))
shefs_ecrps_ss_sset <- readRDS(paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',str_flatten(skill_mods),'_lds=',str_flatten(lds),'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefs-ecrps-ss-sset.rds',sep=''))
shefscc_ecrps_ss_sset <- readRDS(paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',str_flatten(skill_mods),'_lds=',str_flatten(lds),'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefscc-ecrps-ss-sset.rds',sep=''))

source('./src/forecast_verification_functions.R')

#///////////////////////////////////////////////////////////////////////////////////////////////////////////
#Ensemble plots
keysite = which(site_names==site)

hefs_fwd <- hefs_forward[1,keysite,,,]
obs_forward_all_leads_hind <- obs_forward_all_leads_hind[keysite,,]
obs_hefs <- obs[ixx_obs%in%ixx_hefs,keysite]
rm(hefs_dly_out,hefs_forward);gc()

hind_evt_idx <- order(obs_hefs,decreasing=T)[evt_order]
hind_evt_date <- ixx_hefs[hind_evt_idx]

#################syn-HEFS + SWM_nocc ensemble plots################################
#pick event
if(nn_order!=0){
  diff_vec <- abs(obs_hefs[hind_evt_idx] - obs_swm)
  swm_evt_idx <- order(diff_vec,decreasing=F)[nn_order]}
if(nn_order==0|nn_order==99){
  swm_evt_idx <- order(obs_swm,decreasing=T)[evt_order]}
swm_evt_date <- ixx_gen_swm[swm_evt_idx]

shefs_plts = vector('list',length(skill_mods))
show_x = F
ylm_scale = 2.1
labs = letters[1:length(skill_mods)]

hefs_plt_dt <- data.table(x=0:15,hefs=rbind(rep(obs_hefs[hind_evt_idx-lead],dim(hefs_fwd)[1]),hefs_ens_sset[1,,]),
                          obs=c(obs_hefs[hind_evt_idx-lead],obs_forward_all_leads_hind[hind_evt_idx-lead,]))

for(i in 1:length(skill_mods)){
  shefs_plt_dt <- data.table(x=0:15,hefs=rbind(rep(obs_swm[swm_evt_idx-lead],dim(hefs_fwd)[1]),shefs_ens_sset[i,,]),
                             obs=c(obs_fwd_swm[swm_evt_idx-lead,]))
  
  ylm = ylm_scale * max(hefs_plt_dt$obs)
  #convert to long form
  df <- melt(shefs_plt_dt,id='x')
  
  #assign long form df to 'obs' and 'hefs' groups
  df$col <- factor(ifelse(df$variable=="obs",1,2),levels=1:2,labels=c("obs","shefs"))
  
  hefs_ens_plt<-ggplot(df)+theme_minimal()+
    #geom_smooth(mapping=aes(x=x,y=value,group=variable,color=grp,size=grp,alpha=grp),se=F,span=.15)+
    geom_line(mapping=aes(x=x,y=value,group=variable,color=col,linewidth=col,alpha=col))+
    scale_color_manual(values=c('obs'=clrs[[4]],'shefs'=clrs[[9]]), labels=c(TeX('$swm$'),TeX('$sHEFS$')))+
    scale_linewidth_manual(values=c('obs'=1.5,'shefs'=.75), labels=c(TeX('$swm$'),TeX('$sHEFS$')))+
    scale_alpha_manual(values=c('obs'=1,'shefs'=.25), labels=c(TeX('$swm$'),TeX('$sHEFS$')))+
    geom_vline(xintercept = lead,linetype='dotted',linewidth=0.25)+
    annotate('text',x=lead,y=0.95*ylm,label=swm_evt_date,size=4,hjust=-0.25)+
    annotate('text',x=0,y=0.95*ylm,label=paste(labs[i],')',sep=''),size=6,,hjust=0)+
    labs(x='forecast lead (days)',y='flow (kcfs)',title=paste('CRPS_mod:',skill_mods[i]))+
    coord_cartesian(xlim=c(0,15),ylim=c(0,ylm),expand=F)+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(size=1))+
    {if(i==1)
      theme(legend.position = 'inside',
            legend.position.inside = c(.8,.7),
            plot.title = element_text(hjust = 0.5),
            legend.box = 'horizontal',
            legend.key.spacing.y = unit(0.01,'cm'),
            legend.key.spacing.x = unit(0.01,'cm'),
            legend.title=element_blank(),
            legend.text = element_text(size=12))}+
    {if(i!=1)
      theme(legend.position = 'none',
            plot.title = element_text(hjust = 0.5),
            axis.text.y=element_blank(),
            axis.title.y=element_blank())}+
    {if(show_x==F)
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank())}
  
  shefs_plts[[i]] <- hefs_ens_plt}

#shefs_gplot<-marrangeGrob(shefs_plts,nrow=1,ncol=length(skill_mods),top='')
#shefs_gplot

#ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ind_plots/shefs-skillmod_swm_ens-plot_evt=',evt_order,'_nn=',nn_order,'_mods=',str_flatten(skill_mods,collapse='-'),'_samp=',samp_no,'_ld=',lead,'.png',sep=''),shefs_gplot,dpi=320,width=9,height=3,unit='in')


#################syn-HEFS + SWM_cc ensemble plots################################
#pick event
if(nn_order==99){
  swm_evt_idx <- order(obs_swm_cc,decreasing=T)[evt_order]}
swm_evt_date <- ixx_gen_swm[swm_evt_idx]

shefscc_plts = vector('list',length(skill_mods))
show_x = T
labs = letters[(length(skill_mods)*1+1):(length(skill_mods)*2)]

for(i in 1:length(skill_mods)){
  shefs_plt_dt <- data.table(x=0:15,hefs=rbind(rep(obs_swm_cc[swm_evt_idx-lead],dim(hefs_fwd)[1]),shefscc_ens_sset[i,,]),
                             obs=c(obs_fwd_swmcc[swm_evt_idx-lead,]))
  
  ylm = ylm_scale * max(hefs_plt_dt$obs)
  #convert to long form
  df <- melt(shefs_plt_dt,id='x')
  
  #assign long form df to 'obs' and 'hefs' groups
  df$col <- factor(ifelse(df$variable=="obs",1,2),levels=1:2,labels=c("obs","shefs"))
  
  hefs_ens_plt<-ggplot(df)+theme_minimal()+
    #geom_smooth(mapping=aes(x=x,y=value,group=variable,color=grp,size=grp,alpha=grp),se=F,span=.15)+
    geom_line(mapping=aes(x=x,y=value,group=variable,color=col,linewidth=col,alpha=col))+
    scale_color_manual(values=c('obs'=clrs[[2]],'shefs'=clrs[[9]]), labels=c(TeX('$swm_{4C}$'),TeX('$sHEFS$')))+
    scale_linewidth_manual(values=c('obs'=1.5,'shefs'=.75), labels=c(TeX('$swm_{4C}$'),TeX('$sHEFS$')))+
    scale_alpha_manual(values=c('obs'=1,'shefs'=.25), labels=c(TeX('$swm_{4C}$'),TeX('$sHEFS$')))+
    geom_vline(xintercept = lead,linetype='dotted',linewidth=0.25)+
    annotate('text',x=lead,y=0.95*ylm,label=swm_evt_date,size=4,hjust=-0.25)+
    annotate('text',x=0,y=0.95*ylm,label=paste(labs[i],')',sep=''),size=6,hjust=0)+
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
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank())}
  
  shefscc_plts[[i]] <- hefs_ens_plt}

#shefscc_gplot<-marrangeGrob(shefscc_plts,nrow=1,ncol=length(skill_mods),top='')
#shefscc_gplot

#ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ind_plots/shefs-skillmod_swm-cc_ens-plot_evt=',evt_order,'_nn=',nn_order,'_mods=',str_flatten(skill_mods,collapse='-'),'_samp=',samp_no,'_ld=',lead,'.png',sep=''),shefscc_gplot,dpi=320,width=9,height=3,unit='in')

#///////////////////////////////////////////////////////////////////////////////////////////////////////////

#########################################eCRPS boxplot##############################################
ecrps_plts = vector('list',length(lds))
show_x = T
labs = letters[(length(skill_mods)*2+1):(length(skill_mods)*3)]

for(l in 1:length(lds)){
  shefs_mat <- array(NA,c(dim(shefs_ecrps_ss_sset)[2]*dim(shefs_ecrps_ss_sset)[3],length(skill_mods)))
  shefscc_mat <- array(NA,c(dim(shefscc_ecrps_ss_sset)[2]*dim(shefscc_ecrps_ss_sset)[3],length(skill_mods)))
  hefs_mat <- array(NA,c(dim(hefs_ecrps_ss_sset)[3],length(skill_mods)))
  for(i in 1:length(skill_mods)){
    shefs_mat[,i]<-as.vector(shefs_ecrps_ss_sset[i,l,])
    shefscc_mat[,i]<-as.vector(shefscc_ecrps_ss_sset[i,l,])
    hefs_mat[,i]<-hefs_ecrps_ss_sset[i,l,]
  }
  
  ecrps_dt <- data.table(x=1,hefs=hefs_mat,syn=shefs_mat,syncc=shefscc_mat)
  
  df <- melt(ecrps_dt,id='x')
  
  df$x <- rep(rep(c(1:length(skill_mods)),3),each=length(ecrps_dt$x))
  id_vec <- as.vector(df$variable)
  id_vec[str_detect(as.vector(df$variable),'hefs')]<-1
  id_vec[str_detect(as.vector(df$variable),'syn')]<-2
  id_vec[str_detect(as.vector(df$variable),'syncc')]<-3
  df$grp <- factor(x=as.numeric(id_vec),levels=1:3,labels=c("hefs","shefs",'shefscc'))
  
  xlbs = paste(skill_mods)
  
  ecrps_bplt <- ggplot(df)+theme_minimal()+
    geom_boxplot(aes(x=factor(x),y=value,fill=factor(grp),color=factor(grp)),width=.5,outlier.size = .5,outlier.color = clrs[[9]],outlier.alpha = 0.25)+
    scale_fill_manual(name='',values=c('hefs'=str_c(clrs[[9]],'50'),'shefs'=str_c(clrs[[9]],'50'),'shefscc'=str_c(clrs[[9]],'50')), labels=c(TeX('$HEFS$'),TeX('$sHEFS_{swm}$'),TeX('$sHEFS_{swm_{4C}}$')))+
    scale_color_manual(name='',values=c('hefs'=clrs[[1]],'shefs'=clrs[[4]],'shefscc'=clrs[[2]]), labels=c(TeX('$HEFS$'),TeX('$sHEFS_{swm}$'),TeX('$sHEFS_{swm_{4C}}$')))+
    scale_x_discrete(breaks=1:length(skill_mods),labels=xlbs,name=TeX('$CRPS_{mod}$'))+
    scale_y_continuous(name = 'CRPS-SS')+
    coord_cartesian(ylim = c(0,1))+
    labs(title=paste('ld',lds[l]))+
    annotate('text',x=0,y=0.95,label=paste(labs[l],')',sep=''),size=6,hjust=0)+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(linewidth=1),
          plot.title = element_text(hjust=0.5,size=16))+
    {if(l==1) 
      theme(legend.position = 'inside',
            legend.position.inside = c(.8,.35),
            legend.background = element_rect(fill = "#FFFFFF70",color='white'),
            legend.text = element_text(size=12))}+
    {if(l!=1)
      theme(legend.position = 'none',
            axis.text.y=element_blank(),
            axis.title.y=element_blank())}+
    {if(show_x==F)
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank())}
  
  
  ecrps_plts[[l]] <- ecrps_bplt}

#mat<-matrix(1:length(lds),ncol=length(lds))
#ecrps_gplot<-marrangeGrob(ecrps_plts,nrow=1,ncol=length(lds),top='',widths=c(1.15,rep(1,(length(lds)-1))))
#ecrps_gplot

#ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ind_plots/',loc,'-',disp_site,'_ecrps-ss-skillmod-bplots_pcntile=',disp_pcnt,'_lds=',str_flatten(lds,collapse='-'),'_mod=',str_flatten(skill_mods,collapse='-'),'.png',sep=''),ecrps_gplot,dpi=320,width=2.5*length(lds),height=3,unit='in')

layout_mat <- matrix(1:(3*length(lds)),nrow=3,ncol=length(lds),byrow=T)
comb_plot<-marrangeGrob(c(shefs_plts,shefscc_plts,ecrps_plts),nrow=3,ncol=length(lds),top='',layout_matrix = layout_mat, widths = c(1.15,1,1,1),heights = c(1,1,1.15))

ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/forecast_plots/combined_ens-plot_verification-plot_condensed_skillmods_fig5_evt=',evt_order,'_nn=',nn_order,'_lds=',str_flatten(leads,collapse='-'),'_samp=',samp_no,'_mods=',str_flatten(skill_mods,'-'),'.png',sep=''),comb_plot,dpi=320,width=9,height=8.5,unit='in')

##############################################END#############################################################