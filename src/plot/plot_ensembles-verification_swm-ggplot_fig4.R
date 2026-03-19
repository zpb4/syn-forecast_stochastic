

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
library(latex2exp)
#///////////////////////////////////////////////////////////////////////////////////////////////////////////
#Data setup
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

#plot setup
disp_pcnt <- 0.99
n_samp <- 10
save_plots <- F  # T to save plots to png files, F to display in R
clrs<-palette.colors()
#show_col(clrs)

load(paste('out/',loc,'/data_prep_rdata.RData',sep=''))
rm(hefs_forward_cumul,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid,hefs_forward_frac,hefs_inp);gc()

hefs_forward <- readRDS(paste('out/',loc,'/hefs_forward_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'.rds',sep=''))
syn_hefs_forward <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_swm.rds',sep=''))
syn_hefs_forward_cc <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_swm-cc.rds',sep=''))

obs_forward_all_leads_swm <- readRDS(paste('out/',loc,'/obs-fwd_use-bcf=',use_bcf,'_k=',k_swm,'_same-swm=',same_swm,'_cutoff=',yr_cutoff,'_swm.rds',sep=''))
obs_forward_all_leads_swmcc <- readRDS(paste('out/',loc,'/obs-fwd_use-bcf=',use_bcf,'_k=',k_swm,'_same-swm=',same_swm,'_cutoff=',yr_cutoff,'_swm-cc.rds',sep=''))

ixx_gen_swm <- readRDS(paste('out/',loc,'/ixx_gen_swm.rds',sep=''))
ixx_gen <- ixx_obs_forward

shefs_fwd_swm <- syn_hefs_forward[samp_no,1,,,]
shefs_fwd_swmcc <- syn_hefs_forward_cc[samp_no,1,,,]

obs_swm <- obs_forward_all_leads_swm[samp_no,1,,1]
obs_swm_cc <- obs_forward_all_leads_swmcc[samp_no,1,,1]

obs_fwd_swm <- obs_forward_all_leads_swm[samp_no,1,,]
obs_fwd_swmcc <- obs_forward_all_leads_swmcc[samp_no,1,,]

rm(syn_hefs_forward,syn_hefs_forward_cc,obs_forward_all_leads_swm,obs_forward_all_leads_swmcc)
gc()

hefs_ecrps_vec<-readRDS(paste('./out/',loc,'/',site,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_hefs-ecrps-vec.rds',sep=''))
hefs_ecrps_ss_vec<-readRDS(paste('./out/',loc,'/',site,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_hefs-ecrps-ss-vec.rds',sep=''))
hefs_rank_vec<-readRDS(paste('./out/',loc,'/',site,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_hefs-rank-vec.rds',sep=''))

shefs_ecrps_vec<-readRDS(paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefs-ecrps-vec.rds',sep=''))
shefs_ecrps_ss_vec<-readRDS(paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefs-ecrps-ss-vec.rds',sep=''))
shefs_rank_vec<-readRDS(paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefs-rank-vec.rds',sep=''))

shefscc_ecrps_vec<-readRDS(paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefscc-ecrps-vec.rds',sep=''))
shefscc_ecrps_ss_vec<-readRDS(paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefscc-ecrps-ss-vec.rds',sep=''))
shefscc_rank_vec<-readRDS(paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefscc-rank-vec.rds',sep=''))


source('./src/forecast_verification_functions.R')

#///////////////////////////////////////////////////////////////////////////////////////////////////////////
#Ensemble plots
keysite = which(site_names==site)

hefs_fwd <- hefs_forward[1,keysite,,,]
obs_forward_all_leads_hind <- obs_forward_all_leads_hind[keysite,,]
obs_hefs <- obs[ixx_obs%in%ixx_hefs,keysite]


#################HEFS ensemble plots################################
hind_evt_idx <- order(obs_hefs,decreasing=T)[evt_order]
hind_evt_date <- ixx_hefs[hind_evt_idx]

leads = c(1,4,7)
hefs_plts = vector('list',length(leads))
ylm_scale = 2.35
show_x = F
labs = letters[1:length(leads)]

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
    scale_color_manual(values=c('obs'=clrs[[1]],'hefs'=clrs[[9]]), labels=c(TeX('$obs$'),TeX('$HEFS$')))+
    scale_linewidth_manual(values=c('obs'=1.5,'hefs'=.75), labels=c(TeX('$obs$'),TeX('$HEFS$')))+
    scale_alpha_manual(values=c('obs'=1,'hefs'=.35), labels=c(TeX('$obs$'),TeX('$HEFS$')))+
    geom_vline(xintercept = leads[i],linetype='dotted',linewidth=0.25)+
    {if(i==length(leads))
        annotate('text',x=leads[i],y=0.95*ylm,label=hind_evt_date,size=4,hjust=-0.25)}+
    annotate('text',x=0,y=0.95*ylm,label=paste(labs[i],')',sep=''),size=6,,hjust=0)+
    labs(x='forecast lead (days)',y='flow (kcfs)',title=paste('Lead',leads[i]))+
    coord_cartesian(xlim=c(0,15),ylim=c(0,ylm),expand=F)+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(size=1))+
    {if(i==1)
      theme(legend.position = 'inside',
            legend.position.inside = c(.6,.7),
            legend.box = 'horizontal',
            legend.key.spacing.y = unit(0.01,'cm'),
            legend.key.spacing.x = unit(0.01,'cm'),
            legend.title=element_blank(),
            plot.title = element_text(hjust=0.5,size=16),
            legend.text = element_text(size=12))}+
    {if(i!=1)
      theme(legend.position = 'none',
            axis.text.y=element_blank(),
            plot.title = element_text(hjust=0.5,size=16),
            axis.title.y=element_blank())}+
    {if(show_x==F)
      theme(axis.text.x=element_blank(),
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
swm_evt_date <- ixx_gen_swm[swm_evt_idx]

shefs_plts = vector('list',length(leads))
show_x = F
labs = letters[(length(leads)+1):(length(leads)*2)]

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
    scale_color_manual(values=c('obs'=clrs[[4]],'shefs'=clrs[[9]]), labels=c(TeX('$s.hyd_{hist}$'),TeX('$s.HEFS_{hist}$')))+
    scale_linewidth_manual(values=c('obs'=1.5,'shefs'=.75), labels=c(TeX('$s.hyd_{hist}$'),TeX('$s.HEFS_{hist}$')))+
    scale_alpha_manual(values=c('obs'=1,'shefs'=.25), labels=c(TeX('$s.hyd_{hist}$'),TeX('$s.HEFS_{hist}$')))+
    geom_vline(xintercept = leads[i],linetype='dotted',linewidth=0.25)+
    {if(i==length(leads))
      annotate('text',x=leads[i],y=0.95*ylm,label=hind_evt_date,size=4,hjust=-0.25)}+
    annotate('text',x=0,y=0.95*ylm,label=paste(labs[i],')',sep=''),size=6,,hjust=0)+
    labs(x='forecast lead (days)',y='flow (kcfs)')+
    coord_cartesian(xlim=c(0,15),ylim=c(0,ylm),expand=F)+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(size=1))+
    {if(i==1)
      theme(legend.position = 'inside',
            legend.position.inside = c(.65,.7),
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
  
  shefs_plts[[i]] <- hefs_ens_plt}

shefs_gplot<-marrangeGrob(shefs_plts,nrow=1,ncol=length(leads),top='')
#shefs_gplot

ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ind_plots/shefs_swm_ens-plot_evt=',evt_order,'_nn=',nn_order,'_lds=',str_flatten(leads,collapse='-'),'_samp=',samp_no,'_skillmod=',skill_mod,'.png',sep=''),shefs_gplot,dpi=320,width=9,height=3,unit='in')


#################syn-HEFS + SWM_cc ensemble plots################################
#pick event
if(nn_order==99){
  swm_evt_idx <- order(obs_swm_cc,decreasing=T)[evt_order]}
swm_evt_date <- ixx_gen_swm[swm_evt_idx]

shefscc_plts = vector('list',length(leads))
show_x = T
labs = letters[(length(leads)*2+1):(length(leads)*3)]

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
    scale_color_manual(values=c('obs'=clrs[[2]],'shefs'=clrs[[9]]), labels=c(TeX('$s.hyd_{4C}$'),TeX('$s.HEFS_{4C}$')))+
    scale_linewidth_manual(values=c('obs'=1.5,'shefs'=.75), labels=c(TeX('$s.hyd_{4C}$'),TeX('$s.HEFS_{4C}$')))+
    scale_alpha_manual(values=c('obs'=1,'shefs'=.25), labels=c(TeX('$s.hyd_{4C}$'),TeX('$s.HEFS_{4C}$')))+
    geom_vline(xintercept = leads[i],linetype='dotted',linewidth=0.25)+
    {if(i==length(leads))
      annotate('text',x=leads[i],y=0.95*ylm,label=hind_evt_date,size=4,hjust=-0.25)}+
    annotate('text',x=0,y=0.95*ylm,label=paste(labs[i],')',sep=''),size=6,hjust=0)+
    labs(x='forecast lead (days)',y='flow (kcfs)')+
    coord_cartesian(xlim=c(0,15),ylim=c(0,ylm),expand=F)+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(size=1))+
    {if(i==1)
      theme(legend.position = 'inside',
            legend.position.inside = c(.65,.7),
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

shefscc_gplot<-marrangeGrob(shefscc_plts,nrow=1,ncol=length(leads),top='')
#shefscc_gplot

ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ind_plots/shefs_swm-cc_ens-plot_evt=',evt_order,'_nn=',nn_order,'_lds=',str_flatten(leads,collapse='-'),'_samp=',samp_no,'_skillmod=',skill_mod,'.png',sep=''),shefscc_gplot,dpi=320,width=9,height=3,unit='in')

#///////////////////////////////////////////////////////////////////////////////////////////////////////////

#########################################eCRPS-SS boxplot##############################################
lds = leads
ecrps_bplts <- vector('list',length(lds))
labs = letters[(length(leads)*3+1):(length(leads)*4)]

for(i in 1:length(lds)){
  shefs_vec <- as.vector(shefs_ecrps_ss_vec[,,lds[i]])
  shefscc_vec <- as.vector(shefscc_ecrps_ss_vec[,,lds[i]])
  ecrps_dt <- data.table(x=1,hefs=hefs_ecrps_ss_vec[,lds[i]],syn=shefs_vec,syncc=shefscc_vec)
  
  df <- melt(ecrps_dt,id='x')
  
  id_vec <- as.vector(df$variable)
  id_vec[str_detect(as.vector(df$variable),'hefs')]<-1
  id_vec[str_detect(as.vector(df$variable),'syn')]<-2
  id_vec[str_detect(as.vector(df$variable),'syncc')]<-3
  df$grp <- factor(x=as.numeric(id_vec),levels=1:3,labels=c("hefs","shefs","shefscc"))
  
  xlbs = paste('ld',lds[i])
  
  ecrps_bplt <- ggplot(df)+theme_minimal()+
    geom_boxplot(aes(x=factor(x),y=value,fill=factor(grp),col=factor(grp)),width=.5,outlier.size = .5,outlier.color = clrs[[9]],outlier.alpha = 0.25)+
    scale_fill_manual(name='',values=c('hefs'=str_c(clrs[[9]],'50'),'shefs'=str_c(clrs[[9]],'50'),'shefscc'=str_c(clrs[[9]],'50')), labels=c(TeX('$HEFS$'),TeX('$s.HEFS_{hist}$'),TeX('$s.HEFS_{4C}$')))+
    scale_color_manual(name='',values=c('hefs'=clrs[[1]],'shefs'=clrs[[4]],'shefscc'=clrs[[2]]), labels=c(TeX('$HEFS$'),TeX('$s.HEFS_{hist}$'),TeX('$s.HEFS_{4C}$')))+
    scale_x_discrete(breaks=1,labels='',name='')+
    scale_y_continuous(name = 'CRPSS')+
    annotate('text',x=0.35,y=0.95,label=paste(labs[i],')',sep=''),size=6,hjust=0)+
    coord_cartesian(ylim = c(0,1))+
    theme(legend.position = 'inside',
          legend.position.inside = c(0.3,.3),
          legend.background = element_rect(fill = "#FFFFFF90",color='white'),
          axis.text.y=element_text(size=10),
          axis.text.x=element_blank(),
          axis.title=element_text(size=12),
          legend.text = element_text(size=10))+
    {if(i!=1)
      theme(legend.position = 'none')}+
    {if(i!=1)
      theme(axis.text.y=element_blank(),
            axis.title.y=element_blank())}
  
  ecrps_bplts[[i]]<-ecrps_bplt}

mat<-matrix(1:length(lds),ncol=length(lds))
ecrps_gplot<-marrangeGrob(ecrps_bplts,nrow=1,ncol=length(lds),top='',widths=c(1.15,rep(1,(length(lds)-1))))
#ecrps_gplot

ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ind_plots/',loc,'-',disp_site,'_ecrps-ss_bplots_pcntile=',disp_pcnt,'_lds=',str_flatten(lds,collapse='-'),'.png',sep=''),ecrps_gplot,dpi=320,width=2.5*length(lds),height=3,unit='in')

#########################################cumulative rank histogram plot##############################################
#calculations
nsamp = dim(shefs_rank_vec)[1]
show_x = T
labs = letters[(length(leads)*4+1):(length(leads)*5)]

hefs_cumul_frac <- array(NA,c(n_ens+2,length(lds)))
shefs_cumul_frac <- array(NA,c(nsamp,n_ens+2,length(lds)))
shefscc_cumul_frac <- array(NA,c(nsamp,n_ens+2,length(lds)))

for(i in 1:length(lds)){
  hefs_count <- hist(hefs_rank_vec[,lds[i]],breaks=seq(0.5,n_ens+1.5),plot = FALSE)
  hefs_cumul_frac[,i] <- c(0,roll_sum(hefs_count$counts)/length(hefs_rank_vec[,i]))
  for(s in 1:nsamp){
    shefs_count<-hist(shefs_rank_vec[s,,lds[i]],breaks=seq(0.5,n_ens+1.5),plot = FALSE)
    shefs_cumul_frac[s,,i]<-c(0,roll_sum(shefs_count$counts)/length(hefs_rank_vec[,i]))
    shefs_count<-hist(shefscc_rank_vec[s,,lds[i]],breaks=seq(0.5,n_ens+1.5),plot = FALSE)
    shefscc_cumul_frac[s,,i]<-c(0,roll_sum(shefs_count$counts)/length(hefs_rank_vec[,i]))
  }
}

#bounds plot
shefs_rank_upr<-apply(shefs_cumul_frac,c(2,3),max)
shefs_rank_lwr<-apply(shefs_cumul_frac,c(2,3),min)

shefscc_rank_upr<-apply(shefscc_cumul_frac,c(2,3),max)
shefscc_rank_lwr<-apply(shefscc_cumul_frac,c(2,3),min)

#plot swm results
show_x=T
rank_plts <- vector('list',length(lds))

for(i in 1:length(lds)){
  rank_dt <- data.table(x=0:(n_ens+1),syn_upr=shefs_rank_upr[,i],syn_lwr=shefs_rank_lwr[,i],
                        syncc_upr=shefscc_rank_upr[,i],syncc_lwr=shefscc_rank_lwr[,i],hefs=hefs_cumul_frac[,i])
  
  hefs_ens_plt<-ggplot(rank_dt)+theme_minimal()+
    geom_ribbon(mapping=aes(x=x,ymax=syn_upr,ymin=syn_lwr,fill='shefs'),alpha=0.15)+
    geom_ribbon(mapping=aes(x=x,ymax=syncc_upr,ymin=syncc_lwr,fill='shefscc'),alpha=0.15)+
    geom_line(mapping=aes(x=x,y=hefs,color='hefs'),linewidth=1)+
    scale_color_manual(values=c('hefs'=clrs[[9]]), labels=c(TeX('$HEFS$')))+
    scale_fill_manual(values=c('shefs'=clrs[[4]],'shefscc'=clrs[[2]]), labels=c(TeX('$s.HEFS_{hist}$'),TeX('$s.HEFS_{4C}$')))+
    #labs(x='ensemble rank',y='cumulative fraction',caption=paste('ld',lds[i]))+
    labs(x='ensemble rank',y='cumulative fraction')+
    coord_cartesian(xlim=c(0,n_ens+2),ylim=c(-0.01,1),expand=F)+
    annotate('text',x=0,y=0.95,label=paste(labs[i],')',sep=''),size=6,hjust=0)+
    #labs(title=paste(lds[i],'hr'),tag=paste(letters[i],')',sep=''))+
    #labs(title=paste('ld',lds[i]))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(size=.75),
          panel.grid.major.x = element_line(size=.75),
          plot.title = element_text(hjust=0.5,size=16),
          plot.caption = element_text(hjust=0.5,size=16))+
    {if(i==length(lds))
      theme(legend.position = 'inside',
            legend.position.inside = c(.45,.7),
            legend.box = 'vertical',
            legend.key.spacing.y = unit(0.01,'cm'),
            legend.key.spacing.x = unit(0.01,'cm'),
            legend.title=element_blank(),
            legend.text = element_text(size=12))}+
    {if(i!=length(lds))
      theme(legend.position = 'none')}+
    {if(i!=1)
      theme(axis.text.y=element_blank(),
            axis.title.y=element_blank())}+
    {if(show_x==F)
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank())}+
    guides(shape = guide_legend(order=2),col=guide_legend(order=1))
  
  
  rank_plts[[i]] <- hefs_ens_plt}

#rank_plts[[1]]
mat<-matrix(1:length(lds),ncol=length(lds))
rank_gplot<-marrangeGrob(rank_plts,nrow=1,ncol=length(lds),top='',widths=c(1.15,rep(1,(length(lds)-1))))
#rank_gplot

ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ind_plots/',loc,'-',disp_site,'_rankhist-swm-shadeplot_pcntile=',disp_pcnt,'_lds=',str_flatten(lds,collapse='-'),'.png',sep=''),rank_gplot,dpi=320,width=2.75*length(lds),height=3,unit='in')

layout_mat <- matrix(1:(5*length(leads)),nrow=5,ncol=length(leads),byrow=T)
comb_plot<-marrangeGrob(c(hefs_plts,shefs_plts,shefscc_plts,ecrps_bplts,rank_plts),nrow=3,ncol=length(leads),top='',layout_matrix = layout_mat, widths = c(1.15,1,1), heights = c(1.15,1,1,1,1))
ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/forecast_plots/combined_ens-plot_verification-plot_fig4_evt=',evt_order,'_nn=',nn_order,'_lds=',str_flatten(leads,collapse='-'),'_samp=',samp_no,'_skillmod=',skill_mod,'.png',sep=''),comb_plot,dpi=320,width=7,height=11,unit='in')
ggsave(paste('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/manuscript/combined_ens-plot_verification-plot_fig4_evt=',evt_order,'_nn=',nn_order,'_lds=',str_flatten(leads,collapse='-'),'_samp=',samp_no,'_skillmod=',skill_mod,'.png',sep=''),comb_plot,dpi=320,width=7,height=11,unit='in')

########################################################################

