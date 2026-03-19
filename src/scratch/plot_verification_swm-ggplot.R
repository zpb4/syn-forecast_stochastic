

rm(list=ls())
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

setwd('z:/Synthetic-Forecast_Verification/')
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

skill_mod = 0.8
skill_dcy = 0.1
skill_tail = 0.5

#swm data
use_bcf=F
k_swm=106
set.seed(1)

#plot setup
disp_pcnt <- 0.99

ix_swm<-seq(as.Date('2001-01-01'),as.Date('3008-01-08'),'day')
ix2_swm<-as.POSIXlt(ix_swm)

wy_fun<-function(date_vec){
  wy_vec <- date_vec$year
  wy_vec[date_vec$mo%in%c(9,10,11)] <- wy_vec[date_vec$mo%in%c(9,10,11)]+1
  date_vec_wy <- date_vec
  date_vec_wy$year <- wy_vec
  return(date_vec_wy)
}

#plot setup
disp_pcnt <- 0.99
n_samp <- 10
save_plots <- F  # T to save plots to png files, F to display in R
clrs<-palette.colors()

hefs_ecrps_vec<-readRDS(paste('./data/',loc,'-',site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_hefs-ecrps-vec.rds',sep=''))
hefs_rank_vec<-readRDS(paste('./data/',loc,'-',site,'_',cal_val_setup,'_pcnt=',opt_pcnt,'_hefs-rank-vec.rds',sep=''))

shefs_ecrps_vec<-readRDS(paste('./data/',loc,'-',site,'_',cal_val_setup,'_pcnt=',disp_pcnt,'_shefs-ecrps-vec.rds',sep=''))
shefs_rank_vec<-readRDS(paste('./data/',loc,'-',site,'_',cal_val_setup,'_pcnt=',disp_pcnt,'_shefs-rank-vec.rds',sep=''))

shefsref_ecrps_vec<-readRDS(paste('./data/',loc,'-',site,'_',cal_val_setup,'_pcnt=',disp_pcnt,'_shefsref-ecrps-vec.rds',sep=''))
shefsref_rank_vec<-readRDS(paste('./data/',loc,'-',site,'_',cal_val_setup,'_pcnt=',disp_pcnt,'_shefsref-rank-vec.rds',sep=''))

shefscc_ecrps_vec<-readRDS(paste('./data/',loc,'-',site,'_',cal_val_setup,'_pcnt=',disp_pcnt,'_shefscc-ecrps-vec.rds',sep=''))
shefscc_rank_vec<-readRDS(paste('./data/',loc,'-',site,'_',cal_val_setup,'_pcnt=',disp_pcnt,'_shefscc-rank-vec.rds',sep=''))

source('./src/forecast_verification_functions.R')
path = paste('z:/Synthetic-Forecast-v',syn_vers,'-FIRO-DISES/',sep='')
load(paste(path,'out/',loc,'/data_prep_rdata.RData',sep=''))

#///////////////////////////////////////////////////////////////////////////////////////////////////////////

#########################################eCRPS boxplot##############################################
lds = c(1,3,5,7)
ecrps_bplts <- vector('list',length(lds))

for(i in 1:length(lds)){
  #shefs_vec = apply(shefs_ecrps_vec[,,lds],3,as.vector)
  #shefscc_vec = apply(shefscc_ecrps_vec[,,lds],3,as.vector)
  #shefsref_vec = apply(shefsref_ecrps_vec[,,lds],3,as.vector)
  #ecrps_dt <- data.table(x=1,hefs=hefs_ecrps_vec[,lds],syn=shefs_vec,syncc=shefscc_vec,synref=shefsref_vec)
  shefs_vec <- as.vector(shefs_ecrps_vec[,,lds[i]])
  shefscc_vec <- as.vector(shefscc_ecrps_vec[,,lds[i]])
  ecrps_dt <- data.table(x=1,hefs=hefs_ecrps_vec[,lds[i]],syn=shefs_vec,syncc=shefscc_vec)

  df <- melt(ecrps_dt,id='x')

  #df$x <- rep(c(1:4,1:4,1:4,1:4),each=length(ecrps_dt$x))
  id_vec <- as.vector(df$variable)
  id_vec[str_detect(as.vector(df$variable),'hefs')]<-1
  id_vec[str_detect(as.vector(df$variable),'syn')]<-2
  id_vec[str_detect(as.vector(df$variable),'syncc')]<-3
  df$grp <- factor(x=as.numeric(id_vec),levels=1:3,labels=c("hefs","shefs","shefscc"))

  xlbs = paste('ld',lds[i])

  ecrps_bplt <- ggplot(df)+theme_minimal()+
    geom_boxplot(aes(x=factor(x),y=value,fill=factor(grp)),width=.5,outlier.size = .5,outlier.color = clrs[[9]],outlier.alpha = 0.25)+
    scale_fill_manual(name='',values=c('hefs'=clrs[[9]],'shefs'=clrs[[4]],'shefscc'=clrs[[2]]), labels=c(TeX('$HEFS$'),TeX('$sHEFS_{SWM}$'),TeX('$sHEFS_{SWM4C}$')))+
    scale_x_discrete(breaks=1,labels='',name='')+
    scale_y_continuous(name = 'eCRPS')+
    coord_cartesian(ylim = c(0,80))+
    theme(legend.position = 'inside',
          legend.position.inside = c(1,.9),
          axis.text.y=element_text(size=10),
          axis.text.x=element_blank(),
          axis.title=element_text(size=12),
          legend.text = element_text(size=10))+
    {if(i!=1)
      theme(legend.position = 'none',
            axis.text.y=element_blank(),
            axis.title.y=element_blank())}

  ecrps_bplts[[i]]<-ecrps_bplt}

mat<-matrix(1:length(lds),ncol=length(lds))
ecrps_gplot<-marrangeGrob(ecrps_bplts,nrow=1,ncol=length(lds),top='',widths=c(1.15,rep(1,(length(lds)-1))))
ecrps_gplot

ggsave(paste('h:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ind_plots/',loc,'-',disp_site,'_ecrps_bplots_pcntile=',disp_pcnt,'_lds=',str_flatten(lds,collapse='-'),'.png',sep=''),ecrps_gplot,dpi=320,width=2.5*length(lds),height=3,unit='in')

#########################################cumulative rank histogram plot##############################################
#calculations
nsamp = dim(shefs_rank_vec)[1]
refsamp = dim(shefsref_rank_vec)[1]
show_x = T

hefs_cumul_frac <- array(NA,c(n_ens+2,length(lds)))
shefs_cumul_frac <- array(NA,c(nsamp,n_ens+2,length(lds)))
shefscc_cumul_frac <- array(NA,c(nsamp,n_ens+2,length(lds)))
shefsref_cumul_frac <- array(NA,c(refsamp,n_ens+2,length(lds)))

for(i in 1:length(lds)){
  hefs_count <- hist(hefs_rank_vec[,lds[i]],breaks=seq(0.5,n_ens+1.5),plot = FALSE)
  hefs_cumul_frac[,i] <- c(0,roll_sum(hefs_count$counts)/length(hefs_rank_vec[,i]))
  for(s in 1:nsamp){
    shefs_count<-hist(shefs_rank_vec[s,,lds[i]],breaks=seq(0.5,n_ens+1.5),plot = FALSE)
    shefs_cumul_frac[s,,i]<-c(0,roll_sum(shefs_count$counts)/length(hefs_rank_vec[,i]))
    shefs_count<-hist(shefscc_rank_vec[s,,lds[i]],breaks=seq(0.5,n_ens+1.5),plot = FALSE)
    shefscc_cumul_frac[s,,i]<-c(0,roll_sum(shefs_count$counts)/length(hefs_rank_vec[,i]))
  }
  for(s in 1:refsamp){
    shefs_count<-hist(shefsref_rank_vec[s,,lds[i]],breaks=seq(0.5,n_ens+1.5),plot = FALSE)
    shefsref_cumul_frac[s,,i]<-c(0,roll_sum(shefs_count$counts)/length(hefs_rank_vec[,i]))
  }
}

#bounds plot
shefs_rank_upr<-apply(shefs_cumul_frac,c(2,3),max)
shefs_rank_lwr<-apply(shefs_cumul_frac,c(2,3),min)

shefsref_rank_upr<-apply(shefsref_cumul_frac,c(2,3),max)
shefsref_rank_lwr<-apply(shefsref_cumul_frac,c(2,3),min)

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
    scale_fill_manual(values=c('shefs'=clrs[[4]],'shefscc'=clrs[[2]]), labels=c(TeX('$sHEFS_{SWM}$'),TeX('$sHEFS_{SWM4C}$')))+
    labs(x='ensemble rank',y='cumulative fraction',caption=paste('ld',lds[i]))+
    coord_cartesian(xlim=c(0,n_ens+2),ylim=c(-0.01,1),expand=F)+
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
            legend.position.inside = c(.3,.75),
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
rank_gplot

ggsave(paste('h:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ind_plots/',loc,'-',disp_site,'_rankhist-swm-shadeplot_pcntile=',disp_pcnt,'_lds=',str_flatten(lds,collapse='-'),'.png',sep=''),rank_gplot,dpi=320,width=2.75*length(lds),height=3,unit='in')

mat<-matrix(1:(length(lds)*2),ncol=length(lds),byrow=T)
comb_gplot<-marrangeGrob(c(ecrps_bplts,rank_plts),nrow=2,ncol=length(lds),top='',layout_matrix = mat)
ggsave(paste('h:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/',loc,'-',disp_site,'_ecrps-rankhist_pcntile=',disp_pcnt,'_lds=',str_flatten(lds,collapse='-'),'.png',sep=''),comb_gplot,dpi=320,width=9,height=6,unit='in')



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

rank_plts_ln <- vector('list',length(lds))
#Line plot
#plot
show_x = T

for(i in 1:length(lds)){
  rank_dt <- data.table(x=0:(n_ens+1),syn=t(shefs_cumul_frac[,,i]),syncc=t(shefscc_cumul_frac[,,i]),hefs=hefs_cumul_frac[,i])
  
  #convert to long form
  df <- melt(rank_dt,id='x')
  
  #assign long form df to 'obs' and 'hefs' groups
  id_vec <- as.vector(df$variable)
  id_vec[str_detect(as.vector(df$variable),'hefs')]<-1
  id_vec[str_detect(as.vector(df$variable),'syn')]<-2
  id_vec[str_detect(as.vector(df$variable),'syncc')]<-3
  df$grp <- factor(x=as.numeric(id_vec),levels=1:3,labels=c("hefs","shefs","shefscc"))
  
  hefs_ens_plt<-ggplot(df)+theme_minimal()+
    geom_line(mapping=aes(x=x,y=value,group=variable,color=grp,linewidth=grp,alpha=grp))+
    scale_color_manual(values=c('hefs'='gray40','shefs'=clrs[[4]],'shefscc'=clrs[[2]]), labels=c(TeX('$HEFS$'),TeX('$sHEFS$'),TeX('$sHEFS_{4C}$')))+
    scale_linewidth_manual(values=c('hefs'=1.5,'shefs'=.75,'shefscc'=.75))+
    scale_alpha_manual(values=c('hefs'=1,'shefs'=.15,'shefscc'=.15))+
    labs(x='ensemble rank',y='cumulative fraction')+
    coord_cartesian(xlim=c(0,n_ens+2),ylim=c(-0.01,1),expand=F)+
    #labs(title=paste(lds[i],'hr'),tag=paste(letters[i],')',sep=''))+
    labs(title=paste('ld',lds[i]))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(size=1),
          panel.grid.major.x = element_line(size=1),
          plot.title = element_text(hjust=0.5,size=16))+
    {if(i==1)
      theme(legend.position = 'inside',
            legend.position.inside = c(.2,.9),
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
  
  
  rank_plts_ln[[i]] <- hefs_ens_plt}


mat<-matrix(1:length(lds),ncol=length(lds))
rank_ln_gplot<-marrangeGrob(rank_plts_ln,nrow=1,ncol=length(lds),top='',widths=c(1.15,rep(1,(length(lds)-1))))
rank_ln_gplot

ggsave(paste('h:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ind_plots/',loc,'-',disp_site,'_rankhist-lineplot_pcntile=',disp_pcnt,'_lds=',str_flatten(lds,collapse='-'),'.png',sep=''),rank_ln_gplot,dpi=320,width=2.75*length(lds),height=3,unit='in')


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

ref_plts <- vector('list',length(lds))

#plot ref syn-forecast results
show_x = T

for(i in 1:length(lds)){
  rank_dt <- data.table(x=0:(n_ens+1),ref_upr=shefsref_rank_upr[,i],ref_lwr=shefsref_rank_lwr[,i],hefs=hefs_cumul_frac[,i])
  
  hefs_ens_plt<-ggplot(rank_dt)+theme_minimal()+
    geom_ribbon(mapping=aes(x=x,ymax=ref_upr,ymin=ref_lwr,fill='shefs'),alpha=0.15)+
    geom_line(mapping=aes(x=x,y=hefs,color='hefs'),linewidth=1)+
    scale_color_manual(values=c('hefs'=clrs[[9]]), labels=c(TeX('$HEFS$')))+
    scale_fill_manual(values=c('shefs'=clrs[[3]]), labels=c(TeX('$sHEFS_{ref}$')))+
    labs(x='ensemble rank',y='cumulative fraction')+
    coord_cartesian(xlim=c(0,n_ens+2),ylim=c(-0.01,1),expand=F)+
    labs(title=paste('ld',lds[i]))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          panel.grid.major.y = element_line(size=.75),
          panel.grid.major.x = element_line(size=.75),
          plot.title = element_text(hjust=0.5,size=16))+
    {if(i==length(lds))
      theme(legend.position = 'inside',
            legend.position.inside = c(.25,.75),
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
  
  
  ref_plts[[i]] <- hefs_ens_plt}

mat<-matrix(1:length(lds),ncol=length(lds))
ref_gplot<-marrangeGrob(ref_plts,nrow=1,ncol=length(lds),top='',widths=c(1.15,rep(1,(length(lds)-1))))
ref_gplot

ggsave(paste('h:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ind_plots/',loc,'-',disp_site,'_rankhist-ref-shadeplot_pcntile=',disp_pcnt,'_lds=',str_flatten(lds,collapse='-'),'.png',sep=''),ref_gplot,dpi=320,width=2.75*length(lds),height=3,unit='in')


#ggsave(paste('h:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/',loc,'-',disp_site,'_ecrps-rankhist_pcntile=',disp_pcnt,'_lds=',str_flatten(lds,collapse='-'),'.png',sep=''),comb_plot,dpi=320,width=9,height=9,unit='in')

#########################################################END########################################################
