#script to extract subsets to reduce plotting script time
library(stringr)

#setwd('z:/syn-forecast_stochastic//')
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


#skillmod data
skill_mods = c(-0.4,0,0.4,0.8)
skill_dcy = 0.1
skill_tail = 0.5
lead = 5  #plotting lead for ensembles
lds = c(1,3,5,7)  #plotting lead for CRPS

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

load(paste('out/',loc,'/data_prep_rdata.RData',sep=''))
rm(hefs_forward_cumul,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid,hefs_forward_frac,hefs_inp);gc()

ixx_gen_swm <- readRDS(paste('out/',loc,'/ixx_gen_swm.rds',sep=''))
ixx_gen <- ixx_obs_forward

hefs_ens_sset <- array(NA,c(length(skill_mods),leads,n_ens))
shefs_ens_sset <- array(NA,c(length(skill_mods),leads,n_ens))
shefscc_ens_sset <- array(NA,c(length(skill_mods),leads,n_ens))

hefs_ecrps_samp<-readRDS(paste('./out/',loc,'/',site,'_mod=',skill_mods[1],'_dcy=',skill_dcy,'_tail=',skill_tail,'_hefs-ecrps-vec.rds',sep=''))
shefs_ecrps_samp<-readRDS(paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mods[1],'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefs-ecrps-vec.rds',sep=''))

hefs_len <- length(hefs_ecrps_samp[,lds[i]])
shefs_len <- length(as.vector(shefs_ecrps_samp[,,lds[i]]))

hefs_ecrps_sset <- array(NA,c(length(skill_mods),length(lds),hefs_len))
shefs_ecrps_sset <- array(NA,c(length(skill_mods),length(lds),shefs_len))
shefscc_ecrps_sset <- array(NA,c(length(skill_mods),length(lds),shefs_len))

hefs_ecrps_ss_sset <- array(NA,c(length(skill_mods),length(lds),hefs_len))
shefs_ecrps_ss_sset <- array(NA,c(length(skill_mods),length(lds),shefs_len))
shefscc_ecrps_ss_sset <- array(NA,c(length(skill_mods),length(lds),shefs_len))

for(i in 1:length(skill_mods)){
  skill_mod <- skill_mods[i]
  hefs_forward <- readRDS(paste('out/',loc,'/hefs_forward_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'.rds',sep=''))
  syn_hefs_forward <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_swm.rds',sep=''))
  syn_hefs_forward_cc <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_swm-cc.rds',sep=''))

  obs_forward_all_leads_swm <- readRDS(paste('out/',loc,'/obs-fwd_use-bcf=',use_bcf,'_k=',k_swm,'_same-swm=',same_swm,'_cutoff=',yr_cutoff,'_swm.rds',sep=''))
  obs_forward_all_leads_swmcc <- readRDS(paste('out/',loc,'/obs-fwd_use-bcf=',use_bcf,'_k=',k_swm,'_same-swm=',same_swm,'_cutoff=',yr_cutoff,'_swm-cc.rds',sep=''))

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

  keysite = which(site_names==site)
  
  #ensembles at a given lead
  #hefs processing
  hefs_fwd <- hefs_forward[1,keysite,,,]
  obs_hefs <- obs[ixx_obs%in%ixx_hefs,keysite]
  hind_evt_idx <- order(obs_hefs,decreasing=T)[evt_order]
  hind_evt_date <- ixx_hefs[hind_evt_idx]
  
  hefs_out <- t(hefs_fwd[,hind_evt_idx-lead,])
  hefs_ens_sset[i,,] <- hefs_out
  
  #shefs nocc processing
  if(nn_order!=0){
    diff_vec <- abs(obs_hefs[hind_evt_idx] - obs_swm)
    swm_evt_idx <- order(diff_vec,decreasing=F)[nn_order]}
  if(nn_order==0|nn_order==99){
    swm_evt_idx <- order(obs_swm,decreasing=T)[evt_order]}
  
  shefs_out <- t(shefs_fwd_swm[,swm_evt_idx-lead,])
  shefs_ens_sset[i,,] <- shefs_out
  
  #shefs cc processing
  if(nn_order==99){
    swm_evt_idx <- order(obs_swm_cc,decreasing=T)[evt_order]}
  
  shefscc_out <- t(shefs_fwd_swmcc[,swm_evt_idx-lead,])
  shefscc_ens_sset[i,,] <- shefscc_out
  
  #CRPS for a range of lead times
  for(l in 1:length(lds)){
    hefs_vec <- hefs_ecrps_vec[,lds[l]]
    shefs_vec <- as.vector(shefs_ecrps_vec[,,lds[l]])
    shefscc_vec <- as.vector(shefscc_ecrps_vec[,,lds[l]])
    
    hefs_ss_vec <- hefs_ecrps_ss_vec[,lds[l]]
    shefs_ss_vec <- as.vector(shefs_ecrps_ss_vec[,,lds[l]])
    shefscc_ss_vec <- as.vector(shefscc_ecrps_ss_vec[,,lds[l]])
    
    hefs_ecrps_sset[i,l,] <- hefs_vec
    shefs_ecrps_sset[i,l,] <- shefs_vec
    shefscc_ecrps_sset[i,l,] <- shefscc_vec
    
    hefs_ecrps_ss_sset[i,l,] <- hefs_ss_vec
    shefs_ecrps_ss_sset[i,l,] <- shefs_ss_vec
    shefscc_ecrps_ss_sset[i,l,] <- shefscc_ss_vec}
}

saveRDS(hefs_ens_sset,paste('out/',loc,'/hefs_ens-sset_mod=',str_flatten(skill_mods,','),'_dcy=',skill_dcy,'_tail=',skill_tail,'_smp=',samp_no,'.rds',sep=''))
saveRDS(shefs_ens_sset,paste('out/',loc,'/shefs_ens-sset_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',str_flatten(skill_mods,','),'_ld=',lead,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_smp=',samp_no,'_swm.rds',sep=''))
saveRDS(shefscc_ens_sset,paste('out/',loc,'/hefs_ens-sset_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',str_flatten(skill_mods,','),'_ld=',lead,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_smp=',samp_no,'_swm-cc.rds',sep=''))

saveRDS(hefs_ecrps_sset,paste('./out/',loc,'/',site,'_mod=',str_flatten(skill_mods),'_lds=',str_flatten(lds),'_dcy=',skill_dcy,'_tail=',skill_tail,'_hefs-ecrps-sset.rds',sep=''))
saveRDS(shefs_ecrps_sset,paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',str_flatten(skill_mods),'_lds=',str_flatten(lds),'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefs-ecrps-sset.rds',sep=''))
saveRDS(shefscc_ecrps_sset,paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',str_flatten(skill_mods),'_lds=',str_flatten(lds),'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefscc-ecrps-sset.rds',sep=''))

saveRDS(hefs_ecrps_ss_sset,paste('./out/',loc,'/',site,'_mod=',str_flatten(skill_mods),'_lds=',str_flatten(lds),'_dcy=',skill_dcy,'_tail=',skill_tail,'_hefs-ecrps-ss-sset.rds',sep=''))
saveRDS(shefs_ecrps_ss_sset,paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',str_flatten(skill_mods),'_lds=',str_flatten(lds),'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefs-ecrps-ss-sset.rds',sep=''))
saveRDS(shefscc_ecrps_ss_sset,paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',str_flatten(skill_mods),'_lds=',str_flatten(lds),'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefscc-ecrps-ss-sset.rds',sep=''))

rm(list=ls());gc()


######################################################END######################################################
