
args = commandArgs(trailingOnly=TRUE)
print(paste('task #',args[1]))
idx = as.numeric(args[1])

print(paste('calc start',Sys.time()))
#///////////////////////////////////////////////////////////////////////////////////////////////////////////
##idx=4
##setwd('z:/syn-forecast_stochastic/')
#Data setup
syn_vers = 2
loc = 'YRS'
site = 'ORDC1'
cal_val_setup = 'cal'
pcnt_opt = 0.99
obj_pwr = 0
opt_strat = 'ecrps-dts'
same_swm = F
samp_no = 1

skill_mods = c(-0.4,0,0.4,0.8)
skill_mod = skill_mods[idx]
skill_dcy = 0.1
skill_tail = 0.5

#swm data
use_bcf=F
k_swm=106
set.seed(1)
yr_cutoff=500

source('./src/forecast_verification_functions.R')

#plot setup
disp_pcnt <- 0.99

load(paste('./out/',loc,'/data_prep_rdata.RData',sep=''))
rm(hefs_forward,hefs_forward_cumul,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid,hefs_forward_frac,hefs_inp);gc()

hefs_inp <- readRDS(paste('out/',loc,'/hefs_forward_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'.rds',sep=''))
syn_hefs_forward_ref <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_same-swm=',same_swm,'_swm-ref.rds',sep=''))
shefs_ref <- syn_hefs_forward_ref[,1,,,]
rm(syn_hefs_forward_ref)
syn_hefs_forward_swm <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_swm.rds',sep=''))
shefs_swm <- syn_hefs_forward_swm[samp_no,1,,,]
rm(syn_hefs_forward_swm)
syn_hefs_forward_swmcc <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_swm-cc.rds',sep=''))
shefs_swmcc <- syn_hefs_forward_swmcc[samp_no,1,,,]
rm(syn_hefs_forward_swmcc);gc()

obs_forward_all_leads_swm <- readRDS(paste('out/',loc,'/obs-fwd_use-bcf=',use_bcf,'_k=',k_swm,'_same-swm=',same_swm,'_cutoff=',yr_cutoff,'_swm.rds',sep=''))
obs_forward_all_leads_swmcc <- readRDS(paste('out/',loc,'/obs-fwd_use-bcf=',use_bcf,'_k=',k_swm,'_same-swm=',same_swm,'_cutoff=',yr_cutoff,'_swm-cc.rds',sep=''))
obs_fwd_swm <- obs_forward_all_leads_swm[samp_no,1,,]
obs_fwd_swmcc <- obs_forward_all_leads_swmcc[samp_no,1,,]
rm(obs_forward_all_leads_swm,obs_forward_all_leads_swmcc);gc()

ixx_swm <- readRDS(paste('out/',loc,'/ixx_gen_swm.rds',sep=''))
ixx_swm_wy <- wy_fun(ixx_swm)

climo_hefs <- readRDS(paste('out/',loc,'/hefs_climo.rds',sep=''))
climo_shefs <- readRDS(paste('out/',loc,'/synhefs-ref_climo.rds',sep=''))
climo_swm <- readRDS(paste('out/',loc,'/synhefs_climo_swm.rds',sep=''))
climo_swmcc <- readRDS(paste('out/',loc,'/synhefs_climo_swmcc.rds',sep=''))

#///////////////////////////////////////////////////////////////////////////////////////////////////////////
#setup obs to sample from
keysite = which(site_names==site)

hefs_fwd <- hefs_inp[1,keysite,,,]
obs_fwd_hefs <- obs_forward_all_leads_hind[keysite,,]
obs_hefs <- obs[ixx_obs%in%ixx_hefs,keysite]
obs_swm <- obs_fwd_swm[,1]
obs_swm_cc <- obs_fwd_swmcc[,1]

rm(obs_forward_all_leads,obs_forward_all_leads_hind);gc()


###############eCRPS + Rank Histogram####################################
lds<-1:leads  
nevts_hefs <- round((1-disp_pcnt) * length(obs_hefs))

hefs_yrs <- unique(ixx_hefs$year)
swm_yrs <- unique(ixx_swm_wy$year)
samp_vec <- 1:(length(swm_yrs)-length(hefs_yrs))

samps = dim(shefs_ref)[1]
#samp_idx = vector('list',samps)
swm_samps = vector('list',samps)
swmcc_samps = vector('list',samps)
shefs_samps = vector('list',samps)
shefscc_samps = vector('list',samps)
shefs_climo_samps = vector('list',samps)
shefscc_climo_samps = vector('list',samps)

for(i in 1:samps){
  start_yr <- sample(samp_vec,1)
  yrs <- swm_yrs[start_yr:(start_yr+length(hefs_yrs)-1)]
  samp_idx<-(1:length(ixx_swm_wy))[ixx_swm_wy$year%in%yrs]
  swm_samps[[i]]<-obs_swm[samp_idx]
  swmcc_samps[[i]]<-obs_swm_cc[samp_idx]
  shefs_samps[[i]]<-shefs_swm[,samp_idx,]
  shefscc_samps[[i]]<-shefs_swmcc[,samp_idx,]
  shefs_climo_samps[[i]]<-climo_swm[,samp_idx,]
  shefscc_climo_samps[[i]]<-climo_swmcc[,samp_idx,]
}

###eCRPS###
obs_date_loc <- order(obs_hefs,decreasing=TRUE)[1:nevts_hefs]  #index for maximum observation
obs_events <- obs_hefs[obs_date_loc]

hefs_ecrps_vec<-array(NA,c(nevts_hefs,length(lds)))
shefsref_ecrps_vec<-array(NA,c(samps,nevts_hefs,length(lds)))
shefs_ecrps_vec<-array(NA,c(samps,nevts_hefs,length(lds)))
shefscc_ecrps_vec<-array(NA,c(samps,nevts_hefs,length(lds)))

hefs_ecrps_ss_vec<-array(NA,c(nevts_hefs,length(lds)))
shefsref_ecrps_ss_vec<-array(NA,c(samps,nevts_hefs,length(lds)))
shefs_ecrps_ss_vec<-array(NA,c(samps,nevts_hefs,length(lds)))
shefscc_ecrps_ss_vec<-array(NA,c(samps,nevts_hefs,length(lds)))

hefs_rank_vec <- array(NA,c(nevts_hefs,length(lds)))
shefsref_rank_vec <- array(NA,c(samps,nevts_hefs,length(lds)))
shefs_rank_vec <- array(NA,c(samps,nevts_hefs,length(lds)))
shefscc_rank_vec <- array(NA,c(samps,nevts_hefs,length(lds)))

for(ld in 1:length(lds)){
  for(i in 1:nevts_hefs){
    hefs_idx <- obs_date_loc[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
    HEFS <- hefs_fwd[,hefs_idx,lds[ld]]
    CLIMO <- climo_hefs[,hefs_idx,lds[ld]]
    hefs_ecrps_vec[i,ld] <- eCRPS(HEFS,obs_events[i])
    climo_ecrps <- eCRPS(CLIMO,obs_events[i])
    hefs_ecrps_ss_vec[i,ld] <- 1 - (hefs_ecrps_vec[i,ld]/climo_ecrps)
    hefs_rank_vec[i,ld] <- ens_rank(HEFS,obs_events[i])
  }
}

saveRDS(hefs_ecrps_vec,paste('./out/',loc,'/',site,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_hefs-ecrps-vec.rds',sep=''))
saveRDS(hefs_ecrps_ss_vec,paste('./out/',loc,'/',site,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_hefs-ecrps-ss-vec.rds',sep=''))
saveRDS(hefs_rank_vec,paste('./out/',loc,'/',site,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_hefs-rank-vec.rds',sep=''))

#if(idx==1){
  for(s in 1:samps){
    swm_sset <- swm_samps[[s]]
    sset_date_loc <- order(swm_sset,decreasing=TRUE)[1:nevts_hefs]  #index for maximum observation
    sset_events <- swm_sset[sset_date_loc]
    for(ld in 1:length(lds)){
      for(i in 1:nevts_hefs){
        syn_idx <- sset_date_loc[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
        shefs_sset <- shefs_samps[[s]]
        shefs_climo_sset <- shefs_climo_samps[[s]]
        SYN_HEFS <- shefs_sset[,syn_idx,lds[ld]]
        CLIMO <- shefs_climo_sset[,syn_idx,lds[ld]]
        shefs_ecrps_vec[s,i,ld] <- eCRPS(SYN_HEFS,sset_events[i])
        climo_ecrps <- eCRPS(CLIMO,sset_events[i])
        shefs_ecrps_ss_vec[s,i,ld] <- 1 - (shefs_ecrps_vec[s,i,ld]/climo_ecrps)}
    }
    print(s)
  }
  
saveRDS(shefs_ecrps_ss_vec,paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefs-ecrps-ss-vec.rds',sep=''))
saveRDS(shefs_ecrps_vec,paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefs-ecrps-vec.rds',sep=''))

print(paste('swm ecrps end',Sys.time()))
#print(paste('swm ecrps end',Sys.time()))}


#if(idx==2){
  for(s in 1:samps){
    swm_sset <- swm_samps[[s]]
    sset_date_loc <- order(swm_sset,decreasing=TRUE)[1:nevts_hefs]  #index for maximum observation
    sset_events <- swm_sset[sset_date_loc]
    for(ld in 1:length(lds)){
      for(i in 1:nevts_hefs){
        syn_idx <- sset_date_loc[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
        shefs_sset <- shefs_samps[[s]]
        SYN_HEFS <- shefs_sset[,syn_idx,lds[ld]]
        shefs_rank_vec[s,i,ld] <- ens_rank(SYN_HEFS,sset_events[i])}
    }
    print(s)
  }
  saveRDS(shefs_rank_vec,paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefs-rank-vec.rds',sep=''))

print(paste('swm rank end',Sys.time()))
#print(paste('swm rank end',Sys.time()))}

#if(idx==3){
  for(s in 1:samps){
    swm_sset <- swmcc_samps[[s]]
    sset_date_loc <- order(swm_sset,decreasing=TRUE)[1:nevts_hefs]  #index for maximum observation
    sset_events <- swm_sset[sset_date_loc]
    for(ld in 1:length(lds)){
      for(i in 1:nevts_hefs){
        syn_idx <- sset_date_loc[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
        shefs_sset <- shefscc_samps[[s]]
        shefs_climo_sset <- shefscc_climo_samps[[s]]
        SYN_HEFS <- shefs_sset[,syn_idx,lds[ld]]
        CLIMO <- shefs_climo_sset[,syn_idx,lds[ld]]
        shefscc_ecrps_vec[s,i,ld] <- eCRPS(SYN_HEFS,sset_events[i])
        climo_ecrps <- eCRPS(CLIMO,sset_events[i])
        shefscc_ecrps_ss_vec[s,i,ld] <- 1 - (shefscc_ecrps_vec[s,i,ld]/climo_ecrps)}
    }
    print(s)
  }
saveRDS(shefscc_ecrps_vec,paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefscc-ecrps-vec.rds',sep=''))
saveRDS(shefscc_ecrps_ss_vec,paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefscc-ecrps-ss-vec.rds',sep=''))

print(paste('swmcc ecrps end',Sys.time()))
#print(paste('swmcc ecrps end',Sys.time()))}



#if(idx==4){
  for(s in 1:samps){
    swm_sset <- swmcc_samps[[s]]
    sset_date_loc <- order(swm_sset,decreasing=TRUE)[1:nevts_hefs]  #index for maximum observation
    sset_events <- swm_sset[sset_date_loc]
    for(ld in 1:length(lds)){
      for(i in 1:nevts_hefs){
        syn_idx <- sset_date_loc[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
        shefs_sset <- shefscc_samps[[s]]
        SYN_HEFS <- shefs_sset[,syn_idx,lds[ld]]
        shefscc_rank_vec[s,i,ld] <- ens_rank(SYN_HEFS,sset_events[i])}
    }
    #return(shefs_rank)
    print(s)
  }
saveRDS(shefscc_rank_vec,paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefscc-rank-vec.rds',sep=''))

print(paste('swmcc rank end',Sys.time()))
#print(paste('swmcc rank end',Sys.time()))}


#the reference synthetic HEFS were for QCing and not needed to run in primary workflow

#if(idx==5){
  ##for(s in 1:samps){
    ##for(ld in 1:length(lds)){
      ##for(i in 1:nevts_hefs){
        ##syn_idx <- obs_date_loc[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
        ##SYN_HEFS <- shefs_ref[s,,syn_idx,lds[ld]]
        ##CLIMO <- climo_shefs[,syn_idx,lds[ld]]
        ##shefsref_ecrps_vec[s,i,ld] <- eCRPS(SYN_HEFS,obs_events[i])
        ##climo_ecrps <- eCRPS(CLIMO,obs_events[i])
        ##shefsref_ecrps_ss_vec[s,i,ld] <- 1 - (shefsref_ecrps_vec[s,i,ld]/climo_ecrps)}
    ##}
    ##print(s)
  ##}
##saveRDS(shefsref_ecrps_vec,paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefsref-ecrps-vec.rds',sep=''))
##saveRDS(shefsref_ecrps_ss_vec,paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefsref-ecrps-ss-vec.rds',sep=''))
##print(paste('shefs-ref ecrps end',Sys.time()))
#}

#if(idx==6){
  ##for(s in 1:samps){
    ##for(ld in 1:length(lds)){
      ##for(i in 1:nevts_hefs){
        ##syn_idx <- obs_date_loc[i]-lds[ld] #need to back up by lds[ld] because forecasts are in 'forward' format
        ##SYN_HEFS <- shefs_ref[s,,syn_idx,lds[ld]]
        ##shefsref_rank_vec[s,i,ld] <- ens_rank(SYN_HEFS,obs_events[i])}
    ##}
    ##print(s)
  ##}
##saveRDS(shefsref_rank_vec,paste('./out/',loc,'/',site,'_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',site,'_',cal_val_setup,'_same-swm=',same_swm,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_cutoff=',yr_cutoff,'_dpcnt=',disp_pcnt,'_shefsref-rank-vec.rds',sep=''))
##print(paste('shefs-ref rank end',Sys.time()))
#}

print(paste('calc end',Sys.time()))

rm(list=ls());gc()


#####################################################END###################################################