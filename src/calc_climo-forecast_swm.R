
args = commandArgs(trailingOnly=TRUE)
print(paste('task #',args[1]))
idx = as.numeric(args[1])

print(paste('calc start',Sys.time()))

#----------------------------------------
loc = 'YRS'
keysite_name = 'ORDC1'
cal_val_setup = 'cal'
pcnt_opt = 0.99
obj_pwr = 0
opt_strat = 'ecrps-dts'
same_swm = F

use_bcf=F
k_swm=106
yr_cutoff=500

skill_mod = 0
skill_dcy = 0.1
skill_tail = 0.5

load(paste('./out/',loc,'/data_prep_rdata.RData',sep=''))
source('./src/forecast_verification_functions.R')

syn_hefs_forward_ref <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',keysite_name,'_',cal_val_setup,'_same-swm=',same_swm,'_swm-ref.rds',sep=''))
print(dim(syn_hefs_forward_ref))
syn_hefs_forward_swm <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',keysite_name,'_',cal_val_setup,'_same-swm=',same_swm,'_cutoff=',yr_cutoff,'_swm.rds',sep=''))
print(dim(syn_hefs_forward_swm))
syn_hefs_forward_swmcc <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',keysite_name,'_',cal_val_setup,'_same-swm=',same_swm,'_cutoff=',yr_cutoff,'_swm-cc.rds',sep=''))
print(dim(syn_hefs_forward_swmcc))

obs_forward_all_leads_swm <- readRDS(paste('out/',loc,'/obs-fwd_use-bcf=',use_bcf,'_k=',k_swm,'_same-swm=',same_swm,'_cutoff=',yr_cutoff,'_swm.rds',sep=''))
obs_forward_all_leads_swmcc <- readRDS(paste('out/',loc,'/obs-fwd_use-bcf=',use_bcf,'_k=',k_swm,'_same-swm=',same_swm,'_cutoff=',yr_cutoff,'_swm-cc.rds',sep=''))

ixx_gen_swm <- readRDS(paste('out/',loc,'/ixx_gen_swm.rds',sep=''))
ixx_gen <- ixx_obs_forward
n_samp_ref <- dim(syn_hefs_forward_ref)[1]
n_samp_swm <- dim(syn_hefs_forward_swm)[1]

#add a single entry dimension to match synthetic forecasts
hefs_fwd<-array(NA,c(1,dim(hefs_forward)))
hefs_fwd[1,,,,]<-hefs_forward
cur_site <- which(site_names==sites)

if(idx==1){
  climo_farray_hefs <- climo_forecast_swm(ixx_hefs,syn_hefs_forward_ref[1,1,,ixx_obs_forward%in%ixx_hefs,],obs_forward_all_leads[cur_site,ixx_obs_forward%in%ixx_hefs,])
  climo_farray_shefs <- climo_forecast_swm(ixx_obs_forward,syn_hefs_forward_ref[1,1,,,],obs_forward_all_leads[cur_site,,])
  saveRDS(climo_farray_hefs, paste('out/',loc,'/hefs_climo.rds',sep=''))
  saveRDS(climo_farray_shefs, paste('out/',loc,'/synhefs-ref_climo.rds',sep=''))}

if(idx==2){
  ##climo_farray_swm <- array(NA,dim(syn_hefs_forward_swm))
  ##inp <- climo_forecast_swm(ixx_gen_swm,syn_hefs_forward_swm[1,1,,,],obs_forward_all_leads_swm[1,1,,2:(leads+1)])
  ##for(s in 1:n_samp_swm){
    ##inp <- climo_forecast(ixx_gen_swm,syn_hefs_forward_swm[s,1,,,],obs_forward_all_leads_swm[s,1,,2:(leads+1)])
    ##climo_farray_swm[s,1,,,] <- inp}
  climo_farray_swm <- climo_forecast_swm(ixx_gen_swm,syn_hefs_forward_swm[1,1,,,],obs_forward_all_leads_swm[1,1,,2:(leads+1)])
  saveRDS(climo_farray_swm, paste('out/',loc,'/synhefs_climo_swm.rds',sep=''))
  print(paste('swm calc end',Sys.time()))}

if(idx==3){
  ##climo_farray_swmcc <- array(NA,dim(syn_hefs_forward_swmcc))
  ##inp <- climo_forecast_swm(ixx_gen_swm,syn_hefs_forward_swmcc[1,1,,,],obs_forward_all_leads_swmcc[1,1,,2:(leads+1)])
  ##for(s in 1:n_samp_swm){
    ##inp <- climo_forecast(ixx_gen_swm,syn_hefs_forward_swmcc[s,1,,,],obs_forward_all_leads_swmcc[s,1,,2:(leads+1)])
    ##climo_farray_swmcc[s,1,,,] <- inp}
  climo_farray_swmcc <- climo_forecast_swm(ixx_gen_swm,syn_hefs_forward_swmcc[1,1,,,],obs_forward_all_leads_swmcc[1,1,,2:(leads+1)])
  saveRDS(climo_farray_swmcc, paste('out/',loc,'/synhefs_climo_swmcc.rds',sep=''))
  print(paste('swmcc calc end',Sys.time()))}

print(paste('calc end',Sys.time()))

rm(list=ls());gc()

########################################END############################################################