# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 10:38:00 2025

@author: zpb4
"""
import sys
import os
sys.path.insert(0, os.path.abspath('./src'))
import numpy as np
import pandas as pd
import matplotlib as matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import model
import syn_util
from time import localtime, strftime
from datetime import datetime
import matplotlib.dates as mdates
import ensemble_verification_functions as verify
import pickle
import datetime as dt

kcfs_to_tafd = 2.29568411*10**-5 * 86400
max_lds = 15

skill_dcy = 0.1
skill_tail = 0.5

sd_hist = '1990-10-01' 
ed_hist = '2019-08-15'

sd_swm = '2001-01-01'
ed_swm = '3008-01-08'

loc = 'YRS'
site = 'ORDC1'
evt_no = 1
n_evts = 100
sep = 15

tocs_reset = 'none'   # 'tocs' to reset to baseline tocs at beginning of WY, 'firo' to reset to firo pool, 'none' to use continuous storage profile
use_firo_bottom = True
use_firo_top = True
seed = 1
skillmods = ([-0.4,0,0.4,0.8])
yr_cutoff = 500
swm_tst_samp = 2
swmcc_tst_samp = 3

syn_pct = 0.99
syn_setup = 'cal'
syn_optstrat = 'ecrps-dts'
syn_objpwr = 0
syn_path = './'  # path to R synthetic forecast repo for 'r-gen' setting below

use_bcf = False
k_swm = 106
same_swm = False

res_params = pd.read_csv('./data/reservoir-storage-safe-rel.csv')
res_idx = np.where(res_params['Reservoir']==site)
K = res_params['Capacity (TAF)'][res_idx[0]].values[0]
Rmax = res_params['Safe Release (CFS)'][res_idx[0]].values[0] / 1000 * kcfs_to_tafd
ramping_rate_up = res_params['Ramping-up (CFS)'][res_idx[0]].values[0] / 1000 * kcfs_to_tafd 
ramping_rate_down = res_params['Ramping-down (CFS)'][res_idx[0]].values[0] / 1000 * kcfs_to_tafd 
ramping_rate = (ramping_rate_up,ramping_rate_down) 

store_vec = model.store_vec
elev_vec = model.elev_vec
elev_ctrl_vec = model.elev_ctrl_vec
flow_ctrl_vec = model.flow_ctrl_vec
K_ORO = model.K_ORO
firo_space_ORO = model.firo_space_ORO
firo_top_ORO = model.firo_top_ORO

#swm test data
swm_ind = 'swm'
samp_no = swm_tst_samp
Q_full,dowy_full,tocs_full,df_idx_full = syn_util.extract_swm(sd_swm,ed_swm,syn_path,loc=loc,site=site,use_bcf=use_bcf,k_swm=k_swm,same_swm=same_swm,swm_ind=swm_ind,samp_no=samp_no,yr_cutoff=yr_cutoff)
Q_swm_tst = Q_full[:-max_lds]
dowy = dowy_full[:-max_lds]

#swm tst array with multiple skills
swm_tst_arr = np.load('out/%s/%s/swm-test-array_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']
swm_tst_arr_skilltrn = np.load('out/%s/%s/swm-test-array-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']

evt_idx_swm_tst = syn_util.declust_evts_extract(Q_swm_tst,n_evts,sep)

#swm test data
swm_ind = 'swm-cc'
samp_no = swmcc_tst_samp
Q_full,dowy_full,tocs_full,df_idx_full = syn_util.extract_swm(sd_swm,ed_swm,syn_path,loc=loc,site=site,use_bcf=use_bcf,k_swm=k_swm,same_swm=same_swm,swm_ind=swm_ind,samp_no=samp_no,yr_cutoff=yr_cutoff)
Q_swmcc_tst = Q_full[:-max_lds]

evt_idx_swmcc_tst = syn_util.declust_evts_extract(Q_swmcc_tst,n_evts,sep)

#swm tst array with multiple skills
swmcc_tst_arr = np.load('out/%s/%s/swm-test-array_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']
swmcc_tst_arr_skilltrn = np.load('out/%s/%s/swm-test-array-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']

#swm test data
swm_ind = 'swm'
samp_no = swm_tst_samp
evt_idx = syn_util.declust_evts_extract(Q_swm_tst,n_evts,sep)
#cold_idx = np.where((dowy>15) & (dowy<151))
cold_idx = np.where((dowy>349) | (dowy<258))
apr1_idx = np.where(dowy==182)

swm_tst_metrics = np.zeros((14,n_evts,len(skillmods))) 
swm_tst_apr1 = np.zeros((np.shape(apr1_idx)[1],len(skillmods)))

#0. firo prerelease: maximum drawdown prior to flood event
#1. firo overage: maximum overage of the firo pool +/- 15 days of event
#2. firo overage percent: maximum overage of the firo pool +/- 15 days of event expressed as a percentage of flood space [firo_ovg / (K - firo_top)]
#3. maximum release: event based maximum release value +/- 15 days of event
#4. mean release: event based mean release value +/- 15 days of event
#5. max ramping: event based maximum ramping value +/- 15 days of event
#6. mean ramping: event based mean ramping value +/- 15 days of event
#7. post release max drawdown: maximum drawdown occurring within 15 days after the event
#8. post release cumulative loss: cumulative volume below FIRO pool occurring within 15 days after the event
#9. subfiro event storage: event based storage below FIRO pool (S-firo_ovg) +/- 15 days of event
#10. subfiro event storage percent: event based storage below FIRO pool (S-firo_ovg) +/- 15 days of event expressed as a percent of FIRO space [(S-firo_ovg)/(firo_top)]
#11. number of spills: number of spills across entire simulation period
#12. number of spills: volume of spills across entire simulation period
#13. subfiro cold season storage: storage below firo pool in cold season (min firo pool)
#14. Mar 1 storage: 

for k in range(n_evts):
    for j in range(len(skillmods)):         
        swm_tst_metrics[0,k,j] = np.max(swm_tst_arr[j,6,(evt_idx[k]-15):evt_idx[k]+2])
        swm_tst_metrics[1,k,j] = np.max(swm_tst_arr[j,5,(evt_idx[k]-15):(evt_idx[k]+15)])
        swm_tst_metrics[2,k,j] = np.max(swm_tst_arr[j,8,(evt_idx[k]-15):(evt_idx[k]+15)])
        swm_tst_metrics[3,k,j] = np.max(swm_tst_arr[j,1,(evt_idx[k]-15):(evt_idx[k]+15)])
        swm_tst_metrics[4,k,j] = np.mean(swm_tst_arr[j,1,(evt_idx[k]-15):(evt_idx[k]+15)])
        swm_tst_metrics[5,k,j] = np.max(swm_tst_arr[j,9,(evt_idx[k]-15):(evt_idx[k]+15)])
        swm_tst_metrics[6,k,j] = np.mean(swm_tst_arr[j,9,(evt_idx[k]-15):(evt_idx[k]+15)])
        swm_tst_metrics[7,k,j] = np.max(swm_tst_arr[j,6,(evt_idx[k]+1):(evt_idx[k]+15)])
        swm_tst_metrics[8,k,j] = np.sum(swm_tst_arr[j,6,(evt_idx[k]+1):(evt_idx[k]+15)])
        swm_tst_metrics[9,k,j] = np.mean((swm_tst_arr[j,0,(evt_idx[k]-15):(evt_idx[k]+15)] - swm_tst_arr[j,5,(evt_idx[k]-15):(evt_idx[k]+15)]))
        swm_tst_metrics[10,k,j] = np.mean((swm_tst_arr[j,0,(evt_idx[k]-15):(evt_idx[k]+15)] - swm_tst_arr[j,5,(evt_idx[k]-15):(evt_idx[k]+15)])/swm_tst_arr[j,2,(evt_idx[k]-15):(evt_idx[k]+15)]) #mean of sub-FIRO storage normalized by top of FIRO pool
        swm_tst_metrics[11,k,j] = np.count_nonzero(swm_tst_arr[j,3,:])
        swm_tst_metrics[12,k,j] = np.sum(swm_tst_arr[j,3,:])
        swm_tst_metrics[13,k,j] = np.mean(swm_tst_arr[j,0,cold_idx])
        swm_tst_apr1[:,j] = swm_tst_arr[j,0,apr1_idx]

np.savez('out/%s/%s/swm-test-metrics_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff),arr=swm_tst_metrics)
np.savez('out/%s/%s/swm-test-apr1_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff),arr=swm_tst_apr1)

swm_tst_metrics_skilltrn = np.zeros((14,n_evts,len(skillmods))) 
swm_tst_apr1_skilltrn = np.zeros((np.shape(apr1_idx)[1],len(skillmods)))

for k in range(n_evts):
    for j in range(len(skillmods)):         
        swm_tst_metrics_skilltrn[0,k,j] = np.max(swm_tst_arr_skilltrn[j,6,(evt_idx[k]-15):evt_idx[k]+2])
        swm_tst_metrics_skilltrn[1,k,j] = np.max(swm_tst_arr_skilltrn[j,5,(evt_idx[k]-15):(evt_idx[k]+15)])
        swm_tst_metrics_skilltrn[2,k,j] = np.max(swm_tst_arr_skilltrn[j,8,(evt_idx[k]-15):(evt_idx[k]+15)])
        swm_tst_metrics_skilltrn[3,k,j] = np.max(swm_tst_arr_skilltrn[j,1,(evt_idx[k]-15):(evt_idx[k]+15)])
        swm_tst_metrics_skilltrn[4,k,j] = np.mean(swm_tst_arr_skilltrn[j,1,(evt_idx[k]-15):(evt_idx[k]+15)])
        swm_tst_metrics_skilltrn[5,k,j] = np.max(swm_tst_arr_skilltrn[j,9,(evt_idx[k]-15):(evt_idx[k]+15)])
        swm_tst_metrics_skilltrn[6,k,j] = np.mean(swm_tst_arr_skilltrn[j,9,(evt_idx[k]-15):(evt_idx[k]+15)])
        swm_tst_metrics_skilltrn[7,k,j] = np.max(swm_tst_arr_skilltrn[j,6,(evt_idx[k]+1):(evt_idx[k]+15)])
        swm_tst_metrics_skilltrn[8,k,j] = np.sum(swm_tst_arr_skilltrn[j,6,(evt_idx[k]+1):(evt_idx[k]+15)])
        swm_tst_metrics_skilltrn[9,k,j] = np.mean((swm_tst_arr_skilltrn[j,0,(evt_idx[k]-15):(evt_idx[k]+15)] - swm_tst_arr_skilltrn[j,5,(evt_idx[k]-15):(evt_idx[k]+15)]))
        swm_tst_metrics_skilltrn[10,k,j] = np.mean((swm_tst_arr_skilltrn[j,0,(evt_idx[k]-15):(evt_idx[k]+15)] - swm_tst_arr_skilltrn[j,5,(evt_idx[k]-15):(evt_idx[k]+15)])/swm_tst_arr_skilltrn[j,2,(evt_idx[k]-15):(evt_idx[k]+15)]) #mean of sub-FIRO storage normalized by top of FIRO pool
        swm_tst_metrics_skilltrn[11,k,j] = np.count_nonzero(swm_tst_arr_skilltrn[j,3,:])
        swm_tst_metrics_skilltrn[12,k,j] = np.sum(swm_tst_arr_skilltrn[j,3,:])
        swm_tst_metrics_skilltrn[13,k,j] = np.mean(swm_tst_arr_skilltrn[j,0,cold_idx])
        swm_tst_apr1_skilltrn[:,j] = swm_tst_arr_skilltrn[j,0,apr1_idx]
        
np.savez('out/%s/%s/swm-test-metrics-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff),arr=swm_tst_metrics_skilltrn)
np.savez('out/%s/%s/swm-test-apr1-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff),arr=swm_tst_apr1_skilltrn)

#swm-cc test data
swm_ind = 'swm-cc'
samp_no = swmcc_tst_samp
evt_idx = syn_util.declust_evts_extract(Q_swmcc_tst,n_evts,sep)
swmcc_tst_metrics = np.zeros((14,n_evts,len(skillmods))) 
swmcc_tst_apr1 = np.zeros((np.shape(apr1_idx)[1],len(skillmods)))

for k in range(n_evts):
    for j in range(len(skillmods)):         
        swmcc_tst_metrics[0,k,j] = np.max(swmcc_tst_arr[j,6,(evt_idx[k]-15):evt_idx[k]+2])
        swmcc_tst_metrics[1,k,j] = np.max(swmcc_tst_arr[j,5,(evt_idx[k]-15):(evt_idx[k]+15)])
        swmcc_tst_metrics[2,k,j] = np.max(swmcc_tst_arr[j,8,(evt_idx[k]-15):(evt_idx[k]+15)])
        swmcc_tst_metrics[3,k,j] = np.max(swmcc_tst_arr[j,1,(evt_idx[k]-15):(evt_idx[k]+15)])
        swmcc_tst_metrics[4,k,j] = np.mean(swmcc_tst_arr[j,1,(evt_idx[k]-15):(evt_idx[k]+15)])
        swmcc_tst_metrics[5,k,j] = np.max(swmcc_tst_arr[j,9,(evt_idx[k]-15):(evt_idx[k]+15)])
        swmcc_tst_metrics[6,k,j] = np.mean(swmcc_tst_arr[j,9,(evt_idx[k]-15):(evt_idx[k]+15)])
        swmcc_tst_metrics[7,k,j] = np.max(swmcc_tst_arr[j,6,(evt_idx[k]+1):(evt_idx[k]+15)])
        swmcc_tst_metrics[8,k,j] = np.sum(swmcc_tst_arr[j,6,(evt_idx[k]+1):(evt_idx[k]+15)])
        swmcc_tst_metrics[9,k,j] = np.mean((swmcc_tst_arr[j,0,(evt_idx[k]-15):(evt_idx[k]+15)] - swmcc_tst_arr[j,5,(evt_idx[k]-15):(evt_idx[k]+15)]))
        swmcc_tst_metrics[10,k,j] = np.mean((swmcc_tst_arr[j,0,(evt_idx[k]-15):(evt_idx[k]+15)] - swmcc_tst_arr[j,5,(evt_idx[k]-15):(evt_idx[k]+15)])/swmcc_tst_arr[j,2,(evt_idx[k]-15):(evt_idx[k]+15)]) #mean of sub-FIRO storage normalized by top of FIRO pool
        swmcc_tst_metrics[11,k,j] = np.count_nonzero(swmcc_tst_arr[j,3,:])
        swmcc_tst_metrics[12,k,j] = np.sum(swmcc_tst_arr[j,3,:])
        swmcc_tst_metrics[13,k,j] = np.mean(swmcc_tst_arr[j,0,cold_idx])
        swmcc_tst_apr1[:,j] = swmcc_tst_arr[j,0,apr1_idx]

np.savez('out/%s/%s/swm-test-metrics_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff),arr=swmcc_tst_metrics)
np.savez('out/%s/%s/swm-test-apr1_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff),arr=swmcc_tst_apr1)

swmcc_tst_metrics_skilltrn = np.zeros((14,n_evts,len(skillmods))) 
swmcc_tst_apr1_skilltrn = np.zeros((np.shape(apr1_idx)[1],len(skillmods)))

for k in range(n_evts):
    for j in range(len(skillmods)):         
        swmcc_tst_metrics_skilltrn[0,k,j] = np.max(swmcc_tst_arr_skilltrn[j,6,(evt_idx[k]-15):evt_idx[k]+2])
        swmcc_tst_metrics_skilltrn[1,k,j] = np.max(swmcc_tst_arr_skilltrn[j,5,(evt_idx[k]-15):(evt_idx[k]+15)])
        swmcc_tst_metrics_skilltrn[2,k,j] = np.max(swmcc_tst_arr_skilltrn[j,8,(evt_idx[k]-15):(evt_idx[k]+15)])
        swmcc_tst_metrics_skilltrn[3,k,j] = np.max(swmcc_tst_arr_skilltrn[j,1,(evt_idx[k]-15):(evt_idx[k]+15)])
        swmcc_tst_metrics_skilltrn[4,k,j] = np.mean(swmcc_tst_arr_skilltrn[j,1,(evt_idx[k]-15):(evt_idx[k]+15)])
        swmcc_tst_metrics_skilltrn[5,k,j] = np.max(swmcc_tst_arr_skilltrn[j,9,(evt_idx[k]-15):(evt_idx[k]+15)])
        swmcc_tst_metrics_skilltrn[6,k,j] = np.mean(swmcc_tst_arr_skilltrn[j,9,(evt_idx[k]-15):(evt_idx[k]+15)])
        swmcc_tst_metrics_skilltrn[7,k,j] = np.max(swmcc_tst_arr_skilltrn[j,6,(evt_idx[k]+1):(evt_idx[k]+15)])
        swmcc_tst_metrics_skilltrn[8,k,j] = np.sum(swmcc_tst_arr_skilltrn[j,6,(evt_idx[k]+1):(evt_idx[k]+15)])
        swmcc_tst_metrics_skilltrn[9,k,j] = np.mean((swmcc_tst_arr_skilltrn[j,0,(evt_idx[k]-15):(evt_idx[k]+15)] - swmcc_tst_arr_skilltrn[j,5,(evt_idx[k]-15):(evt_idx[k]+15)]))
        swmcc_tst_metrics_skilltrn[10,k,j] = np.mean((swmcc_tst_arr_skilltrn[j,0,(evt_idx[k]-15):(evt_idx[k]+15)] - swmcc_tst_arr_skilltrn[j,5,(evt_idx[k]-15):(evt_idx[k]+15)])/swmcc_tst_arr_skilltrn[j,2,(evt_idx[k]-15):(evt_idx[k]+15)]) #mean of sub-FIRO storage normalized by top of FIRO pool
        swmcc_tst_metrics_skilltrn[11,k,j] = np.count_nonzero(swmcc_tst_arr_skilltrn[j,3,:])
        swmcc_tst_metrics_skilltrn[12,k,j] = np.sum(swmcc_tst_arr_skilltrn[j,3,:])
        swmcc_tst_metrics_skilltrn[13,k,j] = np.mean(swmcc_tst_arr_skilltrn[j,0,cold_idx])
        swmcc_tst_apr1_skilltrn[:,j] = swmcc_tst_arr_skilltrn[j,0,apr1_idx]
        
np.savez('out/%s/%s/swm-test-metrics-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff),arr=swmcc_tst_metrics_skilltrn)
np.savez('out/%s/%s/swm-test-apr1-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff),arr=swmcc_tst_apr1_skilltrn)

##########################################################END#####################################################