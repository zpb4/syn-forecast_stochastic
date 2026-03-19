# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 10:26:24 2025

@author: zpb4
"""

import sys
import os
sys.path.insert(0, os.path.abspath('../Synthetic-Forecast_Verification/src'))
import numpy as np
import pandas as pd
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

rcurve_skillmods = [-0.4,0,0.4,0.8]
hist_idx = np.where(np.array(rcurve_skillmods)==0)[0]
tst_skillmods = [-0.4,0,0.4,0.8]
skill_dcy = 0.1
skill_tail = 0.5
yr_cutoff = 500

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
swm_test_samp = 2
swmcc_test_samp = 3

syn_pct = 0.99
syn_setup = 'cal'
syn_optstrat = 'ecrps-dts'
syn_objpwr = 0
path = './'  # path to R synthetic forecast repo for 'r-gen' setting below

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

#SWM obs data
#swm training data
swm_ind = 'swm'
samp_no = 1
Q_full,dowy_full,tocs_full,df_idx_full= syn_util.extract_swm(sd_swm,ed_swm,path,loc=loc,site=site,use_bcf=use_bcf,k_swm=k_swm,same_swm=same_swm,swm_ind=swm_ind,samp_no=samp_no,yr_cutoff=yr_cutoff)
Q_swm_trn = Q_full[:-max_lds]
dowy = dowy_full[:-max_lds]
tocs = tocs_full[:-max_lds]
df_idx = df_idx_full[:-max_lds]

#swm test data
swm_ind = 'swm'
samp_no = swm_test_samp
Q_full,dowy_full,tocs_full,df_idx_full = syn_util.extract_swm(sd_swm,ed_swm,path,loc=loc,site=site,use_bcf=use_bcf,k_swm=k_swm,same_swm=same_swm,swm_ind=swm_ind,samp_no=samp_no,yr_cutoff=yr_cutoff)
Q_swm_tst = Q_full[:-max_lds]

evt_idx_swm_tst = syn_util.declust_evts_extract(Q_swm_tst,n_evts,sep)

#swm test data
swm_ind = 'swm-cc'
samp_no = swmcc_test_samp
Q_full,dowy_full,tocs_full,df_idx_full = syn_util.extract_swm(sd_swm,ed_swm,path,loc=loc,site=site,use_bcf=use_bcf,k_swm=k_swm,same_swm=same_swm,swm_ind=swm_ind,samp_no=samp_no,yr_cutoff=yr_cutoff)
Q_swmcc_tst = Q_full[:-max_lds]

evt_idx_swmcc_tst = syn_util.declust_evts_extract(Q_swmcc_tst,n_evts,sep)

#load risk curve data
#swm data
swm_ind = 'swm'
samp_no = 1
rcurve_arr_swm = np.zeros((len(rcurve_skillmods),max_lds))
firo_top_swm = np.zeros(len(rcurve_skillmods))
firo_bottom_swm = np.zeros(len(rcurve_skillmods))

for i in range(len(rcurve_skillmods)): 
    skill_mod = rcurve_skillmods[i]
    pars = pickle.load(open('out/%s/%s/param-risk-thresholds_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mod=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_cutoff=%s.pkl'%(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,skill_mod,skill_dcy,skill_tail,swm_ind,samp_no,same_swm,yr_cutoff),'rb'),encoding='latin1')
    #pars = pickle.load(open('data/%s/%s/param-risk-thresholds_tocs-reset=%s_fixed=%s_use-firo-bottom=%s_seed-%s_mod=%s_dcy=%s_tail=%s.pkl'%(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,skill_mod,skill_dcy,skill_tail),'rb'),encoding='latin1')
    risk_curve = syn_util.create_param_risk_curve((pars['lo'],pars['hi'],pars['pwr'],pars['no_risk'],pars['all_risk']),lds=max_lds)
    firo_top = pars['firo_top']
    firo_bottom = pars['firo_bottom']
    rcurve_arr_swm[i] = risk_curve
    firo_top_swm[i] = firo_top
    firo_bottom_swm[i] = firo_bottom


#swm-cc data
swm_ind = 'swm-cc'
samp_no = 1
#CRPSS = 0, swm trained
rcurve_arr_swmcc = np.zeros((len(rcurve_skillmods),max_lds))
firo_top_swmcc = np.zeros(len(rcurve_skillmods))
firo_bottom_swmcc = np.zeros(len(rcurve_skillmods))

for i in range(len(rcurve_skillmods)): 
    skill_mod = rcurve_skillmods[i]
    pars = pickle.load(open('out/%s/%s/param-risk-thresholds_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mod=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_cutoff=%s.pkl'%(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,skill_mod,skill_dcy,skill_tail,swm_ind,samp_no,same_swm,yr_cutoff),'rb'),encoding='latin1')
    #pars = pickle.load(open('data/%s/%s/param-risk-thresholds_tocs-reset=%s_fixed=%s_use-firo-bottom=%s_seed-%s_mod=%s_dcy=%s_tail=%s.pkl'%(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,skill_mod,skill_dcy,skill_tail),'rb'),encoding='latin1')
    risk_curve = syn_util.create_param_risk_curve((pars['lo'],pars['hi'],pars['pwr'],pars['no_risk'],pars['all_risk']),lds=max_lds)
    firo_top = pars['firo_top']
    firo_bottom = pars['firo_bottom']
    rcurve_arr_swmcc[i] = risk_curve
    firo_top_swmcc[i] = firo_top
    firo_bottom_swmcc[i] = firo_bottom
    
#swm tst array with multiple skills
swm_ind = 'swm'
samp_no = swm_test_samp

swm_tst_arr = np.zeros((len(tst_skillmods),10,np.shape(Q_swm_trn)[0]))

for i in range(len(tst_skillmods)):
    skill_mod = tst_skillmods[i]
    risk_curve = rcurve_arr_swm[hist_idx,:][0]
    Qf = syn_util.extract_skillmod_swm(sd_swm,ed_swm,samp_no,path,loc,site,syn_pct,syn_objpwr,syn_optstrat,syn_setup,same_swm,K,skill_mod,skill_dcy,skill_tail,swm_ind,yr_cutoff)
    Qf = Qf[:,:,:max_lds] # just use 14 lead days
    ne = np.shape(Qf)[1]
    nl = max_lds    
    Qf_summed = np.cumsum(Qf, axis=2)
    Qf_summed_sorted = np.sort(Qf_summed, axis = 1)
    ix = ((1 - risk_curve) * (ne)).astype(np.int32)-1
    S, R, firot, firob, spill, Q_cp, rel_leads, firo_ovg, ddown = model.simulate_nonjit(firo_pool_top=firo_top_swm[hist_idx][0], firo_pool_bottom=firo_bottom_swm[hist_idx][0], ix=ix, Q=Q_swm_tst, Qf=Qf_summed_sorted, dowy=dowy, tocs=tocs, K = K, Rmax = Rmax, ramping_rate = ramping_rate, store_vec=store_vec, elev_vec=elev_vec, elev_ctrl_vec=elev_ctrl_vec, flow_ctrl_vec=flow_ctrl_vec,policy='firo', tocs_reset='none',fixed_top=use_firo_top, fixed_bottom=use_firo_bottom)
    swm_tst_arr[i,0,:] = S
    swm_tst_arr[i,1,:] = R
    swm_tst_arr[i,2,:] = firot
    swm_tst_arr[i,3,:] = spill
    swm_tst_arr[i,4,:] = rel_leads
    swm_tst_arr[i,5,:] = firo_ovg
    swm_tst_arr[i,6,:] = ddown
    swm_tst_arr[i,7,:] = firob
    swm_tst_arr[i,8,:] = firo_ovg / (K_ORO - firot) #firo overage expressed as percent of flood space encroachment
    swm_tst_arr[i,9,:] = np.append(np.abs(np.diff(R)),0) #ramping rates

np.savez('out/%s/%s/swm-test-array_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(tst_skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff),arr=swm_tst_arr)

swm_tst_arr_skilltrn = np.zeros((len(tst_skillmods),10,np.shape(Q_swm_trn)[0]))

for i in range(len(tst_skillmods)):
    skill_mod = tst_skillmods[i]
    risk_curve = rcurve_arr_swm[i,:]
    Qf = syn_util.extract_skillmod_swm(sd_swm,ed_swm,samp_no,path,loc,site,syn_pct,syn_objpwr,syn_optstrat,syn_setup,same_swm,K,skill_mod,skill_dcy,skill_tail,swm_ind,yr_cutoff)
    Qf = Qf[:,:,:max_lds] # just use 14 lead days
    ne = np.shape(Qf)[1]
    nl = max_lds    
    Qf_summed = np.cumsum(Qf, axis=2)
    Qf_summed_sorted = np.sort(Qf_summed, axis = 1)
    ix = ((1 - risk_curve) * (ne)).astype(np.int32)-1
    S, R, firot, firob, spill, Q_cp, rel_leads, firo_ovg, ddown = model.simulate_nonjit(firo_pool_top=firo_top_swm[i], firo_pool_bottom=firo_bottom_swm[i], ix=ix, Q=Q_swm_tst, Qf=Qf_summed_sorted, dowy=dowy, tocs=tocs, K = K, Rmax = Rmax, ramping_rate = ramping_rate, store_vec=store_vec, elev_vec=elev_vec, elev_ctrl_vec=elev_ctrl_vec, flow_ctrl_vec=flow_ctrl_vec,policy='firo', tocs_reset='none',fixed_top=use_firo_top, fixed_bottom=use_firo_bottom)
    swm_tst_arr_skilltrn[i,0,:] = S
    swm_tst_arr_skilltrn[i,1,:] = R
    swm_tst_arr_skilltrn[i,2,:] = firot
    swm_tst_arr_skilltrn[i,3,:] = spill
    swm_tst_arr_skilltrn[i,4,:] = rel_leads
    swm_tst_arr_skilltrn[i,5,:] = firo_ovg
    swm_tst_arr_skilltrn[i,6,:] = ddown
    swm_tst_arr_skilltrn[i,7,:] = firob
    swm_tst_arr_skilltrn[i,8,:] = firo_ovg / (K_ORO - firot) #firo overage expressed as percent of flood space encroachment
    swm_tst_arr_skilltrn[i,9,:] = np.append(np.abs(np.diff(R)),0) #ramping rates

np.savez('out/%s/%s/swm-test-array-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(tst_skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff),arr=swm_tst_arr_skilltrn)

#swm tst array with multiple skills
swm_ind = 'swm-cc'
samp_no = swmcc_test_samp

swmcc_tst_arr = np.zeros((len(tst_skillmods),10,np.shape(Q_swm_trn)[0]))

for i in range(len(tst_skillmods)):
    skill_mod = tst_skillmods[i]
    risk_curve = rcurve_arr_swm[hist_idx,:][0]
    Qf = syn_util.extract_skillmod_swm(sd_swm,ed_swm,samp_no,path,loc,site,syn_pct,syn_objpwr,syn_optstrat,syn_setup,same_swm,K,skill_mod,skill_dcy,skill_tail,swm_ind,yr_cutoff)
    Qf = Qf[:,:,:max_lds] # just use 14 lead days
    ne = np.shape(Qf)[1]
    nl = max_lds    
    Qf_summed = np.cumsum(Qf, axis=2)
    Qf_summed_sorted = np.sort(Qf_summed, axis = 1)
    ix = ((1 - risk_curve) * (ne)).astype(np.int32)-1
    S, R, firot, firob, spill, Q_cp, rel_leads, firo_ovg, ddown = model.simulate_nonjit(firo_pool_top=firo_top_swm[hist_idx][0], firo_pool_bottom=firo_bottom_swm[hist_idx][0], ix=ix, Q=Q_swmcc_tst, Qf=Qf_summed_sorted, dowy=dowy, tocs=tocs, K = K, Rmax = Rmax, ramping_rate = ramping_rate, store_vec=store_vec, elev_vec=elev_vec, elev_ctrl_vec=elev_ctrl_vec, flow_ctrl_vec=flow_ctrl_vec,policy='firo', tocs_reset='none',fixed_top=use_firo_top, fixed_bottom=use_firo_bottom)
    swmcc_tst_arr[i,0,:] = S
    swmcc_tst_arr[i,1,:] = R
    swmcc_tst_arr[i,2,:] = firot
    swmcc_tst_arr[i,3,:] = spill
    swmcc_tst_arr[i,4,:] = rel_leads
    swmcc_tst_arr[i,5,:] = firo_ovg
    swmcc_tst_arr[i,6,:] = ddown
    swmcc_tst_arr[i,7,:] = firob
    swmcc_tst_arr[i,8,:] = firo_ovg / (K_ORO - firot) #firo overage expressed as percent of flood space encroachment
    swmcc_tst_arr[i,9,:] = np.append(np.abs(np.diff(R)),0) #ramping rates

np.savez('out/%s/%s/swm-test-array_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(tst_skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff),arr=swmcc_tst_arr)

swmcc_tst_arr_skilltrn = np.zeros((len(tst_skillmods),10,np.shape(Q_swm_trn)[0]))

for i in range(len(tst_skillmods)):
    skill_mod = tst_skillmods[i]
    risk_curve = rcurve_arr_swmcc[i,:]
    Qf = syn_util.extract_skillmod_swm(sd_swm,ed_swm,samp_no,path,loc,site,syn_pct,syn_objpwr,syn_optstrat,syn_setup,same_swm,K,skill_mod,skill_dcy,skill_tail,swm_ind,yr_cutoff)
    Qf = Qf[:,:,:max_lds] # just use 14 lead days
    ne = np.shape(Qf)[1]
    nl = max_lds    
    Qf_summed = np.cumsum(Qf, axis=2)
    Qf_summed_sorted = np.sort(Qf_summed, axis = 1)
    ix = ((1 - risk_curve) * (ne)).astype(np.int32)-1
    S, R, firot, firob, spill, Q_cp, rel_leads, firo_ovg, ddown = model.simulate_nonjit(firo_pool_top=firo_top_swmcc[i], firo_pool_bottom=firo_bottom_swmcc[i], ix=ix, Q=Q_swmcc_tst, Qf=Qf_summed_sorted, dowy=dowy, tocs=tocs, K = K, Rmax = Rmax, ramping_rate = ramping_rate, store_vec=store_vec, elev_vec=elev_vec, elev_ctrl_vec=elev_ctrl_vec, flow_ctrl_vec=flow_ctrl_vec,policy='firo', tocs_reset='none',fixed_top=use_firo_top, fixed_bottom=use_firo_bottom)
    swmcc_tst_arr_skilltrn[i,0,:] = S
    swmcc_tst_arr_skilltrn[i,1,:] = R
    swmcc_tst_arr_skilltrn[i,2,:] = firot
    swmcc_tst_arr_skilltrn[i,3,:] = spill
    swmcc_tst_arr_skilltrn[i,4,:] = rel_leads
    swmcc_tst_arr_skilltrn[i,5,:] = firo_ovg
    swmcc_tst_arr_skilltrn[i,6,:] = ddown
    swmcc_tst_arr_skilltrn[i,7,:] = firob
    swmcc_tst_arr_skilltrn[i,8,:] = firo_ovg / (K_ORO - firot) #firo overage expressed as percent of flood space encroachment
    swmcc_tst_arr_skilltrn[i,9,:] = np.append(np.abs(np.diff(R)),0) #ramping rates

np.savez('out/%s/%s/swm-test-array-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(tst_skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff),arr=swmcc_tst_arr_skilltrn)


########################################################END#####################################################
