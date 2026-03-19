# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 11:26:54 2024

@author: zpb4
"""
import sys
import os
sys.path.insert(0, os.path.abspath('../Synthetic-Forecast_Verification/src'))
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
import matplotlib.colors as colors
import ensemble_verification_functions as verify
import pickle
import datetime as dt
import pyextremes

col_cb = sns.color_palette('colorblind')
sns.palplot(col_cb)  #plot coloblind palette for comparison
colv000 = sns.color_palette('Greens',10)
colv0 = sns.color_palette('Blues',10)
colv00 = sns.color_palette('Purples',10)
colv1 = sns.color_palette('PuRd',10)
colv2 = sns.color_palette('YlOrBr',10)
sns.palplot(colv0)
sns.palplot(colv00)
sns.palplot(colv000)
sns.palplot(colv1)
sns.palplot(colv2)

swm_skn4 = col_cb[2]
swm_sk0 = col_cb[9]  
swm_skmod = colv1[4]  
swm_skhi = colv1[7]  

swm_sk0_skn4 = colv000[3]
swm_sk0_sk0 = colv0[3]
swm_sk0_skmod = colv0[6]
swm_sk0_skhi = colv0[9]

swmcc_sk0_skn4 = colv000[3]
swmcc_sk0_sk0 = colv00[3]
swmcc_sk0_skmod = colv00[6]
swmcc_sk0_skhi = colv00[9]

swmcc_skn4 = col_cb[2]
swmcc_sk0 = colv2[4]  
swmcc_skmod = colv2[7]  
swmcc_skhi = colv1[9]  

cb_grn = col_cb[2]
cb_org = col_cb[3]
cb_brwn = col_cb[5]
cb_gry = col_cb[7]

#sns.palplot((col_cb[4],cv1,col_cb[3],cv2)) #base colors pretty close to 'colorblind' palette
kcfs_to_tafd = 2.29568411*10**-5 * 86400
max_lds = 15

skillmods = [-0.4,0,0.4,0.8]
skill_dcy = 0.1
skill_tail = 0.5
yr_cutoff = 500
swm_tst_samp = 2
swmcc_tst_samp = 3

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
use_firo_bottom = False
use_firo_top = False
seed = 41

syn_pct = 0.99
syn_setup = 'cal'
syn_optstrat = 'ecrps-dts'
syn_objpwr = 0
syn_path = '../Synthetic-Forecast-v2-FIRO-DISES'  # path to R synthetic forecast repo for 'r-gen' setting below

use_bcf = False
k_swm = 106
same_swm = False

kcfs_to_tafd = 2.29568411*10**-5 * 86400

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

#SWM obs data
#swm training data
swm_ind = 'swm'
samp_no = 1
Q_full,dowy_full,tocs_full,df_idx_full= syn_util.extract_swm(sd_swm,ed_swm,'../Synthetic-Forecast-v2-FIRO-DISES',loc=loc,site=site,use_bcf=use_bcf,k_swm=k_swm,same_swm=same_swm,swm_ind=swm_ind,samp_no=samp_no,yr_cutoff=yr_cutoff)
Q_swm_trn = Q_full[:-max_lds]
dowy = dowy_full[:-max_lds]
tocs = tocs_full[:-max_lds]
df_idx = df_idx_full[:-max_lds]

#swm test data
swm_ind = 'swm'
samp_no = swm_tst_samp
Q_full,dowy_full,tocs_full,df_idx_full = syn_util.extract_swm(sd_swm,ed_swm,'../Synthetic-Forecast-v2-FIRO-DISES',loc=loc,site=site,use_bcf=use_bcf,k_swm=k_swm,same_swm=same_swm,swm_ind=swm_ind,samp_no=samp_no,yr_cutoff=yr_cutoff)
Q_swm_tst = Q_full[:-max_lds]

evt_idx_swm_tst = syn_util.declust_evts_extract(Q_swm_tst,n_evts,sep)

#swm test data
swm_ind = 'swm-cc'
samp_no = swmcc_tst_samp
Q_full,dowy_full,tocs_full,df_idx_full = syn_util.extract_swm(sd_swm,ed_swm,'../Synthetic-Forecast-v2-FIRO-DISES',loc=loc,site=site,use_bcf=use_bcf,k_swm=k_swm,same_swm=same_swm,swm_ind=swm_ind,samp_no=samp_no,yr_cutoff=yr_cutoff)
Q_swmcc_tst = Q_full[:-max_lds]

evt_idx_swmcc_tst = syn_util.declust_evts_extract(Q_swmcc_tst,n_evts,sep)

#load risk curve data
#swm data
swm_ind = 'swm'
samp_no = 1
rcurve_arr_swm = np.zeros((len(skillmods),max_lds))
firo_top_swm = np.zeros(len(skillmods))
firo_bottom_swm = np.zeros(len(skillmods))

for i in range(len(skillmods)): 
    skill_mod = skillmods[i]
    pars = pickle.load(open('data/%s/%s/param-risk-thresholds_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mod=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_cutoff=%s.pkl'%(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,skill_mod,skill_dcy,skill_tail,swm_ind,samp_no,same_swm,yr_cutoff),'rb'),encoding='latin1')
    #pars = pickle.load(open('data/%s/%s/param-risk-thresholds_tocs-reset=%s_fixed=%s_use-firo-bottom=%s_seed-%s_mod=%s_dcy=%s_tail=%s.pkl'%(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,skill_mod,skill_dcy,skill_tail),'rb'),encoding='latin1')
    risk_curve = syn_util.create_param_risk_curve((pars['lo'],pars['hi'],pars['pwr'],pars['no_risk'],pars['all_risk']),lds=max_lds)
    firo_top = pars['firo_top']
    firo_bottom = pars['firo_bottom']
    rcurve_arr_swm[i] = risk_curve
    firo_top_swm = firo_top
    firo_bottom_swm = firo_bottom


#swm-cc data
swm_ind = 'swm-cc'
samp_no = 1
#CRPSS = 0, swm trained
rcurve_arr_swmcc = np.zeros((len(skillmods),max_lds))
firo_top_swmcc = np.zeros(len(skillmods))
firo_bottom_swmcc = np.zeros(len(skillmods))

for i in range(len(skillmods)): 
    skill_mod = skillmods[i]
    pars = pickle.load(open('data/%s/%s/param-risk-thresholds_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mod=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_cutoff=%s.pkl'%(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,skill_mod,skill_dcy,skill_tail,swm_ind,samp_no,same_swm,yr_cutoff),'rb'),encoding='latin1')
    #pars = pickle.load(open('data/%s/%s/param-risk-thresholds_tocs-reset=%s_fixed=%s_use-firo-bottom=%s_seed-%s_mod=%s_dcy=%s_tail=%s.pkl'%(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,skill_mod,skill_dcy,skill_tail),'rb'),encoding='latin1')
    risk_curve = syn_util.create_param_risk_curve((pars['lo'],pars['hi'],pars['pwr'],pars['no_risk'],pars['all_risk']),lds=max_lds)
    firo_top = pars['firo_top']
    firo_bottom = pars['firo_bottom']
    rcurve_arr_swmcc[i] = risk_curve
    firo_top_swmcc = firo_top
    firo_bottom_swmcc = firo_bottom

#swm tst array with multiple skills
swm_ind = 'swm'
samp_no = swm_tst_samp

swm_tst_arr = np.load('data/%s/%s/swm-test-array_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']
swm_tst_arr_skilltrn = np.load('data/%s/%s/swm-test-array-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']

#swm tst array with multiple skills
swm_ind = 'swm-cc'
samp_no = swmcc_tst_samp

swmcc_tst_arr = np.load('data/%s/%s/swm-test-array_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']
swmcc_tst_arr_skilltrn = np.load('data/%s/%s/swm-test-array-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']


#............................................................................................................
#peaks over threshold processing for comparison
th = 0.99
pot_dys = 7     # window for POT analysis; maximum distance between peaks in a cluster over the threshold
pk_thresh = 0  # TAF/d limit for 

Q_inp_swm = pd.Series(Q_swm_tst,index=df_idx)
Q_thresh_swm = np.sort(Q_inp_swm)[round(th*len(Q_inp_swm))]

Q_inp_swmcc = pd.Series(Q_swmcc_tst,index=df_idx)
Q_thresh_swmcc = np.sort(Q_inp_swmcc)[round(th*len(Q_inp_swmcc))]

pot_swm = pyextremes.get_extremes(Q_inp_swm, 'POT',threshold=Q_thresh_swm,r=str(pot_dys*24)+'h')
pot_swm_np = pd.Series.to_numpy(pot_swm)

pot_dur_swm = np.zeros_like(pot_swm_np)
avg_mag_swm = np.zeros_like(pot_swm_np)
pks_swm = np.zeros_like(pot_swm_np)
spills_swm_pot = np.zeros((len(skillmods),len(pot_swm_np)))
spillvol_swm_pot = np.zeros((len(skillmods),len(pot_swm_np)))
Rmx_swm_pot = np.zeros((len(skillmods),len(pot_swm_np)))
Rmn_swm_pot = np.zeros((len(skillmods),len(pot_swm_np)))
fovg_swm_pot = np.zeros((len(skillmods),len(pot_swm_np)))

for i in range(len(pot_dur_swm)):
    evt_idx = np.where(Q_swm_tst == pot_swm_np[i])[0]
    sset_idx =np.arange(len(Q_swm_tst))[max((evt_idx-sep),0)[0]:min((evt_idx+sep),len(Q_swm_tst))[0]]
    sset = Q_swm_tst[sset_idx]
    sset_pks = np.concat((([0]),np.diff(sset)))
    sset_pks_ind = np.zeros_like(sset_pks)
    sset_pks_ind[sset_pks>0] = 1
    sset_pks_ind[sset_pks<0] = -1
    # select peaks where a) rising and falling limb are above the 'pk_thresh' value and b) peaks are above the specified threshold
    sset_pks_selec = np.where((np.concat((np.diff(sset_pks_ind),([0])))==-2) & (sset_pks>pk_thresh) & (np.concat((np.diff(sset),([0])))<(-pk_thresh)) & (sset>Q_thresh_swm))
    pks_swm[i] = np.shape(sset_pks_selec)[1]
    pks_vec = np.zeros_like(sset)
    pks_vec[sset>Q_thresh_swm] = 1
    pks_vec[sset<=Q_thresh_swm] = 0
    pks_ind = np.diff(pks_vec)
    pre_vec = np.flip(pks_ind[:sep])
    post_vec = pks_ind[sep:]
    pot_dur_swm[i] = np.min(np.where(pre_vec!=0)) + np.min(np.where(post_vec!=0)) + 1
    pre_idx = sep - np.min(np.where(pre_vec!=0))
    pst_idx = sep + np.min(np.where(post_vec!=0))
    avg_mag_swm[i] = np.sum(sset[pre_idx:(pst_idx+1)]) / pot_dur_swm[i]
    #pks_swm[i] = min(np.sum(pks_ind>0),np.sum(pks_ind<0))
    for k in range(len(skillmods)):
        spills_swm_pot[k,i] = np.count_nonzero(swm_tst_arr[k,3,sset_idx])
        spillvol_swm_pot[k,i] = np.sum(swm_tst_arr[k,3,sset_idx])
        Rmx_swm_pot[k,i] = np.max(swm_tst_arr[k,1,sset_idx])
        Rmn_swm_pot[k,i] = np.mean(swm_tst_arr[k,1,sset_idx])
        fovg_swm_pot[k,i] = np.max(swm_tst_arr[k,8,sset_idx])

pot_swmcc = pyextremes.get_extremes(Q_inp_swmcc, 'POT',threshold=Q_thresh_swmcc,r=str(sep*24)+'h')
pot_swmcc_np = pd.Series.to_numpy(pot_swmcc)

pot_dur_swmcc = np.zeros_like(pot_swmcc_np)
avg_mag_swmcc = np.zeros_like(pot_swmcc_np)
pks_swmcc = np.zeros_like(pot_swmcc_np)
spills_swmcc_pot = np.zeros((len(skillmods),len(pot_swmcc_np)))
spillvol_swmcc_pot = np.zeros((len(skillmods),len(pot_swmcc_np)))
Rmx_swmcc_pot = np.zeros((len(skillmods),len(pot_swmcc_np)))
Rmn_swmcc_pot = np.zeros((len(skillmods),len(pot_swmcc_np)))
fovg_swmcc_pot = np.zeros((len(skillmods),len(pot_swmcc_np)))

for i in range(len(pot_dur_swmcc)):
    evt_idx = np.where(Q_swmcc_tst == pot_swmcc_np[i])[0]
    sset_idx =np.arange(len(Q_swmcc_tst))[max((evt_idx-sep),0)[0]:min((evt_idx+sep),len(Q_swmcc_tst))[0]]
    sset = Q_swmcc_tst[sset_idx]
    sset_pks = np.concat((([0]),np.diff(sset)))
    sset_pks_ind = np.zeros_like(sset_pks)
    sset_pks_ind[sset_pks>0] = 1
    sset_pks_ind[sset_pks<0] = -1
    # select peaks where a) rising and falling limb are above the 'pk_thresh' value and b) peaks are above the specified threshold
    sset_pks_selec = np.where((np.concat((np.diff(sset_pks_ind),([0])))==-2) & (sset_pks>pk_thresh) & (np.concat((np.diff(sset),([0])))<(-pk_thresh)) & (sset>Q_thresh_swm))
    pks_swmcc[i] = np.shape(sset_pks_selec)[1]
    pks_vec = np.zeros_like(sset)
    pks_vec[sset>Q_thresh_swm] = 1
    pks_vec[sset<=Q_thresh_swm] = 0
    pks_ind = np.diff(pks_vec)
    pre_vec = np.flip(pks_ind[:sep])
    post_vec = pks_ind[sep:]
    pot_dur_swmcc[i] = np.min(np.where(pre_vec!=0)) + np.min(np.where(post_vec!=0)) + 1
    pre_idx = sep - np.min(np.where(pre_vec!=0))
    pst_idx = sep + np.min(np.where(post_vec!=0))
    avg_mag_swmcc[i] = np.sum(sset[pre_idx:(pst_idx+1)]) / pot_dur_swm[i]
    #pks_swm[i] = min(np.sum(pks_ind>0),np.sum(pks_ind<0))
    for k in range(len(skillmods)):
        spills_swmcc_pot[k,i] = np.count_nonzero(swmcc_tst_arr[k,3,sset_idx])
        spillvol_swmcc_pot[k,i] = np.sum(swmcc_tst_arr[k,3,sset_idx])
        Rmx_swmcc_pot[k,i] = np.max(swmcc_tst_arr[k,1,sset_idx])
        Rmn_swmcc_pot[k,i] = np.mean(swmcc_tst_arr[k,1,sset_idx])
        fovg_swmcc_pot[k,i] = np.max(swmcc_tst_arr[k,8,sset_idx])

#compare magnitude-duration distributions of swm & swm-cc
plt.scatter(pot_dur_swm,pot_swm_np,color=cb_gry,alpha=0.25)
plt.scatter(pot_dur_swmcc+0.4,pot_swmcc_np,color=cb_org,alpha=0.25)
plt.xlim([1,19])
plt.xticks(ticks=np.arange(21)+0.2,labels=np.arange(21))
plt.xlabel('Event duration (days > 99%)')
plt.ylabel('Peak flow (TAF/d)')
plt.ylim([0,1200])
plt.text(15,1000,'$n=%s$' %(str(len(pot_dur_swm))),color=cb_gry,fontsize='large')
plt.text(15,900,'$n=%s$' %(str(len(pot_dur_swmcc))),color=cb_org,fontsize='large')
plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)
plt.show()

#compare number of peaks between swm & swm-cc
swm_hist = np.histogram(pks_swm,bins=np.int64(np.max(pks_swm)-1))
swmcc_hist = np.histogram(pks_swmcc,bins=np.int64(np.max(pks_swmcc)-1))
swm_histx = swm_hist[1][:-1]-0.2
swm_histy = swm_hist[0]/np.sum(swm_hist[0])
swmcc_histx = swmcc_hist[1][:-1]+0.2
swmcc_histy = swmcc_hist[0]/np.sum(swmcc_hist[0])
plt.bar(swm_histx,swm_histy,color=cb_gry,width=0.35)
plt.bar(swmcc_histx,swmcc_histy,color=cb_org,width=0.35)
plt.legend(['$swm$','$swm_{4C}$'],loc='upper right',fontsize='large',title_fontsize='large',frameon=False)
plt.xlabel('Number of peaks > 99%')
plt.ylabel('Density')
plt.text(4.5,0.45,'$n=%s$' %(str(len(pot_dur_swm))),color=cb_gry,fontsize='large')
plt.text(4.5,0.4,'$n=%s$' %(str(len(pot_dur_swmcc))),color=cb_org,fontsize='large')
plt.text(4,0.3,'$Q_{99}=%s\ TAF/d$' %(str(round(Q_thresh_swm,1))),color=cb_gry,fontsize='large')
plt.text(4,0.25,'$Q_{99}=%s\ TAF/d$' %(str(round(Q_thresh_swmcc,1))),color=cb_org,fontsize='large')
plt.show()

#--------------------------------------------------------------------------------------------------
#map event duration & magnitude & max release
bounds = np.array([0,50,100,150,200,250,300])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

fig = plt.figure(layout='constrained',figsize=(12,6))
gs0 = fig.add_gridspec(2,5,width_ratios=(1,1,1,1,0.1))
ax1 = fig.add_subplot(gs0[0,0])

rmx_img = np.zeros((len(Rmx_swm_pot[0,:]),1))
rmx_img[:,0] = Rmx_swm_pot[0,:]

ax1.scatter(pot_dur_swm,pot_swm_np,c=Rmx_swm_pot[0,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,19])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(19)+0.2,labels=np.arange(19))
#ax1.set_xlabel('Event duration (days > 99%)')
ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'a)',fontsize='xx-large',fontweight='bold')
ax1.text(6,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(6,1000,'$CRPS-SS_{test}=%s_{swm}$' %(skillmods[0]),fontsize='medium')
ax1.xaxis.set_ticklabels([])
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.0_{swm}$')
#cax = ax1.imshow(Rmx_swm0tst_pot,cmap=plt.cm.coolwarm,norm=norm)
#fig.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical',ax=ax1)
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[0,1])

ax1.scatter(pot_dur_swm,pot_swm_np,c=Rmx_swm_pot[1,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,19])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(19)+0.2,labels=np.arange(19))
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
#ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'b)',fontsize='xx-large',fontweight='bold')
ax1.text(6,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(6,1000,'$CRPS-SS_{test}=%s_{swm}$' %(skillmods[1]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[0,2])

ax1.scatter(pot_dur_swm,pot_swm_np,c=Rmx_swm_pot[2,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,19])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(19)+0.2,labels=np.arange(19))
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
#ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'c)',fontsize='xx-large',fontweight='bold')
ax1.text(6,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(6,1000,'$CRPS-SS_{test}=%s_{swm}$' %(skillmods[2]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[0,3])

ax1.scatter(pot_dur_swm,pot_swm_np,c=Rmx_swm_pot[3,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,19])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(19)+0.2,labels=np.arange(19))
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
#ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'d)',fontsize='xx-large',fontweight='bold')
ax1.text(6,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(6,1000,'$CRPS-SS_{test}=%s_{swm}$' %(skillmods[3]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,0])

ax1.scatter(pot_dur_swmcc,pot_swmcc_np,c=Rmx_swmcc_pot[0,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,19])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(19)+0.2,labels=np.arange(19))
ax1.set_xlabel('Event duration (days > 99%)')
ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'e)',fontsize='xx-large',fontweight='bold')
ax1.text(6,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(6,1000,'$CRPS-SS_{test}=%s_{swm4C}$' %(skillmods[0]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.0_{swm4C}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,1])
ax1.scatter(pot_dur_swmcc,pot_swmcc_np,c=Rmx_swmcc_pot[1,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,19])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(19)+0.2,labels=np.arange(19))
ax1.yaxis.set_ticklabels([])
ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'f)',fontsize='xx-large',fontweight='bold')
ax1.text(6,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(6,1000,'$CRPS-SS_{test}=%s_{swm4C}$' %(skillmods[1]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm4C}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,2])
ax1.scatter(pot_dur_swmcc,pot_swmcc_np,c=Rmx_swmcc_pot[2,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,19])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(19)+0.2,labels=np.arange(19))
ax1.yaxis.set_ticklabels([])
ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'g)',fontsize='xx-large',fontweight='bold')
ax1.text(6,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(6,1000,'$CRPS-SS_{test}=%s_{swm4C}$' %(skillmods[2]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm4C}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,3])
ax1.scatter(pot_dur_swmcc,pot_swmcc_np,c=Rmx_swmcc_pot[3,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,19])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(19)+0.2,labels=np.arange(19))
ax1.yaxis.set_ticklabels([])
ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'h)',fontsize='xx-large',fontweight='bold')
ax1.text(6,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(6,1000,'$CRPS-SS_{test}=%s_{swm4C}$' %(skillmods[3]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm4C}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,4])
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
im1 = ax1.imshow(rmx_img,cmap=plt.cm.coolwarm,norm=norm,alpha=0.25)
cbar_ax = fig.add_subplot(gs0[:,4])
cbar = fig.colorbar(im1, cax=cbar_ax)
cbar.set_label('$R_{max}\ (TAF/d)$',fontsize='large')

fig.savefig('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ops_plots/2x4_magnitude-duration-Rmax_use-firo-top=%s_use-firo-bottom=%s_seed-%s_evt=%s_mods=%s_dcy=%s_tail=%s_swm=%s_same-swm=%s_cutoff=%s.png' %(use_firo_top,use_firo_bottom,seed,evt_no,str(skillmods),skill_dcy,skill_tail,swm_ind,same_swm,yr_cutoff),dpi=300,bbox_inches='tight')


#--------------------------------------------------------------------------------------------------
#map # of peaks & magnitude & max release
bounds = np.array([0,50,100,150,200,250,300])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

fig = plt.figure(layout='constrained',figsize=(12.5,6))
gs0 = fig.add_gridspec(2,5,width_ratios=(1,1,1,1,0.1))
ax1 = fig.add_subplot(gs0[0,0])

rmx_img = np.zeros((len(Rmx_swm_pot[0,:]),1))
rmx_img[:,0] = Rmx_swm_pot[0,:]

ax1.scatter(pks_swm,pot_swm_np,c=Rmx_swm_pot[0,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,10])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(10)+0.2,labels=np.arange(10))
#ax1.set_xlabel('Event duration (days > 99%)')
ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'a)',fontsize='xx-large',fontweight='bold')
ax1.text(3,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(3,1000,'$CRPS-SS_{test}=%s_{swm}$' %(skillmods[0]),fontsize='medium')
ax1.xaxis.set_ticklabels([])
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.0_{swm}$')
#cax = ax1.imshow(Rmx_swm0tst_pot,cmap=plt.cm.coolwarm,norm=norm)
#fig.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical',ax=ax1)
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[0,1])

ax1.scatter(pks_swm,pot_swm_np,c=Rmx_swm_pot[1,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,10])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(10)+0.2,labels=np.arange(10))
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
#ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'b)',fontsize='xx-large',fontweight='bold')
ax1.text(3,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(3,1000,'$CRPS-SS_{test}=%s_{swm}$' %(skillmods[1]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[0,2])

ax1.scatter(pks_swm,pot_swm_np,c=Rmx_swm_pot[2,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,10])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(10)+0.2,labels=np.arange(10))
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
#ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'c)',fontsize='xx-large',fontweight='bold')
ax1.text(3,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(3,1000,'$CRPS-SS_{test}=%s_{swm}$' %(skillmods[2]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[0,3])

ax1.scatter(pks_swm,pot_swm_np,c=Rmx_swm_pot[3,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,10])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(10)+0.2,labels=np.arange(10))
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
#ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'d)',fontsize='xx-large',fontweight='bold')
ax1.text(3,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(3,1000,'$CRPS-SS_{test}=%s_{swm}$' %(skillmods[3]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,0])

ax1.scatter(pks_swmcc,pot_swmcc_np,c=Rmx_swmcc_pot[0,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,10])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(10)+0.2,labels=np.arange(10))
ax1.set_xlabel('Number of peaks > 99%')
ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'e)',fontsize='xx-large',fontweight='bold')
ax1.text(3,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(3,1000,'$CRPS-SS_{test}=%s_{swm4C}$' %(skillmods[0]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.0_{swm4C}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,1])
ax1.scatter(pks_swmcc,pot_swmcc_np,c=Rmx_swmcc_pot[1,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,10])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(10)+0.2,labels=np.arange(10))
ax1.yaxis.set_ticklabels([])
ax1.set_xlabel('Number of peaks > 99%')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'f)',fontsize='xx-large',fontweight='bold')
ax1.text(3,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(3,1000,'$CRPS-SS_{test}=%s_{swm4C}$' %(skillmods[1]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm4C}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,2])
ax1.scatter(pks_swmcc,pot_swmcc_np,c=Rmx_swmcc_pot[2,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,10])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(10)+0.2,labels=np.arange(10))
ax1.yaxis.set_ticklabels([])
ax1.set_xlabel('Number of peaks > 99%')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'g)',fontsize='xx-large',fontweight='bold')
ax1.text(3,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(3,1000,'$CRPS-SS_{test}=%s_{swm4C}$' %(skillmods[2]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm4C}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,3])
ax1.scatter(pks_swmcc,pot_swmcc_np,c=Rmx_swmcc_pot[3,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,10])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(10)+0.2,labels=np.arange(10))
ax1.yaxis.set_ticklabels([])
ax1.set_xlabel('Number of peaks > 99%')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'h)',fontsize='xx-large',fontweight='bold')
ax1.text(3,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(3,1000,'$CRPS-SS_{test}=%s_{swm4C}$' %(skillmods[3]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm4C}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,4])
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
im1 = ax1.imshow(rmx_img,cmap=plt.cm.coolwarm,norm=norm,alpha=0.25)
cbar_ax = fig.add_subplot(gs0[:,4])
cbar = fig.colorbar(im1, cax=cbar_ax)
cbar.set_label('$R_{max}\ (TAF/d)$',fontsize='large')

fig.savefig('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ops_plots/2x4_magnitude-nopeaks-Rmax_use-firo-top=%s_use-firo-bottom=%s_seed-%s_evt=%s_mods=%s_dcy=%s_tail=%s_swm=%s_same-swm=%s_cutoff=%s.png' %(use_firo_top,use_firo_bottom,seed,evt_no,str(skillmods),skill_dcy,skill_tail,swm_ind,same_swm,yr_cutoff),dpi=300,bbox_inches='tight')

#--------------------------------------------------------------------------------------------------
#map event duration & magnitude & FIRO pool overage
#bounds = np.array([0,100,200,300,400,500,600])
#bounds = np.array([0,50,100,150,200,250,300,350])
bounds = np.array([0,0.5,0.6,0.7,0.8,0.9,1.0])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

fig = plt.figure(layout='constrained',figsize=(12,6))
gs0 = fig.add_gridspec(2,5,width_ratios=(1,1,1,1,0.1))
ax1 = fig.add_subplot(gs0[0,0])

rmx_img = np.zeros((len(fovg_swm_pot[0,:]),1))
rmx_img[:,0] = fovg_swm_pot[0,:]

ax1.scatter(pot_dur_swm,pot_swm_np,c=fovg_swm_pot[0,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,19])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(19)+0.2,labels=np.arange(19))
#ax1.set_xlabel('Event duration (days > 99%)')
ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'a)',fontsize='xx-large',fontweight='bold')
ax1.text(6,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(6,1000,'$CRPS-SS_{test}=%s_{swm}$' %(skillmods[0]),fontsize='medium')
ax1.xaxis.set_ticklabels([])
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.0_{swm}$')
#cax = ax1.imshow(Rmx_swm0tst_pot,cmap=plt.cm.coolwarm,norm=norm)
#fig.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical',ax=ax1)
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[0,1])

ax1.scatter(pot_dur_swm,pot_swm_np,c=fovg_swm_pot[1,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,19])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(19)+0.2,labels=np.arange(19))
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
#ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'b)',fontsize='xx-large',fontweight='bold')
ax1.text(6,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(6,1000,'$CRPS-SS_{test}=%s_{swm}$' %(skillmods[1]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[0,2])

ax1.scatter(pot_dur_swm,pot_swm_np,c=fovg_swm_pot[2,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,19])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(19)+0.2,labels=np.arange(19))
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
#ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'c)',fontsize='xx-large',fontweight='bold')
ax1.text(6,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(6,1000,'$CRPS-SS_{test}=%s_{swm}$' %(skillmods[2]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[0,3])

ax1.scatter(pot_dur_swm,pot_swm_np,c=fovg_swm_pot[3,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,19])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(19)+0.2,labels=np.arange(19))
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
#ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'d)',fontsize='xx-large',fontweight='bold')
ax1.text(6,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(6,1000,'$CRPS-SS_{test}=%s_{swm}$' %(skillmods[3]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,0])

ax1.scatter(pot_dur_swmcc,pot_swmcc_np,c=fovg_swmcc_pot[0,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,19])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(19)+0.2,labels=np.arange(19))
ax1.set_xlabel('Event duration (days > 99%)')
ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'e)',fontsize='xx-large',fontweight='bold')
ax1.text(6,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(6,1000,'$CRPS-SS_{test}=%s_{swm4C}$' %(skillmods[0]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.0_{swm4C}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,1])
ax1.scatter(pot_dur_swmcc,pot_swmcc_np,c=fovg_swmcc_pot[1,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,19])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(19)+0.2,labels=np.arange(19))
ax1.yaxis.set_ticklabels([])
ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'f)',fontsize='xx-large',fontweight='bold')
ax1.text(6,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(6,1000,'$CRPS-SS_{test}=%s_{swm4C}$' %(skillmods[0]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm4C}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,2])
ax1.scatter(pot_dur_swmcc,pot_swmcc_np,c=fovg_swmcc_pot[2,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,19])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(19)+0.2,labels=np.arange(19))
ax1.yaxis.set_ticklabels([])
ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'g)',fontsize='xx-large',fontweight='bold')
ax1.text(6,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(6,1000,'$CRPS-SS_{test}=%s_{swm4C}$' %(skillmods[2]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm4C}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)


ax1 = fig.add_subplot(gs0[1,3])
ax1.scatter(pot_dur_swmcc,pot_swmcc_np,c=fovg_swmcc_pot[3,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,19])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(19)+0.2,labels=np.arange(19))
ax1.yaxis.set_ticklabels([])
ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'h)',fontsize='xx-large',fontweight='bold')
ax1.text(6,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(6,1000,'$CRPS-SS_{test}=%s_{swm4C}$' %(skillmods[3]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm4C}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)


ax1 = fig.add_subplot(gs0[1,4])
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
im1 = ax1.imshow(rmx_img,cmap=plt.cm.coolwarm,norm=norm,alpha=0.25)
cbar_ax = fig.add_subplot(gs0[:,4])
cbar = fig.colorbar(im1, cax=cbar_ax)
cbar.set_label('$FIRO\ pool\ overage \ (TAF)$',fontsize='large')

fig.savefig('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ops_plots/2x4_magnitude-duration-fovg_use-firo-top=%s_use-firo-bottom=%s_seed-%s_evt=%s_mods=%s_dcy=%s_tail=%s_swm=%s_same-swm=%s_cutoff=%s.png' %(use_firo_top,use_firo_bottom,seed,evt_no,str(skillmods),skill_dcy,skill_tail,swm_ind,same_swm,yr_cutoff),dpi=300,bbox_inches='tight')

#--------------------------------------------------------------------------------------------------
#map event peaks & magnitude & FIRO pool overage
#bounds = np.array([0,100,200,300,400,500,600])
#bounds = np.array([0,50,100,150,200,250,300,350])
bounds = np.array([0,0.5,0.6,0.7,0.8,0.9,1.0])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

fig = plt.figure(layout='constrained',figsize=(12.5,6))
gs0 = fig.add_gridspec(2,5,width_ratios=(1,1,1,1,0.1))
ax1 = fig.add_subplot(gs0[0,0])

rmx_img = np.zeros((len(fovg_swm_pot[0,:]),1))
rmx_img[:,0] = fovg_swm_pot[0,:]

ax1.scatter(pks_swm,pot_swm_np,c=fovg_swm_pot[0,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,10])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(10)+0.2,labels=np.arange(10))
#ax1.set_xlabel('Event duration (days > 99%)')
ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'a)',fontsize='xx-large',fontweight='bold')
ax1.text(3,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(3,1000,'$CRPS-SS_{test}=%s_{swm}$' %(skillmods[0]),fontsize='medium')
ax1.xaxis.set_ticklabels([])
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.0_{swm}$')
#cax = ax1.imshow(Rmx_swm0tst_pot,cmap=plt.cm.coolwarm,norm=norm)
#fig.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical',ax=ax1)
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[0,1])

ax1.scatter(pks_swm,pot_swm_np,c=fovg_swm_pot[1,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,10])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(10)+0.2,labels=np.arange(10))
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
#ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'b)',fontsize='xx-large',fontweight='bold')
ax1.text(3,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(3,1000,'$CRPS-SS_{test}=%s_{swm}$' %(skillmods[1]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[0,2])

ax1.scatter(pks_swm,pot_swm_np,c=fovg_swm_pot[2,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,10])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(10)+0.2,labels=np.arange(10))
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
#ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'c)',fontsize='xx-large',fontweight='bold')
ax1.text(3,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(3,1000,'$CRPS-SS_{test}=%s_{swm}$' %(skillmods[2]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[0,3])

ax1.scatter(pks_swm,pot_swm_np,c=fovg_swm_pot[3,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,10])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(10)+0.2,labels=np.arange(10))
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
#ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'d)',fontsize='xx-large',fontweight='bold')
ax1.text(3,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(3,1000,'$CRPS-SS_{test}=%s_{swm}$' %(skillmods[3]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,0])

ax1.scatter(pks_swmcc,pot_swmcc_np,c=fovg_swmcc_pot[0,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,10])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(10)+0.2,labels=np.arange(10))
ax1.set_xlabel('Number of peaks > 99%')
ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'e)',fontsize='xx-large',fontweight='bold')
ax1.text(3,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(3,1000,'$CRPS-SS_{test}=%s_{swm4C}$' %(skillmods[0]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.0_{swm4C}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,1])
ax1.scatter(pks_swmcc,pot_swmcc_np,c=fovg_swmcc_pot[0,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,10])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(10)+0.2,labels=np.arange(10))
ax1.yaxis.set_ticklabels([])
ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'f)',fontsize='xx-large',fontweight='bold')
ax1.text(3,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(3,1000,'$CRPS-SS_{test}=%s_{swm4C}$' %(skillmods[1]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm4C}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,2])
ax1.scatter(pks_swmcc,pot_swmcc_np,c=fovg_swmcc_pot[2,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,10])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(10)+0.2,labels=np.arange(10))
ax1.yaxis.set_ticklabels([])
ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'g)',fontsize='xx-large',fontweight='bold')
ax1.text(3,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(3,1000,'$CRPS-SS_{test}=%s_{swm4C}$' %(skillmods[2]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm4C}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,1])
ax1.scatter(pks_swmcc,pot_swmcc_np,c=fovg_swmcc_pot[3,:],cmap=plt.cm.coolwarm,alpha=0.25,norm=norm)
ax1.set_xlim([1,10])
ax1.set_ylim([0,1200])
ax1.set_xticks(ticks=np.arange(10)+0.2,labels=np.arange(10))
ax1.yaxis.set_ticklabels([])
ax1.set_xlabel('Event duration (days > 99%)')
#ax1.set_ylabel('Peak flow (TAF/d)')
ax1.text(1,1100,'h)',fontsize='xx-large',fontweight='bold')
ax1.text(3,1100,'$CRPS-SS_{train}=0.0_{swm}$',fontsize='medium')
ax1.text(3,1000,'$CRPS-SS_{test}=%s_{swm4C}$' %(skillmods[3]),fontsize='medium')
#ax1.set_title('$CRPS-SS_{train}=0.0_{swm};\ CRPS-SS_{test}=0.8_{swm4C}$')
#plt.colorbar(label='$R_{max}\ (TAF/d)$',orientation='vertical')
#plt.legend(['$swm$','$swm_{4C}$'],loc='upper left',fontsize='large',title_fontsize='large',frameon=False)

ax1 = fig.add_subplot(gs0[1,4])
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
im1 = ax1.imshow(rmx_img,cmap=plt.cm.coolwarm,norm=norm,alpha=0.25)
cbar_ax = fig.add_subplot(gs0[:,4])
cbar = fig.colorbar(im1, cax=cbar_ax)
cbar.set_label('$FIRO\ pool\ overage \ (TAF)$',fontsize='large')

fig.savefig('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ops_plots/2x4_magnitude-nopeaks-fovg_use-firo-top=%s_use-firo-bottom=%s_seed-%s_evt=%s_mods=%s_dcy=%s_tail=%s_swm=%s_same-swm=%s_cutoff=%s.png' %(use_firo_top,use_firo_bottom,seed,evt_no,str(skillmods),skill_dcy,skill_tail,swm_ind,same_swm,yr_cutoff),dpi=300,bbox_inches='tight')

#--------------------------------------end--------------------------------------
