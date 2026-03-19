# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 11:26:54 2024

@author: zpb4
"""
import sys
import os
sys.path.insert(0, os.path.abspath('../Synthetic-Forecast_Verification/src/scratch'))
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


kcfs_to_tafd = 2.29568411*10**-5 * 86400
max_lds = 15

skillmods = skillmods = ([-0.4,0,0.4,0.8])
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
firo_space_ORO = model.firo_space_ORO
firo_top_ORO = model.firo_top_ORO

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

#load data
#load FVA firo pool runs
use_firo_bottom = True
use_firo_top = True
#swm data
swm_ind = 'swm'
samp_no = swm_tst_samp

swm_tst_metrics = np.load('data/%s/%s/swm-test-metrics_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']
swm_tst_metrics_skilltrn = np.load('data/%s/%s/swm-test-metrics-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']
swm_tst_apr1 = np.load('data/%s/%s/swm-test-apr1_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']
swm_tst_apr1_skilltrn = np.load('data/%s/%s/swm-test-apr1-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']

#swmcc data
swm_ind = 'swm-cc'
samp_no = swmcc_tst_samp

swmcc_tst_metrics = np.load('data/%s/%s/swm-test-metrics_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']
swmcc_tst_metrics_skilltrn = np.load('data/%s/%s/swm-test-metrics-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']
swmcc_tst_apr1 = np.load('data/%s/%s/swm-test-apr1_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']
swmcc_tst_apr1_skilltrn = np.load('data/%s/%s/swm-test-apr1-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']


#load flexible FIRO pool runs
use_firo_bottom = False
use_firo_top = False
#swm data
swm_ind = 'swm'
samp_no = swm_tst_samp

swm_tstfpool_metrics = np.load('data/%s/%s/swm-test-metrics_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']
swm_tstfpool_metrics_skilltrn = np.load('data/%s/%s/swm-test-metrics-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']
swm_tstfpool_apr1 = np.load('data/%s/%s/swm-test-apr1_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']
swm_tstfpool_apr1_skilltrn = np.load('data/%s/%s/swm-test-apr1-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']

#swmcc data
swm_ind = 'swm-cc'
samp_no = swmcc_tst_samp

swmcc_tstfpool_metrics = np.load('data/%s/%s/swm-test-metrics_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']
swmcc_tstfpool_metrics_skilltrn = np.load('data/%s/%s/swm-test-metrics-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']
swmcc_tstfpool_apr1 = np.load('data/%s/%s/swm-test-apr1_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']
swmcc_tstfpool_apr1_skilltrn = np.load('data/%s/%s/swm-test-apr1-skilltrain_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_nevts=%s_cutoff=%s.npz' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,str(skillmods),skill_dcy,skill_tail,swm_ind,samp_no,same_swm,n_evts,yr_cutoff))['arr']

#/////////////////////////////////////////////////////////////////////////////////////////////////////////////
#.............................................................................................................
#Plot 1: Storage metric w/ and w/o FIRO pool flexibility and w/ and w/o CC
#Event based metrics
#panel 1
swm_ind = 'swm'
nb = len(skillmods)

ylow = 0.75
ylm = ([ylow*K_ORO,K_ORO+(1-ylow)*K_ORO*0.05])
xlm = ([0.5,nb+0.5])
upperq = 0.9
lowerq = 0.1
use_err = True
ref_line = model.top_of_varflood_ORO

#store_inp1 = swm_tst_metrics[13,:,:]
#store_inp2 = swm_tst_metrics_skilltrn[13,:,:] 

store_inp1 = swm_tst_apr1
store_inp2 = swm_tst_apr1_skilltrn

x1 = np.arange(nb)+0.8
x2 = np.arange(nb)+1.2

y1 = np.mean(store_inp1,axis=0)
y2 = np.mean(store_inp2,axis=0)

yerr1 = np.zeros((2,nb))
yerr1[1,:] = np.quantile(store_inp1,upperq,axis=0) - y1
yerr1[0,:] = y1 - np.quantile(store_inp1,lowerq,axis=0)

yerr2 = np.zeros((2,nb))
yerr2[1,:] = np.quantile(store_inp2,upperq,axis=0) - y2
yerr2[0,:] = y2 - np.quantile(store_inp2,lowerq,axis=0)

fig = plt.figure(layout='constrained',figsize=(8,8))
gs0 = fig.add_gridspec(2,2)
ax1 = fig.add_subplot(gs0[0])

#plt.rcParams['figure.figsize'] = [4,3]
if use_err == True:
    ax1.bar(x1,y1,color=swm_sk0,width=0.35,yerr=yerr1)
    ax1.bar(x2,y2,color=([swm_skn4,swm_sk0,swm_skmod,swm_skhi]),width=0.35,yerr=yerr2)
if use_err == False:
    ax1.bar(x1,y1,color=swm_sk0,width=0.35)
    ax1.bar(x2,y2,color=([swm_skn4,swm_sk0,swm_skmod,swm_skhi]),width=0.35)
ax1.set_title('Fixed FIRO pool',fontsize='large')
#ax1.set_xlabel('$CRPSS-SS_{test-swm}$',fontsize='large')
ax1.set_ylim(ylm)
ax1.set_xlim(xlm)
ax1.axhline(K_ORO,color=cb_gry,linestyle='--',linewidth=1)
#l1, = ax1.plot(0,0,color=cb_gry,linestyle='--',linewidth=1)
#l2, = ax1.plot(0,0,color=cb_grn,linestyle='--',linewidth=1)
#ax1.legend([l1,l2],['K','WCM Apr 1 limit'],loc='upper left',fontsize='small',frameon=False)
ax1.set_ylabel('Storage (TAF)',fontsize='large')
#ax1.axhline(ref_line,color=cb_grn,linestyle='--',linewidth=1)
#ax1.set_xticks(ticks=(x1+0.2))
#ax1.xaxis.set_ticklabels(['0','0.4','0.8'])
ax1.xaxis.set_ticklabels([])
ax1.tick_params(axis='both',which='major',labelsize='large')

#panel 2
ax1 = fig.add_subplot(gs0[1])

store_inp1 = swm_tstfpool_apr1
store_inp2 = swm_tstfpool_apr1_skilltrn

x1 = np.arange(nb)+0.8
x2 = np.arange(nb)+1.2

y1 = np.mean(store_inp1,axis=0)
y2 = np.mean(store_inp2,axis=0)

yerr1 = np.zeros((2,nb))
yerr1[1,:] = np.quantile(store_inp1,upperq,axis=0) - y1
yerr1[0,:] = y1 - np.quantile(store_inp1,lowerq,axis=0)

yerr2 = np.zeros((2,nb))
yerr2[1,:] = np.quantile(store_inp2,upperq,axis=0) - y2
yerr2[0,:] = y2 - np.quantile(store_inp2,lowerq,axis=0)

#plt.rcParams['figure.figsize'] = [4,3]
if use_err == True:
    ax1.bar(x1,y1,color=swm_sk0,width=0.35,yerr=yerr1)
    ax1.bar(x2,y2,color=([swm_skn4,swm_sk0,swm_skmod,swm_skhi]),width=0.35,yerr=yerr2)
if use_err == False:
    ax1.bar(x1,y1,color=swm_sk0,width=0.35)
    ax1.bar(x2,y2,color=([swm_skn4,swm_sk0,swm_skmod,swm_skhi]),width=0.35)
ax1.set_title('Flex FIRO pool',fontsize='large')
#ax1.set_xlabel('$CRPSS-SS_{test-swm}$',fontsize='large')
ax1.set_ylim(ylm)
ax1.set_xlim(xlm)
ax1.axhline(K_ORO,color=cb_gry,linestyle='--',linewidth=1)
#ax1.set_ylabel('Storage (TAF)',fontsize='large')
#ax1.axhline(ref_line,color=cb_grn,linestyle='--',linewidth=1)
ax1.yaxis.set_ticklabels([])
ax1.xaxis.set_ticklabels([])
#ax1.set_xticks(ticks=(x1+0.2))
#ax1.xaxis.set_ticklabels(['0','0.4','0.8'])
ax1.tick_params(axis='both',which='major',labelsize='large')

l1, = ax1.bar(x2[0],y2[0],color=swm_skn4,width=0.35)
l2, = ax1.bar(x1[1],y1[1],color=swm_sk0,width=0.35)
l3, = ax1.bar(x2[2],y2[2],color=swm_skmod,width=0.35)
l4, = ax1.bar(x2[3],y2[3],color=swm_skhi,width=0.35)
ax1.legend([l1,l2,l3,l4],['$%s_{swm}$' %(skillmods[0]),'$%s_{swm}$' %(skillmods[1]),'$%s_{swm}$' %(skillmods[2]),'$%s_{swm}$' %(skillmods[3])],loc='center left',fontsize='large',bbox_to_anchor=(1.0, 0.5),title='$CRPSS-SS_{train}$',title_fontsize='large',frameon=False)

#panel 3
ax1 = fig.add_subplot(gs0[2])

store_inp1 = swmcc_tst_apr1
store_inp2 = swmcc_tst_apr1_skilltrn

x1 = np.arange(nb)+0.8
x2 = np.arange(nb)+1.2

y1 = np.mean(store_inp1,axis=0)
y2 = np.mean(store_inp2,axis=0)

yerr1 = np.zeros((2,nb))
yerr1[1,:] = np.quantile(store_inp1,upperq,axis=0) - y1
yerr1[0,:] = y1 - np.quantile(store_inp1,lowerq,axis=0)

yerr2 = np.zeros((2,nb))
yerr2[1,:] = np.quantile(store_inp2,upperq,axis=0) - y2
yerr2[0,:] = y2 - np.quantile(store_inp2,lowerq,axis=0)

#plt.rcParams['figure.figsize'] = [4,3]
if use_err == True:
    ax1.bar(x1,y1,color=swm_sk0,width=0.35,yerr=yerr1)
    ax1.bar(x2,y2,color=([swmcc_skn4,swmcc_sk0,swmcc_skmod,swmcc_skhi]),width=0.35,yerr=yerr2)
if use_err == False:
    ax1.bar(x1,y1,color=swm_sk0,width=0.35)
    ax1.bar(x2,y2,color=([swmcc_skn4,swmcc_sk0,swmcc_skmod,swmcc_skhi]),width=0.35)
#ax1.set_title('Flex FIRO pool',fontsize='large')
ax1.set_xlabel('$CRPSS-SS_{test}$',fontsize='large')
ax1.set_ylim(ylm)
ax1.set_xlim(xlm)
ax1.axhline(K_ORO,color=cb_gry,linestyle='--',linewidth=1)
ax1.set_ylabel('Storage (TAF)',fontsize='large')
#ax1.axhline(ref_line,color=cb_grn,linestyle='--',linewidth=1)
#ax1.yaxis.set_ticklabels([])
#ax1.xaxis.set_ticklabels([])
ax1.set_xticks(ticks=(x1+0.2))
ax1.xaxis.set_ticklabels(['-0.4','0','0.4','0.8'])
ax1.tick_params(axis='both',which='major',labelsize='large')

#panel 4
ax1 = fig.add_subplot(gs0[3])

store_inp1 = swmcc_tstfpool_apr1
store_inp2 = swmcc_tstfpool_apr1_skilltrn

x1 = np.arange(nb)+0.8
x2 = np.arange(nb)+1.2

y1 = np.mean(store_inp1,axis=0)
y2 = np.mean(store_inp2,axis=0)

yerr1 = np.zeros((2,nb))
yerr1[1,:] = np.quantile(store_inp1,upperq,axis=0) - y1
yerr1[0,:] = y1 - np.quantile(store_inp1,lowerq,axis=0)

yerr2 = np.zeros((2,nb))
yerr2[1,:] = np.quantile(store_inp2,upperq,axis=0) - y2
yerr2[0,:] = y2 - np.quantile(store_inp2,lowerq,axis=0)

#plt.rcParams['figure.figsize'] = [4,3]
if use_err == True:
    ax1.bar(x1,y1,color=swm_sk0,width=0.35,yerr=yerr1)
    ax1.bar(x2,y2,color=([swmcc_skn4,swmcc_sk0,swmcc_skmod,swmcc_skhi]),width=0.35,yerr=yerr2)
if use_err == False:
    ax1.bar(x1,y1,color=swm_sk0,width=0.35)
    ax1.bar(x2,y2,color=([swmcc_skn4,swmcc_sk0,swmcc_skmod,swmcc_skhi]),width=0.35)
#ax1.set_title('Flex FIRO pool',fontsize='large')
ax1.set_xlabel('$CRPSS-SS_{test}$',fontsize='large')
ax1.set_ylim(ylm)
ax1.set_xlim(xlm)
ax1.axhline(K_ORO,color=cb_gry,linestyle='--',linewidth=1)
#ax1.set_ylabel('Storage (TAF)',fontsize='large')
#ax1.axhline(ref_line,color=cb_grn,linestyle='--',linewidth=1)
ax1.yaxis.set_ticklabels([])
#ax1.xaxis.set_ticklabels([])
ax1.set_xticks(ticks=(x1+0.2))
ax1.xaxis.set_ticklabels(['0.4','0','0.4','0.8'])
ax1.tick_params(axis='both',which='major',labelsize='large')

l1, = ax1.bar(x1[0],y1[0],color=swm_sk0,width=0.35)
l2, = ax1.bar(x2[0],y2[0],color=swmcc_skn4,width=0.35)
l3, = ax1.bar(x2[1],y2[1],color=swmcc_sk0,width=0.35)
l4, = ax1.bar(x2[2],y2[2],color=swmcc_skmod,width=0.35)
l5, = ax1.bar(x2[3],y2[3],color=swmcc_skhi,width=0.35)
ax1.legend([l1,l2,l3,l4,l5],['$%s_{swm}$' %(skillmods[1]),'$%s_{swm4c}$' %(skillmods[0]),'$%s_{swm4c}$' %(skillmods[1]),'$%s_{swm4c}$' %(skillmods[2]),'$%s_{swm4c}$' %(skillmods[3])],loc='center left',fontsize='large',bbox_to_anchor=(1.0, 0.5),title='$CRPSS-SS_{train}$',title_fontsize='large',frameon=False)

syn_util.fig_title(fig,'April 1 Storage',loc=(0.45,1.),fontsize='xx-large',fontweight='bold',rotation=0,ha='center',va='center')
syn_util.fig_title(fig,'Historical',loc=(-0.02,0.725),fontsize='xx-large',fontweight='bold',rotation=90,ha='center',va='center')
syn_util.fig_title(fig,'4C Warming',loc=(-0.02,0.275),fontsize='xx-large',fontweight='bold',rotation=90,ha='center',va='center')

fig.savefig('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ops_plots/2x2_aggregate-storage-comparison_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_same-swm=%s_nevts=%s_cutoff=%s.png' %(use_firo_top,use_firo_bottom,seed,skillmods,skill_dcy,skill_tail,swm_ind,same_swm,n_evts,yr_cutoff),dpi=300,bbox_inches='tight')

#.............................................................................................................
#Plot 2: Spill volume metric w/ and w/o FIRO pool flexibility and w/ and w/o CC
#Event based metrics
#panel 1
swm_ind = 'swm'
nb = len(skillmods)

ylm1 = ([0,450])
ylm2 = ([0,7000])
xlm = ([0.5,nb+0.5])
upperq = 0.9
lowerq = 0.1
use_err = False

#store_inp1 = swm_tst_metrics[13,:,:]
#store_inp2 = swm_tst_metrics_skilltrn[13,:,:] 

spill_inp1 = swm_tst_metrics[12,:,:]
spill_inp2 = swm_tst_metrics_skilltrn[12,:,:] 

x1 = np.arange(nb)+0.8
x2 = np.arange(nb)+1.2

y1 = np.mean(spill_inp1,axis=0)
y2 = np.mean(spill_inp2,axis=0)

fig = plt.figure(layout='constrained',figsize=(8,8))
gs0 = fig.add_gridspec(2,2)
ax1 = fig.add_subplot(gs0[0])

#plt.rcParams['figure.figsize'] = [4,3]
ax1.bar(x1,y1,color=swm_sk0,width=0.35)
ax1.bar(x2,y2,color=([swm_skn4,swm_sk0,swm_skmod,swm_skhi]),width=0.35)
ax1.set_title('Fixed FIRO pool',fontsize='large')
#ax1.set_xlabel('$CRPSS-SS_{test-swm}$',fontsize='large')
ax1.set_ylim(ylm1)
ax1.set_xlim(xlm)
ax1.set_ylabel('Spill Volume (TAF)',fontsize='large')
ax1.axhline(ref_line,color=cb_grn,linestyle='--',linewidth=1)
#ax1.set_xticks(ticks=(x1+0.2))
#ax1.xaxis.set_ticklabels(['0','0.4','0.8'])
ax1.xaxis.set_ticklabels([])
ax1.tick_params(axis='both',which='major',labelsize='large')

#panel 2
ax1 = fig.add_subplot(gs0[1])

spill_inp1 = swm_tstfpool_metrics[12,:,:]
spill_inp2 = swm_tstfpool_metrics_skilltrn[12,:,:]

y1 = np.mean(spill_inp1,axis=0)
y2 = np.mean(spill_inp2,axis=0)

#plt.rcParams['figure.figsize'] = [4,3]
ax1.bar(x1,y1,color=swm_sk0,width=0.35)
ax1.bar(x2,y2,color=([swm_skn4,swm_sk0,swm_skmod,swm_skhi]),width=0.35)
ax1.set_title('Flex FIRO pool',fontsize='large')
#ax1.set_xlabel('$CRPSS-SS_{test-swm}$',fontsize='large')
ax1.set_ylim(ylm1)
ax1.set_xlim(xlm)
#ax1.set_ylabel('Storage (TAF)',fontsize='large')
ax1.axhline(ref_line,color=cb_grn,linestyle='--',linewidth=1)
ax1.yaxis.set_ticklabels([])
ax1.xaxis.set_ticklabels([])
#ax1.set_xticks(ticks=(x1+0.2))
#ax1.xaxis.set_ticklabels(['0','0.4','0.8'])
ax1.tick_params(axis='both',which='major',labelsize='large')
l1, = ax1.bar(x2[0],y2[0],color=swm_skn4,width=0.35)
l2, = ax1.bar(x1[0],y1[0],color=swm_sk0,width=0.35)
l3, = ax1.bar(x2[1],y2[1],color=swm_skmod,width=0.35)
l4, = ax1.bar(x2[2],y2[2],color=swm_skhi,width=0.35)
ax1.legend([l1,l2,l3,l4],['$%s_{swm}$' %(skillmods[0]),'$%s_{swm}$' %(skillmods[1]),'$%s_{swm}$' %(skillmods[2]),'$%s_{swm}$' %(skillmods[3])],loc='center left',fontsize='large',bbox_to_anchor=(1.0, 0.5),title='$CRPSS-SS_{train}$',title_fontsize='large',frameon=False)

#panel 3
ax1 = fig.add_subplot(gs0[2])

spill_inp1 = swmcc_tst_metrics[12,:,:]
spill_inp2 = swmcc_tst_metrics_skilltrn[12,:,:]

x1 = np.arange(nb)+0.8
x2 = np.arange(nb)+1.2

y1 = np.mean(spill_inp1,axis=0)
y2 = np.mean(spill_inp2,axis=0)

#plt.rcParams['figure.figsize'] = [4,3]
ax1.bar(x1,y1,color=swm_sk0,width=0.35)
ax1.bar(x2,y2,color=([swmcc_skn4,swmcc_sk0,swmcc_skmod,swmcc_skhi]),width=0.35)
#ax1.set_title('Flex FIRO pool',fontsize='large')
ax1.set_xlabel('$CRPSS-SS_{test}$',fontsize='large')
ax1.set_ylim(ylm2)
ax1.set_xlim(xlm)
ax1.set_ylabel('Spill Volume (TAF)',fontsize='large')
#ax1.yaxis.set_ticklabels([])
#ax1.xaxis.set_ticklabels([])
ax1.set_xticks(ticks=(x1+0.2))
ax1.xaxis.set_ticklabels(['-0.4','0','0.4','0.8'])
ax1.tick_params(axis='both',which='major',labelsize='large')

#panel 4
ax1 = fig.add_subplot(gs0[3])

spill_inp1 = swmcc_tstfpool_metrics[12,:,:]
spill_inp2 = swmcc_tstfpool_metrics_skilltrn[12,:,:]

x1 = np.arange(nb)+0.8
x2 = np.arange(nb)+1.2

y1 = np.mean(spill_inp1,axis=0)
y2 = np.mean(spill_inp2,axis=0)

#plt.rcParams['figure.figsize'] = [4,3]
ax1.bar(x1,y1,color=swm_sk0,width=0.35)
ax1.bar(x2,y2,color=([swmcc_skn4,swmcc_sk0,swmcc_skmod,swmcc_skhi]),width=0.35)
#ax1.set_title('Flex FIRO pool',fontsize='large')
ax1.set_xlabel('$CRPSS-SS_{test}$',fontsize='large')
ax1.set_ylim(ylm2)
ax1.set_xlim(xlm)
#ax1.set_ylabel('Storage (TAF)',fontsize='large')
#ax1.axhline(ref_line,color=cb_grn,linestyle='--',linewidth=1)
ax1.yaxis.set_ticklabels([])
#ax1.xaxis.set_ticklabels([])
ax1.set_xticks(ticks=(x1+0.2))
ax1.xaxis.set_ticklabels(['-0.4','0','0.4','0.8'])
ax1.tick_params(axis='both',which='major',labelsize='large')
l1, = ax1.bar(x1[0],y1[0],color=swm_sk0,width=0.35)
l2, = ax1.bar(x2[0],y2[0],color=swmcc_skn4,width=0.35)
l3, = ax1.bar(x2[1],y2[1],color=swmcc_sk0,width=0.35)
l4, = ax1.bar(x2[2],y2[2],color=swmcc_skmod,width=0.35)
l5, = ax1.bar(x2[3],y2[3],color=swmcc_skhi,width=0.35)
ax1.legend([l1,l2,l3,l4,l5],['$%s_{swm}$' %(skillmods[1]),'$%s_{swm4c}$' %(skillmods[0]),'$%s_{swm4c}$' %(skillmods[1]),'$%s_{swm4c}$' %(skillmods[2]),'$%s_{swm4c}$' %(skillmods[2])],loc='center left',fontsize='large',bbox_to_anchor=(1.0, 0.5),title='$CRPSS-SS_{train}$',title_fontsize='large',frameon=False)
syn_util.fig_title(fig,'Spill Volume',loc=(0.45,1.),fontsize='xx-large',fontweight='bold',rotation=0,ha='center',va='center')
syn_util.fig_title(fig,'Historical',loc=(-0.02,0.725),fontsize='xx-large',fontweight='bold',rotation=90,ha='center',va='center')
syn_util.fig_title(fig,'4C Warming',loc=(-0.02,0.275),fontsize='xx-large',fontweight='bold',rotation=90,ha='center',va='center')

fig.savefig('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ops_plots/2x2_aggregate-spillvol-comparison_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_same-swm=%s_nevts=%s_cutoff=%s.png' %(use_firo_top,use_firo_bottom,seed,skillmods,skill_dcy,skill_tail,swm_ind,same_swm,n_evts,yr_cutoff),dpi=300,bbox_inches='tight')

#.............................................................................................................
#Plot 3: # of spills metric w/ and w/o FIRO pool flexibility and w/ and w/o CC
#Event based metrics
#panel 1
swm_ind = 'swm'
nb = len(skillmods)

ylm1 = ([0,6])
ylm2 = ([0,60])
xlm = ([0.5,nb+0.5])
upperq = 0.9
lowerq = 0.1
use_err = False

#store_inp1 = swm_tst_metrics[13,:,:]
#store_inp2 = swm_tst_metrics_skilltrn[13,:,:] 

spill_inp1 = swm_tst_metrics[11,:,:]
spill_inp2 = swm_tst_metrics_skilltrn[11,:,:] 

x1 = np.arange(nb)+0.8
x2 = np.arange(nb)+1.2

y1 = np.mean(spill_inp1,axis=0)
y2 = np.mean(spill_inp2,axis=0)

fig = plt.figure(layout='constrained',figsize=(8,8))
gs0 = fig.add_gridspec(2,2)
ax1 = fig.add_subplot(gs0[0])

#plt.rcParams['figure.figsize'] = [4,3]
ax1.bar(x1,y1,color=swm_sk0,width=0.35)
ax1.bar(x2,y2,color=([swm_skn4,swm_sk0,swm_skmod,swm_skhi]),width=0.35)
ax1.set_title('Fixed FIRO pool',fontsize='large')
#ax1.set_xlabel('$CRPSS-SS_{test-swm}$',fontsize='large')
ax1.set_ylim(ylm1)
ax1.set_xlim(xlm)
ax1.set_ylabel('No. of Spills',fontsize='large')
ax1.axhline(ref_line,color=cb_grn,linestyle='--',linewidth=1)
#ax1.set_xticks(ticks=(x1+0.2))
#ax1.xaxis.set_ticklabels(['0','0.4','0.8'])
ax1.xaxis.set_ticklabels([])
ax1.tick_params(axis='both',which='major',labelsize='large')

#panel 2
ax1 = fig.add_subplot(gs0[1])

spill_inp1 = swm_tstfpool_metrics[11,:,:]
spill_inp2 = swm_tstfpool_metrics_skilltrn[11,:,:]

y1 = np.mean(spill_inp1,axis=0)
y2 = np.mean(spill_inp2,axis=0)

#plt.rcParams['figure.figsize'] = [4,3]
ax1.bar(x1,y1,color=swm_sk0,width=0.35)
ax1.bar(x2,y2,color=([swm_skn4,swm_sk0,swm_skmod,swm_skhi]),width=0.35)
ax1.set_title('Flex FIRO pool',fontsize='large')
#ax1.set_xlabel('$CRPSS-SS_{test-swm}$',fontsize='large')
ax1.set_ylim(ylm1)
ax1.set_xlim(xlm)
#ax1.set_ylabel('Storage (TAF)',fontsize='large')
ax1.axhline(ref_line,color=cb_grn,linestyle='--',linewidth=1)
ax1.yaxis.set_ticklabels([])
ax1.xaxis.set_ticklabels([])
#ax1.set_xticks(ticks=(x1+0.2))
#ax1.xaxis.set_ticklabels(['0','0.4','0.8'])
ax1.tick_params(axis='both',which='major',labelsize='large')

l1, = ax1.bar(x2[0],y2[0],color=swm_skn4,width=0.35)
l2, = ax1.bar(x1[0],y1[0],color=swm_sk0,width=0.35)
l3, = ax1.bar(x2[1],y2[1],color=swm_skmod,width=0.35)
l4, = ax1.bar(x2[2],y2[2],color=swm_skhi,width=0.35)
ax1.legend([l1,l2,l3,l4],['$%s_{swm}$' %(skillmods[0]),'$%s_{swm}$' %(skillmods[1]),'$%s_{swm}$' %(skillmods[2]),'$%s_{swm}$' %(skillmods[3])],loc='center left',fontsize='large',bbox_to_anchor=(1.0, 0.5),title='$CRPSS-SS_{train}$',title_fontsize='large',frameon=False)

#panel 3
ax1 = fig.add_subplot(gs0[2])

spill_inp1 = swmcc_tst_metrics[11,:,:]
spill_inp2 = swmcc_tst_metrics_skilltrn[11,:,:]

x1 = np.arange(nb)+0.8
x2 = np.arange(nb)+1.2

y1 = np.mean(spill_inp1,axis=0)
y2 = np.mean(spill_inp2,axis=0)

#plt.rcParams['figure.figsize'] = [4,3]
ax1.bar(x1,y1,color=swm_sk0,width=0.35)
ax1.bar(x2,y2,color=([swmcc_skn4,swmcc_sk0,swmcc_skmod,swmcc_skhi]),width=0.35)
#ax1.set_title('Flex FIRO pool',fontsize='large')
ax1.set_xlabel('$CRPSS-SS_{test}$',fontsize='large')
ax1.set_ylim(ylm2)
ax1.set_xlim(xlm)
ax1.set_ylabel('No. of Spills',fontsize='large')
#ax1.yaxis.set_ticklabels([])
#ax1.xaxis.set_ticklabels([])
ax1.set_xticks(ticks=(x1+0.2))
ax1.xaxis.set_ticklabels(['-0.4','0','0.4','0.8'])
ax1.tick_params(axis='both',which='major',labelsize='large')

#panel 4
ax1 = fig.add_subplot(gs0[3])

spill_inp1 = swmcc_tstfpool_metrics[11,:,:]
spill_inp2 = swmcc_tstfpool_metrics_skilltrn[11,:,:]

x1 = np.arange(nb)+0.8
x2 = np.arange(nb)+1.2

y1 = np.mean(spill_inp1,axis=0)
y2 = np.mean(spill_inp2,axis=0)

#plt.rcParams['figure.figsize'] = [4,3]
ax1.bar(x1,y1,color=swm_sk0,width=0.35)
ax1.bar(x2,y2,color=([swmcc_skn4,swmcc_sk0,swmcc_skmod,swmcc_skhi]),width=0.35)
#ax1.set_title('Flex FIRO pool',fontsize='large')
ax1.set_xlabel('$CRPSS-SS_{test}$',fontsize='large')
ax1.set_ylim(ylm2)
ax1.set_xlim(xlm)
#ax1.set_ylabel('Storage (TAF)',fontsize='large')
#ax1.axhline(ref_line,color=cb_grn,linestyle='--',linewidth=1)
ax1.yaxis.set_ticklabels([])
#ax1.xaxis.set_ticklabels([])
ax1.set_xticks(ticks=(x1+0.2))
ax1.xaxis.set_ticklabels(['-0.4','0','0.4','0.8'])
ax1.tick_params(axis='both',which='major',labelsize='large')

l1, = ax1.bar(x1[0],y1[0],color=swm_sk0,width=0.35)
l2, = ax1.bar(x2[0],y2[0],color=swmcc_skn4,width=0.35)
l3, = ax1.bar(x2[1],y2[1],color=swmcc_sk0,width=0.35)
l4, = ax1.bar(x2[2],y2[2],color=swmcc_skmod,width=0.35)
l5, = ax1.bar(x2[3],y2[3],color=swmcc_skhi,width=0.35)
ax1.legend([l1,l2,l3,l4,l5],['$%s_{swm}$' %(skillmods[1]),'$%s_{swm4c}$' %(skillmods[0]),'$%s_{swm4c}$' %(skillmods[1]),'$%s_{swm4c}$' %(skillmods[2]),'$%s_{swm4c}$' %(skillmods[2])],loc='center left',fontsize='large',bbox_to_anchor=(1.0, 0.5),title='$CRPSS-SS_{train}$',title_fontsize='large',frameon=False)

syn_util.fig_title(fig,'No. of Spills',loc=(0.45,1.),fontsize='xx-large',fontweight='bold',rotation=0,ha='center',va='center')
syn_util.fig_title(fig,'Historical',loc=(-0.02,0.725),fontsize='xx-large',fontweight='bold',rotation=90,ha='center',va='center')
syn_util.fig_title(fig,'4C Warming',loc=(-0.02,0.275),fontsize='xx-large',fontweight='bold',rotation=90,ha='center',va='center')

fig.savefig('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ops_plots/2x2_aggregate-spillno-comparison_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mods=%s_dcy=%s_tail=%s_swm=%s_same-swm=%s_nevts=%s_cutoff=%s.png' %(use_firo_top,use_firo_bottom,seed,skillmods,skill_dcy,skill_tail,swm_ind,same_swm,n_evts,yr_cutoff),dpi=300,bbox_inches='tight')

#--------------------------------------end--------------------------------------
