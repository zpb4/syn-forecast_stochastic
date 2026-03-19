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
import ensemble_verification_functions as verify
import pickle
import datetime as dt

col_cb = sns.color_palette('colorblind')
sns.palplot(col_cb)  #plot coloblind palette for comparison
colv0 = sns.color_palette('Blues',10)
colv000 = sns.color_palette('Greens',10)
colv00 = sns.color_palette('Purples',10)
colv1 = sns.color_palette('PuRd',10)
colv2 = sns.color_palette('YlOrBr',10)
sns.palplot(colv0)
sns.palplot(colv00)
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
"""
swmcc_sk0_skn4 = colv000[3]
swmcc_sk0_sk0 = colv00[3]
swmcc_sk0_skmod = colv00[6]
swmcc_sk0_skhi = colv00[9]
"""
swmcc_sk0_skn4 = colv000[3]
swmcc_sk0_sk0 = colv0[3]
swmcc_sk0_skmod = colv0[6]
swmcc_sk0_skhi = colv0[9]

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
swm_evt_no = 1
swmcc_evt_no = 1
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


#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#plots
swm_ind = 'swm'
#1. 3 x 2 plot of data subsets as defined below
sns.set_theme()
sns.set_style('ticks')
sns.set_context('paper')

R99_swmn4tst = np.sort(swm_tst_arr[0,1,:])[int(0.99*len(swm_tst_arr[0,1,:]))]
R99_swm0tst = np.sort(swm_tst_arr[1,1,:])[int(0.99*len(swm_tst_arr[0,1,:]))]
R99_swm4tst = np.sort(swm_tst_arr[2,1,:])[int(0.99*len(swm_tst_arr[0,1,:]))]
R99_swm8tst = np.sort(swm_tst_arr[3,1,:])[int(0.99*len(swm_tst_arr[0,1,:]))]

R99_swmn4tst_skilltrn = np.sort(swm_tst_arr_skilltrn[0,1,:])[int(0.99*len(swm_tst_arr[0,1,:]))]
R99_swm0tst_skilltrn = np.sort(swm_tst_arr_skilltrn[1,1,:])[int(0.99*len(swm_tst_arr[0,1,:]))]
R99_swm4tst_skilltrn = np.sort(swm_tst_arr_skilltrn[2,1,:])[int(0.99*len(swm_tst_arr[0,1,:]))]
R99_swm8tst_skilltrn = np.sort(swm_tst_arr_skilltrn[3,1,:])[int(0.99*len(swm_tst_arr[0,1,:]))]


firo_top = swm_tst_arr[0,2,:] 
firo_bot = swm_tst_arr[0,7,:] 

#f, ((a1,a2,a7,a8),(a3,a4,a9,a10),(a5,a6,a11,a12)) = plt.subplots(3,4,figsize=(12,7),sharex=True,layout='constrained')
f = plt.figure(layout='constrained',figsize=(12,8))
#gs0 = fig.add_gridspec(2,2,width_ratios=[1,2])
#ax1 = fig.add_subplot(gs0[0])
gs0 = f.add_gridspec(1,2,width_ratios=(1.25,1))
gs1 = gs0[0].subgridspec(2,1)
gs2 = gs0[1].subgridspec(3,1)
#gs3 = gs0[3].subgridspec(3,2)

a0 = f.add_subplot(gs1[0])

#f, (a0,a1) = plt.subplots(1,2,width_ratios=[1,2])
#plt.subplot(1,2,1)
a0.plot(np.arange(max_lds)+1,rcurve_arr_swm[1,:],linewidth=2,color=swm_sk0)
#a0.plot(np.arange(max_lds)+1,rcurve_arr_swm[1,:],linewidth=2,color=swm_skmod)
#a0.plot(np.arange(max_lds)+1,rcurve_arr_swm[2,:],linewidth=2,color=swm_skhi)
a0.set_xlabel('Lead time (days)',fontsize='large')
a0.set_ylabel('Risk threshold',fontsize='large')
#a0.set_title('Risk curve',fontsize='x-large',fontweight='bold')
a0.set_ylim([0,1.01])
a0.set_xlim([1,max_lds])
#a0.text(6,0.1,'FIRO pool: '+str(round(pars['firo_pool'],3)) + ' (%s TAF)' %(round((pars['firo_pool']+0.5)*K_scale)),fontsize='x-large')
a0.text(1.5,.9,'a)',fontsize='xx-large',fontweight='bold')
a0.tick_params(axis='both',which='major',labelsize='large')
a0.legend(['$CRPSS-SS_{train}=%s_{swm}$' %(skillmods[1])],loc='lower right',fontsize='large',frameon=False)

evt_idx = evt_idx_swm_tst[swm_evt_no-1]
#evt_idx = np.where(swm_tst_arr[0,3,:] == np.max(swm_tst_arr[0,3,:]))[0][0]
#max_spill = np.argmax(swm_tst_arr[0,3,np.where(swm_tst_arr[0,3,:]>0)[0]])
#evt_idx = np.where(swm_tst_arr[0,3,:]>0)[0][max_spill]
#evt_idx = np.where(swm_tst_arr[0,3,:]>0)[0][3]

#max_idx = int(np.where(Q_swm_tst == np.sort(Q_swm_tst)[-evt_no])[0])
evt_date = df_idx.values[evt_idx].strftime("%Y-%m-%d")
evt_mo = evt_date[5:7]
if np.int64(evt_mo) >=8 and np.int64(evt_mo) <=12:
    evt_ref = '1984'+evt_date[4:]
else:
    evt_ref = '1985'+evt_date[4:]
ref_yr_st = '1984-08-01'
ref_yr_ed = '1985-06-30'
st_idx = len(pd.date_range(ref_yr_st,evt_ref,freq='D'))
ed_idx = len(pd.date_range(evt_ref,ref_yr_ed,freq='D'))
#st_sl_date = df_idx.values[evt_idx-365].strftime("%Y-%m-%d")
#ed_sl_date = df_idx.values[evt_idx+365].strftime("%Y-%m-%d")

dt_idx = np.arange(evt_idx-st_idx,evt_idx+ed_idx)

spill_swmn4tst = np.where(swm_tst_arr[0,3,dt_idx]>0)
spill_swm0tst = np.where(swm_tst_arr[1,3,dt_idx]>0)
spill_swm4tst = np.where(swm_tst_arr[2,3,dt_idx]>0)
spill_swm8tst = np.where(swm_tst_arr[3,3,dt_idx]>0)

spills_swmn4tst = np.count_nonzero(swm_tst_arr[0,3,:])
spills_swm0tst = np.count_nonzero(swm_tst_arr[1,3,:])
spills_swm4tst = np.count_nonzero(swm_tst_arr[2,3,:])
spills_swm8tst = np.count_nonzero(swm_tst_arr[3,3,:])

spills_swmn4tst_skilltrn = np.count_nonzero(swm_tst_arr_skilltrn[0,3,:])
spills_swm0tst_skilltrn = np.count_nonzero(swm_tst_arr_skilltrn[1,3,:])
spills_swm4tst_skilltrn = np.count_nonzero(swm_tst_arr_skilltrn[2,3,:])
spills_swm8tst_skilltrn = np.count_nonzero(swm_tst_arr_skilltrn[3,3,:])

a00 = f.add_subplot(gs1[1])

dt_format=mdates.DateFormatter('%m/%d')
a00.xaxis.set_major_formatter(dt_format)

#1b. plot baseline ops vs hefs 
l1, = a00.plot(df_idx[dt_idx], tocs[dt_idx]/1000, c=cb_gry,alpha=0.5,linewidth=1)
l2, = a00.plot(df_idx[dt_idx], firo_top[dt_idx]/1000, c = cb_brwn,alpha=0.5,linewidth=1)
l7, = a00.plot(df_idx[dt_idx], firo_bot[dt_idx]/1000, c = cb_brwn,alpha=0.5,linewidth=1)
l3, = a00.plot(df_idx[dt_idx], swm_tst_arr[0,0,dt_idx]/1000, c=swm_sk0_skn4,alpha=0.6,linewidth=1.5)
l4, = a00.plot(df_idx[dt_idx], swm_tst_arr[1,0,dt_idx]/1000, c=swm_sk0_sk0,alpha=0.6,linewidth=1.5)
l5, = a00.plot(df_idx[dt_idx], swm_tst_arr[2,0,dt_idx]/1000, c=swm_sk0_skmod,alpha=0.6,linewidth=1.5)
l6, = a00.plot(df_idx[dt_idx], swm_tst_arr[3,0,dt_idx]/1000, c=swm_sk0_skhi,alpha=0.6,linewidth=1.5)
a00.legend([l3,l4,l5,l6,l1,l2],['$S_{%s}$' %(skillmods[0]),'$S_{%s}$' %(skillmods[1]),'$S_{%s}$' %(skillmods[2]),'$S_{%s}$' %(skillmods[3]),'TOCS','FIRO'],loc='lower left',fontsize='large',frameon=False)
a00.set_ylabel('Storage (MAF)',fontsize='large')
a00.set_xlabel('Date',fontsize='large')
a00.set_ylim([0.65*K/1000, K*1.1/1000])
a00.axhline(3.524,linewidth=1,color='gray',linestyle='--',alpha=0.5)
a00.text(df_idx[dt_idx[180]-150],K*1.07/1000,r'$\overline{S}:$'+str(round(np.mean(swm_tst_arr[0,0,:])/1000,3)) +' MAF',color=swm_sk0_skn4,fontsize='large')
a00.text(df_idx[dt_idx[180]-150],K*1.05/1000,r'$\overline{S}:$'+str(round(np.mean(swm_tst_arr[1,0,:])/1000,3)) +' MAF',color=swm_sk0_sk0,fontsize='large')
a00.text(df_idx[dt_idx[180]-150],K*1.03/1000,r'$\overline{S}:$'+str(round(np.mean(swm_tst_arr[2,0,:])/1000,3)) +' MAF',color=swm_sk0_skmod,fontsize='large')
a00.text(df_idx[dt_idx[180]-150],K*1.01/1000,r'$\overline{S}:$'+str(round(np.mean(swm_tst_arr[3,0,:])/1000,3)) +' MAF',color=swm_sk0_skhi,fontsize='large')
a00.text(df_idx[dt_idx[180]-30],K*1.07/1000,r'$R_{99}:$'+str(round(R99_swmn4tst,1))+' TAF/d',color=swm_sk0_skn4,fontsize='large')
a00.text(df_idx[dt_idx[180]-30],K*1.05/1000,r'$R_{99}:$'+str(round(R99_swm0tst,1))+' TAF/d',color=swm_sk0_sk0,fontsize='large')
a00.text(df_idx[dt_idx[180]-30],K*1.03/1000,r'$R_{99}:$'+str(round(R99_swm4tst,1))+' TAF/d',color=swm_sk0_skmod,fontsize='large')
a00.text(df_idx[dt_idx[180]-30],K*1.01/1000,r'$R_{99}:$'+str(round(R99_swm8tst,1))+' TAF/d',color=swm_sk0_skhi,fontsize='large')
a00.text(df_idx[dt_idx[180]+90],K*1.07/1000,r'$Spills:$'+str(round(spills_swmn4tst)),color=swm_sk0_skn4,fontsize='large')
a00.text(df_idx[dt_idx[180]+90],K*1.05/1000,r'$Spills:$'+str(round(spills_swm0tst)),color=swm_sk0_sk0,fontsize='large')
a00.text(df_idx[dt_idx[180]+90],K*1.03/1000,r'$Spills:$'+str(round(spills_swm4tst)),color=swm_sk0_skmod,fontsize='large')
a00.text(df_idx[dt_idx[180]+90],K*1.01/1000,r'$Spills:$'+str(round(spills_swm8tst)),color=swm_sk0_skhi,fontsize='large')
a00.tick_params(axis='both',which='major',labelsize='large')
a00.text(df_idx[dt_idx[180]-185],K*1.05/1000,'b)',fontsize='xx-large',fontweight='bold')

evt = df_idx[evt_idx].strftime('%Y-%m-%d')
pad = 15
#evt_idx=df_idx.get_loc(evt)
st,ed = df_idx[evt_idx-pad],df_idx[evt_idx+pad]
start,end = str(st)[0:10],str(ed)[0:10]
dt_idx = np.arange(pad*2+1) + (evt_idx-pad)
dt_format=mdates.DateFormatter('%m/%d')

a1 = f.add_subplot(gs2[0])

ymn = 2300 #0.95*min(np.min(sim_data1[:,dt_idx,0]),np.min(sim_data2[:,dt_idx,0]))
ymx = 1.2* (K-ymn)+ymn

l1, = a1.plot(df_idx[dt_idx], firo_top[dt_idx]/1000, c = cb_brwn, linewidth=1,alpha=0.5)
a1.plot(df_idx[dt_idx], firo_bot[dt_idx]/1000, c = cb_brwn, linewidth=1,alpha=0.5)
l4, = a1.plot(df_idx[dt_idx], tocs[dt_idx]/1000, c = 'gray', linewidth=1,alpha=0.5)
a1.axhline(3.524,linewidth=1,color='gray',linestyle='--',alpha=0.5)
leg1 = a1.legend([l1,l4],['FIRO','TOCS'],bbox_transform=a1.transAxes,loc=(.75,.6),fontsize='large',frameon=False)
a1.add_artist(leg1)
l1, = a1.plot(df_idx[dt_idx], swm_tst_arr[0,0,dt_idx]/1000, c=swm_sk0_skn4,alpha=0.75)
l2, = a1.plot(df_idx[dt_idx], swm_tst_arr[1,0,dt_idx]/1000, c=swm_sk0,alpha=0.75)
l3, = a1.plot(df_idx[dt_idx], swm_tst_arr[2,0,dt_idx]/1000, c=swm_sk0_skmod,alpha=0.75)
l4, = a1.plot(df_idx[dt_idx], swm_tst_arr[3,0,dt_idx]/1000, c=swm_sk0_skhi,alpha=0.75)
a1.legend([l1,l2,l3,l4],['$S_{%s}$' %(skillmods[0]),'$S_{%s}$' %(skillmods[1]),'$S_{%s}$' %(skillmods[2]),'$S_{%s}$' %(skillmods[3])],bbox_transform=a1.transAxes,loc=(.01,.5),fontsize='large',frameon=False)
a1.axvline(df_idx[evt_idx],linewidth=0.5,color='black',linestyle='--')
a1.text(df_idx[evt_idx],(1.05*(K-ymn)+ymn)/1000,evt,fontsize='large')
a1.text(df_idx[evt_idx-pad+1],(.1*(ymx-ymn)+ymn)/1000,'c)',fontsize='xx-large',fontweight='bold')
a1.set_ylabel('Storage (MAF)',fontsize='large')
a1.set_ylim([ymn/1000, ymx/1000])
a1.set_xlim([st, ed])
a1.xaxis.set_major_formatter(dt_format)
#plt.gcf().autofmt_xdate()
a1.xaxis.set_ticklabels([])
a1.tick_params(axis='both',which='major',labelsize='large')

a3 = f.add_subplot(gs2[1])

l1, = a3.plot(df_idx[dt_idx], swm_tst_arr[0,1,dt_idx] / kcfs_to_tafd,c=swm_sk0_skn4,alpha=0.75)
l2, = a3.plot(df_idx[dt_idx], swm_tst_arr[1,1,dt_idx] / kcfs_to_tafd,c=swm_sk0_sk0,alpha=0.75)
l3, = a3.plot(df_idx[dt_idx], swm_tst_arr[2,1,dt_idx] / kcfs_to_tafd, c=swm_sk0_skmod,alpha=0.75)
l4, = a3.plot(df_idx[dt_idx], swm_tst_arr[3,1,dt_idx] / kcfs_to_tafd, c=swm_sk0_skhi,alpha=0.75)
l5, = a3.plot(df_idx[dt_idx], Q_swm_tst[dt_idx] / kcfs_to_tafd,c='black',linewidth=3,alpha=0.25)
a3.legend([l1,l2,l3,l4,l5],['$R_{%s}$' %(skillmods[0]),'$R_{%s}$' %(skillmods[1]),'$R_{%s}$' %(skillmods[2]),'$R_{%s}$' %(skillmods[3]),'$Q$'],bbox_transform=a3.transAxes,loc=(0.8,.2),fontsize='large',frameon=False)
a3.axhline(Rmax / kcfs_to_tafd, color='red')
#a3.axhline(R_thresh, color=chefs,linewidth=0.5,alpha=0.5)
a3.axvline(df_idx[evt_idx],linewidth=0.5,color='black',linestyle='--')
#a3.text(df_idx[evt_idx],1.15 * Rmax_scale / kcfs_to_tafd,evt,fontsize='small')
a3.text(df_idx[evt_idx-int(0.75*pad)],1.1*Rmax / kcfs_to_tafd,'$Peak = %s kcfs$' %(str(round(Q_swm_tst[evt_idx]/kcfs_to_tafd,1))),color='black',fontsize='large')
a3.text(df_idx[evt_idx+int(pad/2)],1.025*Rmax / kcfs_to_tafd,'$R_{max}$',color='red',fontsize='large')
a3.set_ylabel('Streamflow (kcfs)',fontsize='large')
a3.text(df_idx[evt_idx-pad+1],1.1 * Rmax / kcfs_to_tafd,'d)',fontsize='xx-large',fontweight='bold')
a3.set_ylim([0, 1.25 * Rmax / kcfs_to_tafd])
a3.set_xlim([st, ed])
a3.xaxis.set_major_formatter(dt_format)
#plt.gcf().autofmt_xdate()
a3.xaxis.set_ticklabels([])
a3.tick_params(axis='both',which='major',labelsize='large')

a5 = f.add_subplot(gs2[2])

l1, = a5.plot(df_idx[dt_idx], swm_tst_arr[0,4,dt_idx], c=swm_sk0_skn4,alpha=0.75)
l2, = a5.plot(df_idx[dt_idx], swm_tst_arr[1,4,dt_idx], c=swm_sk0_sk0,alpha=0.75)
l3, = a5.plot(df_idx[dt_idx], swm_tst_arr[2,4,dt_idx], c=swm_sk0_skmod,alpha=0.75)
l4, = a5.plot(df_idx[dt_idx], swm_tst_arr[3,4,dt_idx], c=swm_sk0_skhi,alpha=0.75)
#a5.legend([l1,l2,l3],['$R^{ld}_{CRPSS%s}$' %(skillmods[0]),'$R^{ld}_{CRPSS%s}$' %(skillmods[1]),'$R^{ld}_{CRPSS%s}$' %(skillmods[2])],loc='upper right',fontsize='large',frameon=False)
a5.axvline(df_idx[evt_idx],linewidth=0.5,color='black',linestyle='--')
#a5.text(df_idx[evt_idx-int(pad*0.7)],14,evt,fontsize='small')
a5.set_ylabel('Release lead time (days)',fontsize='large')
a5.text(df_idx[evt_idx-pad+1],13,'e)',fontsize='xx-large',fontweight='bold')
a5.set_ylim([0, 15.1])
a5.set_xlim([st, ed])
a5.tick_params(axis='both',which='major',labelsize='large')

a5.xaxis.set_major_formatter(dt_format)
plt.setp(plt.xticks()[1],rotation=45,ha='right')

f.savefig('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ops_plots/3x2_risk-curve-ts-S-Qcp-Rlead_noskilltrn_use-firo-top=%s_use-firo-bottom=%s_seed-%s_evt=%s_mods=%s_dcy=%s_tail=%s_swm=%s_same-swm=%s_cutoff=%s.png' %(use_firo_top,use_firo_bottom,seed,swm_evt_no,str(skillmods),skill_dcy,skill_tail,swm_ind,same_swm,yr_cutoff),dpi=300,bbox_inches='tight')

#sys.modules[__name__].__dict__.clear()


#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#plots
swm_ind = 'swm-cc'
#1. 3 x 2 plot of data subsets as defined below
sns.set_theme()
sns.set_style('ticks')
sns.set_context('paper')

R99_swmn4tst = np.sort(swmcc_tst_arr[0,1,:])[int(0.99*len(swmcc_tst_arr[0,1,:]))]
R99_swm0tst = np.sort(swmcc_tst_arr[1,1,:])[int(0.99*len(swmcc_tst_arr[0,1,:]))]
R99_swm4tst = np.sort(swmcc_tst_arr[2,1,:])[int(0.99*len(swmcc_tst_arr[0,1,:]))]
R99_swm8tst = np.sort(swmcc_tst_arr[3,1,:])[int(0.99*len(swmcc_tst_arr[0,1,:]))]

R99_swmn4tst_skilltrn = np.sort(swmcc_tst_arr_skilltrn[0,1,:])[int(0.99*len(swmcc_tst_arr[0,1,:]))]
R99_swm0tst_skilltrn = np.sort(swmcc_tst_arr_skilltrn[1,1,:])[int(0.99*len(swmcc_tst_arr[0,1,:]))]
R99_swm4tst_skilltrn = np.sort(swmcc_tst_arr_skilltrn[2,1,:])[int(0.99*len(swmcc_tst_arr[0,1,:]))]
R99_swm8tst_skilltrn = np.sort(swmcc_tst_arr_skilltrn[3,1,:])[int(0.99*len(swmcc_tst_arr[0,1,:]))]

#f, ((a1,a2,a7,a8),(a3,a4,a9,a10),(a5,a6,a11,a12)) = plt.subplots(3,4,figsize=(12,7),sharex=True,layout='constrained')
f = plt.figure(layout='constrained',figsize=(12,8))
#gs0 = fig.add_gridspec(2,2,width_ratios=[1,2])
#ax1 = fig.add_subplot(gs0[0])
gs0 = f.add_gridspec(1,2,width_ratios=(1.25,1))
gs1 = gs0[0].subgridspec(2,1)
gs2 = gs0[1].subgridspec(3,1)
#gs3 = gs0[3].subgridspec(3,2)

a0 = f.add_subplot(gs1[0])

#f, (a0,a1) = plt.subplots(1,2,width_ratios=[1,2])
#plt.subplot(1,2,1)
a0.plot(np.arange(max_lds)+1,rcurve_arr_swm[1,:],linewidth=2,color=swm_sk0,alpha=0.5)
#a0.plot(np.arange(max_lds)+1,rcurve_arr_swmcc[0,:],linewidth=2,color=chefs)
#a0.plot(np.arange(max_lds)+1,rcurve_arr_swmcc[1,:],linewidth=2,color=cv1)
#a0.plot(np.arange(max_lds)+1,rcurve_arr_swmcc[2,:],linewidth=2,color=cv2)
a0.set_xlabel('Lead time (days)',fontsize='large')
a0.set_ylabel('Risk threshold',fontsize='large')
#a0.set_title('Risk curve',fontsize='x-large',fontweight='bold')
a0.set_ylim([0,1.01])
a0.set_xlim([1,max_lds])
#a0.text(6,0.1,'FIRO pool: '+str(round(pars['firo_pool'],3)) + ' (%s TAF)' %(round((pars['firo_pool']+0.5)*K_scale)),fontsize='x-large')
a0.text(1.5,.9,'a)',fontsize='xx-large',fontweight='bold')
a0.tick_params(axis='both',which='major',labelsize='large')
a0.legend(['$CRPSS-SS_{train}=%s_{swm}$' %(skillmods[1])],loc='lower right',fontsize='large',frameon=False)
#a0.legend(['$CRPSS_{train-swm}=%s$' %(skillmods[0]),'$CC-CRPSS%s$' %(skillmods[0]),'$CC-CRPSS%s$'%(skillmods[1]),'$CC-CRPSS%s$'%(skillmods[2])],loc='lower right',fontsize='large',frameon=False)

evt_idx = evt_idx_swmcc_tst[swmcc_evt_no-1]

evt_date = df_idx.values[evt_idx].strftime("%Y-%m-%d")
evt_mo = evt_date[5:7]
if np.int64(evt_mo) >=8 and np.int64(evt_mo) <=12:
    evt_ref = '1984'+evt_date[4:]
else:
    evt_ref = '1985'+evt_date[4:]
ref_yr_st = '1984-08-01'
ref_yr_ed = '1985-06-30'
st_idx = len(pd.date_range(ref_yr_st,evt_ref,freq='D'))
ed_idx = len(pd.date_range(evt_ref,ref_yr_ed,freq='D'))

dt_idx = np.arange(evt_idx-st_idx,evt_idx+ed_idx)

spill_swmn4tst = np.where(swmcc_tst_arr[0,3,dt_idx]>0)
spill_swm0tst = np.where(swmcc_tst_arr[1,3,dt_idx]>0)
spill_swm4tst = np.where(swmcc_tst_arr[2,3,dt_idx]>0)
spill_swm8tst = np.where(swmcc_tst_arr[3,3,dt_idx]>0)

spills_swmn4tst = np.count_nonzero(swmcc_tst_arr[0,3,:])
spills_swm0tst = np.count_nonzero(swmcc_tst_arr[1,3,:])
spills_swm4tst = np.count_nonzero(swmcc_tst_arr[2,3,:])
spills_swm8tst = np.count_nonzero(swmcc_tst_arr[3,3,:])

spills_swmn4tst_skilltrn = np.count_nonzero(swmcc_tst_arr_skilltrn[0,3,:])
spills_swm0tst_skilltrn = np.count_nonzero(swmcc_tst_arr_skilltrn[1,3,:])
spills_swm4tst_skilltrn = np.count_nonzero(swmcc_tst_arr_skilltrn[2,3,:])
spills_swm8tst_skilltrn = np.count_nonzero(swmcc_tst_arr_skilltrn[3,3,:])

a00 = f.add_subplot(gs1[1])

dt_format=mdates.DateFormatter('%m/%d')
a00.xaxis.set_major_formatter(dt_format)

#1b. plot baseline ops vs hefs 
l1, = a00.plot(df_idx[dt_idx], tocs[dt_idx]/1000, c=cb_gry,alpha=0.5,linewidth=1)
l2, = a00.plot(df_idx[dt_idx], firo_top[dt_idx]/1000, c = cb_brwn,alpha=0.5,linewidth=1)
l7, = a00.plot(df_idx[dt_idx], firo_bot[dt_idx]/1000, c = cb_brwn,alpha=0.5,linewidth=1)
l3, = a00.plot(df_idx[dt_idx], swmcc_tst_arr[0,0,dt_idx]/1000, c=swmcc_sk0_skn4,alpha=0.6,linewidth=1.5)
l4, = a00.plot(df_idx[dt_idx], swmcc_tst_arr[1,0,dt_idx]/1000, c=swmcc_sk0_sk0,alpha=0.6,linewidth=1.5)
l5, = a00.plot(df_idx[dt_idx], swmcc_tst_arr[2,0,dt_idx]/1000, c=swmcc_sk0_skmod,alpha=0.6,linewidth=1.5)
l6, = a00.plot(df_idx[dt_idx], swmcc_tst_arr[3,0,dt_idx]/1000, c=swmcc_sk0_skhi,alpha=0.6,linewidth=1.5)
a00.legend([l3,l4,l5,l6,l1,l2],['$S_{%s}$' %(skillmods[0]),'$S_{%s}$' %(skillmods[1]),'$S_{%s}$' %(skillmods[2]),'$S_{%s}$' %(skillmods[3]),'TOCS','FIRO'],loc='lower left',fontsize='large',frameon=False)
a00.set_ylabel('Storage (MAF)',fontsize='large')
a00.set_xlabel('Date',fontsize='large')
a00.set_ylim([0.65*K/1000, K*1.1/1000])
a00.axhline(3.524,linewidth=1,color='gray',linestyle='--',alpha=0.5)
a00.text(df_idx[dt_idx[180]-150],K*1.07/1000,r'$\overline{S}:$'+str(round(np.mean(swmcc_tst_arr[0,0,:])/1000,3)) +' MAF',color=swmcc_sk0_skn4,fontsize='large')
a00.text(df_idx[dt_idx[180]-150],K*1.05/1000,r'$\overline{S}:$'+str(round(np.mean(swmcc_tst_arr[1,0,:])/1000,3)) +' MAF',color=swmcc_sk0_sk0,fontsize='large')
a00.text(df_idx[dt_idx[180]-150],K*1.03/1000,r'$\overline{S}:$'+str(round(np.mean(swmcc_tst_arr[2,0,:])/1000,3)) +' MAF',color=swmcc_sk0_skmod,fontsize='large')
a00.text(df_idx[dt_idx[180]-150],K*1.01/1000,r'$\overline{S}:$'+str(round(np.mean(swmcc_tst_arr[3,0,:])/1000,3)) +' MAF',color=swmcc_sk0_skhi,fontsize='large')
a00.text(df_idx[dt_idx[180]-30],K*1.07/1000,r'$R_{99}:$'+str(round(R99_swmn4tst,1))+' TAF/d',color=swmcc_sk0_skn4,fontsize='large')
a00.text(df_idx[dt_idx[180]-30],K*1.05/1000,r'$R_{99}:$'+str(round(R99_swm0tst,1))+' TAF/d',color=swmcc_sk0_sk0,fontsize='large')
a00.text(df_idx[dt_idx[180]-30],K*1.03/1000,r'$R_{99}:$'+str(round(R99_swm4tst,1))+' TAF/d',color=swmcc_sk0_skmod,fontsize='large')
a00.text(df_idx[dt_idx[180]-30],K*1.01/1000,r'$R_{99}:$'+str(round(R99_swm8tst,1))+' TAF/d',color=swmcc_sk0_skhi,fontsize='large')
a00.text(df_idx[dt_idx[180]+90],K*1.07/1000,r'$Spills:$'+str(round(spills_swmn4tst)),color=swmcc_sk0_skn4,fontsize='large')
a00.text(df_idx[dt_idx[180]+90],K*1.05/1000,r'$Spills:$'+str(round(spills_swm0tst)),color=swmcc_sk0_sk0,fontsize='large')
a00.text(df_idx[dt_idx[180]+90],K*1.03/1000,r'$Spills:$'+str(round(spills_swm4tst)),color=swmcc_sk0_skmod,fontsize='large')
a00.text(df_idx[dt_idx[180]+90],K*1.01/1000,r'$Spills:$'+str(round(spills_swm8tst)),color=swmcc_sk0_skhi,fontsize='large')
a00.text(df_idx[dt_idx[180]-185],K*1.05/1000,'b)',fontsize='xx-large',fontweight='bold')

evt = df_idx[evt_idx].strftime('%Y-%m-%d')
pad = 15
#evt_idx=df_idx.get_loc(evt)
st,ed = df_idx[evt_idx-pad],df_idx[evt_idx+pad]
start,end = str(st)[0:10],str(ed)[0:10]
dt_idx = np.arange(pad*2+1) + (evt_idx-pad)
dt_format=mdates.DateFormatter('%m/%d')

a1 = f.add_subplot(gs2[0])

ymn = 2300 #0.95*min(np.min(sim_data1[:,dt_idx,0]),np.min(sim_data2[:,dt_idx,0]))
ymx = 1.2* (K-ymn)+ymn

l1, = a1.plot(df_idx[dt_idx], firo_top[dt_idx]/1000, c = cb_brwn, linewidth=1,alpha=0.5)
a1.plot(df_idx[dt_idx], firo_bot[dt_idx]/1000, c = cb_brwn, linewidth=1,alpha=0.5)
l4, = a1.plot(df_idx[dt_idx], tocs[dt_idx]/1000, c = 'gray', linewidth=1,alpha=0.5)
a1.axhline(3.524,linewidth=1,color='gray',linestyle='--',alpha=0.5)
leg1 = a1.legend([l1,l4],['FIRO','TOCS'],bbox_transform=a1.transAxes,loc=(.75,.6),fontsize='large',frameon=False)
a1.add_artist(leg1)
l1, = a1.plot(df_idx[dt_idx], swmcc_tst_arr[0,0,dt_idx]/1000, c=swmcc_sk0_skn4,alpha=0.75)
l2, = a1.plot(df_idx[dt_idx], swmcc_tst_arr[1,0,dt_idx]/1000, c=swmcc_sk0_sk0,alpha=0.75)
l3, = a1.plot(df_idx[dt_idx], swmcc_tst_arr[2,0,dt_idx]/1000, c=swmcc_sk0_skmod,alpha=0.75)
l4, = a1.plot(df_idx[dt_idx], swmcc_tst_arr[3,0,dt_idx]/1000, c=swmcc_sk0_skhi,alpha=0.75)
a1.legend([l1,l2,l3,l4],['$S_{%s}$' %(skillmods[0]),'$S_{%s}$' %(skillmods[1]),'$S_{%s}$' %(skillmods[2]),'$S_{%s}$' %(skillmods[3])],bbox_transform=a1.transAxes,loc=(.01,.5),fontsize='large',frameon=False)
a1.axvline(df_idx[evt_idx],linewidth=0.5,color='black',linestyle='--')
a1.text(df_idx[evt_idx],(1.05*(K-ymn)+ymn)/1000,evt,fontsize='large')
a1.set_ylabel('Storage (MAF)',fontsize='large')
a1.text(df_idx[evt_idx-pad+1],(.1*(ymx-ymn)+ymn)/1000,'c)',fontsize='xx-large',fontweight='bold')
a1.set_ylim([ymn/1000, ymx/1000])
a1.set_xlim([st, ed])
a1.xaxis.set_major_formatter(dt_format)
#plt.gcf().autofmt_xdate()
a1.xaxis.set_ticklabels([])
a1.tick_params(axis='both',which='major',labelsize='large')

a3 = f.add_subplot(gs2[1])

l1, = a3.plot(df_idx[dt_idx], swmcc_tst_arr[0,1,dt_idx] / kcfs_to_tafd,c=swmcc_sk0_skn4,alpha=0.75)
l2, = a3.plot(df_idx[dt_idx], swmcc_tst_arr[1,1,dt_idx] / kcfs_to_tafd,c=swmcc_sk0_sk0,alpha=0.75)
l3, = a3.plot(df_idx[dt_idx], swmcc_tst_arr[2,1,dt_idx] / kcfs_to_tafd, c=swmcc_sk0_skmod,alpha=0.75)
l4, = a3.plot(df_idx[dt_idx], swmcc_tst_arr[3,1,dt_idx] / kcfs_to_tafd, c=swmcc_sk0_skhi,alpha=0.75)
l5, = a3.plot(df_idx[dt_idx], Q_swmcc_tst[dt_idx] / kcfs_to_tafd,c=cb_org,linewidth=3,alpha=0.25)
a3.legend([l1,l2,l3,l4,l5],['$R_{%s}$' %(skillmods[0]),'$R_{%s}$' %(skillmods[1]),'$R_{%s}$' %(skillmods[2]),'$R_{%s}$' %(skillmods[3]),'$Q$'],bbox_transform=a3.transAxes,loc=(0.8,.2),fontsize='large',frameon=False)
a3.axhline(Rmax / kcfs_to_tafd, color='red')
#a3.axhline(R_thresh, color=chefs,linewidth=0.5,alpha=0.5)
a3.axvline(df_idx[evt_idx],linewidth=0.5,color='black',linestyle='--')
#a3.text(df_idx[evt_idx],1.15 * Rmax_scale / kcfs_to_tafd,evt,fontsize='small')
a3.text(df_idx[evt_idx-int(0.75*pad)],1.1*Rmax / kcfs_to_tafd,'$Peak = %s kcfs$' %(str(round(Q_swmcc_tst[evt_idx]/kcfs_to_tafd,1))),color=cb_org,fontsize='large')
a3.text(df_idx[evt_idx+int(pad/2)],1.025*Rmax / kcfs_to_tafd,'$R_{max}$',color='red',fontsize='large')
a3.set_ylabel('Streamflow (kcfs)',fontsize='large')
a3.text(df_idx[evt_idx-pad+1],1.1 * Rmax / kcfs_to_tafd,'d)',fontsize='xx-large',fontweight='bold')
a3.set_ylim([0, 1.25 * Rmax / kcfs_to_tafd])
a3.set_xlim([st, ed])
a3.xaxis.set_major_formatter(dt_format)
#plt.gcf().autofmt_xdate()
a3.xaxis.set_ticklabels([])
a3.tick_params(axis='both',which='major',labelsize='large')

a5 = f.add_subplot(gs2[2])

l1, = a5.plot(df_idx[dt_idx], swmcc_tst_arr[0,4,dt_idx], c=swmcc_sk0_skn4,alpha=0.75)
l2, = a5.plot(df_idx[dt_idx], swmcc_tst_arr[1,4,dt_idx], c=swmcc_sk0_sk0,alpha=0.75)
l3, = a5.plot(df_idx[dt_idx], swmcc_tst_arr[2,4,dt_idx], c=swmcc_sk0_skmod,alpha=0.75)
l4, = a5.plot(df_idx[dt_idx], swmcc_tst_arr[3,4,dt_idx], c=swmcc_sk0_skhi,alpha=0.75)
#a5.legend([l1,l2,l3],['$R^{ld}_{CRPSS%s}$' %(skillmods[0]),'$R^{ld}_{CRPSS%s}$' %(skillmods[1]),'$R^{ld}_{CRPSS%s}$' %(skillmods[2])],loc='upper right',fontsize='large',frameon=False)
a5.axvline(df_idx[evt_idx],linewidth=0.5,color='black',linestyle='--')
#a5.text(df_idx[evt_idx-int(pad*0.7)],14,evt,fontsize='small')
a5.set_ylabel('Release lead time (days)',fontsize='large')
a5.text(df_idx[evt_idx-pad+1],13,'e)',fontsize='xx-large',fontweight='bold')
a5.set_ylim([0, 15.1])
a5.set_xlim([st, ed])
a5.tick_params(axis='both',which='major',labelsize='large')

a5.xaxis.set_major_formatter(dt_format)
plt.setp(plt.xticks()[1],rotation=45,ha='right')

f.savefig('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/ops_plots/3x2_risk-curve-ts-S-Qcp-Rlead_noskilltrn_use-firo-top=%s_use-firo-bottom=%s_seed-%s_evt=%s_mods=%s_dcy=%s_tail=%s_swm=%s_same-swm=%s_cutoff=%s.png' %(use_firo_top,use_firo_bottom,seed,swmcc_evt_no,str(skillmods),skill_dcy,skill_tail,swm_ind,same_swm,yr_cutoff),dpi=300,bbox_inches='tight')

#sys.modules[__name__].__dict__.clear()

#--------------------------------------end--------------------------------------
