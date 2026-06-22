# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 09:43:22 2026

@author: zpb4
"""
import sys
import os
sys.path.insert(0, os.path.abspath('./src'))
import numpy as np
import pandas as pd
import matplotlib as matplotlib
import matplotlib.pyplot as plt

hypso = pd.read_csv('./data/inp_Oroville_Hypso.csv')
store_vec = np.array(hypso['storage'] / 1e6)
elev_vec = np.array(hypso['elev'])

CtrlOutlet = pd.read_csv('./data/inp_Oroville_CtrlOutlet.csv')
elev_ctrl_vec = np.array(CtrlOutlet['elev']) 
flow_ctrl_vec = np.array(CtrlOutlet['flow'])/1000  #change cfs to kcfs

fig, ax1 = plt.subplots()

l1, = ax1.plot(flow_ctrl_vec,elev_ctrl_vec,linewidth=2)
l2, = ax1.plot(0,0,color='gray',linestyle='--')
l3, = ax1.plot(0,0,color='gray')
l4, = ax1.plot(0,0,color='brown',linestyle='--')
l5, = ax1.plot(0,0,color='red')

ax1.axhline(901,color='gray',linestyle='--')
ax1.axhline(848.5,color='gray')
ax1.axhline(859.5,color='brown',linestyle='--')
ax1.axhline(835,color='brown',linestyle='--')
ax1.axvline(150,color='red')

ax1.set_xlabel('Flow (kcfs)')
ax1.set_ylabel('Elevation (ft)')
ax1.set_ylim([750,910])
ax1.set_xlim([0,300])

def elev_to_store(elev):
    return np.interp(elev,elev_vec,store_vec)

def store_to_elev(store):
    return np.interp(store,store_vec,elev_vec)
    
secax = ax1.secondary_yaxis('right',functions=(elev_to_store, store_to_elev))
secax.set_ylabel('Storage (MAF)')

ax1.legend([l1,l2,l3,l4,l5],['Release Limit','Top of Gross Pool','Top of Conservation','Winter FIRO Space','Max Release'],loc='lower right',fontsize='large',frameon=False)

fig.savefig('e:/Projects/FIRO/firo_syn-forecast_stochastic-CC/figs/manuscript/ORO-hypsometry.png' ,dpi=300,bbox_inches='tight')


#####################################################END################################################################################