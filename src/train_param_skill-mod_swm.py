
import sys
import os
sys.path.insert(0, os.path.abspath('./src'))
import numpy as np
import pandas as pd
import xarray as xr
import model
import syn_util
from scipy.optimize import differential_evolution as DE
from time import localtime, strftime
from numba import njit
import pickle
from datetime import datetime

now=datetime.now()
print('opt start',now.strftime("%H:%M:%S"))

idx = int(sys.argv[1])-1

sd = '1990-10-01' 
ed = '2019-08-15'

sd_swm = '2001-01-01'
ed_swm = '3008-01-08'

loc = 'YRS'
site = 'ORDC1'
vers_out = 'cal'
max_lds = 15
skill_mods = ([-0.4,0,0.4,0.8])
skill_mod = skill_mods[idx]
skill_dcy = 0.1
skill_tail = 0.5

tocs_reset = 'none'   # 'tocs' to reset to baseline tocs at beginning of WY, 'firo' to reset to firo pool, 'none' to use continuous storage profile
use_firo_top = False    #if T, use FVA published top, if F, model optimizes the FIRO top along with risk curve
use_firo_bottom = False  #if T, use FVA published bottom, if F, no restriction besides elev dependent release
seed = 1
opt_pct = 0.9

syn_pct = 0.99
syn_setup = 'cal'
syn_optstrat = 'ecrps-dts'
syn_objpwr = 0
syn_path = './'  # path to R synthetic forecast repo for 'r-gen' setting below

#swm settings
use_bcf = False
k_swm = 106
same_swm = False
swm_ind = 'swm-cc'
samp_no = 1
yr_cutoff = 500

#DE optimizer settings
maxiter = 1000
pol = False
tol = 1e-3  
wkrs = 1

kcfs_to_tafd = 2.29568411*10**-5 * 86400

Q_full,dowy_full,tocs_full,df_idx_full = syn_util.extract_swm(sd,ed,syn_path,loc=loc,site=site,use_bcf=use_bcf,k_swm=k_swm,same_swm=same_swm,swm_ind=swm_ind,samp_no=samp_no,yr_cutoff=yr_cutoff)
Q = Q_full[:-max_lds]
dowy = dowy_full[:-max_lds]
tocs = tocs_full[:-max_lds]
df_idx = df_idx_full[:-max_lds]

#print(np.max(Q)/kcfs_to_tafd)

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
#fixed_firo_bottom_ORO = model.firo_bottomfrac_ORO
#fixed_firo_top_ORO = model.firo_frac_ORO
fixed_firo_top_ORO = ((model.firo_top_ORO-model.TOCS_ORO)/(model.K_ORO - model.TOCS_ORO))
fixed_firo_bottom_ORO = ((model.TOCS_ORO - model.firo_bottom_ORO)/(model.TOCS_ORO))

Qf = syn_util.extract_skillmod_swm(sd,ed,samp_no,syn_path,loc,site,syn_pct,syn_objpwr,syn_optstrat,syn_setup,same_swm,K,skill_mod,skill_dcy,skill_tail,swm_ind,yr_cutoff)
Qf = Qf[:,:,:max_lds] # just use 14 lead days
ne = np.shape(Qf)[1]
nl = max_lds

# sum and sort the forecasts    
Qf_summed = np.cumsum(Qf, axis=2)
Qf_summed_sorted = np.sort(Qf_summed, axis = 1)

@njit
def create_param_risk_curve(x,lds=max_lds):
    lo = x[0]
    hi = x[1]
    pwr = x[2]
    no_risk = int(x[3])
    all_risk = int(x[4])
    if all_risk > 0:
        hi = 1

    if no_risk > 0:
        lo = 0

    ld_sset = np.arange(lds-(no_risk+all_risk)+2)+1
    ld_sset2 = np.arange(lds-(no_risk+all_risk)+1)+1
    ld_sset3 = np.arange(lds-(no_risk+all_risk))+1
    if pwr != 0:
        risk_curve_sset = (np.exp(pwr*ld_sset)-np.exp(pwr)) / (np.exp(pwr*ld_sset[len(ld_sset)-1])-np.exp(pwr))
        risk_curve_sset2 = (np.exp(pwr*ld_sset2)-np.exp(pwr)) / (np.exp(pwr*ld_sset2[len(ld_sset2)-1])-np.exp(pwr))
        risk_curve_sset3 = (np.exp(pwr*ld_sset3)-np.exp(pwr)) / (np.exp(pwr*ld_sset3[len(ld_sset3)-1])-np.exp(pwr))
    if pwr == 0:
        risk_curve_sset = np.interp(ld_sset,np.array([1,ld_sset[-1]]),np.array([0,1]))
        risk_curve_sset2 = np.interp(ld_sset2,np.array([1,ld_sset2[-1]]),np.array([0,1]))
        risk_curve_sset3 = np.interp(ld_sset3,np.array([1,ld_sset3[-1]]),np.array([0,1]))
    if no_risk > 0 and all_risk > 0:
        risk_curve = np.concatenate((np.zeros(no_risk), risk_curve_sset[1:-1],np.ones(all_risk)))
    if no_risk > 0 and all_risk == 0:
        risk_curve = np.concatenate((np.zeros(no_risk), risk_curve_sset2[1:])) * hi
    if no_risk == 0 and all_risk > 0:
        risk_curve_sset = (risk_curve_sset - risk_curve_sset[1]) / (1-risk_curve_sset[1]) * (1-lo)
        risk_curve = np.concatenate((risk_curve_sset+lo, np.ones(all_risk)))[1:-1]
    if no_risk == 0 and all_risk == 0:
        risk_curve_sset = (risk_curve_sset3 - risk_curve_sset3[0]) / (1-risk_curve_sset3[0]) * (1-lo) * hi
        risk_curve = (risk_curve_sset+lo)[:]

    return risk_curve
 
if use_firo_bottom==False and use_firo_top==False:   
    @njit
    def opt_wrapper(x):  # x is an array of the decision variables, where x[0] is the size of firo_pool from 0-0.25
        #global Q,Q_MSG,Qf,Qf_MSG,dowy,tocs,df_idx
        risk_thresholds = create_param_risk_curve(x[2:])
        ix = ((1 - risk_thresholds) * (ne)).astype(np.int32)-1
        S, R, firo_top, firo_bottom, spill, Q_cp, rel_leads, firo_ovg, ddown = model.simulate(firo_pool_top=x[0], firo_pool_bottom=x[1], ix=ix, Q=Q, Qf=Qf_summed_sorted, dowy=dowy, tocs=tocs, K = K, Rmax = Rmax, ramping_rate = ramping_rate, store_vec=store_vec, elev_vec=elev_vec, elev_ctrl_vec=elev_ctrl_vec, flow_ctrl_vec=flow_ctrl_vec, policy='firo', tocs_reset='none', fixed_top=use_firo_top, fixed_bottom=use_firo_bottom)

        #obj = model.objective(S, firo_ovg, Q_cp, K, Rmax, spill, pcnt=opt_pct)
        obj = model.objective(S, firo_ovg, ddown, Q_cp, K, Rmax, spill, dowy, firo_top)
        print(obj)
  
        return obj

    # decision variable bounds: firo_top, firo_bottom, lo, hi, pwr, no_risk, all_risk
    #bounds = [(0.1,1.)] + [(0.25,1.)] + [(0.,0.99),(0.,1.),(-2.,2.),(0.,6.),(0.,6.)] 
    bounds = [(0.,1.)] + [(0.,0.5)] + [(0.,0.99),(0.,1.),(-2.,2.),(0.,6.),(0.,6.)] 
    
    def opt_result(x):
        opt = DE(opt_wrapper, bounds = bounds, disp = False, maxiter = maxiter, tol = tol, polish = pol, init='sobol', workers=wkrs, seed=x)
        
        return opt.x,opt.fun

    opt_out,opt_fun = opt_result(seed)
    
elif use_firo_bottom==True and use_firo_top==False:   
    @njit
    def opt_wrapper(x):  # x is an array of the decision variables, where x[0] is the size of firo_pool from 0-0.25
        #global Q,Q_MSG,Qf,Qf_MSG,dowy,tocs,df_idx
        risk_thresholds = create_param_risk_curve(x[1:])
        ix = ((1 - risk_thresholds) * (ne)).astype(np.int32)-1
        S, R, firo_top, firo_bottom, spill, Q_cp, rel_leads, firo_ovg, ddown = model.simulate(firo_pool_top=x[0], firo_pool_bottom=fixed_firo_bottom_ORO, ix=ix, Q=Q, Qf=Qf_summed_sorted, dowy=dowy, tocs=tocs, K = K, Rmax = Rmax, ramping_rate = ramping_rate, store_vec=store_vec, elev_vec=elev_vec, elev_ctrl_vec=elev_ctrl_vec, flow_ctrl_vec=flow_ctrl_vec, policy='firo', tocs_reset='none', fixed_top=use_firo_top, fixed_bottom=use_firo_bottom)

        #obj = model.objective(S, firo_ovg, Q_cp, K, Rmax, spill, pcnt=opt_pct)
        obj = model.objective(S, firo_ovg, ddown, Q_cp, K, Rmax, spill, dowy, firo_top)
        print(obj)
  
        return obj

    # decision variable bounds: firo_top, firo_bottom, lo, hi, pwr, no_risk, all_risk
    #bounds = [(0.1,1.)] + [(0.,0.99),(0.,1.),(-2.,2.),(0.,6.),(0.,6.)] 
    bounds = [(0.,1.)] + [(0.,0.99),(0.,1.),(-2.,2.),(0.,6.),(0.,6.)] 
    
    def opt_result(x):
        opt = DE(opt_wrapper, bounds = bounds, disp = False, maxiter = maxiter, tol = tol, polish = pol, init='sobol', workers=wkrs, seed=x)
        
        return opt.x,opt.fun

    opt_out,opt_fun = opt_result(seed)  
    opt_out = np.insert(opt_out,1,fixed_firo_bottom_ORO)
    
elif use_firo_bottom==False and use_firo_top==True:   
    @njit
    def opt_wrapper(x):  # x is an array of the decision variables, where x[0] is the size of firo_pool from 0-0.25
        #global Q,Q_MSG,Qf,Qf_MSG,dowy,tocs,df_idx
        risk_thresholds = create_param_risk_curve(x[1:])
        ix = ((1 - risk_thresholds) * (ne)).astype(np.int32)-1
        S, R, firo_top, firo_bottom, spill, Q_cp, rel_leads, firo_ovg, ddown = model.simulate(firo_pool_top=fixed_firo_top_ORO, firo_pool_bottom=x[0], ix=ix, Q=Q, Qf=Qf_summed_sorted, dowy=dowy, tocs=tocs, K = K, Rmax = Rmax, ramping_rate = ramping_rate, store_vec=store_vec, elev_vec=elev_vec, elev_ctrl_vec=elev_ctrl_vec, flow_ctrl_vec=flow_ctrl_vec, policy='firo', tocs_reset='none', fixed_top=use_firo_top, fixed_bottom=use_firo_bottom)

        #obj = model.objective(S, firo_ovg, Q_cp, K, Rmax, spill, pcnt=opt_pct)
        obj = model.objective(S, firo_ovg, ddown, Q_cp, K, Rmax, spill, dowy, firo_top)
        print(obj)
  
        return obj

    # decision variable bounds: firo_top, firo_bottom, lo, hi, pwr, no_risk, all_risk
    #bounds = [(0.25,fixed_firo_bottom_ORO)] + [(0.,0.99),(0.,1.),(-2.,2.),(0.,6.),(0.,6.)] 
    bounds = [(0.,0.5)] + [(0.,0.99),(0.,1.),(-2.,2.),(0.,6.),(0.,6.)] 
    
    def opt_result(x):
        opt = DE(opt_wrapper, bounds = bounds, disp = False, maxiter = maxiter, tol = tol, polish = pol, init='sobol', workers=wkrs, seed=x)
        
        return opt.x,opt.fun

    opt_out,opt_fun = opt_result(seed)  
    opt_out = np.append(fixed_firo_top_ORO,opt_out)
    
elif use_firo_bottom==True and use_firo_top==True:   
    @njit
    def opt_wrapper(x):  # x is an array of the decision variables, where x[0] is the size of firo_pool from 0-0.25
        #global Q,Q_MSG,Qf,Qf_MSG,dowy,tocs,df_idx
        risk_thresholds = create_param_risk_curve(x)
        ix = ((1 - risk_thresholds) * (ne)).astype(np.int32)-1
        S, R, firo_top, firo_bottom, spill, Q_cp, rel_leads, firo_ovg, ddown = model.simulate(firo_pool_top=fixed_firo_top_ORO, firo_pool_bottom=fixed_firo_bottom_ORO, ix=ix, Q=Q, Qf=Qf_summed_sorted, dowy=dowy, tocs=tocs, K = K, Rmax = Rmax, ramping_rate = ramping_rate, store_vec=store_vec, elev_vec=elev_vec, elev_ctrl_vec=elev_ctrl_vec, flow_ctrl_vec=flow_ctrl_vec, policy='firo', tocs_reset='none', fixed_top=use_firo_top, fixed_bottom=use_firo_bottom)

        #obj = model.objective(S, firo_ovg, Q_cp, K, Rmax, spill, pcnt=opt_pct)
        obj = model.objective(S, firo_ovg, ddown, Q_cp, K, Rmax, spill, dowy, firo_top)
        print(obj)
  
        return obj

    # decision variable bounds: firo_top, firo_bottom, lo, hi, pwr, no_risk, all_risk
    bounds = [(0.,0.99),(0.,1.),(-2.,2.),(0.,6.),(0.,6.)] 
    
    def opt_result(x):
        opt = DE(opt_wrapper, bounds = bounds, disp = False, maxiter = maxiter, tol = tol, polish = pol, init='sobol', workers=wkrs, seed=x)
        
        return opt.x,opt.fun

    opt_out,opt_fun = opt_result(seed)  
    opt_out = np.append(([fixed_firo_top_ORO,fixed_firo_bottom_ORO]),opt_out)

opt_params = opt_out[2:]
opt_params[3] = int(opt_params[3])
opt_params[4] = int(opt_params[4])
if opt_params[3] > 0:
    opt_params[0] = 0
if opt_params[4] > 0:
    opt_params[1] = 1
    
#print FIRO pool sizes and parameter values
"""
if use_firo_top==True:
    firo_upr = ((fixed_firo_bottom_ORO*(1-fixed_firo_top_ORO))+fixed_firo_top_ORO)*K
    firo_lwr = opt_out[1]*K
elif use_firo_top==False:
    firo_upr = ((opt_out[1]*(1-opt_out[0]))+opt_out[0])*K
    firo_lwr = opt_out[1]*K
"""

firo_lwr = model.TOCS_ORO - (opt_out[1] * model.TOCS_ORO)
firo_upr = model.TOCS_ORO + (opt_out[0] * (K - model.TOCS_ORO))

#risk_thresholds = np.zeros(max_lds)
risk_thresholds = create_param_risk_curve(opt_out[2:])
ix = ((1 - risk_thresholds) * (ne)).astype(np.int32)-1
S, R, firo_top, firo_bottom, spill, Q_cp, rel_leads, firo_ovg, ddown = model.simulate_nonjit(firo_pool_top=opt_out[0],firo_pool_bottom=opt_out[1], ix=ix, Q=Q, Qf=Qf_summed_sorted, dowy=dowy, tocs=tocs, K = K, Rmax = Rmax, ramping_rate = ramping_rate, store_vec=store_vec, elev_vec=elev_vec, elev_ctrl_vec=elev_ctrl_vec, flow_ctrl_vec=flow_ctrl_vec, policy='firo', tocs_reset='none', fixed_top=use_firo_top, fixed_bottom=use_firo_bottom)
obj = model.objective(S, firo_ovg, ddown, Q_cp, K, Rmax, spill, dowy, firo_top)


f = open("./opt_perf_logs/firot=%s_firob=%s/%s-perf-skillmod-%s_firot=%s_firob=%s_cutoff=%s_seed=%s.txt" %(use_firo_top,use_firo_bottom,swm_ind,skill_mod,use_firo_top,use_firo_bottom,yr_cutoff,seed), "w")
print('Skillmod:' + str(skill_mod),file=f)
print('obj='+str(obj)+' spills='+str(np.count_nonzero(spill))+' Rmax='+str(np.max(R))+' fovgmax='+str(np.max(firo_ovg))+' ddownmax='+str(np.max(ddown)),file=f)
print('Firo bottom:' +str(firo_lwr),file=f)
print('Firo top:' + str(firo_upr),file=f)
print('Opt params:'+' lo='+str(round(opt_params[0],3))+' hi='+str(round(opt_params[1],3))+' pwr='+str(round(opt_params[2],3))+' no_risk='+str(int(opt_params[3]))+' all_risk='+str(int(opt_params[4])),file=f)
f.close()

params = {'obj': opt_fun, 'firo_top': opt_out[0], 'firo_bottom': opt_out[1], 'lo': opt_params[0],'hi': opt_params[1], 'pwr': opt_params[2], 'no_risk': opt_params[3], 'all_risk': opt_params[4]}
pickle.dump(params,open('out/%s/%s/param-risk-thresholds_tocs-reset=%s_use-firo-top=%s_use-firo-bottom=%s_seed-%s_mod=%s_dcy=%s_tail=%s_swm=%s_samp=%s_same-swm=%s_cutoff=%s.pkl' %(loc,site,tocs_reset,use_firo_top,use_firo_bottom,seed,skill_mod,skill_dcy,skill_tail,swm_ind,samp_no,same_swm,yr_cutoff),'wb'))

now=datetime.now()
print('opt end',now.strftime("%H:%M:%S"))



########################################################END#########################################################