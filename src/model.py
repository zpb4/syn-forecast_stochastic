import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
from util import water_day
from numba import njit
import calendar

kcfs_to_tafd = 2.29568411*10**-5 * 86400

site = 'ORDC1'

K_ORO = 3524 # TAF
TOCS_ORO = 2788
firo_top_ORO = 2938
firo_bottom_ORO = 2612
top_of_varflood_ORO = 3163
bottom_of_varflood_ORO = 2958
varflood_diff_ORO = top_of_varflood_ORO - bottom_of_varflood_ORO
firo_space_ORO = firo_top_ORO - firo_bottom_ORO
tocs_frac_ORO = TOCS_ORO/K_ORO
firo_topfrac_ORO = firo_top_ORO/K_ORO
firo_bottomfrac_ORO = firo_bottom_ORO/K_ORO
firo_frac_ORO = (firo_top_ORO - firo_bottom_ORO) / (K_ORO - firo_bottom_ORO)

hypso = pd.read_csv('./data/inp_Oroville_Hypso.csv')
store_vec = np.array(hypso['storage'] / 1000)
elev_vec = np.array(hypso['elev'])

CtrlOutlet = pd.read_csv('./data/inp_Oroville_CtrlOutlet.csv')
elev_ctrl_vec = np.array(CtrlOutlet['elev']) 
flow_ctrl_vec = np.array(CtrlOutlet['flow'])

Spillway = pd.read_csv('./data/inp_Oroville_Spillway.csv')
elev_spill_vec = np.array(Spillway['elev']) 
flow_spill_vec = np.array(Spillway['flow'])

res_params = pd.read_csv('./data/reservoir-storage-safe-rel.csv')
res_idx = np.where(res_params['Reservoir']==site)
Rmax = res_params['Safe Release (CFS)'][res_idx[0]].values[0] / 1000 * kcfs_to_tafd
ramping_rate_up = res_params['Ramping-up (CFS)'][res_idx[0]].values[0] / 1000 * kcfs_to_tafd 
ramping_rate_down = res_params['Ramping-down (CFS)'][res_idx[0]].values[0] / 1000 * kcfs_to_tafd 

def extract(sd,ed,forecast_type,syn_sample,Rsyn_path,syn_vers,forecast_param,loc,site,opt_pcnt,objpwr,optstrat,gen_setup,K):

    if forecast_type=='hefs':
        path = '%s/out/%s/Qf-%s.nc' % (Rsyn_path,loc,forecast_type)
    elif forecast_type=='syn' and syn_vers=='v1':
        path = '%s/out/%s/Qf-%s_%s_%s.nc' % (Rsyn_path,loc,forecast_type+forecast_param,site,gen_setup)
    elif forecast_type=='syn' and syn_vers=='v2':
        path = '%s/out/%s/Qf-%s_pcnt=%s_objpwr=%s_optstrat=%s_%s_%s.nc' % (Rsyn_path,loc,forecast_type,opt_pcnt,objpwr,optstrat,site,gen_setup)
        
    da = xr.open_dataset(path)[forecast_type]
    df = pd.read_csv('%s/data/%s/observed_flows.csv' %(Rsyn_path,loc), index_col=0, parse_dates=True)[sd:ed]

    df = df * kcfs_to_tafd
    Q = df[site].values 

    dowy = np.array([water_day(d,calendar.isleap(d.year)) for d in df.index])
    tocs = get_tocs(dowy,K)
    
    ado = {'ADOC1':0}
    nhg = {'MSGC1L':0,'NHGC1':1}
    lam = {'HOPC1L':0,'LAMC1':1,'UKAC1':2}
    yrs = {'MRYC1L':0,'NBBC1':1,'ORDC1':2}
    locs = {'ADO':ado,'NHG':nhg,'LAM':lam,'YRS':yrs}
    
    site_id = locs[loc][site]
    
    # (ensemble: 4, site: 2, date: 15326, trace: 42, lead: 15)
    if forecast_type == 'hefs':
        Qf = da.sel(ensemble=0, site=site_id, date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)
    if forecast_type == 'syn':
        Qf = da.sel(ensemble=int(syn_sample[5:])-1, site=site_id, date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)
    
    #recommend not presorting ensembles because it will mix and match ensemble members
    #Qf.sort(axis = 1)
    #Qf_MSG.sort(axis = 1)
    df_idx = df.index
    
    return Q,Qf,dowy,tocs,df_idx

def extract_skillmod(sd,ed,forecast_type,syn_sample,Rsyn_path,syn_vers,loc,site,opt_pcnt,objpwr,optstrat,gen_setup,K,skill_mod,skill_dcy,skill_tail):

    if forecast_type=='hefs':
        path = '%s/out/%s/Qf-%s_mod=%s_dcy=%s_tail=%s.nc' % (Rsyn_path,loc,forecast_type,skill_mod,skill_dcy,skill_tail)
    elif forecast_type=='syn' and syn_vers=='v2':
        path = '%s/out/%s/Qf-%s_pcnt=%s_objpwr=%s_optstrat=%s_%s_%s_mod=%s_dcy=%s_tail=%s.nc' % (Rsyn_path,loc,forecast_type,opt_pcnt,objpwr,optstrat,site,gen_setup,skill_mod,skill_dcy,skill_tail)
        
    da = xr.open_dataset(path)[forecast_type]
    df = pd.read_csv('%s/data/%s/observed_flows.csv' %(Rsyn_path,loc), index_col=0, parse_dates=True)[sd:ed]

    df = df * kcfs_to_tafd
    Q = df[site].values 

    dowy = np.array([water_day(d,calendar.isleap(d.year)) for d in df.index])
    tocs = get_tocs(dowy,K)
    
    ado = {'ADOC1':0}
    nhg = {'MSGC1L':0,'NHGC1':1}
    lam = {'HOPC1L':0,'LAMC1':1,'UKAC1':2}
    yrs = {'MRYC1L':0,'NBBC1':1,'ORDC1':2}
    locs = {'ADO':ado,'NHG':nhg,'LAM':lam,'YRS':yrs}
    
    site_id = locs[loc][site]
    
    # (ensemble: 4, site: 2, date: 15326, trace: 42, lead: 15)
    if forecast_type == 'hefs':
        Qf = da.sel(ensemble=0, site=site_id, date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)
    if forecast_type == 'syn':
        Qf = da.sel(ensemble=int(syn_sample[5:])-1, site=site_id, date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)
    
    #recommend not presorting ensembles because it will mix and match ensemble members
    #Qf.sort(axis = 1)
    #Qf_MSG.sort(axis = 1)
    df_idx = df.index
    
    return Q,Qf,dowy,tocs,df_idx

def extract_scale(sd,ed,forecast_type,syn_sample,Rsyn_path,syn_vers,forecast_param,loc,site,opt_pcnt,gen_setup,K,scale_site,event_no,rtn_period):

    if forecast_type=='hefs':
        path = '%s/out/%s/Qf-%s_scaled_%s_evt=%s_rtn=%s.nc' % (Rsyn_path,loc,forecast_type,scale_site,event_no,rtn_period)
    elif forecast_type=='syn' and syn_vers=='v1':
        path = '%s/out/%s/Qf-%s_%s_%s.nc' % (Rsyn_path,loc,forecast_type+forecast_param,site,gen_setup)
    elif forecast_type=='syn' and syn_vers=='v2':
        path = '%s/out/%s/Qf-%s_pcnt=%s_%s_%s_scaled_%s_evt=%s_rtn=%s.nc' % (Rsyn_path,loc,forecast_type,opt_pcnt,site,gen_setup,scale_site,event_no,rtn_period)
        
    da = xr.open_dataset(path)[forecast_type]
    df = pd.read_csv('%s/data/%s/observed_flows_scaled_%s_evt=%s_rtn=%s.csv' %(Rsyn_path,loc,scale_site,event_no,rtn_period), index_col=0, parse_dates=True)[sd:ed]

    df = df * kcfs_to_tafd
    Q = df[site].values 

    dowy = np.array([water_day(d,calendar.isleap(d.year)) for d in df.index])
    tocs = get_tocs(dowy,K)
    
    ado = {'ADOC1':0}
    nhg = {'MSGC1L':0,'NHGC1':1}
    lam = {'HOPC1L':0,'LAMC1':1,'UKAC1':2}
    yrs = {'MRYC1L':0,'NBBC1':1,'ORDC1':2}
    locs = {'ADO':ado,'NHG':nhg,'LAM':lam,'YRS':yrs}
    
    site_id = locs[loc][site]
    
    # (ensemble: 4, site: 2, date: 15326, trace: 42, lead: 15)
    if forecast_type == 'hefs':
        Qf = da.sel(ensemble=0, site=site_id, date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)
    if forecast_type == 'syn':
        Qf = da.sel(ensemble=int(syn_sample[5:])-1, site=site_id, date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)
    
    #recommend not presorting ensembles because it will mix and match ensemble members
    #Qf.sort(axis = 1)
    #Qf_MSG.sort(axis = 1)
    df_idx = df.index
    
    return Q,Qf,dowy,tocs,df_idx

def extract86(sd,ed,Rsyn_path,loc,site,K):
    path = '%s/out/%s/Qf-hefs86.nc' % (Rsyn_path,loc)
    da = xr.open_dataset(path)['hefs']
    df = pd.read_csv('%s/data/%s/observed_flows.csv' %(Rsyn_path,loc), index_col=0, parse_dates=True)[sd:ed]

    df = df * kcfs_to_tafd
    Q = df[site].values 
    
    ado = {'ADOC1':0}
    nhg = {'MSGC1L':0,'NHGC1':1}
    lam = {'HOPC1L':0,'LAMC1':1,'UKAC1':2}
    yrs = {'MRYC1L':0,'NBBC1':1,'ORDC1':2}
    locs = {'ADO':ado,'NHG':nhg,'LAM':lam,'YRS':yrs}
    
    site_id = locs[loc][site]

    dowy = np.array([water_day(d,calendar.isleap(d.year)) for d in df.index])
    tocs = get_tocs(dowy,K)
    
    Qf = da.sel(ensemble=0, site=site_id, date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)

    df_idx = df.index
    
    return Q,Qf,dowy,tocs,df_idx

# helper functions
@njit
def get_tocs(d,K):
    tp = np.array([0, 15, 182, 258, 349, 364, 365], dtype=np.float64)
    tocs_inc = (1-tocs_frac_ORO)/30
    sp = np.array([(1-(tocs_inc*16))*K, (1-(tocs_inc*30))*K, (1-(tocs_inc*30))*K, K, K, (1-(tocs_inc*15))*K, (1-(tocs_inc*15))*K], dtype=np.float64)
    return np.interp(d, tp, sp)
"""
@njit
def firo_top(d, firo_pool_top,firo_pool_bottom,K):
  tp = np.array([0, 15, 151, 152, 258, 349, 364, 365], dtype=np.float64)
  firo_inc = (1-(firo_pool_bottom+((1-firo_pool_bottom)*firo_pool_top)))/30
  sp = np.array([(1-(firo_inc*16))*K, (1-(firo_inc*30))*K, (1-(firo_inc*30))*K, min((1-(firo_inc*30))*K + varflood_diff_ORO,K_ORO), K, K, (1-(firo_inc*15))*K, (1-(firo_inc*15))*K], dtype=np.float64)
  return np.interp(d, tp, sp)

@njit
def firo_lwr(d, firo_pool_bottom,K):
  if firo_pool_bottom <= tocs_frac_ORO:
      tp = np.array([0, 15, 151, 182, 258, 349, 364, 365], dtype=np.float64)
      firo_inc = ((1-firo_pool_bottom))/30
      sp = np.array([(1-(firo_inc*16))*K, (1-(firo_inc*30))*K, (1-(firo_inc*30))*K, TOCS_ORO, K, K, (1-(firo_inc*15))*K, (1-(firo_inc*15))*K], dtype=np.float64)
  elif firo_pool_bottom > tocs_frac_ORO:
      tp = np.array([0, 15, 151, 258, 349, 364, 365], dtype=np.float64)
      firo_inc = ((1-firo_pool_bottom))/30
      sp = np.array([(1-(firo_inc*16))*K, (1-(firo_inc*30))*K, (1-(firo_inc*30))*K, K, K, (1-(firo_inc*15))*K, (1-(firo_inc*15))*K], dtype=np.float64)
  return np.interp(d, tp, sp)
"""
"""
#plots for firo bounds and TOCS
dates = pd.date_range('1983-08-01','1984-06-30',freq='D')
dates2 = pd.date_range('1984-08-01','1985-06-30',freq='D')
dowy = np.array([water_day(d,calendar.isleap(d.year)) for d in dates])
tocs = get_tocs(dowy,K_ORO)

plt.plot(dates,tocs,color='black')
plt.plot(dates,firo_top(dowy,firo_frac_ORO,firo_bottomfrac_ORO,K_ORO),color='blue',linestyle='--')
plt.plot(dates,firo_lwr(dowy,firo_bottomfrac_ORO,K_ORO),color='orange',linestyle='--')
plt.show()

plt.plot(dates,np.interp(tocs,store_vec,elev_vec),color='black')
plt.plot(dates,np.interp(firo_top(dowy,firo_frac_ORO,firo_bottomfrac_ORO,K_ORO),store_vec,elev_vec),color='blue',linestyle='--')
plt.plot(dates,np.interp(firo_lwr(dowy,firo_bottomfrac_ORO,K_ORO),store_vec,elev_vec),color='orange',linestyle='--')
plt.show()

plt.plot(dates,tocs,color='black')
plt.plot(dates,firo_top(dowy,1,0.9,K_ORO),color='blue',linestyle='--')
plt.plot(dates,firo_lwr(dowy,0.9,K_ORO),color='orange',linestyle='--')
plt.show()
"""

@njit
def firo_top(d, firo_pool_top,K,tocs_low):
  tp = np.array([0, 15, 151, 152, 258, 349, 364, 365], dtype=np.float64)
  firo_inc = ((1-firo_pool_top)*(K-tocs_low))/30
  sp = np.array([(K - firo_inc*16), (K - firo_inc*30), (K - firo_inc*30), min((K - firo_inc*30) + varflood_diff_ORO,K), K, K, (K - firo_inc*15), (K - firo_inc*15)], dtype=np.float64)
  return np.interp(d, tp, sp)

@njit
def firo_lwr(d, firo_pool_bottom,K,tocs_low):
    tp = np.array([0, 15, 151, 182, 258, 349, 364, 365], dtype=np.float64)
    firo_inc = (firo_pool_bottom * tocs_low)/30
    sp = np.array([(tocs_low - firo_inc*16), (tocs_low - firo_inc*30), (tocs_low - firo_inc*30), tocs_low, K, K, (tocs_low - firo_inc*15), (tocs_low - firo_inc*15)], dtype=np.float64)
 
    return np.interp(d, tp, sp)

"""
#plots for firo bounds and TOCS
dates = pd.date_range('1983-08-01','1984-06-30',freq='D')
dates2 = pd.date_range('1984-08-01','1985-06-30',freq='D')
dowy = np.array([water_day(d,calendar.isleap(d.year)) for d in dates])
tocs = get_tocs(dowy,K_ORO)

plt.plot(dates,tocs,color='black')
plt.plot(dates,firo_top(dowy,((firo_top_ORO-TOCS_ORO)/(K_ORO - TOCS_ORO)),K_ORO,TOCS_ORO),color='blue',linestyle='--')
plt.plot(dates,firo_lwr(dowy,((TOCS_ORO - firo_bottom_ORO)/(TOCS_ORO)),K_ORO,TOCS_ORO),color='orange',linestyle='--')
plt.show()

plt.plot(dates,tocs,color='black')
plt.plot(dates,firo_top(dowy,.5,K_ORO,TOCS_ORO),color='blue',linestyle='--')
plt.plot(dates,firo_lwr(dowy,.2,K_ORO,TOCS_ORO),color='orange',linestyle='--')
plt.show()

plt.plot(dates,np.interp(tocs,store_vec,elev_vec),color='black')
plt.plot(dates,np.interp(firo_top(dowy,firo_frac_ORO,firo_bottomfrac_ORO,K_ORO),store_vec,elev_vec),color='blue',linestyle='--')
plt.plot(dates,np.interp(firo_lwr(dowy,firo_bottomfrac_ORO,K_ORO),store_vec,elev_vec),color='orange',linestyle='--')
plt.show()

plt.plot(dates,tocs,color='black')
plt.plot(dates,firo_top(dowy,1,0.9,K_ORO),color='blue',linestyle='--')
plt.plot(dates,firo_lwr(dowy,0.9,K_ORO),color='orange',linestyle='--')
plt.show()
"""
@njit
def create_risk_curve(x):
  # convert 0-1 values to non-decreasing risk curve
    x_copy = np.copy(x)
    for i in range(1, len(x_copy)):
        x_copy[i] = x_copy[i-1] + (1 - x_copy[i-1]) * x_copy[i]
    return x_copy

@njit
def clip(x, l, u):
    return max(min(x, u), l)

@njit
def ctrl_rel_limit(S,store_vec,elev_vec,elev_ctrl_vec,flow_ctrl_vec):
    S_elev = np.interp(S,store_vec,elev_vec)
    ctrl_flow_limit = np.interp(S_elev,elev_ctrl_vec,flow_ctrl_vec) / 1000 * kcfs_to_tafd
    
    return ctrl_flow_limit

ctrl_rel_limit(3000,store_vec,elev_vec,elev_ctrl_vec,flow_ctrl_vec)

@njit
def daily_opt(S0, fmax, Q, Qf_summed_sorted, ix):
    ne,nl = Qf_summed_sorted.shape
    rel_ld = 0
    Qf_quantiles = np.zeros(nl)
    # iterate over lead times
    for l in range(nl):
        if ix[l] != -1:
            Qf_quantiles[l] = Qf_summed_sorted[ix[l],l]
        elif ix[l] == -1:
            Qf_quantiles[l] = 0
    
    #a risk threshold of near 100% means you accept all risk and don't react to forecasts at all
    ##Qf_quantiles[ix==0] = 0
    # vectorized computation of future storage
    Sf_q = S0 + Qf_quantiles
    
    # vectorized computation to find where future storage exceeds fmax and calculate the required release
    releases = np.maximum(Sf_q - fmax,0) / (np.arange(nl) + 1)
    
    # find the max release and the corresponding lead time
    R = np.max(releases)
    if R>0:
        rel_ld = np.argmax(releases)+1
    
    return R, rel_ld


@njit
def simulate(firo_pool_top, firo_pool_bottom, ix, Q, Qf, dowy, tocs, K, Rmax, ramping_rate, store_vec, elev_vec, elev_ctrl_vec, flow_ctrl_vec, policy='baseline', tocs_reset='tocs', summer_rel_rule='firo',fixed_top=False,fixed_bottom=False):
  lds = np.shape(Qf)[2]
  T = len(Q)
  S = np.full(T+1, np.nan)
  R = np.full(T, np.nan)
  Q_cp = np.full(T, np.nan) # downstream flow at the control point
  Q_sim = np.zeros(T+1)
  Q_sim[:T] = Q
  spill = np.zeros(T)
  firo_ovg = np.zeros(T)
  ddown = np.zeros(T)
  rel_leads = np.zeros(T)
  tocs = get_tocs(dowy,K)
  
  ##manual change to use either flat tocs or precip-based tocs
  #tocs = get_tocs(dowy) # storage guide curve that corresponds to conservative >=12 in antecedent precip
  #tocs = tocs            # use precip-based tocs
  """
  firo_btm = firo_lwr(dowy,firo_pool_bottom,K)
  firo_tp = firo_top(dowy,firo_pool_top,firo_pool_bottom,K) # firo pool guide curve
  if fixed_top == True:
      firo_tp = firo_top(dowy,firo_frac_ORO,firo_bottomfrac_ORO,K) # firo pool guide curve
  if fixed_bottom == True:
      firo_btm = firo_lwr(dowy,firo_bottomfrac_ORO,K) # firo pool guide curve
  """
  firo_btm = firo_lwr(dowy,firo_pool_bottom,K,np.min(tocs))
  firo_tp = firo_top(dowy,firo_pool_top,K,np.min(tocs)) # firo pool guide curve
  if fixed_top == True:
      firo_tp = firo_top(dowy,((firo_top_ORO-TOCS_ORO)/(K_ORO - TOCS_ORO)),K,np.min(tocs)) # firo pool guide curve
  if fixed_bottom == True:
      firo_btm = firo_lwr(dowy,((TOCS_ORO - firo_bottom_ORO)/(TOCS_ORO)),K,np.min(tocs)) # firo pool guide curve
  #S[0] = np.unique(tocs[dowy==7])[0] #start at baseline winter TOCS for all simulation configurations
  #S[0] = K #start full
  S[0] = firo_tp[0] # start at firo pool limit
  for t in range(T):
    if dowy[t] == 0 and tocs_reset == 'tocs':
        S[t] = np.unique(tocs[dowy==7])[0] # start every WY at winter TOCS
      
    elif dowy[t] == 0 and tocs_reset == 'firo':
        S[t] = np.unique(firo_tp[dowy==7])[0] # start every WY at winter FIRO pool
    
    elif dowy[t] == 0 and tocs_reset == 'none':
        S[t] = S[t]
      
    if policy == 'baseline':
        R[t] = clip(S[t] + Q_sim[t+1] - tocs[t], 0, Rmax) #ensure releases will not exceeed downstream Rmax

    elif policy == 'firo' and summer_rel_rule == 'firo':
      R[t],rel_leads[t] = daily_opt(S[t], np.min(firo_tp[t:(t+lds)]), Q_sim[t+1], Qf[t,:,:], ix) # added the daily firo pool value
      if R[t] < (S[t] + Q_sim[t+1] - K): # added to try to fix spill issue 1/10/24
          R[t] = S[t] + Q_sim[t+1] - K # override release calculated from daily_opt to ensure no spill when res is full
      if R[t] > Rmax: # added 4/24/24
          R[t] = Rmax
    
    elif policy == 'firo' and summer_rel_rule == 'baseline':
      #if in summer pool, execute baseline operations
      if dowy[t]>258 and dowy[t]<350:
          R[t] = clip(S[t] + Q_sim[t+1] - tocs[t], 0, Rmax) #ensure releases will not exceeed downstream Rmax
          
      #o/w execute the firo rule
      else:
          R[t],rel_leads[t] = daily_opt(S[t], np.min(firo_tp[t:(t+lds)]), Q_sim[t+1], Qf[t,:,:], ix) # added the daily firo pool value
          if R[t] < (S[t] + Q_sim[t+1] - K): # added to try to fix spill issue 1/10/24
              R[t] = S[t] + Q_sim[t+1] - K # override release calculated from daily_opt to ensure no spill when res is full
          if R[t] > Rmax: 
              R[t] = Rmax
              
    if R[t] > S[t] + Q_sim[t+1]: # release limited by water available
        R[t] = S[t] + Q_sim[t+1]
    
    firo_ddown_limit = firo_btm[t]
    if S[t] + Q_sim[t+1] - R[t] < firo_ddown_limit: # can't draw down below FIRO space for flood release
        R[t] = max(S[t] + Q_sim[t+1] - firo_ddown_limit,0) #release the amount that bottoms out in the FIRO space
      
    if R[t] - R[t-1] > ramping_rate[0]: #release increases limited by ramp-up constraints
        R[t] = R[t-1] + ramping_rate[0]
     
    if R[t-1] - R[t] > ramping_rate[1]: #release decreases limited by ramp-down constraints
        R[t] = R[t-1] - ramping_rate[1]
    
    rel_limit = ctrl_rel_limit(S[t],store_vec,elev_vec,elev_ctrl_vec,flow_ctrl_vec)
    if R[t] > rel_limit:
        rel_diff = max((Q_sim[t+1] - rel_limit)/2,0)  # take 1/2 inflow volume minus current release limit (~mean); only when positive (inflows exceeding release capacity)
        R[t] = min(R[t],ctrl_rel_limit(S[t]+rel_diff,store_vec,elev_vec,elev_ctrl_vec,flow_ctrl_vec)) #reset R_t to current storage + 1/2 inflow
    
    if S[t] + Q_sim[t+1] - R[t] > K: # spill
        spill[t] = S[t] + Q_sim[t+1] - R[t] - K
    
    S[t+1] = S[t] + Q_sim[t+1] - R[t] - spill[t]
    firo_ovg[t] = max(S[t]-firo_tp[t],0)
    ddown[t] = max(firo_tp[t]-S[t],0)
    Q_cp[t] = R[t] # downstream flow at the control point

  return S[:T],R,firo_tp,firo_btm,spill,Q_cp,rel_leads,firo_ovg,ddown

def simulate_nonjit(firo_pool_top, firo_pool_bottom, ix, Q, Qf, dowy, tocs, K, Rmax, ramping_rate, store_vec, elev_vec, elev_ctrl_vec, flow_ctrl_vec, policy='baseline', tocs_reset='tocs', summer_rel_rule='firo',fixed_top=False,fixed_bottom=False):
  lds = np.shape(Qf)[2]
  T = len(Q)
  S = np.full(T+1, np.nan)
  R = np.full(T, np.nan)
  Q_cp = np.full(T, np.nan) # downstream flow at the control point
  Q_sim = np.zeros(T+1)
  Q_sim[:T] = Q
  spill = np.zeros(T)
  firo_ovg = np.zeros(T)
  ddown = np.zeros(T)
  rel_leads = np.zeros(T)
  tocs = get_tocs(dowy,K)
  
  ##manual change to use either flat tocs or precip-based tocs
  #tocs = get_tocs(dowy) # storage guide curve that corresponds to conservative >=12 in antecedent precip
  #tocs = tocs            # use precip-based tocs
  
  """
  firo_btm = firo_lwr(dowy,firo_pool_bottom,K)
  firo_tp = firo_top(dowy,firo_pool_top,firo_pool_bottom,K) # firo pool guide curve
  if fixed_top == True:
      firo_tp = firo_top(dowy,firo_frac_ORO,firo_bottomfrac_ORO,K) # firo pool guide curve
  if fixed_bottom == True:
      firo_btm = firo_lwr(dowy,firo_bottomfrac_ORO,K) # firo pool guide curve
  """
  firo_btm = firo_lwr(dowy,firo_pool_bottom,K,np.min(tocs))
  firo_tp = firo_top(dowy,firo_pool_top,K,np.min(tocs)) # firo pool guide curve
  if fixed_top == True:
      firo_tp = firo_top(dowy,((firo_top_ORO-TOCS_ORO)/(K_ORO - TOCS_ORO)),K,np.min(tocs)) # firo pool guide curve
  if fixed_bottom == True:
      firo_btm = firo_lwr(dowy,((TOCS_ORO - firo_bottom_ORO)/(TOCS_ORO)),K,np.min(tocs)) # firo pool guide curve
  #S[0] = np.unique(tocs[dowy==7])[0] #start at baseline winter TOCS for all simulation configurations
  #S[0] = K #start full
  S[0] = firo_tp[0] # start at firo pool limit
  for t in range(T):
    if dowy[t] == 0 and tocs_reset == 'tocs':
          S[t] = np.unique(tocs[dowy==7])[0] # start every WY at winter TOCS
      
    elif dowy[t] == 0 and tocs_reset == 'firo':
      S[t] = np.unique(firo_tp[dowy==7])[0] # start every WY at winter FIRO pool
    
    elif dowy[t] == 0 and tocs_reset == 'none':
      S[t] = S[t]
      
    if policy == 'baseline':
      R[t] = clip(S[t] + Q_sim[t+1] - tocs[t], 0, Rmax) #ensure releases will not exceeed downstream Rmax

    elif policy == 'firo' and summer_rel_rule == 'firo':
      R[t],rel_leads[t] = daily_opt(S[t], np.min(firo_tp[t:(t+lds)]), Q_sim[t+1], Qf[t,:,:], ix) # added the daily firo pool value
      if R[t] < (S[t] + Q_sim[t+1] - K): # added to try to fix spill issue 1/10/24
          R[t] = S[t] + Q_sim[t+1] - K # override release calculated from daily_opt to ensure no spill when res is full
      if R[t] > Rmax: # added 4/24/24
          R[t] = Rmax
    
    elif policy == 'firo' and summer_rel_rule == 'baseline':
      #if in summer pool, execute baseline operations
      if dowy[t]>258 and dowy[t]<350:
          R[t] = clip(S[t] + Q_sim[t+1] - tocs[t], 0, Rmax) #ensure releases will not exceeed downstream Rmax
          
      #o/w execute the firo rule
      else:
          R[t],rel_leads[t] = daily_opt(S[t], np.min(firo_tp[t:(t+lds)]), Q_sim[t+1], Qf[t,:,:], ix) # added the daily firo pool value
          if R[t] < (S[t] + Q_sim[t+1] - K): # added to try to fix spill issue 1/10/24
              R[t] = S[t] + Q_sim[t+1] - K # override release calculated from daily_opt to ensure no spill when res is full
          if R[t] > Rmax: 
              R[t] = Rmax
              
    if R[t] > S[t] + Q_sim[t+1]: # release limited by water available
        R[t] = S[t] + Q_sim[t+1]
    
    firo_ddown_limit = firo_btm[t]
    if S[t] + Q_sim[t+1] - R[t] < firo_ddown_limit: # can't draw down below FIRO space for flood release
        R[t] = max(S[t] + Q_sim[t+1] - firo_ddown_limit,0) #release the amount that bottoms out in the FIRO space
    
    if R[t] - R[t-1] > ramping_rate[0]: #release increases limited by ramp-up constraints
        R[t] = R[t-1] + ramping_rate[0]
     
    if R[t-1] - R[t] > ramping_rate[1]: #release decreases limited by ramp-down constraints
        R[t] = R[t-1] - ramping_rate[1]
      
    rel_limit = ctrl_rel_limit(S[t],store_vec,elev_vec,elev_ctrl_vec,flow_ctrl_vec)
    if R[t] > rel_limit:
        rel_diff = max((Q_sim[t+1] - rel_limit)/2,0)  # take 1/2 inflow volume minus current release limit (~mean); only when positive (inflows exceeding release capacity)
        R[t] = min(R[t],ctrl_rel_limit(S[t]+rel_diff,store_vec,elev_vec,elev_ctrl_vec,flow_ctrl_vec)) #reset R_t to current storage + 1/2 inflow
    
    if S[t] + Q_sim[t+1] - R[t] > K: # spill
        spill[t] = S[t] + Q_sim[t+1] - R[t] - K
    
    S[t+1] = S[t] + Q_sim[t+1] - R[t] - spill[t]
    firo_ovg[t] = max(S[t]-firo_tp[t],0)
    ddown[t] = max(firo_tp[t]-S[t],0)
    Q_cp[t] = R[t] # downstream flow at the control point

  return S[:T],R,firo_tp,firo_btm,spill,Q_cp,rel_leads,firo_ovg,ddown


@njit
def objective(S,firo_ovg,ddown,R,K,Rmax,spill,dowy,firo_top):
    #idx = np.arange(len(S))
    idx = np.where((dowy > 0) & (dowy <258))
    obj = 0
    obj -= (S[idx]-firo_ovg[idx]).mean()/K  #mean storage below FIRO pool normalized by K (maximize storage excluding FIRO pool overages)
    ##obj += np.mean((R[idx]/Rmax)**2) * 10**-2
    ramping = np.diff(R[idx])
    obj += np.mean((np.abs(ramping[ramping>0]/ramping_rate_up))**2) * 10**-2
    obj += np.mean((np.abs(ramping[ramping<0]/ramping_rate_down))**2) * 10**-2
    norm_firo_ovg = firo_ovg[idx]/(K_ORO-firo_top[idx])
    norm_firo_ovg[np.isnan(norm_firo_ovg)]=0
    obj += np.mean((norm_firo_ovg)**2) * 10**-1 #normalize FIRO overages by flood space (space btwn FIRO top and K), square and sum (minimize FIRO pool overages)
    #obj += np.sum((ddown[idx]/(K_ORO))**2) #normalize FIRO overages by flood pool, square and sum 
    #obj += np.sum(spill**2)         # square and sum spill volume (minimize spill volume)
    obj += np.sum(spill>0) * 10

    return obj

