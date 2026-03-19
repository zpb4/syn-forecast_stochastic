import numpy as np 
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import scipy.stats as sstat
import xarray as xr
from util import water_day
#import cvxpy as cvx
from numba import njit
import calendar
import datetime as dt
import cftime

# using cvxpy for the inner loop prevents numba compiling

kcfs_to_tafd = 2.29568411*10**-5 * 86400
K_ORO = 3524 # TAF
TOCS_ORO = 2788
firo_ORO = 2938
firo_bottom = 2612
firo_space = firo_ORO - firo_bottom
tocs_frac = TOCS_ORO/K_ORO
firo_frac = firo_ORO/K_ORO
ORO_fixed_pool = (firo_ORO - TOCS_ORO) / (K_ORO - TOCS_ORO)

def extract_obs(sd,ed,Rsyn_path,loc,site):
    
    df = pd.read_csv('%s/data/%s/observed_flows.csv' %(Rsyn_path,loc), index_col=0, parse_dates=True)[sd:ed]

    df = df * kcfs_to_tafd
    df_idx = df.index
    Q = df[site].values 

    dowy = np.array([water_day(d,calendar.isleap(d.year)) for d in df.index])
    
    return Q,dowy,df_idx

def extract_swm(sd,ed,Rsyn_path,loc,site,use_bcf,k_swm,same_swm,swm_ind,samp_no,yr_cutoff):

    if use_bcf == True:
      use_bcf = 'TRUE'
    if use_bcf == False:
      use_bcf = 'FALSE'
    if same_swm == True:
      same_swm = 'TRUE'
    if same_swm == False:
      same_swm = 'FALSE'
    
    df = pd.read_csv('%s/out/%s/obs_use-bcf=%s_k=%s_same-swm=%s_cutoff=%s_%s.csv' %(Rsyn_path,loc,use_bcf,k_swm,same_swm,yr_cutoff,swm_ind))
    Q = df['SWM '+str(samp_no)].values * kcfs_to_tafd
    df_idx = df['Date'].apply(lambda x: dt.datetime.strptime(x,'%Y-%m-%d') if type(x)==str else pd.NaT)

    dowy = np.array([water_day(d,calendar.isleap(d.year)) for d in df_idx])
    tocs = get_tocs(dowy,K_ORO)

    return Q,dowy,tocs,df_idx

def extract_obs_scale(sd,ed,Rsyn_path,loc,site,scale_site,event_no,rtn_period):
    
    df = pd.read_csv('%s/data/%s/observed_flows_scaled_%s_evt=%s_rtn=%s.csv' %(Rsyn_path,loc,scale_site,event_no,rtn_period), index_col=0, parse_dates=True)[sd:ed]

    df = df * kcfs_to_tafd
    df_idx = df.index
    Q = df[site].values 

    dowy = np.array([water_day(d,calendar.isleap(d.year)) for d in df.index])
    
    return Q,dowy,df_idx
    
def extract(sd,ed,forecast_type,syn_sample,Rsyn_path,syn_vers,forecast_param,loc,site,opt_pcnt,gen_setup):

    if forecast_type=='hefs':
        path = '%s/out/%s/Qf-%s.nc' % (Rsyn_path,loc,forecast_type)
    elif forecast_type=='syn' and syn_vers=='v1':
        path = '%s/out/%s/Qf-%s_%s_%s.nc' % (Rsyn_path,loc,forecast_type+forecast_param,site,gen_setup)
    elif forecast_type=='syn' and syn_vers=='v2':
        path = '%s/out/%s/Qf-%s_pcnt=%s_%s_%s.nc' % (Rsyn_path,loc,forecast_type,opt_pcnt,site,gen_setup)
        
    da = xr.open_dataset(path)[forecast_type]
    df = pd.read_csv('%s/data/%s/observed_flows.csv' %(Rsyn_path,loc), index_col=0, parse_dates=True)[sd:ed]

    df = df * kcfs_to_tafd
    Q = df[site].values 

    dowy = np.array([water_day(d,calendar.isleap(d.year)) for d in df.index])
    tocs = get_tocs(dowy)
    
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

def extract_skillmod_swm(sd,ed,syn_sample,Rsyn_path,loc,site,opt_pcnt,objpwr,optstrat,gen_setup,same_swm,K,skill_mod,skill_dcy,skill_tail,swm_ind,yr_cutoff):
    if same_swm == True:
      same_swm = 'TRUE'
    if same_swm == False:
      same_swm = 'FALSE'

    path = '%s/out/%s/Qf-syn_pcnt=%s_objpwr=%s_optstrat=%s_%s_%s_same-swm=%s_mod=%s_dcy=%s_tail=%s_cutoff=%s_%s.nc' % (Rsyn_path,loc,opt_pcnt,objpwr,optstrat,site,gen_setup,same_swm,skill_mod,skill_dcy,skill_tail,yr_cutoff,swm_ind)
        
    da = xr.open_dataset(path)['syn']
    
    Qf = da.sel(ensemble=int(syn_sample)-1, site=0).values * kcfs_to_tafd # np.array (time, trace, lead)
    
    return Qf

def extract_scale(sd,ed,forecast_type,syn_sample,Rsyn_path,syn_vers,forecast_param,loc,site,opt_pcnt,gen_setup,scale_site,event_no,rtn_period):

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
    tocs = get_tocs(dowy)
    
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

def extract86(sd,ed,Rsyn_path,loc,site):
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
    tocs = get_tocs(dowy)
    
    Qf = da.sel(ensemble=0, site=site_id, date=df.index).values * kcfs_to_tafd # np.array (time, trace, lead)

    df_idx = df.index
    
    return Q,Qf,dowy,tocs,df_idx


def get_tocs(d,K):
    tp = np.array([0, 15, 182, 258, 349, 364, 365], dtype=np.float64)
    tocs_inc = (1-tocs_frac)/30
    sp = np.array([(1-(tocs_inc*16))*K, (1-(tocs_inc*30))*K, (1-(tocs_inc*30))*K, K, K, (1-(tocs_inc*15))*K, (1-(tocs_inc*15))*K], dtype=np.float64)
    return np.interp(d, tp, sp)

def create_param_risk_curve(x,lds):
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
        #risk_curve_sset = risk_curve_sset3 * (hi - lo) + lo
        risk_curve = (risk_curve_sset+lo)[:]
        #risk_curve = risk_curve_sset

    return risk_curve

def create_risk_curve(x):
	# convert 0-1 values to non-decreasing risk curve
    x_copy = np.copy(x)
    for i in range(1, len(x_copy)):
        x_copy[i] = x_copy[i-1] + (1 - x_copy[i-1]) * x_copy[i]
    return x_copy



"""
lo = 0.2
hi = 0.5
pwr = 1
no_risk = 0
all_risk = 0
ls = 14

x = ([lo,hi,pwr,no_risk,all_risk])

rc = create_param_risk_curve(x, ls)
rc
"""

def declust_evts_extract(Q,n_evts,sep):
    rnk_data = sstat.rankdata(-Q)
    srt_rnks_idx = np.argsort(rnk_data)[0:n_evts*10]

    vec = np.zeros(len(srt_rnks_idx))
    for i in range(len(srt_rnks_idx)):
        vec[i]=np.where(rnk_data==np.min(rnk_data[max((srt_rnks_idx[i]-sep),0):min(len(Q),(srt_rnks_idx[i]+sep))]))[0][0]
    declust_evts=np.unique(vec)
    sel_idx=np.argsort(-Q[np.int64(declust_evts)])
    evt_idx = np.int64(declust_evts[sel_idx])[0:n_evts]
    
    return evt_idx

def fig_title(
    fig: matplotlib.figure.Figure, txt: str, loc=(0.5,0.98), fontdict=None, **kwargs
):
    """Alternative to fig.suptitle that behaves like ax.set_title.
    DO NOT use with suptitle.

    See also:
    https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_title.html
    https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.suptitle.html
    https://stackoverflow.com/a/77063164/8954109
    """
    if fontdict is not None:
        kwargs = {**fontdict, **kwargs}
    if "fontsize" not in kwargs and "size" not in kwargs:
        kwargs["fontsize"] = plt.rcParams["axes.titlesize"]

    if "fontweight" not in kwargs and "weight" not in kwargs:
        kwargs["fontweight"] = plt.rcParams["figure.titleweight"]

    if "verticalalignment" not in kwargs:
        kwargs["verticalalignment"] = "top"
    if "horizontalalignment" not in kwargs:
        kwargs["horizontalalignment"] = "center"

    # Tell the layout engine that our text is using space at the top of the figure
    # so that tight_layout does not break.
    # Is there a more direct way to do this?
    fig.suptitle(" ")
    text = fig.text(loc[0], loc[1], txt, transform=fig.transFigure, in_layout=True, **kwargs)

    return text