#!/usr/bin/env python
# coding: utf-8

# # Evaluating the results from PDAF

# In[1]:


import numpy as np
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
import xarray as xr
from xmitgcm import open_mdsdataset
import glob
import cmocean.cm as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#==========
# settings
#==========
rdir = '/home/hajsong/PDAF/tutorial_global_oce_latlon_v2/run_pdaf/'
# rdir = '/home/hajsong/PDAF/tutorial_global_oce_latlon_v1.16_omi_hj/run_pdaf'
nens = 4
obsfreq = 10
files = sorted(glob.glob(rdir+'ETAN_analysis.*.meta'))
tidx=[]
for it in range(len(files)):
    tidx.append(int(files[it][-12:-5]))
nt = len(tidx)

#=================================
# create dictionary for variables
#=================================
extra_variables = dict(
    THETA_analysis = dict(dims=['k','j','i'], 
                          attrs=dict(standard_name='T_analysis',
                                     long_name='Theta after data assimilation', 
                                     units='degC')),
    THETA_forecast = dict(dims=['k','j','i'], 
                          attrs=dict(standard_name='T_forecast',
                                     long_name='Theta before data assimilation', 
                                     units='degC')),
    SALT_analysis = dict(dims=['k','j','i'], 
                         attrs=dict(standard_name='S_analysis',
                                    long_name='Salt after data assimilation', 
                                    units='psu')),
    SALT_forecast = dict(dims=['k','j','i'], 
                         attrs=dict(standard_name='S_forecast',
                                    long_name='Salt before data assimilation', 
                                    units='psu')),
    ETAN_analysis = dict(dims=['j','i'], 
                         attrs=dict(standard_name='eta_analysis',
                                    long_name='SSH after data assimilation', 
                                    units='m')),
    ETAN_forecast = dict(dims=['j','i'], 
                         attrs=dict(standard_name='eta_forecast',
                                    long_name='SSH before data assimilation', 
                                    units='m'))
)


for n in range(1,nens+1):
    trname = 'THETA_%03d_analysis' % n
    extra_variables[trname] = dict(dims=['k','j','i'], 
                                   attrs=dict(standard_name='%s' % trname,
                                              long_name="ensemble member for %s" % trname,
                                              units="degC")
                                  )
for n in range(1,nens+1):
    trname = 'THETA_%03d_forecast' % n
    extra_variables[trname] = dict(dims=['k','j','i'], 
                                   attrs=dict(standard_name='%s' % trname,
                                              long_name="ensemble member for %s" % trname,
                                              units="degC")
                                  )
for n in range(1,nens+1):
    trname = 'SALT_%03d_analysis' % n
    extra_variables[trname] = dict(dims=['k','j','i'], 
                                   attrs=dict(standard_name='%s' % trname,
                                              long_name="ensemble member for %s" % trname,
                                              units="psu")
                                  )
for n in range(1,nens+1):
    trname = 'SALT_%03d_forecast' % n
    extra_variables[trname] = dict(dims=['k','j','i'], 
                                   attrs=dict(standard_name='%s' % trname,
                                              long_name="ensemble member for %s" % trname,
                                              units="psu")
                                  )
for n in range(1,nens+1):
    trname = 'ETAN_%03d_analysis' % n
    extra_variables[trname] = dict(dims=['k','j','i'], 
                                   attrs=dict(standard_name='%s' % trname,
                                              long_name="ensemble member for %s" % trname,
                                              units="m")
                                  )
for n in range(1,nens+1):
    trname = 'ETAN_%03d_forecast' % n
    extra_variables[trname] = dict(dims=['k','j','i'], 
                                   attrs=dict(standard_name='%s' % trname,
                                              long_name="ensemble member for %s" % trname,
                                              units="m")
                                  )


#=====================
#  Read model output
#=====================
ds = open_mdsdataset(rdir, rdir, iters=tidx, delta_t=86400,
                     prefix=extra_variables.keys(), extra_variables=extra_variables)    # takes too long without "iters=None"


#====================
#  Read observations
#====================
# with open(rdir+'/Eta_obs.0000000030.data', 'rb') as f:
#     obs = np.fromfile(f, '>f4')
# obs = obs.reshape(40,90)
obs = 0.5


#==========================
#  Evaluating performance
#==========================
#----------
# 1. RMSE
#----------
fmean = ds.ETAN_forecast.where(ds.hFacC[0]>0)
amean = ds.ETAN_analysis.where(ds.hFacC[0]>0)

rms_eta_f = np.sqrt(np.nanmean((fmean.data - obs).reshape(nt,-1)**2, axis=1))
rms_eta_a = np.sqrt(np.nanmean((amean.data - obs).reshape(nt,-1)**2, axis=1))

rms_f = np.round(np.asarray(rms_eta_f), decimals=4)
rms_a = np.round(np.asarray(rms_eta_a), decimals=4)

tstep = [str(tidx[it]) for it in range(nt)]
print('==================')
print(' [ RMS residual ]')
print('==================')
print('timestep:  '+ " ".join([f'{ x:<8}' for x in tstep]))
print('----------'*nt)
print('forecast:  '+ " ".join([f'{ x:.6f}' for x in rms_f]))
print('analysis:  '+ " ".join([f'{ x:.6f}' for x in rms_a]))
print(' ')

#------------
# 2. Spread
#------------

s_a = xr.zeros_like(ds.ETAN_analysis)
s_f = xr.zeros_like(ds.ETAN_forecast)

for ie in range(nens):
    aname = 'ETAN_%03d_analysis' % (ie+1)
    fname = 'ETAN_%03d_forecast' % (ie+1)
    avar = ds[aname].where(ds.hFacC[0]>0)
    fvar = ds[fname].where(ds.hFacC[0]>0)
    s_a += (amean-avar)**2
    s_f += (fmean-fvar)**2
    
s_a = s_a / (nens-1)
s_f = s_f / (nens-1)


sp_eta_f = np.zeros(len(ds.time))
sp_eta_a = np.zeros(len(ds.time))
sp_eta_f = [np.nanmean(s_a[i]) for i in range(len(ds.time))]
sp_eta_a = [np.nanmean(s_f[i]) for i in range(len(ds.time))]

print('=====================')
print(' [ ensemble spread ]')
print('=====================')
print('timestep:  '+ " ".join([f'{ x:<8}' for x in tstep]))
print('----------'*nt)
print('forecast:  '+ " ".join([f'{ x:.6f}' for x in sp_eta_f]))
print('analysis:  '+ " ".join([f'{ x:.6f}' for x in sp_eta_a]))
print(' ')


#--------------------------
# 3. Check extreme values
#--------------------------
# (a) temperature
var_max_f = np.zeros(nt)
var_min_f = np.zeros(nt)
var_max_a = np.zeros(nt)
var_min_a = np.zeros(nt)

for ie in range(nens):
    aname = 'THETA_%03d_analysis' % (ie+1)
    fname = 'THETA_%03d_forecast' % (ie+1)
    avar = ds[aname].where(ds.hFacC>0)
    fvar = ds[fname].where(ds.hFacC>0)
    max_a = np.nanmax(np.asarray(avar).reshape(nt,-1), axis=1)
    min_a = np.nanmin(np.asarray(avar).reshape(nt,-1), axis=1)
    max_f = np.nanmax(np.asarray(fvar).reshape(nt,-1), axis=1)
    min_f = np.nanmin(np.asarray(fvar).reshape(nt,-1), axis=1)
    var_max_f = np.maximum(var_max_f, max_f)
    var_min_f = np.minimum(var_min_f, min_f)
    var_max_a = np.maximum(var_max_a, max_a)
    var_min_a = np.minimum(var_min_a, min_a)

print('======================')
print(' [ extreme vales, T ]')
print('======================')
print('timestep:  '+ " ".join([f'{ x:<8}' for x in tstep]))
print('----------'*nt)
print('forecast:  '+ " ".join([f'{ x:.5f}' for x in var_max_f]))
print('        :  '+ " ".join([f'{ x:.5f}' for x in var_min_f]))
print('analysis:  '+ " ".join([f'{ x:.5f}' for x in var_max_a]))
print('        :  '+ " ".join([f'{ x:.5f}' for x in var_min_a]))
print(' ')

# (b) salinity
var_max_f = np.ones(nt)*35
var_min_f = np.ones(nt)*35
var_max_a = np.ones(nt)*35
var_min_a = np.ones(nt)*35

for ie in range(nens):
    aname = 'SALT_%03d_analysis' % (ie+1)
    fname = 'SALT_%03d_forecast' % (ie+1)
    avar = ds[aname].where(ds.hFacC>0)
    fvar = ds[fname].where(ds.hFacC>0)
    max_a = np.nanmax(np.asarray(avar).reshape(nt,-1), axis=1)
    min_a = np.nanmin(np.asarray(avar).reshape(nt,-1), axis=1)
    max_f = np.nanmax(np.asarray(fvar).reshape(nt,-1), axis=1)
    min_f = np.nanmin(np.asarray(fvar).reshape(nt,-1), axis=1)
    var_max_f = np.maximum(var_max_f, max_f)
    var_min_f = np.minimum(var_min_f, min_f)
    var_max_a = np.maximum(var_max_a, max_a)
    var_min_a = np.minimum(var_min_a, min_a)

print('======================')
print(' [ extreme vales, S ]')
print('======================')
print('timestep:  '+ " ".join([f'{ x:<8}' for x in tstep]))
print('----------'*nt)
print('forecast:  '+ " ".join([f'{ x:.5f}' for x in var_max_f]))
print('        :  '+ " ".join([f'{ x:.5f}' for x in var_min_f]))
print('analysis:  '+ " ".join([f'{ x:.5f}' for x in var_max_a]))
print('        :  '+ " ".join([f'{ x:.5f}' for x in var_min_a]))
print(' ')
