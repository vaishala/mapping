#! /home/mas001/msh/python_environments/base_env/latest/bin/python

'''
This script takes an input of new pairs of HR-corrected PM2.5 estimate from GEMMACH and MODIS matched locations
This script creates a gridded map of a certain chosen grid size 
The script then calculates the average of all the points over the 2 month range that fall within a single grid
'''
from pathlib import Path  
import numpy as np
import netCDF4 as nc
import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pandas as pd
from scipy.interpolate import Rbf
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import xarray as xr
f1 = '/home/vat000/vaishala_test/data_and_plots/MODIS_3km_GM_Paul_2_5km_File_Pairs_15month/naps_gm_modis_pairs_2018_ALT.nc'
ds1 = nc.Dataset(f1)
lat=ds1['lat_gemmach'][:]
lon=ds1['lon_gemmach'][:]
pm_est=ds1["new_pm25_estimate"][:]

pm_naps=ds1["pm25_naps"][:]
gm_pm25=ds1['gemmach_mass'][:]
gm_aod=ds1['gemmach_aod'][:]
modis_aod=ds1['modis_aod'][:]
gm_ratio=gm_pm25/gm_aod

lonmin=np.amin(lon)
lonmax=np.amax(lon)
latmin=np.amin(lat)
latmax=np.amax(lat)



gridsize = 0.1 #1/10 of a degree

x=np.arange(int(lonmin-1), int(lonmax+1), gridsize) 
y=np.arange(int(latmin-1), int(latmax+1), gridsize)

df_final=pd.DataFrame({})
print(np.shape(x))
print(np.shape(y))

for i in range(np.shape(x)[0]-1):
    for j in range(np.shape(y)[0]-1):
        if np.shape(np.where((y[j]<lat)&(lat<y[j+1]) & (x[i]<lon)&(lon<x[i+1])))!=(1, 0):
            ind=np.where((y[j]<lat)&(lat<y[j+1]) & (x[i]<lon)&(lon<x[i+1])) 
            pm_est_avg=np.average(pm_est[ind])
            gm_ratio_avg=np.average(gm_ratio[ind])
            pm_naps_avg=np.average(pm_naps[ind])
            gm_aod_avg=np.average(gm_aod[ind])
            gm_pm25_avg=np.average(gm_pm25[ind])
            modis_aod_avg=np.average(modis_aod[ind])

            #Change the column names here
            df_pair=pd.DataFrame({"lon":(x[i]+x[i+1])/2,"lat":(y[j+1]+y[j])/2,"pm_est_avg":pm_est_avg,"num_obs":np.shape(ind)[1],"gm_ratio_avg":gm_ratio_avg,"pm_naps_avg":pm_naps_avg,"gm_aod_avg":gm_aod_avg,"gm_pm25_avg":gm_pm25_avg,"modis_aod_avg":modis_aod_avg},index=[0])
            df_final=pd.concat([df_final,df_pair],ignore_index=True)
            print("index",i,j)

ds = xr.Dataset.from_dataframe(df_final)
filepath = Path('/home/vat000/vaishala_test/data_and_plots/pm25_product_data/Paul_Gridded_Data/2018_modis_gm_naps_10th_ofdeg_gridded_data.nc')  #Change filename based on subject
ds.to_netcdf(filepath)
