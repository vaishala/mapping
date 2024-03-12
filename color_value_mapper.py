#! /home/mas001/msh/python_environments/base_env/latest/bin/python

'''
This script takes an input of an nc file with a variable we want to plot on a map and generates a basemap with the variable values represented by the colorbar
'''




#Different Basemaps based on region you want to plot
'''






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
f1 = '/home/vat000/vaishala_test/data_and_plots/MODIS_GEM-MACH_2km_File_Pairs/MODIS_GEMMACH_PM25estimate_pairs_new_with_second_lowest_mass.nc'
ds1 = nc.Dataset(f1)


lon=ds1["lon_modis"][:]
lat=ds1["lat_modis"][:]

mass_gemmach_secondlowest=ds1['mass_gemmach_secondlowest'][:]
print(np.amin(mass_gemmach_secondlowest))
mass_gemmach_lowest=ds1['mass_gemmach'][:]
print(np.amin(mass_gemmach_lowest))

sum=mass_gemmach_lowest+mass_gemmach_secondlowest

mindiff=-20
maxdiff=-2
difference2=difference[np.where((difference>mindiff)&(difference<maxdiff))]
lat2=lat[np.where((difference>mindiff)&(difference<maxdiff))]
lon2=lon[np.where((difference>mindiff)&(difference<maxdiff))]


fig=plt.figure()
map = Basemap(projection='lcc', resolution='h',lat_1=61,lat_2=48,lat_0=55,lon_0=-110, width=3E6, height=2E6)
#map.scatter(lon,lat,latlon=True,c=difference,s=0.05,marker='d',edgecolors='k', linewidths=0.0005)#zorder=10, linewidths=0.05, marker='D')#,cmap='Spectral_r',s=0.001,edgecolors='face',vmin=0)#, linewidths=0.05)

map.scatter(lon2,lat2,latlon=True,c=difference2,s=20,marker='d',edgecolors='k', linewidths=0.0005)#zorder=10, linewidths=0.05, marker='D')#,cmap='Spectral_r',s=0.001,edgecolors='face',vmin=0)#, linewidths=0.05)
map.drawcountries(color='k',linewidth=0.3)
map.drawcoastlines(color='k',linewidth=0.1)
map.drawparallels(np.arange(-90,90,5),labels=[1,0,0,0],color='k',linewidth=0.5,dashes=[1,1])
map.drawmeridians(np.arange(-180,180,10),labels=[0,0,0,1],color='k',linewidth=0.5,dashes=[1,1])
plt.colorbar(label=u'Total PM2.5 Difference: Lowest Layer-Second Lowest Layer (\u03bcg m^-3)')
output_path = '/home/vat000/vaishala_test/data_and_plots/plots/mass/layer_difference_pm25_min.png'
plt.savefig(output_path,dpi=1200)
