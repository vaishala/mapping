'''
The script 
'''
from pathlib import Path  
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

filepath = Path('/home/vat000/vaishala_test/data_and_plots/pm25_product_data/Paul_Gridded_Data/2017_modis_gm_naps_10th_ofdeg_gridded_data.nc')  #Change filename based on subject
# filepath = Path('/home/vat000/vaishala_test/data_and_plots/pm25_product_data/Paul_Gridded_Data/2018_modis_gm_naps_10th_ofdeg_gridded_data.nc')  #Change filename based on subject

ds1 = nc.Dataset(filepath)

lon=ds1["lon"][:]
lat=ds1["lat"][:]

pm_est_avg=ds1["pm_est_avg"][:]
num_obs=ds1["num_obs"][:]
gm_ratio_avg=ds1["gm_ratio_avg"][:]
naps_pm_avg=ds1["pm_naps_avg"][:]
modis_aod_avg=ds1["modis_aod_avg"][:]
fig=plt.figure()
map = Basemap(projection='lcc', resolution='h',lat_1=61,lat_2=48,lat_0=55,lon_0=-110, width=3E6, height=2E6)
map.drawcountries(color='k',linewidth=0.5)
map.drawstates(color='k',linewidth=0.1)
map.drawcoastlines(color='k',linewidth=0.1)
map.drawparallels(np.arange(-90,90,5),labels=[1,0,0,0],color='k',linewidth=0.5,dashes=[1,1])
map.drawmeridians(np.arange(-180,180,10),labels=[0,0,0,1],color='k',linewidth=0.5,dashes=[1,1])

#Modify these 4 lines
ind=np.where((pm_est_avg<1000)[0]) & (num_obs>=20))
subject=pm_est_avg #could be number of observatoins or any exported variables
map.scatter(lon[ind],lat[ind],latlon=True,c=subject[ind],s=1,marker='s',edgecolors='k', linewidths=0.001,cmap='Spectral_r')

plt.colorbar(label=u'PM2.5 Estimate 2017(\u03bcg m^-3) ')#'GEM-MACH Ratio')#HR-Corrected PM2.5 Estimate (\u03bcg m^-3)') #change this as needed
output_path = '/home/vat000/vaishala_test/data_and_plots/plots/Paul_GM_Run_2017-2018/Oil_Sands_2017-2018/PAUL_MODIS_NAPS/2017_modis_gm_naps_gridded_pm25_est_map.png' #Change this based on subject

plt.savefig(output_path,dpi=1200)
