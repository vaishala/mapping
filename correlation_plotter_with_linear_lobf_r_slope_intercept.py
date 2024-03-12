#! /home/mas001/msh/python_environments/base_env/latest/bin/python

'''
This script takes an input of pairs of PM2.5 estimate from GEMMACH and MODIS matched to a NAPS PM2.5 surface observation
This script creates a correlation plot between the pairs
'''

import numpy as np
import netCDF4 as nc
import pandas as pd
from pathlib import Path  
import datetime
import scipy.stats as stats
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


filepath = '/home/vat000/vaishala_test/data_and_plots/MODIS_3km_GM_Paul_2_5km_File_Pairs_15month/MODIS_GEMMACH_PM25estimate_pairs_2017_ALT.nc'
filepath2 = '/home/vat000/vaishala_test/data_and_plots/MODIS_3km_GM_Paul_2_5km_File_Pairs_15month/MODIS_GEMMACH_PM25estimate_pairs_2018_ALT.nc'
ds1 = nc.Dataset(filepath)
ds2 = nc.Dataset(filepath2)


modis_aod_2017=ds1['aod_modis'][:]
modis_aod_2018=ds2['aod_modis'][:]
gm_aod_2017=ds1['aod_gemmach'][:]
gm_aod_2018=ds2['aod_gemmach'][:]

X=np.concatenate((modis_aod_2017,modis_aod_2018))
print(np.shape(X))
Y=np.concatenate((gm_aod_2017,gm_aod_2018))

df=pd.DataFrame({'x':X,'y':Y})
fig=plt.figure()
ax=fig.add_subplot(111)
plt.scatter(df.x,df.y,s=0.5)
plt.xlabel(u'MODIS AOD')
plt.ylabel(u'Paul AOD')

r,p=stats.pearsonr(df.x,df.y)
theta = np.polyfit(df.x, df.y, 1)
y_line = theta[1] + theta[0] * df.x
m=theta[0]
b=theta[1]
plt.plot(df.x, y_line, 'r')
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax.text(0.65,0.6,'slope='+str(np.around(theta[0],decimals=5))+'\noffset='+str(np.around(theta[1],decimals=5)),transform=fig.transFigure,fontsize=10,bbox=props,weight='bold')
plt.annotate('r = {:.2f}'.format(r), xy=(0.75, 0.5), xycoords='axes fraction',bbox=props,weight='bold')


plt.grid(True,alpha=0.5)
output_path = '/home/vat000/vaishala_test/data_and_plots/plots/Paul_GM_Run_2017-2018/Oil_Sands_2017-2018/PAULXMODIS/modis_aod_vs_paul_aod_correlation.png'
plt.savefig(output_path)
