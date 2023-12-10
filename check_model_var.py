import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import glob
import numpy.ma as ma
from netCDF4 import Dataset
import postpro_func
from mpl_toolkits.basemap import Basemap

data ='/home/mhaghigh/run_new_version_icon/NWP_LAM_DOM01_20201222T180000Z_0018.nc'
data_con='/home/mhaghigh/run_new_version_icon/remap_NWP_LAM_DOM01_20201222T180000Z_0018.nc'
nc = Dataset(data,'r')
var = nc.variables['tqc'][0, :, :]*1000
nc_con = Dataset(data_con,'r')
var_con = nc_con.variables['tqc'][0, :, :]*1000
fig, ((ax0, ax1)) = plt.subplots(1, 2, figsize=(17, 12))
titel = 'Test Run'
fs_titel = 20
fig.suptitle(titel, fontsize=fs_titel)
lwp_bound = np.arange(0, 300, 1)
print(var)
print(np.max(var))
limit = np.zeros(4)
limit[0] = 14
limit[1] = 26
limit[2] = -174
limit[3] = -147
gridSize = 0.02
cbar_label = [ '$\mathrm{cm^{-3}}$', '']
titel_var = ['diffNd', 'Nd_rccn', 'Nd']
postpro_func.visulize_model(ax0, var, lwp_bound, cbar_label[0],titel_var[0], limit)
postpro_func.visulize_model(ax1, var_con, lwp_bound, cbar_label[0],titel_var[0], limit)
plt.savefig('testlwp_model.png')
plt.show()