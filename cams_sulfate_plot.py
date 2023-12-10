# this code indicates  CAMS reanalyse sulfate aerosls concentration as e geographical distribution
import os

import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
fs_title = 20
# set the cordinate
lat_max = 90
lat_min = -90
lon_max = 180
lon_min = -180
increment = 0.7
# set the coridinate in regional domain close to Kila volcano
map_lat_max = 35
map_lat_min = 5
map_lon_max = -140
map_lon_min = -180


def find_ij(reg_lat_max, reg_lat_min, reg_lon_max, reg_lon_min):
    i_start = int((reg_lon_min ) / increment)
    i_end = int((reg_lon_max ) / increment)
    j_end = int((lat_max -reg_lat_min) / increment)
    j_start = int((lat_max - reg_lat_max) / increment)
    return j_start, j_end, i_start, i_end
reg_lon_min = 180
reg_lon_max = 220
reg_lat_min = 5
reg_lat_max = 35
j_start, j_end, i_start, i_end = find_ij(reg_lat_max, reg_lat_min, reg_lon_max, reg_lon_min)


# read in CAMS reanalyse data:
input_name = '/home/mhaghigh/kilaue/kila_cams/kilauea_cams_monthly/holo_ml.nc'
input_cams = Dataset(input_name, 'r')
lat = input_cams.variables['lat'][j_start:j_end]
lon = input_cams.variables['lon'][i_start:i_end]
sulfate = input_cams.variables['aermr11'][:, :, j_start:j_end, i_start:i_end]
sulfate_mean = np.mean(sulfate, axis=1)


#set the date of data
start = np.datetime64('2014-09-01')
stop = np.datetime64('2014-09-08')
time_day = np.arange(start, stop, dtype='datetime64[D]')
#thereshold for plotting
vmin = 0.
vmax  = 1.4e-8
fig, ax = plt.subplots(2, 2, figsize=(17, 12))
#plt.tight_layout()
#plot the data
def plot(i,j,day):
    n_time = day*8
    map = Basemap(ax = ax[i,j] ,projection='cyl', llcrnrlat=map_lat_min, urcrnrlat=map_lat_max, llcrnrlon=map_lon_min,
              urcrnrlon=map_lon_max, resolution='c')
    map.drawcoastlines()
    ax[i, j].set_title(time_day[day], fontsize=fs_title)
    l1 = map.imshow(sulfate_mean[n_time,:,:], cmap = 'binary',vmin = vmin, vmax = vmax)
    #plt.colorbar(l1,ax=ax[i,j])
    map.drawparallels(np.arange(0., 85., 5.), linewidth=1.2, labels=[1, 0, 0, 0], color='grey', zorder=2, latmax=90,
                  )
    map.drawmeridians(np.arange(-180., 180.,5.), linewidth=1.2, labels=[0, 0, 0, 1], color='grey', zorder=3, latmax=90,
                 )
    return l1


#plot(0,0,1)
#plot(0,1,4)
#plot(1,0,10)
#l1 =plot(1,1,20)
plot(0,0,1)
plot(0,1,2)
plot(1,0,4)
l1 =plot(1,1,5)
#add cbar for all subplots
cax = fig.add_axes([0.15, 0.06, 0.7, 0.02])
cbar = fig.colorbar(l1, cax=cax,orientation = 'horizontal')
cbar.set_label('Sulphate Aerosol Mixing Ratio $kg/kg$', fontsize=fs_title)

#save the output
plt.savefig('geo_dist_sulfate_CAMS_2014.pdf')
plt.show()
