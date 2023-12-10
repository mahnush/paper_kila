import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

n_day = 7
fs_label = 15
# set the min_max coridinate in regional domain close to Kila volcano
lat_max = 90
lat_min = -90
lon_max = 180
lon_min = -180
increment = 0.7


def find_ij(reg_lat_max, reg_lat_min, reg_lon_max, reg_lon_min):
    i_start = int((reg_lon_min ) / increment)
    i_end = int((reg_lon_max ) / increment)
    j_end = int((lat_max -reg_lat_min) / increment)
    j_start = int((lat_max - reg_lat_max) / increment)
    return j_start, j_end, i_start, i_end
#set the cordinate of reigon of volcano
reg_lon_min = 180
reg_lon_max = 220
reg_lat_min = 5
reg_lat_max = 35
j_start, j_end, i_start, i_end = find_ij(reg_lat_max, reg_lat_min, reg_lon_max, reg_lon_min)

#to read in data:
# read in CAMS reanalyse data:
input_name = '/work/bb1093/b380900/pycharm/kila_cams_icon.nc'
input_cams = Dataset(input_name, 'r')
time = input_cams.variables['time'][:]
lat = input_cams.variables['lat'][j_start:j_end]
lon = input_cams.variables['lon'][i_start:i_end]
sulfate = input_cams.variables['aermr11'][:, :, j_start:j_end, i_start:i_end]

#set the date of data
start = np.datetime64('2014-09-01')
stop = np.datetime64('2014-09-08')
time_day = np.arange(start, stop, dtype='datetime64[D]')
print(time_day)

#compute the mean value over reigon for Sulfate

sulfate_mean_vertical = np.nanmean(sulfate,axis=1)
sulfate_mean_lat = np.nanmean(sulfate_mean_vertical,axis=1)
sulfate_mean_lon = np.nanmean(sulfate_mean_lat,axis=1)

sulfate_daily_mean = np.zeros(int(len(time)/8.))

#compute the daily mean in a way that mean every 3 hour data for determin the daily mean
n_start = 0
for n in range(n_day):
    n_end= n_start +8
    temp = sulfate_mean_lon[n_start:n_end]
    sulfate_daily_mean[n] = np.nanmean(temp, axis=0)
    n_start = n_start+8



#plot the time series
fig, ax = plt.subplots(1, 1, figsize=(17, 12))
ax.plot(time_day,sulfate_daily_mean, color ='C1', label ='daily-mean')
ax.set_xticks(time_day)
ax.set_xlabel('Time', fontsize=fs_label)
ax.set_ylabel('Sulphate Aerosol Mixing Ratio $kg/kg$', fontsize=fs_label)
ax.legend()
#ax.tick_params(axis='x',labelsize=fs_label)
#ax.tick_params(axis='y',labelsize=fs_label)
plt.xticks(rotation = 45)

#save
plt.savefig('time_series_2014.pdf')
plt.show()
