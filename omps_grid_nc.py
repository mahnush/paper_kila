import numpy as np
from netCDF4 import Dataset
import h5py
import MODIS_functions
import postpro_func
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
# change the domain 8.8.2022
limit = np.zeros(4)
#limit[0] = 15
#limit[1] = 25
#limit[2] = -169
#limit[3] = -149
#limit[0] = 10
#limit[1] = 30
#limit[2] = -175
#limit[3] = -145
limit[0] = 6
limit[1] = 34
limit[2] = -178
limit[3] = -142
gridSize = 1
scientific_data = 'SCIENCE_DATA'
so2_PBL = 'ColumnAmountSO2_PBL'
so2_TRL = 'ColumnAmountSO2_TRL'
so2_TRM = 'ColumnAmountSO2_TRM'
so2_STL = 'ColumnAmountSO2_STL'
geolocation_data = 'GEOLOCATION_DATA'
latitute = 'Latitude'
longitude = 'Longitude'


def read_data( first_data_set, second_data_set):
    hf_in = h5py.File(FILE_NAME, 'r')
    # print(list(hf_in.keys()))
    # print(list(hf_in['SCIENCE_DATA']))
    read_var = hf_in[first_data_set][second_data_set]
    temp_var = read_var[:]
    temp_var = temp_var.flatten()
    return temp_var


file_path = '/home/mhaghigh/kilaue/omps_kila/24DEC/'
fileList = open(file_path + 'file_list.txt', 'r+')
allLat = []
allLon = []
allso2_PBL = []
allso2_TRL = []
allso2_TRM = []
allso2_STL = []

for FILE_NAME in fileList:
    FILE_NAME = FILE_NAME.strip()
    FILE_NAME = file_path + FILE_NAME
    print(FILE_NAME)

    if len(allso2_PBL) == 0:
        allso2_PBL = read_data(scientific_data, so2_PBL)
        allso2_TRL = read_data(scientific_data, so2_TRL)
        allso2_TRM = read_data(scientific_data, so2_TRM)
        allso2_STL = read_data(scientific_data, so2_STL)
        allLat = read_data(geolocation_data, latitute)
        allLon = read_data(geolocation_data, longitude)
    elif len(allso2_PBL) > 0 :
        allso2_PBL = np.concatenate((allso2_PBL, read_data(scientific_data, so2_PBL)), axis = 0)
        allso2_TRL = np.concatenate((allso2_TRL, read_data(scientific_data, so2_TRL)), axis = 0)
        allso2_TRM = np.concatenate((allso2_TRM, read_data(scientific_data, so2_TRM)), axis = 0)
        allso2_STL = np.concatenate((allso2_STL, read_data(scientific_data, so2_STL)), axis = 0)
        allLat = np.concatenate((allLat, read_data(geolocation_data, latitute)), axis = 0)
        allLon = np.concatenate((allLon, read_data(geolocation_data, longitude)), axis = 0)

so2_PBL_grid = MODIS_functions.grid(limit, gridSize, allso2_PBL, allLat, allLon)
so2_TRL_grid = MODIS_functions.grid(limit, gridSize, allso2_TRL, allLat, allLon)
so2_TRM_grid = MODIS_functions.grid(limit, gridSize, allso2_TRM, allLat, allLon)
so2_STL_grid = MODIS_functions.grid(limit, gridSize, allso2_STL, allLat, allLon)
lat_grid, lon_grid, nlon, nlat = MODIS_functions.grid_coordinate(limit, gridSize)
print(np.shape(so2_PBL_grid))
so2_PBL_grid_smooth = savgol_filter(so2_PBL_grid, 5, 2, mode='nearest')
so2_PBL_grid_smooth[:, 24:29] = 0
fs_titel = 20
#to remove some  noise part:
#so2_PBL_grid_smooth[0:10, :] = 0
#so2_PBL_grid_smooth[27:31, :] = 0
titel = '24 DEC 2020'
bounds = np.arange(0, 10, 1)
cbar_label = 'So2 amount (DU)'
title = ['Boundary Layer', 'Lower Troposphere(3km)', 'Middle Troposphere (8km)', 'Stratosphere(18km)']
fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2, 2, figsize = (17, 12))
fig.suptitle(titel, fontsize = fs_titel)
postpro_func.visulize_sat_revised(ax0, so2_PBL_grid_smooth, lat_grid, lon_grid, bounds, cbar_label, title[0], limit)
postpro_func.visulize_sat_revised(ax1, so2_TRL_grid, lat_grid, lon_grid, bounds, cbar_label, title[1], limit)
postpro_func.visulize_sat_revised(ax2, so2_TRM_grid, lat_grid, lon_grid, bounds, cbar_label, title[2], limit)
postpro_func.visulize_sat_revised(ax3, so2_STL_grid, lat_grid, lon_grid, bounds, cbar_label, title[3], limit)
plt.savefig('icon_smooth_so2_concentration_24DEC.png')
plt.savefig('./outputs/so2_concentration_27DEC_new_domain.pdf')
# to write OMPS data on netcdf file
ncout = Dataset('/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_24DEC_icon_domain_smooth.nc', mode = "w",
                format = 'NETCDF4_CLASSIC')
# nlon = 52
# nlat = 41
print(nlon, nlat)

ncout.createDimension('lat', nlat)
ncout.createDimension('lon', nlon)
lat_o = ncout.createVariable('lat', np.float32, ('lon', 'lat'))
lon_o = ncout.createVariable('lon', np.float32, ('lon', 'lat'))
so2_PBL_o = ncout.createVariable('so2_PBL', np.float32, ('lon', 'lat'))
so2_TRL_o = ncout.createVariable('so2_TRL', np.float32, ('lon', 'lat'))
so2_TRM_o = ncout.createVariable('so2_TRM', np.float32, ('lon', 'lat'))
so2_STL_o = ncout.createVariable('so2_STL', np.float32, ('lon', 'lat'))

so2_PBL_o[:] = so2_PBL_grid_smooth[:]
so2_TRL_o[:] = so2_TRL_grid[:]
so2_TRM_o[:] = so2_TRM_grid[:]
so2_STL_o[:] = so2_STL_grid[:]
lat_o[:] = lat_grid[:]
lon_o[:] = lon_grid[:]

plt.show()
