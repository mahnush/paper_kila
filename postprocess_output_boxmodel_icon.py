import numpy as np
from netCDF4 import Dataset

# ---open netcdf file

ifileS = '/projekt5/climate/mhaghigh/output_CAM/sajede/remap_output_15032011.nc'
nc = Dataset(ifileS, 'r')
# ---open output file
ncout = Dataset('/projekt5/climate/mhaghigh/output_CAM/sajede/output_icon_test.nc', mode = "w",
                format = 'NETCDF4_CLASSIC')
ncout.description = 'ICON_CCN_file'

# ---read inputs
time = nc.variables['time'][:]
lev = nc.variables['lev'][:]
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]

CCN_01 = nc.variables['CCN_01'][:, ::-1, :, :]
CCN_02 = nc.variables['CCN_02'][:, ::-1, :, :]
CCN_03 = nc.variables['CCN_03'][:, ::-1, :, :]
CCN_04 = nc.variables['CCN_04'][:, ::-1, :, :]
CCN_05 = nc.variables['CCN_05'][:, ::-1, :, :]
CCN_06 = nc.variables['CCN_06'][:, ::-1, :, :]
CCN_07 = nc.variables['CCN_07'][:, ::-1, :, :]
CCN_08 = nc.variables['CCN_08'][:, ::-1, :, :]
CCN_09 = nc.variables['CCN_09'][:, ::-1, :, :]
CCN_10 = nc.variables['CCN_10'][:, ::-1, :, :]
z = nc.variables['z'][:, ::-1, ::-1, :]

# ---define vertical updraft
w = np.array([0.01, 0.0278, 0.0774, 0.215, 0.599, 1.67, 4.64, 12.9, 35.9, 100])

# ---get_dimentions
ntime = len(time)
nlev = len(lev)
nlat = len(lat)
nlon = len(lon)
nw = len(w)

# ---define local variable
nccn = np.zeros((ntime, nw, nlev, nlat, nlon))
z_avg = np.zeros((nlev, nlat, nlon))
# ---average hight over time
for it in range(ntime):
    z_avg[:, :, :] = z_avg[:, :, :] + z[it, :, :, :]
z_avg[:, :, :] = z_avg[:, :, :] / ntime
# ---create N_CCn file
nccn[:, 0, :, :, :] = CCN_01[:, :, :, :]
nccn[:, 1, :, :, :] = CCN_02[:, :, :, :]
nccn[:, 2, :, :, :] = CCN_03[:, :, :, :]
nccn[:, 3, :, :, :] = CCN_04[:, :, :, :]
nccn[:, 4, :, :, :] = CCN_05[:, :, :, :]
nccn[:, 5, :, :, :] = CCN_06[:, :, :, :]
nccn[:, 6, :, :, :] = CCN_07[:, :, :, :]
nccn[:, 7, :, :, :] = CCN_08[:, :, :, :]
nccn[:, 8, :, :, :] = CCN_09[:, :, :, :]
nccn[:, 9, :, :, :] = CCN_10[:, :, :, :]

# ---------------------
ncout.createDimension('lon', nlon)
ncout.createDimension('lat', nlat)
ncout.createDimension('lev', nlev)
ncout.createDimension('w', nw)
ncout.createDimension('time', ntime)
lon_o = ncout.createVariable('lon', np.float32, ('lon',))
lat_o = ncout.createVariable('lat', np.float32, ('lat',))
lev_o = ncout.createVariable('lev', np.float32, ('lev',))
W = ncout.createVariable('w', np.float32, ('w',))
time_o = ncout.createVariable('time', np.float32, ('time',))
N_CCN = ncout.createVariable('N_CCN', np.float32, ('time', 'w', 'lev', 'lat', 'lon'))
LH = ncout.createVariable('LH', np.float32, ('lev', 'lat', 'lon'))

# ------------------------------
lon_o.long_name = 'longitude'
lat_o.long_name = 'lattitude'
time_o.long_name = 'time'
N_CCN.long_name = 'CCN number concentration'
W.long_name = 'vertical velocity'
LH.long_name = 'top height of grid cell'
lev_o.long_name = 'level number'
lat_o.units = 'degrees_north'
lon_o.units = 'degrees_east'
time_o.units = 'hours since 2014-02-10T00:00:00'
LH.units = 'meters'
W.units = 'meters per second'
N_CCN.units = 'particles per cubic meter'

# ----------
lon_o[:] = lon[:]
lat_o[:] = lat[:]
lev_o[:] = lev[:]
W[:] = w[:]
time_o[:] = time[:]
LH[:] = z_avg[:]
N_CCN[:] = nccn[:]
ncout.close
