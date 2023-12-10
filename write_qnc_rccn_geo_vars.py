import numpy as np
from netCDF4 import Dataset

ncout = Dataset('/home/b/b380900/icon_nwp_newversion/icon/experiments/icon_lam_lim_kila_1dom/ncfiles_select_var/vars_3dim_per_22.nc',mode="w",format='NETCDF4_CLASSIC')
#ipath='/work/bb1093/b380900/holouraun_nccn/output_nccn_alpha/new_run_ctt/output_uns/remap_files/'
ipath= '/home/b/b380900/icon_nwp_newversion/icon/experiments/icon_lam_lim_kila_1dom/'
fileList=open('/home/b/b380900/icon_nwp_newversion/icon/experiments/icon_lam_lim_kila_1dom/ncfiles_select_var/file_list_22dec.txt','r+')
#fileList=open('/work/bb1093/b380900/holouraun_nccn/output_nccn_alpha/new_run_ctt/output_uns/remap_files/file_list_1sep.txt','r+')
nlon = 1801
nlat = 1401
ns = 1
nlev = 75
#ns=24
qnc = np.zeros((ns, nlev, nlat, nlon))
rccn = np.zeros((ns, nlev, nlat, nlon))
geo = np.zeros((ns, nlev, nlat, nlon))

ih=0

for FILE_NAME in fileList:
    FILE_NAME=FILE_NAME.strip()
    FILE_NAME=ipath+FILE_NAME
    print(FILE_NAME)
    nc_dw = Dataset(FILE_NAME,'r')
    qnc[ih, :, :, :] = nc_dw.variables['qnc'][0, :, :, :]
    rccn[ih, :, :, :] = nc_dw.variables['rccn'][0, 0:75, :, :]
    geo[ih, :, :, :] = nc_dw.variables['geopot'][0, :, :, :]
    ih+= 1
    print(ih)
lon=nc_dw.variables['lon'][:]
lat=nc_dw.variables['lat'][:]
limit = np.zeros(4)
limit[0] = 6
limit[1] = 34
limit[2] = -178
limit[3] = -142

ncout.createDimension('lon',nlon)
ncout.createDimension('lat',nlat)
ncout.createDimension('time',ns)
ncout.createDimension('lev', nlev)
lon_o = ncout.createVariable('lon', np.float32, ('lon',))
lat_o = ncout.createVariable('lat', np.float32, ('lat',))
re_dw_mean_o = ncout.createVariable('qnc', np.float32, ('time', 'lev', 'lat','lon'))
tau_dw_mean_o = ncout.createVariable('rccn', np.float32, ('time', 'lev', 'lat','lon'))
lwp_dw_mean_o = ncout.createVariable('geo', np.float32, ('time', 'lev', 'lat', 'lon'))
lon_o[:] = lon[:]
lat_o[:] = lat[:]


re_dw_mean_o[:] = qnc[:]
tau_dw_mean_o[:] = rccn[:]
lwp_dw_mean_o[:] = geo[:]
