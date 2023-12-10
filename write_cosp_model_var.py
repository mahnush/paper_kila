#his codereads the cosp variables from model run and write it in a nc file for each day.

import numpy as np
from netCDF4 import Dataset

ncout = Dataset('/home/b/b380900/icon_nwp_newversion/icon/experiments/icon_lam_lim_kila_1dom/ncfiles_select_var/con_22dec.nc',mode="w",format='NETCDF4_CLASSIC')
#ipath='/work/bb1093/b380900/holouraun_nccn/output_nccn_alpha/new_run_ctt/output_uns/remap_files/'
ipath= '/home/b/b380900/icon_nwp_newversion/icon/experiments/icon_lam_lim_kila_1dom_con/'
fileList=open('/home/b/b380900/icon_nwp_newversion/icon/experiments/icon_lam_lim_kila_1dom/ncfiles_select_var/file_list_22dec.txt','r+')
#fileList=open('/work/bb1093/b380900/holouraun_nccn/output_nccn_alpha/new_run_ctt/output_uns/remap_files/file_list_1sep.txt','r+')
nlon=1801
nlat=1401
ns=2
#ns=24
re_dw=np.zeros((ns,nlat,nlon))
tau_dw=np.zeros((ns,nlat,nlon))
lwp_dw=np.zeros((ns,nlat,nlon))
nd_dw=np.zeros((ns,nlat,nlon))
ctt_dw=np.zeros((ns,nlat,nlon))
ih=0

for FILE_NAME in fileList:
    FILE_NAME=FILE_NAME.strip()
    FILE_NAME=ipath+FILE_NAME
    print(FILE_NAME)
    nc_dw = Dataset(FILE_NAME,'r')
    re_dw[ih,:,:] = nc_dw.variables['modis_Cloud_Particle_Size_Water_Mean'][0,:,:]
    tau_dw[ih,:,:] =nc_dw.variables['modis_Optical_Thickness_Water_Mean'][0,:,:]
    lwp_dw[ih,:,:] = nc_dw.variables['modis_Liquid_Water_Path_Mean'][0,:,:]
    ih+=1
    print(ih)
re_mask = np.ma.masked_where(re_dw==0, re_dw)
nd_dw=1.37e-5*(tau_dw**0.5)*(re_mask**-2.5)*1e-6

lon=nc_dw.variables['lon'][:]
lat=nc_dw.variables['lat'][:]
limit = np.zeros(4)
limit[0] = 6
limit[1] = 34
limit[2] = -178
limit[3] = -142
#fig, (ax0) = plt.subplots(1, 1, figsize=(17, 12))

#postpro_func.visulize_model(ax0,nd,bounds= np.arange(0,1000,10), cbar_label='', titel_figure='',
#                            map_limit=limit)
#plt.savefig('nd_control_23DEC.png')
#plt.show()
ns = 2
ncout.createDimension('lon',nlon)
ncout.createDimension('lat',nlat)
ncout.createDimension('time',ns)
lon_o = ncout.createVariable('lon',np.float32,('lon',))
lat_o= ncout.createVariable('lat',np.float32,('lat',))
re_dw_mean_o = ncout.createVariable('re_dw',np.float32,('time','lat','lon'))
tau_dw_mean_o = ncout.createVariable('tau_dw',np.float32,('time','lat','lon'))
lwp_dw_mean_o= ncout.createVariable('lwp_dw',np.float32,('time', 'lat','lon'))
qnc_dw_mean_o = ncout.createVariable('nd_dw',np.float32,('time', 'lat','lon'))
lon_o[:]=lon[:]
lat_o[:]=lat[:]
print(np.shape(re_dw_mean_o))
re_dw_mean_o[:] = re_dw[:]
tau_dw_mean_o[:] = tau_dw[:]
lwp_dw_mean_o[:] = lwp_dw[:]
qnc_dw_mean_o[:] = nd_dw[:]


