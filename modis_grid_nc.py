from pyhdf.SD import SD, SDC
import numpy as np
from netCDF4 import Dataset
import MODIS_functions
import numpy.ma as ma

# import modis_postprocess

# --set the limit for Lat and Lon and grid size
limit = np.zeros(4)
limit[0] = 6
limit[1] = 34
limit[2] = -178
limit[3] = -142
#limit[0] = 15
#limit[1] = 25
#limit[2] = -169
#limit[3] = -149
gridSize = 0.02
gridSize_5km = 0.05

# --opening L2 and geo file list
fileList = open('/home/mhaghigh/kilaue/modis_kil/22dec/file_list.txt', 'r+')
fileList_geo = open('/home/mhaghigh/kilaue/modis_kil/22dec/g_file_list.txt', 'r+')
geo_file = [line for line in fileList_geo.readlines()]

latgrid, longrid, nlon, nlat= MODIS_functions.grid_coordinate(limit, gridSize)
latgrid_5km, longrid_5km, nlon_5km, nlat_5km = MODIS_functions.grid_coordinate(limit, gridSize_5km)
# --defining the array to put all lat and lon and data
allLat = []
allLon = []
allreff = []
alltau = []
allctp = []
alllwp = []
allphase = []
allcf = []
allLat_5km = []
allLon_5km = []
ipath = '/home/mhaghigh/kilaue/modis_kil/22dec/'
# name of variables
reff = 'Cloud_Effective_Radius'
tau = 'Cloud_Optical_Thickness'
ctp = 'cloud_top_temperature_1km'
optical_phase = 'Cloud_Phase_Optical_Properties'
lwp = 'Cloud_Water_Path'
cf = 'Cloud_Fraction'
# --reading the data on data set
k = 0

for FILE_NAME in fileList:
    FILE_NAME = FILE_NAME.strip()
    file_gn = geo_file[k]
    file_gn = file_gn.strip()
    # print(k)
    FILE_NAME = ipath + FILE_NAME
    print(FILE_NAME)
    file_gn = ipath + file_gn
    print(file_gn)
    # set the modis file L2 & geo
    file_i = SD(FILE_NAME, SDC.READ)
    file_g = SD(file_gn, SDC.READ)
    if len(MODIS_functions.read_cloud_var(file_i, reff)) == 0:
        allLat, allLon = MODIS_functions.read_coordinate(file_g)
        allreff = MODIS_functions.read_cloud_var(file_i, reff)
        alltau = MODIS_functions.read_cloud_var(file_i, tau)
        allctp = MODIS_functions.read_cloud_var(file_i, ctp)
        allphase = MODIS_functions.read_cloud_var(file_i, optical_phase)
        alllwp = MODIS_functions.read_cloud_var(file_i, lwp)
        allcf = MODIS_functions.read_cloud_var(file_i, cf)
        allLat_5km, allLon_5km = MODIS_functions.read_coordinate(file_i)
    elif len(MODIS_functions.read_cloud_var(file_i, reff)) > 0:
        allLat = np.concatenate((allLat, MODIS_functions.read_coordinate(file_g)[0]), axis=0)
        allLon = np.concatenate((allLon, MODIS_functions.read_coordinate(file_g)[1]), axis=0)
        allreff = np.concatenate((allreff, MODIS_functions.read_cloud_var(file_i, reff)), axis=0)
        alltau = np.concatenate((alltau, MODIS_functions.read_cloud_var(file_i, tau)), axis=0)
        allctp = np.concatenate((allctp, MODIS_functions.read_cloud_var(file_i, ctp)), axis=0)
        allphase = np.concatenate((allphase, MODIS_functions.read_cloud_var(file_i, optical_phase)), axis=0)
        alllwp = np.concatenate((alllwp, MODIS_functions.read_cloud_var(file_i, lwp)), axis=0)
        allcf = np.concatenate((allcf, MODIS_functions.read_cloud_var(file_i, cf)), axis=0)
        allLat_5km = np.concatenate((allLat_5km, MODIS_functions.read_coordinate(file_i)[0]), axis=0)
        allLon_5km = np.concatenate((allLon_5km, MODIS_functions.read_coordinate(file_i)[1]), axis=0)
    k = k + 1
print(allcf)
#cf_grid = MODIS_functions.grid(limit, gridSize_5km, allcf, allLat_5km, allLon_5km)
cf_grid = MODIS_functions.grid(limit, gridSize, allcf, allLat_5km, allLon_5km)
print(np.min(cf_grid))
reff_grid = MODIS_functions.grid(limit, gridSize, allreff, allLat, allLon)
lwp_grid = MODIS_functions.grid(limit, gridSize, alllwp, allLat, allLon)
phase_grid = MODIS_functions.grid(limit, gridSize, allphase, allLat, allLon)
tau_grid = MODIS_functions.grid(limit, gridSize, alltau, allLat, allLon)
ctp_grid = MODIS_functions.grid(limit, gridSize, allctp, allLat, allLon)

# mask ice phase clouds
reff_grid = ma.masked_where(phase_grid == 3, reff_grid)
lwp_grid = ma.masked_where(phase_grid == 3, lwp_grid)
tau_grid = ma.masked_where(phase_grid == 3, tau_grid)
reff_grid = ma.masked_where(reff_grid < 4, reff_grid)
reff_grid_test = reff_grid[ctp_grid > 273]
print(np.max(reff_grid_test))

# to compute Nd
nd_grid = (1.37e-5 * (tau_grid ** 0.5) * ((reff_grid * 1e-6) ** -2.5)) * 1e-6
nd_grid_mask_cf = ma.masked_where(cf_grid< 0.9, nd_grid)
# plt.savefig('nd_warm.pdf')
# plt.show()
#regrid_cloud_fraction


ncout = Dataset('/home/mhaghigh/kilaue/modis_kil/KV_ncfile/test2_modis_kila_2km_22Dec.nc', mode="w", format='NETCDF4_CLASSIC')

ncout.createDimension('lat', nlat)
ncout.createDimension('lon', nlon)
ncout.createDimension('lat_5km', nlat_5km)
ncout.createDimension('lon_5km', nlon_5km)
lat_o = ncout.createVariable('lat', np.float32, ('lon', 'lat'))
lon_o = ncout.createVariable('lon', np.float32, ('lon', 'lat'))
lat_o_5km = ncout.createVariable('lat_5km', np.float32, ('lon_5km', 'lat_5km'))
lon_o_5km = ncout.createVariable('lon_5km', np.float32, ('lon_5km', 'lat_5km'))
re_o = ncout.createVariable('re', np.float32, ('lon', 'lat'))
tau_o = ncout.createVariable('tau', np.float32, ('lon', 'lat'))
lwp_o = ncout.createVariable('lwp', np.float32, ('lon', 'lat'))
nd_o = ncout.createVariable('nd', np.float32, ('lon', 'lat'))
#cf_o = ncout.createVariable('cf', np.float32, ('lon_5km', 'lat_5km'))
cf_o = ncout.createVariable('cf', np.float32, ('lon', 'lat'))
re_o[:] = reff_grid[:]
tau_o[:] = tau_grid[:]
lwp_o[:] = lwp_grid[:]
nd_o[:] = nd_grid_mask_cf[:]
cf_o[:] = cf_grid[:]
lat_o[:] = latgrid[:]
lon_o[:] = longrid[:]
lat_o_5km[:] = latgrid_5km[:]
lon_o_5km[:] = longrid_5km[:]
