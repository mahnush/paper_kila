#this code did not work for me

from pyhdf.SD import SD, SDC
import numpy as np
from netCDF4 import Dataset
import grid_func, h5py
import numpy.ma as ma
import visu_MODIS_MODEL
import matplotlib.pyplot as plt
import modis_grid_nc
#--set the limit for Lat and Lon and grid size
limit=np.zeros(4)
limit[0]=6
limit[1]=34
limit[2]=-178
limit[3]=-142
gridSize=0.02

fileList=open('/home/mhaghigh/kilaue/modis_kil/24dec/file_list.txt','r+')
fileList_geo=open('/home/mhaghigh/kilaue/modis_kil/24dec/g_file_list.txt','r+')
geo_file = [line for line in fileList_geo.readlines()]

#--defining the array to put all lat and lon and data
allLat = []
allLon = []
allreff = []
alltau = []
allctp = []
alllwp = []
allphase = []

ipath='/home/mhaghigh/kilaue/modis_kil/24dec/'
# name of variables
reff = 'Cloud_Effective_Radius'
tau = 'Cloud_Optical_Thickness'
ctp = 'cloud_top_temperature_1km'
optical_phase ='Cloud_Phase_Optical_Properties'
lwp = 'Cloud_Water_Path'
#--reading the data on data set
k=0

for FILE_NAME in fileList:
    FILE_NAME = FILE_NAME.strip()
    file_gn = geo_file[k]
    file_gn = file_gn.strip()
    #print(k)
    FILE_NAME = ipath+FILE_NAME
    print(FILE_NAME)
    file_gn = ipath+file_gn
    print(file_gn)
    #set the modis file L2 & geo
    file_i = SD(FILE_NAME, SDC.READ)
    file_g = SD(file_gn, SDC.READ)
    if len(modis_grid_nc.read_cloud_var(reff)) == 0:
        allLat, allLon = modis_grid_nc.read_coordinate()
        allreff = modis_grid_nc.read_cloud_var(reff)
        alltau = modis_grid_nc.read_cloud_var(tau)
        allctp = modis_grid_nc.read_cloud_var(ctp)
        allphase = modis_grid_nc.read_cloud_var(optical_phase)
        alllwp = modis_grid_nc. read_cloud_var(lwp)
    elif len(modis_grid_nc.read_cloud_var(reff)) > 0:
       allLat = np.concatenate((allLat, modis_grid_nc.read_coordinate()[0]), axis=0)
       allLon = np.concatenate((allLon, modis_grid_nc.read_coordinate()[1]), axis=0)
       allreff = np.concatenate((allreff, modis_grid_nc.read_cloud_var(reff )), axis=0)
       alltau = np.concatenate((alltau, modis_grid_nc.read_cloud_var(tau)), axis=0)
       allctp = np.concatenate((allctp, modis_grid_nc.read_cloud_var(ctp)), axis=0)
       allphase = np.concatenate((allphase, modis_grid_nc.read_cloud_var(optical_phase)), axis=0)
       alllwp = np.concatenate((alllwp, modis_grid_nc.read_cloud_var(lwp)), axis=0)
    k=k+1

bins = {}
bins['lat'] = np.arange(6, 34.0, 0.02)
bins['lon'] = np.arange(-178, -141.98, 0.02)

allreff[allreff<0.] = 0.
reff_test = np.histogram2d(allLon, allLat, bins=[bins['lon'], bins['lat']], weights=allreff)[0]
coordinate_test = np.histogram2d(allLon, allLat,  bins=[bins['lon'], bins['lat']])[0]
#reff_test = ma.masked_invalid(reff_test)
var_reff = reff_test/coordinate_test
#----
x = np.arange(-178, -142, 0.02)
y = np.arange(6, 34, 0.02)
xv, yv = np.meshgrid(y, x)
print(np.shape(yv))
re_bound = np.arange(0,32,2)
cbar_label = '$\mathrm{{\mu}m}$'
titel_figure = 'gg'
#var_reff = ma.filled(var_reff,fill_value=0)
#var_reff = ma.masked_invalid(var_reff)
fig, (ax0) = plt.subplots(1, 1, figsize=(17, 12))
print(np.shape(var_reff))
print(var_reff)
#var_reff =var_reff[0:1800, 0:1400]
#print(np.shape(var_reff))
#modis_postprocess.visulize(ax0, var_reff,yv, yv, re_bound ,cbar_label,
#             titel_figure, limit)

#plt.savefig('ss.png')
#plt.show()