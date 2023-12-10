from pyhdf.SD import SD, SDC
import numpy as np
import MODIS_functions
import matplotlib.pyplot as plt
import postpro_func

limit=np.zeros(4)
limit[0]=47.5
limit[1]= 54.5
limit[2]= 4.5
limit[3]= 14.5
gridSize=0.01


#--opening L1 and geo file list
fileList=open('/home/mhaghigh/test_jes/file_list.txt','r+')
fileList_geo=open('/home/mhaghigh/test_jes/g_file_list.txt','r+')
geo_file = [line for line in fileList_geo.readlines()]

latgrid ,longrid = MODIS_functions.grid_coordinate(limit,gridSize)
#--defining the array to put all lat and lon and data
allLat = []
allLon = []


ipath = '/home/mhaghigh/test_jes/'
fileList=open(ipath + 'file_list.txt','r+')
fileList_geo=open(ipath + 'g_file_list.txt','r+')

k=0
#var = 'EV_250_Aggr1km_RefSB'
var = 'EV_1KM_Emissive'
lat = 'Latitude'
lon = 'Longitude'

for FILE_NAME in fileList:
    FILE_NAME = FILE_NAME.strip()
    file_gn = geo_file[k]
    file_gn = file_gn.strip()
    #print(k)
    FILE_NAME = ipath+FILE_NAME
    print(FILE_NAME)
    file_gn = ipath+file_gn
    print(file_gn)
    #set the modis file L1 & geo
    file_i = SD(FILE_NAME, SDC.READ)
    file_g = SD(file_gn, SDC.READ)
    allLat, allLon = MODIS_functions.read_coordinate(file_g)
    allvar = MODIS_functions.read_level1_var(file_i,var)




var_grid = MODIS_functions.grid(limit, gridSize, allvar, allLat, allLon )
print(np.shape(var_grid))
print(np.amax(var_grid))
print(np.amin(var_grid))
print(var_grid)
gr_lat , gr_lon = MODIS_functions.grid_coordinate(limit,gridSize)


bounds = np.arange(0, 10, 0.5)
cbar_label = ''
titel = ''
fig, (ax0) = plt.subplots(1, 1, figsize=(17, 12))
postpro_func.visulize_sat_revised(ax0, var_grid, gr_lat, gr_lon, bounds, cbar_label, titel, limit)
plt.savefig('test2.png')
plt.show()