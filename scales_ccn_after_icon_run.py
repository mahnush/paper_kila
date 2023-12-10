import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import postpro_func
import numpy.ma as ma
from scipy import signal
# read in the control NCCN file
data_path = '/home/mhaghigh/kilaue/2020-12-22_ICON_test.nc'




# open the SO2 file
omps_path_22Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_22DEC_nccn_domain_smooth.nc'
omps_path_23Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_23DEC_nccn_domain.nc'
omps_path_24Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_24DEC_nccn_domain.nc'
omps_path_25Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_25DEC_nccn_domain.nc'
omps_path_26Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_26DEC_nccn_domain.nc'
omps_path_27Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_27DEC_nccn_domain.nc'
dataset_so2 = [omps_path_22Dec, omps_path_23Dec, omps_path_24Dec, omps_path_25Dec, omps_path_26Dec, omps_path_27Dec]
so2_scale = (np.zeros((241, 480)))+1

def so2_scale_factor(omps_path) :
    so2 = postpro_func.read_nc(omps_path, 'so2_PBL')
    lat_so2 = postpro_func.read_nc(omps_path, 'lat')
    lon_so2 = postpro_func.read_nc(omps_path, 'lon')
    so2 = np.transpose(so2)
    lon_so2 = np.transpose(lon_so2)
    lat_so2 = np.transpose(lat_so2)
    so2 = np.ma.filled(so2, fill_value = 0.0)
    so2_back = so2[so2 <= 1.]
    so2_plume = so2[so2 > 1.0]
    #so2[so2<=1] = 0.0
    mean_plume = np.ma.mean(so2_plume)
    mean_back = np.ma.mean(so2_back)
    scale_number = mean_plume/mean_back

    return so2, lat_so2, lon_so2, scale_number
def get_control_ccn():
    nccn = postpro_func.read_nc(data_path, 'N_CCN')
    height = postpro_func.read_nc(data_path, 'LH')
    w = postpro_func.read_nc(data_path, 'w')
    lon = postpro_func.read_nc(data_path, 'lon')
    lat = postpro_func.read_nc(data_path, 'lat')
    lev = postpro_func.read_nc(data_path, 'lev')
    time = postpro_func.read_nc(data_path, 'time')
    return nccn, w, height, lat, lon, lev, time
def interpolate_scal_factor(omps_path):
    import scipy.interpolate as sci
    ori_so2 = so2_scale_factor(omps_path)[0]
    ori_lat_so2 = so2_scale_factor(omps_path)[1]
    ori_lon_so2 = so2_scale_factor(omps_path)[2]
    ori_lat_so2 = ori_lat_so2[:, 0]
    ori_lon_so2 = ori_lon_so2[0, :]
    lat_fine = np.arange(10,30,0.75)
    lon_fine = np.arange(-175,-145,0.75)
    f = sci.RectBivariateSpline(ori_lat_so2, ori_lon_so2, ori_so2 )
    scale_interp = f(lat_fine, lon_fine)
    return scale_interp
def visualize(row, j, updraft, panels, avg_ccn):
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib as mpl
    from matplotlib.colors import LinearSegmentedColormap

    font_size = 30
    cmap = [(0.0, 0.0, 0.0)] + [(cm.jet(i)) for i in range( 1 , 256 )]
    cmap = mpl.colors.ListedColormap( cmap )
    bounds = np.arange(0, 1000, 100)
    norm = mpl.colors.BoundaryNorm( bounds , cmap.N )
    cmap_2 = LinearSegmentedColormap.from_list( "" , ["lightskyblue" , "steelblue" , "green" , "yellowgreen" ,
                                                      "yellow" , "gold" , "red" , "firebrick" , "darkred"] )

    m = Basemap( ax = axs[row , j] , projection = 'cyl' , llcrnrlat = -90 , urcrnrlat = 90 , llcrnrlon = -180 ,
                 urcrnrlon = 180 , resolution = 'c' )
    m.drawcoastlines( )
    #m.drawmeridians( np.arange( -180. , 180. , 5. ) , linewidth = 1.2 , labels = [0 , 0 , 0 , 1] , fontsize = 20 ,
    #                 color = 'black' , zorder = 3 , latmax = 90 ,
    #                 )
    #m.drawparallels( np.arange( 0. , 85. , 5. ) , linewidth = 1.2 , labels = [1 , 0 , 0 , 0] , fontsize = 20 ,
    #                 color = 'black' , zorder = 3 , latmax = 90 ,
    #                 )

    axs[row , j].set_title( updraft , fontsize = font_size , loc = 'center' )
    axs[row , j].set_title( panels , loc = 'left' , fontsize = font_size )

    axs[row , j] = m.imshow( avg_ccn , cmap = cmap_2 , norm = norm )
    cax = fig.add_axes( [0.15 , 0.05 , 0.69 , 0.02] )  # left, bottom, width,height
    cbar_bounds = bounds
    cbar_ticks = cbar_bounds
    cbar = fig.colorbar( axs[row , j] , cax = cax , norm = norm , boundaries = cbar_bounds ,
                         ticks = cbar_ticks , orientation = 'horizontal' )
    cbar.ax.tick_params( labelsize = font_size )
    cbar.set_label( 'CCN ($\mathrm{cm^{-3}}$)' , fontsize = font_size )

    return
nt = np.size(get_control_ccn()[6])
nw = np.size(get_control_ccn()[1])
nlev = np.size(get_control_ccn()[5])
nlat = 27
nlon = 40
scaled_ccn = np.zeros((nt, nw, nlev, nlat, nlon))
scaled_ccn_f = np.zeros((nt, nw, nlev, nlat, nlon))
scaled_fac =  np.zeros((nt, nw, nlev, nlat, nlon))+1
for t, date in enumerate(dataset_so2):
#    nccn_control = get_control_ccn()[0]
#    nccn_control = nccn_control[t, :, :, :, :]

    so2_interpolated = interpolate_scal_factor(date)

    scale_number = round(so2_scale_factor(date)[3])

 #   scaled_ccn[t, :, :, :, :] = np.where(so2_interpolated < 1, nccn_control, nccn_control*scale_number)
    scaled_fac[t, :, 0:17, :, :] = np.where(so2_interpolated < 1.5, scaled_fac[t, :, 0:17, :, :], scaled_fac[t, :, 0:17, :, :]*scale_number)
    scaled_fac[t, :, 17:30, :, :] = np.where(so2_interpolated < 1.5, scaled_fac[t, :, 17:30, :, :], scaled_fac[t, :, 17:30, :, :]*2.0)
    scaled_fac[t, :, 31:60, :, :] = np.where(so2_interpolated < 1.5, scaled_fac[t, :, 31:60, :, :], scaled_fac[t, :, 31:60, :, :]*2.0)
nccn_petrurbed = np.zeros((8, 10, 60, 241, 480))
so2_scale = (np.zeros((6, 10, 60, 241, 480)))+1
for t, date in enumerate(dataset_so2):
    so2_scale[t, :, :, 133:160, 6:46] = scaled_fac[t, :, :, :, :]
nccn_control = get_control_ccn()[0]
print(np.shape(nccn_control))

nccn_petrurbed[:, :, :, :, :] = nccn_control[:, :, :, :, :]*so2_scale[0, :, :, :, :]
#nccn_petrurbed=signal.savgol_filter(nccn_petrurbed,53,3)
w = get_control_ccn()[1]
z_avg = get_control_ccn()[2]
lat = get_control_ccn()[3]
lon = get_control_ccn()[4]
lev = get_control_ccn()[5]
time = get_control_ccn()[6]
nt = np.size(time)
nlev = np.size(lev)
nlon = np.size(lon)
nlat = np.size(lat)
nw = np.size(w)

ncout = Dataset('/home/mhaghigh/kilaue/kila_nccn/per_OUTPUT_Kila_ICON_after_run.nc', mode="w",
                format='NETCDF4_CLASSIC')
ncout.description = 'ICON_CCN_file'
ncout.createDimension('lon', nlon)
ncout.createDimension('lat', nlat)
ncout.createDimension('lev', nlev)
ncout.createDimension('w', nw)
ncout.createDimension('time', nt)
lon_o = ncout.createVariable('lon', np.float32, ('lon',))
lat_o = ncout.createVariable('lat', np.float32, ('lat',))
lev_o = ncout.createVariable('lev', np.float32, ('lev',))
W = ncout.createVariable('w', np.float32, ('w',))
time_o = ncout.createVariable('time', np.float32, ('time',))
N_CCN = ncout.createVariable('N_CCN', np.float32, ('time', 'w', 'lev', 'lat', 'lon'))
LH = ncout.createVariable('LH', np.float32, ('lev', 'lat', 'lon'))
lon_o.long_name = 'longitude'
lat_o.long_name = 'lattitude'
time_o.long_name = 'time'
N_CCN.long_name = 'CCN number concentration'
W.long_name = 'vertical velocity'
LH.long_name = 'top height of grid cell'
lev_o.long_name = 'level number'
lat_o.units = 'degrees_north'
lon_o.units = 'degrees_east'
time_o.units = 'hours since 2020-12-22T00:00:00'
LH.units = 'meters'
W.units = 'meters per second'
N_CCN.units = 'particles per cubic meter'


lon_o[:] = lon[:]
lat_o[:] = lat[:]
lev_o[:] = lev[:]
W[:] = w[:]
time_o[:] = time[:]
LH[:] = z_avg[:]
N_CCN[:] = nccn_petrurbed[:]
#nccn_perturbed_lev_mean_22dec = np.mean(nccn_petrurbed[3, :, :, :], axis= 0)*1e-6
#print(nccn_perturbed_lev_mean_22dec)
#nccn_control = get_control_ccn()[0]
#nccn_control_lev_mean_22dec = np.mean(nccn_control[ 3, :, :, :], axis= 0)*1e-6
##panels_1 = ['perturbed', 'control', 'c', 'd']
#panel_1 = ['w = 0.215m/s', 'w = 4.64 m/s']
#fig , axs = plt.subplots(2, 2, figsize = (30, 20))
#axs.flatten()
#visualize(0, 0, panel_1[0], panels_1[0], nccn_perturbed_lev_mean_22dec)
#visualize(0, 1, panel_1[0], panels_1[1], nccn_control_lev_mean_22dec)
#import cartopy.crs as ccrs
#lat = get_control_ccn()[3]
#lon = get_control_ccn()[4]
#ax = plt.axes(projection=ccrs.PlateCarree())
#cnplot = ax.contourf(lon, lat,nccn_perturbed_lev_mean_22dec ,
#                             vmin=0,
##                             vmax=1000,
#                             zorder=0,
#                             transform=ccrs.PlateCarree())
#cbar = plt.colorbar(cnplot, orientation='horizontal', pad=0.1, shrink=0.8)
#ax1 = plt.axes(projection=ccrs.PlateCarree())
#cnplot2 = ax1.contourf(lon, lat,  nccn_control_lev_mean_22dec,
#                             vmin=0,
#                             vmax=1000,
#                             zorder=0,
 #                            transform=ccrs.PlateCarree())
#cbar = plt.colorbar(cnplot2, orientation='horizontal', pad=0.1, shrink=0.8)
#nccn_perturbed_lev_mean_22dec.plot.contourf(ax=ax, vmin=0, vmax=1000, transform=ccrs.PlateCarree(), cbar_kwargs={'shrink': 0.7})
#ax.coastlines()
#nccn_perturbed_lev_mean_22dec = np.mean(scaled_ccn[1, 6, :, :, :], axis= 0)*1e-6
#nccn_control = get_control_ccn()[0]
#nccn_control_lev_mean_22dec = np.mean(nccn_control[1, 6, :, :, :], axis= 0)*1e-6
#visualize(1, 0, panel_1[1], panels_1[0], nccn_perturbed_lev_mean_22dec)
#visualize(1, 1, panel_1[1], panels_1[1], nccn_control_lev_mean_22dec)
#plt.savefig('test_per_con_nccn_23Dec.png')
#fig, (axs0) = plt.subplots(1, 1)
#nccn_perturbed_lev_mean_22dec = nccn_perturbed_lev_mean_22dec.flatten()
#nccn_control_lev_mean_22dec = nccn_control_lev_mean_22dec.flatten()
#numbin = np.arange(0, 20000, 100)
#def weight(var):
    #weight = (1 + np.zeros(len(var))) / len(var)
#    weight = np.zeros_like(var) + 1. / (var.size)
#    return weight

#axs0.hist(nccn_perturbed_lev_mean_22dec,bins = numbin,weights = weight(nccn_perturbed_lev_mean_22dec),histtype= 'step' ,log= 'TRUE')
#axs0.hist(nccn_control_lev_mean_22dec, bins = numbin, weights = weight(nccn_control_lev_mean_22dec), histtype= 'step', log = 'TRUE')
#plt.savefig('nccn_23dec')
plt.show()


