# this code creates vertical profile of nccn inside and outside of plume for con and per files.
import numpy as np
import postpro_func
import matplotlib.pyplot as plt



data_path_control = '/home/mhaghigh/kilaue/kila_nccn/sel_icon_remap_nccn_kila.nc'
#data_path_per = '/home/mhaghigh/kilaue/kila_nccn/OUTPUT_Kila_ICON_scaled_scaled_vertically_exponetial_test.nc'
data_path_per = '/home/mhaghigh/kilaue/kila_nccn/OUTPUT_Kila_ICON_scaled_scaled_vertically_exponetial_final.nc'
# open the SO2
omps_path_22Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_22DEC_nccn_domain_smooth.nc'
omps_path_23Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_23DEC_nccn_domain_smooth.nc'
omps_path_24Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_24DEC_nccn_domain_smooth.nc'
omps_path_25Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_25DEC_nccn_domain_smooth.nc'


def get_ccn(data_path):
    nccn = postpro_func.read_nc(data_path, 'N_CCN') * 1e-6
    height = postpro_func.read_nc(data_path, 'LH')
    w = postpro_func.read_nc(data_path, 'w')
    lon = postpro_func.read_nc(data_path, 'lon')
    lat = postpro_func.read_nc(data_path, 'lat')
    lev = postpro_func.read_nc(data_path, 'lev')
    time = postpro_func.read_nc(data_path, 'time')
    return nccn, w, height, lat, lon, lev, time

def get_specific_time_w( data_path, date, w ) :
    nccn = get_ccn(data_path)[0]
    ts_1 = date*8
    ts_2 = ts_1 + 8
    nccn_daily = np.mean(nccn[ts_1:ts_2, w, :, :, :], axis = 0)
    w_values = get_ccn(data_path)[1]
    return nccn_daily, w_values


def interpolate_omps(omps_path, data_path):
    import scipy.interpolate as sci
    ori_omps = postpro_func.read_nc(omps_path, 'so2_PBL')
    ori_lat_so2 = postpro_func.read_nc(omps_path, 'lat')
    ori_lon_so2 = postpro_func.read_nc(omps_path, 'lon')
    ori_omps = np.transpose(ori_omps)
    ori_lat_so2 = np.transpose(ori_lat_so2)
    ori_lon_so2 = np.transpose(ori_lon_so2)
    ori_omps = np.ma.filled(ori_omps, fill_value = 1)
    ori_lat_so2 = ori_lat_so2[:, 0]
    ori_lon_so2 = ori_lon_so2[0, :]
    lat_fine = get_ccn(data_path)[3]
    lon_fine = get_ccn(data_path)[4]
    f = sci.RectBivariateSpline(ori_lat_so2, ori_lon_so2, ori_omps)
    scale_interp = f(lat_fine, lon_fine)
    return scale_interp


def CCN_in_out_plume(omps_path, data_path, ccn_mean):
    avg_ccn_inside = np.zeros((60))
    avg_ccn_outside = np.zeros((60))
    so2 = interpolate_omps(omps_path, data_path)
    for ik in range(60):
        temp_con_inside = ccn_mean[ik, :, :][so2 > 1.0]
        avg_ccn_inside[ik] = np.mean(temp_con_inside, axis = 0)
        temp_con_outside = ccn_mean[ik, :, :][so2 <= 1.0]
        avg_ccn_outside[ik] = np.mean(temp_con_outside, axis = 0)
    return avg_ccn_inside, avg_ccn_outside


def get_height(data_path):
    height = get_ccn(data_path)[2]
    height_lat_mean = np.mean(height, axis = 1)
    height_mean = np.mean(height_lat_mean, axis = 1)
    height_axis = np.round(height_mean) * 1e-3
    return height_axis


def vertical_profile_CCN(omps_path, data_path_per, data_path_con, date, w, ax):
    ccn_mean_per = get_specific_time_w(data_path_per, date, w)[0]
    ccn_in_per = CCN_in_out_plume(omps_path, data_path_per, ccn_mean_per)[0]
    ccn_out_per = CCN_in_out_plume(omps_path, data_path_per, ccn_mean_per)[1]
    ccn_mean_con = get_specific_time_w(data_path_con, date, w)[0]
    ccn_in_con = CCN_in_out_plume(omps_path, data_path_con, ccn_mean_con)[0]
    ccn_out_con = CCN_in_out_plume(omps_path, data_path_con, ccn_mean_con)[1]
    height = get_height(data_path_con)
    w_values = get_specific_time_w(data_path_per, date, w)[1]
    ccn_in_per = ccn_in_per[0 :30]
    ccn_out_per = ccn_out_per[0 :30]
    ccn_in_con = ccn_in_con[0 :30]
    ccn_out_con = ccn_out_con[0 :30]
    height = height[0 :30]
    line_width = 2
    w = 'w = ' + str(w_values[w]) + ' m/s'
    ax.plot(ccn_in_per, height, linewidth = line_width, label = 'in_volcano')
    ax.plot(ccn_out_per, height, linewidth = line_width, label = 'out_volcano')
    ax.plot(ccn_in_con, height, linewidth = line_width, label = 'in_no_volcano')
    ax.plot(ccn_out_con, height, linewidth = line_width, label = 'out_no_volcano')
    ax.legend(fontsize = fs)
    ax.tick_params(axis = 'x', labelsize = fs)  # to Set Matplotlib Tick Labels Font Size
    ax.tick_params(axis = 'y', labelsize = fs)
    ax.set_xlabel("CCN ($\mathrm{cm^{-3}}$)", fontsize = fs)
    ax.set_ylabel("Height (km)", fontsize = fs)
    #ax.annotate(w, xy= (10, 10), fontsize = fs)
    ax.set_title(w, fontsize = fs, loc =  'left')

    return


def visualize(row, j, updraft, panels, avg_ccn):
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib as mpl
    from matplotlib.colors import LinearSegmentedColormap

    font_size = 20
    cmap = [(0.0, 0.0, 0.0)] + [(cm.jet(i)) for i in range( 1 , 256 )]
    cmap = mpl.colors.ListedColormap( cmap )
    bounds = np.arange(0, 2200, 200)
    norm = mpl.colors.BoundaryNorm( bounds , cmap.N )
    cmap_2 = LinearSegmentedColormap.from_list( "" , ["lightskyblue" , "steelblue" , "green" , "yellowgreen" ,
                                                      "yellow" , "gold" , "red" , "firebrick" , "darkred"] )

    font_size_map = 10
    m = Basemap( ax = axs[row , j] , projection = 'cyl' , llcrnrlat = 0 , urcrnrlat = 40 , llcrnrlon = -180 ,
                 urcrnrlon = -140 , resolution = 'c' )
    m.drawcoastlines( )
    m.drawmeridians( np.arange( -180. , 180. , 10. ) , linewidth = 1.2 , labels = [0 , 0 , 0 , 1] , fontsize = font_size_map,
                     color = 'black' , zorder = 3 , latmax = 90,
                     )
    m.drawparallels(np.arange(0., 85., 10.), linewidth = 1.2 , labels = [1 , 0 , 0 , 0] , fontsize = font_size_map ,
                     color = 'black' , zorder = 3 , latmax = 90,
                     )

    axs[row , j].set_title( updraft , fontsize = font_size , loc = 'center', fontweight="bold")
    axs[row, j].set_title(panels, loc = 'left', fontsize = font_size)

    axs[row , j] = m.imshow( avg_ccn , cmap = cmap_2, norm = norm )
    cax = fig.add_axes([0.15, 0.05, 0.69, 0.02])  # left, bottom, width,height
    cbar_bounds = bounds
    cbar_ticks = cbar_bounds
    cbar = fig.colorbar(axs[row, j], cax = cax, norm = norm, boundaries = cbar_bounds,
                         ticks = cbar_ticks, orientation = 'horizontal')
    cbar.ax.tick_params( labelsize = font_size )
    cbar.set_label( 'CCN ($\mathrm{cm^{-3}}$)' , fontsize = font_size )
    return
fig, axes = plt.subplots(nrows = 5, ncols = 2, figsize = (15, 16))
i = 0
fs = 15
num_step = 0
first_timestep = num_step*8
second_timestep = first_timestep + 8
for ax in axes.flatten():
    # the number 0,1,2,.. indicates the index of omps_data which encounters for the day of simulation
    vertical_profile_CCN(omps_path_22Dec, data_path_per, data_path_control, num_step, i, ax)
    i = i + 1
# plt.tight_layout()
#fig.suptitle(' 22 Dec 2020', fontsize = fs)
# plt.subplots_adjust(left=0.1,bottom=0.2,right=0.9, top=0.8, wspace=0.1, hspace=0.5)
plt.subplots_adjust(top = 0.95, bottom = 0.05)

plt.savefig('ccn_vertical.png')
plt.savefig('ccn_vertical.pdf')
plt.show()
plt.close()
panels_1 = ['perturbed', 'control', 'c', 'd']
panel_1 = ['w = 0.215m/s' , 'w = 4.64 m/s']
nccn_perturbed = get_ccn(data_path_per)[0]
nccn_perturbed_time_mean = np.mean(nccn_perturbed[first_timestep:second_timestep, 3, :, :, :], axis= 0)
nccn_perturbed_lev_mean = np.mean(nccn_perturbed_time_mean[:, :, :], axis= 0)
nccn_control = get_ccn(data_path_control)[0]
nccn_control_time_mean = np.mean(nccn_control[first_timestep:second_timestep, 3, :, :, :], axis= 0)
nccn_control_lev_mean = np.mean(nccn_control_time_mean[:, :, :], axis= 0)

fig , axs = plt.subplots(2, 2, figsize = (12, 11))
visualize(0, 0, panel_1[0], panels_1[0], nccn_perturbed_lev_mean)
visualize(0, 1, panel_1[0], panels_1[1], nccn_control_lev_mean)

nccn_perturbed_time_mean = np.mean(nccn_perturbed[first_timestep:second_timestep, 6, :, :, :], axis= 0)
nccn_perturbed_lev_mean = np.mean(nccn_perturbed_time_mean[:, :, :], axis= 0)
nccn_control_time_mean = np.mean(nccn_control[first_timestep:second_timestep, 6, :, :, :], axis= 0)
nccn_control_lev_mean = np.mean(nccn_control_time_mean[:, :, :], axis= 0)
visualize(1, 0, panel_1[1], panels_1[0], nccn_perturbed_lev_mean)
visualize(1, 1, panel_1[1], panels_1[1], nccn_control_lev_mean)
#fig.suptitle(' 22 Dec 2020', fontsize = fs)
#plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.subplots_adjust(top = 0.95, bottom = 0.1,left = 0.125, right = 0.9, wspace = 0.0002)
#plt.tight_layout()
plt.savefig('ccn_geo.pdf')
plt.savefig('ccn_geo.png')
plt.show()
