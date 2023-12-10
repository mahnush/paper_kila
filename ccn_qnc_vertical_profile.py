import matplotlib.pyplot as plt

import postpro_func
import numpy as np

# open the SO2 file
omps_path_22Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_22DEC_icon_domain_smooth.nc'
omps_path_23Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_23DEC_nccn_domain_smooth.nc'

model_per_22dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/ncfiles_select_var/vars_3dim_per_22.nc'
model_per_23dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/ncfiles_select_var/per_23dec.nc'
model_per_24dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/ncfiles_select_var/per_24dec.nc'
#model_per_25dec = '/home/mhaghigh/kilaue/model_kila/ncfiles_onetime/per_25dec.nc'

model_con_22dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/ncfiles_select_var/vars_3dim_con_22.nc'
model_con_23dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/ncfiles_select_var/con_23dec.nc'
model_con_24dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/ncfiles_select_var/con_24dec.nc'
var_rccn = 'rccn'
var_qnc = 'qnc'
def read_nc(file_name, var_name):
    from netCDF4 import Dataset
    nc = Dataset(file_name, 'r')
    var = nc.variables[var_name][0, :, :, :]
    return var
def read_nc_cordinate(file_name):
    from netCDF4 import Dataset
    file = file_name
    nc = Dataset(file, 'r')
    lat = nc.variables['lat'][:]
    lon =  nc.variables['lon'][:]
    return lat, lon
def get_height(data_path):
    var_name = 'geo'
    g = 9.8
    height = read_nc(data_path, var_name)
    height_lat_mean = np.mean(height, axis = 1)
    height_mean = np.mean(height_lat_mean, axis = 1)
    height_mean = height_mean/g
    height_axis = np.round(height_mean) * 1e-3
    return height_axis


def resd_test():
    ccn_per = read_nc(model_per_22dec, var_rccn)
    ccn_con = read_nc(model_con_22dec, var_rccn)
    qnc_per = read_nc(model_per_22dec, var_qnc)
    qnc_con = read_nc(model_con_22dec, var_qnc)
    ccn_per_mask = np.ma.masked_where(ccn_per <= 0, ccn_per)
    #ccn_per_mask = ccn_per
    ccn_con_mask = np.ma.masked_where(ccn_con <= 0, ccn_con)
    #ccn_con_mask = ccn_con
    qnc_per_mask = np.ma.masked_where(qnc_per <= 0, qnc_per)
    qnc_con_mask = np.ma.masked_where(qnc_con <= 0, qnc_con)
    return ccn_per_mask,ccn_con_mask, qnc_per_mask, qnc_con_mask

def mean_lev(var):
    avg_var = np.zeros((45))
    for ik in range(30,75):
        print(ik)
        iz = ik-30
        temp_var = var[ik, :, :].flatten()
        temp_var = temp_var.compressed()
        avg_var[iz] = np.mean(temp_var, axis = 0)
    return avg_var
height = get_height(model_con_22dec)
height = height[30:75]
def plot_vertically():
    fig, ax = plt.subplots(figsize= (20, 15))
    line_width = 2
    ccn_per_mask, ccn_con_mask, qnc_per_mask, qnc_con_mask = resd_test()
    avg_ccn_per_mask = mean_lev(ccn_per_mask)
    ax.plot(avg_ccn_per_mask, height, linewidth = line_width, label = 'ccn_volcano')
    avg_ccn_con_mask = mean_lev(ccn_con_mask)
    ax.plot(avg_ccn_con_mask , height, linewidth = line_width, label = 'ccn_no_volcano')
    return


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
    lat_fine = read_nc_cordinate(model_con_22dec)[0]
    lon_fine = read_nc_cordinate(model_con_22dec)[1]
    f = sci.RectBivariateSpline(ori_lat_so2, ori_lon_so2, ori_omps)
    scale_interp = f(lat_fine, lon_fine)
    return scale_interp

def var_in_out_plume(omps_path, data_path, var):
    nlev = 75
    avg_var_inside = np.zeros((nlev))
    avg_var_outside = np.zeros((nlev))
    so2 = interpolate_omps(omps_path, data_path)
    for ik in range(nlev):
        temp_con_inside = var[ik, :, :][so2 > 1.0]
        temp_con_inside = temp_con_inside.compressed()
        avg_var_inside[ik] = np.mean(temp_con_inside, axis = 0)
        temp_con_outside = var[ik, :, :][so2 < 1.0]
        temp_con_outside = temp_con_outside.compressed()
        avg_var_outside[ik] = np.mean(temp_con_outside, axis = 0)

    return avg_var_inside, avg_var_outside




def vertical_profile_CCN(omps_path, data_path_per, data_path_con, ax):
    ccn_mean_per = read_nc(data_path_per, var_rccn)
    ccn_mean_per = ccn_mean_per*1e-6
    ccn_mean_per_mask = np.ma.masked_where(ccn_mean_per <=0, ccn_mean_per)
    ccn_mean_con = read_nc(data_path_con, var_rccn)
    ccn_mean_con = ccn_mean_con*1e-6
    ccn_mean_con_mask = np.ma.masked_where(ccn_mean_con <=0, ccn_mean_con)
    qnc_mean_per = read_nc(data_path_per, var_qnc)
    qnc_mean_per = qnc_mean_per*1e-6
    qnc_mean_per_mask = np.ma.masked_where(qnc_mean_per <=0, qnc_mean_per)
    qnc_mean_con = read_nc(data_path_con, var_qnc)
    qnc_mean_con = qnc_mean_con*1e-6
    qnc_mean_con_mask = np.ma.masked_where(qnc_mean_con <=0, qnc_mean_con)
    ccn_in_per = var_in_out_plume(omps_path, data_path_per, ccn_mean_per_mask)[0]
    ccn_out_per = var_in_out_plume(omps_path, data_path_per, ccn_mean_per_mask)[1]
    ccn_in_con = var_in_out_plume(omps_path, data_path_con, ccn_mean_con_mask)[0]
    ccn_out_con = var_in_out_plume(omps_path, data_path_con, ccn_mean_con_mask)[1]
    qnc_in_per = var_in_out_plume(omps_path, data_path_per, qnc_mean_per_mask)[0]
    qnc_out_per = var_in_out_plume(omps_path, data_path_per, qnc_mean_per_mask)[1]
    qnc_in_con = var_in_out_plume(omps_path, data_path_con, qnc_mean_con_mask)[0]
    qnc_out_con = var_in_out_plume(omps_path, data_path_con, qnc_mean_con_mask)[1]
    height = get_height(data_path_con)
    line_width = 2
    fs = 20
    ax.plot(ccn_in_per, height, linewidth = line_width, label = 'ccn_in_volcano')
    #ax.plot(ccn_out_per, height, linewidth = line_width, label = 'ccn_out_volcano')
    ax.plot(ccn_in_con, height, linewidth = line_width, label = 'ccn_in_no_volcano')
    ax.plot(ccn_out_con, height, linewidth = line_width, label = 'ccn_out_no_volcano')
    ax.plot(qnc_in_per, height, linewidth = line_width, label = 'qnc_in_volcano')
    ax.plot(qnc_out_per, height, linewidth = line_width, label = 'qnc_out_volcano')
    ax.plot(qnc_in_con, height, linewidth = line_width, label = 'qnc_in_no_volcano')
    ax.plot(qnc_out_con, height, linewidth = line_width, label = 'qnc_out_no_volcano')
    ax.legend(fontsize = fs)
    ax.tick_params(axis = 'x', labelsize = fs)  # to Set Matplotlib Tick Labels Font Size
    ax.tick_params(axis = 'y', labelsize = fs)
    ax.set_xlabel(" ($\mathrm{cm^{-3}}$)", fontsize = fs)
    ax.set_ylabel("Height (km)", fontsize = fs)
    #ax.annotate(w, xy= (10, 10), fontsize = fs)
    return


#fig, axs0 = plt.subplots(figsize= (20, 15))
#vertical_profile_CCN(omps_path_22Dec, model_per_22dec, model_con_22dec,axs0)
plot_vertically()
plt.show()