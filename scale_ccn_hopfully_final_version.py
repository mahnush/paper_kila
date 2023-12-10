import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import postpro_func
import numpy.ma as ma
from scipy.signal import savgol_filter
#this code is written to scale ccn in a way that scales the nccn vertically exponentially reducing
# read in the control NCCN file
data_path = '/home/mhaghigh/kilaue/kila_nccn/sel_icon_remap_nccn_kila.nc'

# open the SO2 file
omps_path_22Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_22DEC_nccn_domain_smooth.nc'
omps_path_23Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_23DEC_nccn_domain_smooth.nc'
omps_path_24Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_24DEC_nccn_domain_smooth.nc'
omps_path_25Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_25DEC_nccn_domain_smooth.nc'
omps_path_26Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_26DEC_nccn_domain_smooth.nc'

dataset_so2 = [omps_path_22Dec, omps_path_23Dec, omps_path_24Dec, omps_path_25Dec, omps_path_26Dec]
#read in OPMS files in a domain taht they are -- lat:10,30 lon = -175,-145--
def so2_data_read_in(omps_path) :
    so2 = postpro_func.read_nc(omps_path, 'so2_PBL')
    lat_so2 = postpro_func.read_nc(omps_path, 'lat')
    lon_so2 = postpro_func.read_nc(omps_path, 'lon')

    so2 = np.transpose(so2)
    lon_so2 = np.transpose(lon_so2)
    lat_so2 = np.transpose(lat_so2)
    lat_so2 = lat_so2[:, 0]
    lon_so2 = lon_so2[0, :]
    so2 = np.ma.filled(so2, fill_value = 0.0)
    return so2, lat_so2, lon_so2


def get_pressure():
    import math
    nlev = 60
    p = np.zeros((nlev))
    scale_lev = np.zeros((nlev))
    hyam = [10, 29.2126712799072, 51.0365734100342, 79.6423835754395,
        115.060134887695, 157.533828735352, 207.681701660156, 266.637420654297,
        336.233856201172, 419.295028686523, 520.134567260742, 644.434539794922,
        798.439300537109, 989.247619628906, 1225.65466308594, 1518.55743408203,
        1881.45709228516, 2331.08129882812, 2888.15515136719, 3578.35656738281,
        4433.5, 5462.36401367188, 6662.32543945312, 8035.84252929688,
        9570.59033203125, 11226.7866210938, 12926.3857421875, 14577.5654296875,
        16099.6401367188, 17432.3291015625, 18536.439453125, 19391.40234375,
        19988.6572265625, 20326.0341796875, 20407.171875, 20240.94140625,
        19840.8662109375, 19224.5400390625, 18413.0537109375, 17430.4130859375,
        16302.9580078125, 15058.7856445312, 13727.1713867188, 12337.9887695312,
        10921.1298828125, 9505.9287109375, 8120.57983398438, 6791.55908203125,
        5543.04663085938, 4396.34582519531, 3369.30493164062, 2475.73815917969,
        1724.84619140625, 1120.63717651367, 661.347671508789, 338.863739013672,
        138.141567230225, 36.6284935474396, 3.68387150764465, 0 ]
    hybm = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        3.79117482225411e-05, 0.000268609241175, 0.00113827553286683,
        0.00344813754782081, 0.0081120147369802, 0.0159103930927813,
        0.0273995194584131, 0.0429057851433754, 0.0626121200621128,
        0.0866042636334896, 0.114848602563143, 0.147203415632248,
        0.183430127799511, 0.223204538226128, 0.266128048300743,
        0.311738923192024, 0.359523519873619, 0.408927544951439,
        0.459367260336876, 0.510240748524666, 0.560939162969589,
        0.610857933759689, 0.659408032894135, 0.706027209758759,
        0.750191211700439, 0.79142501950264, 0.829314172267914,
        0.863515913486481, 0.893770396709442, 0.919912099838257,
        0.941880911588669, 0.959733366966248, 0.973653972148895,
        0.983966410160065, 0.991144776344299, 0.995824784040451, 0.998815059661865]
    p_s = np.zeros((nlev))+1.e5
    p[:] = hybm[:] * p_s[:] +  hyam[:]
    p = p[::-1]
    return (p)


#get the control nccn file in whole glob and 3 time step per day.


def get_control_ccn():
    nccn = postpro_func.read_nc(data_path, 'N_CCN')
    height = postpro_func.read_nc(data_path, 'LH')
    w = postpro_func.read_nc(data_path, 'w')
    lon = postpro_func.read_nc(data_path, 'lon')
    lat = postpro_func.read_nc(data_path, 'lat')
    lev = postpro_func.read_nc(data_path, 'lev')
    time = postpro_func.read_nc(data_path, 'time')
    return nccn, w, height, lat, lon, lev, time


#function that actually interpolates var in coarser lat and lon in to the fine lat and lon


def interpolate_function(coarse_lat, coarse_lon, fine_lat, fine_lon, var):
    import scipy.interpolate as sci
    f = sci.RectBivariateSpline(coarse_lat, coarse_lon, var )
    var_interp = f(fine_lat, fine_lon)
    return var_interp

def interpolate_smooth_so2(omps_path):
    lat_so2_domain_nccn_reso = np.arange(10, 30, 0.75)
    lon_so2_domain_nccn_reso = np.arange(-175, -145, 0.75)
    so2, lat_so2, lon_so2 = so2_data_read_in(omps_path)
    so2_interpolated = interpolate_function(lat_so2, lon_so2, lat_so2_domain_nccn_reso,lon_so2_domain_nccn_reso, so2)
    so2_interpolated_smooth = savgol_filter(so2_interpolated, 5, 2, mode='nearest')
    so2_interpolated_smooth[25:31, :] = 1.
    nlat = np.size(lat_so2_domain_nccn_reso)
    nlon = np.size(lon_so2_domain_nccn_reso)
    return so2_interpolated_smooth, nlat , nlon


def compute_scale_field(omps_path):
    so2 = interpolate_smooth_so2(omps_path)[0]
    nlat_so2 = interpolate_smooth_so2(omps_path)[1]
    nlon_so2 = interpolate_smooth_so2(omps_path)[2]

    so2_back = so2[so2 <= 1.]
    so2_plume = so2[so2 > 1.0]
    mean_plume = np.ma.mean(so2_plume)
    mean_back = np.ma.mean(so2_back)
    scale_fac = np.zeros((nlat_so2, nlon_so2))
    scale_fac_num = np.round(mean_plume/1.0)
    scale_fac = scale_fac + scale_fac_num
    scale_fac[so2 <= 1.0] = 1.0
    return scale_fac


def compute_scale_linear(omps_path):
    p_s = 1e5
    p_2 = np.zeros((27, 40, 60))
    scale_test = np.zeros((27, 40, 60))
    p = get_pressure()
    so2_scale = compute_scale_field(omps_path)
    b = (p_s/so2_scale - p[59]) * (1/(p_s-p[59]))
    m = (1-b)/p_s
    for k in range(60):
        p_2[:,:,k] = m[:,:]*p[k] + b[:,:]
        scale_test[:,:,k] = p_2[:,:,k]* so2_scale[:,:]

    return scale_test


def compute_scale_exponetial(omps_path):
    import math
    p_s = 1e5
    scale_lev = np.zeros((60))
    scale_test = np.zeros((27, 40, 60))+1
    p = get_pressure()
    so2_scale = compute_scale_field(omps_path)
    for k in range(nlev-1):
        exponent = (1-(p_s/p[k]))
        scale_lev[k] = (math.exp(exponent))
        scale_treshold = np.max(so2_scale)*scale_lev[k]
       #print(scale_treshold)
        if scale_treshold >= 1:
          scale_test[:, :, k] = so2_scale[:, :] * scale_lev[k]
       #scale_test[:, :, k] = so2_scale[:, :] * scale_lev[k]
        scale_test[so2_scale<=1] = 1.
    return scale_test


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
scaled_fac_global = np.zeros((nt, nw, nlev, nlat, nlon))+1
end_time = 8
s_time = 0


for t, date in enumerate(dataset_so2):
    s_time = t*8
    for hour in range(s_time, end_time):
        print(hour)
        scal_factor_so2 = compute_scale_field(date)
        #print(np.min(scal_factor_so2))
        #print(np.max(scal_factor_so2))

        for windex in range(nw):
            for k in range(nlev):
                scaled_fac_global[hour, windex, k,  13:40, 6:46] = compute_scale_exponetial(date)[:, :, k]
    end_time = end_time + 8

nccn_control = get_control_ccn()[0]
nccn_perturbed = scaled_fac_global * nccn_control


ncout = Dataset('/home/mhaghigh/kilaue/kila_nccn/OUTPUT_Kila_ICON_scaled_scaled_vertically_exponetial_test.nc', mode="w",
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
N_CCN[:] = nccn_perturbed[:]
z_lon = np.mean(z_avg, axis = 1)
z_mean = np.mean(z_lon,axis = 1)

