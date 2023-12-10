import matplotlib.pyplot as plt
import postpro_func
import numpy as np
import numpy.ma as ma

model_per_22dec = '/home/mhaghigh/kilaue/model_kila/ncfiles_onetime/per_22dec.nc'
model_con_22dec = '/home/mhaghigh/kilaue/model_kila/ncfiles_onetime/con_22dec.nc'
modis_22dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/modis_kila_22Dec.nc'
control_var = postpro_func.read_nc(model_con_22dec,'re_dw')
per_var = postpro_func.read_nc(model_per_22dec,'re_dw')
modis_var = postpro_func.read_nc(modis_22dec,'re')
modis_var = np.transpose(modis_var)
lat_modis = postpro_func.read_nc(modis_22dec,'lat')
lat_model = postpro_func.read_nc(model_per_22dec,'lat')
print(lat_modis[:,450:951])
print(lat_model[450:951])
threashod = 4
def get_var_ready(var):
    var= var
    var = var.flatten()
    var =  ma.masked_where(var<=threashod, var)
    var = var.compressed()
    return var
control_var = get_var_ready(control_var)
per_var = get_var_ready(per_var)
modis_var = get_var_ready(modis_var)
print(np.shape(modis_var))
fig, (axs0) = plt.subplots(1, 1)


def weight(var):
    #weight = (1 + np.zeros(len(var))) / len(var)
    weight = np.zeros_like(var) + 1. / (var.size)
    return weight
numbin = np.arange(0,60,1)
line_width = 3
axs0.hist(control_var, bins=numbin, weights=weight(control_var),  histtype='step',
          linewidth=line_width, color = 'red', label='no_Volcano ')

axs0.hist(per_var, bins=numbin, weights=weight(per_var),  histtype='step',
          linewidth=line_width, color = 'blue', label='Volcano ')
axs0.hist(modis_var, bins=numbin, weights=weight(modis_var),  histtype='step',
          linewidth=line_width, color = 'black', label='modis ')
axs0.legend(loc='upper right')
plt.savefig('re_whole.png')
plt.show()