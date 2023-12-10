#this code creates a geogrophical disribution of difference between control and perturbed runs.

import matplotlib.pyplot as plt
import numpy as np
import postpro_func
model_per_22dec = '/home/mhaghigh/kilaue/model_kila/ncfiles_onetime/per_22dec.nc'
model_per_23dec = '/home/mhaghigh/kilaue/model_kila/ncfiles_onetime/per_23dec.nc'
model_per_24dec = '/home/mhaghigh/kilaue/model_kila/ncfiles_onetime/per_24dec.nc'
model_per_25dec = '/home/mhaghigh/kilaue/model_kila/ncfiles_onetime/per_25dec.nc'
model_per = [model_per_22dec, model_per_23dec, model_per_24dec, model_per_25dec]
model_con_22dec = '/home/mhaghigh/kilaue/model_kila/ncfiles_onetime/con_22dec.nc'
model_con_23dec = '/home/mhaghigh/kilaue/model_kila/ncfiles_onetime/con_23dec.nc'
model_con_24dec = '/home/mhaghigh/kilaue/model_kila/ncfiles_onetime/con_24dec.nc'
model_con_25dec = '/home/mhaghigh/kilaue/model_kila/ncfiles_onetime/con_25dec.nc'
model_con = [model_con_22dec, model_con_23dec, model_con_24dec, model_con_25dec]
var_name_model = ['nd_dw','lwp_dw', 'tau_dw', 're_dw']
titel_var = ['Nd', 'LWP', 'tau', 'reff']
titel_day = ['22 Dec ', ' 23 Dec', '24 Dec', '25 Dec']
cbar_label = ['$\mathrm{cm^{-3}}$', '$\mathrm{g\,m^{-2}}$', '','$\mathrm{{\mu}m}$']
fs_label = 30
fs_titel = 30
nd_bound = np.arange(-200,204,4)
re_bound = np.arange(1,32,2)
lwp_bound = np.arange(-200,204,4)
tau_bound = np.arange(1,110,10)
limit = np.zeros(4)
limit[0] = 6
limit[1] = 34
limit[2] = -178
limit[3] = -142
gridSize = 0.02
i = 0
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(30, 20))

for ax in axes.flatten():
    control_var = postpro_func.read_nc(model_con[i], var_name_model[0])
    perturbed_var = postpro_func.read_nc(model_per[i], var_name_model[0])
    diff_per_con = perturbed_var - control_var
    postpro_func.visulize_model_diff(ax, diff_per_con, nd_bound, cbar_label[0], titel_day[i], limit)
    i =i+1
fig.suptitle(' Nd' , fontsize=fs_titel)
plt.tight_layout()
plt.savefig('different_Nd.png')
plt.show()