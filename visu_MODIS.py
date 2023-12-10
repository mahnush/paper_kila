#this code creates a geogrophical distribusion  MODIS
#different variables can be chossen but here the Nd and LWP are plotted.
import matplotlib.pyplot as plt
import numpy as np
import postpro_func
#necessary inputs
file_name  = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/test2_modis_kila_2km_22Dec.nc'
var_name = ['nd','lwp', 'cf', 're']
var_name_model = ['nd_dw','lwp_dw', 'tau_dw', 're_dw']
titel_var = ['Nd', 'LWP', 'cf', 'reff']
titel_kind = [' (MODIS)', ' (ICON)']
cbar_label = ['$\mathrm{cm^{-3}}$', '$\mathrm{g\,m^{-2}}$', '','$\mathrm{{\mu}m}$']
fs_label = 20
fs_titel = 20
nd_bound = np.arange(0, 201, 1)
re_bound = np.arange(0, 32, 2)
lwp_bound = np.arange(0, 201, 2)
tau_bound = np.arange(0, 100, 4)
limit = np.zeros(4)
#limit[0] = 6
#limit[1] = 34
#limit[2] = -178
#limit[3] = -142
#gridSize = 0.02
limit[0] = 6
limit[1] = 34
limit[2] = -178
limit[3] = -142
gridSize = 0.02
titel = '24 DEC 2020'

print(np.min(postpro_func.read_nc(file_name, var_name[2])*100))
fig, ((ax0,ax1), (ax2,ax3)) = plt.subplots(2, 2, figsize=(17, 12))
fig.suptitle(titel, fontsize=fs_titel)
postpro_func.visulize_sat_revised(ax0, postpro_func.read_nc(file_name, var_name[3]),
                                  postpro_func.read_nc(file_name, 'lat'), postpro_func.read_nc(file_name, 'lon'),
                                  re_bound, cbar_label[3], titel_var[3] + titel_kind[0], limit)
#2end subplot
postpro_func.visulize_sat_revised(ax1, postpro_func.read_nc(file_name, var_name[2]) * 100,
                                  postpro_func.read_nc(file_name, 'lat'), postpro_func.read_nc(file_name, 'lon'),
                                  tau_bound, cbar_label[2], titel_var[2] + titel_kind[0], limit)

#3rd subplot
postpro_func.visulize_sat_revised(ax2, postpro_func.read_nc(file_name, var_name[0]),
                                  postpro_func.read_nc(file_name, 'lat'), postpro_func.read_nc(file_name, 'lon'),
                                  nd_bound, cbar_label[0], titel_var[0] + titel_kind[0], limit)
#4th subplot
postpro_func.visulize_sat_revised(ax3, postpro_func.read_nc(file_name, var_name[1]),
                                  postpro_func.read_nc(file_name, 'lat'), postpro_func.read_nc(file_name, 'lon'),
                                  lwp_bound, cbar_label[1], titel_var[1] + titel_kind[0], limit)
plt.tight_layout()
plt.savefig('ggmodis_kila_2km_23Dec.png')
#plt.savefig('28Dec.pdf')
plt.show()