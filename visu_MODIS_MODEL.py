#this code creates a geogrophical distribusion between model and MODIS
#different variables can be chossen but here the Nd and LWP are plotted.
import matplotlib.pyplot as plt
import numpy as np
import postpro_func
#necessary inputs
file_name = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/modis_kila_2km_22Dec.nc'
file_name_model = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/ncfiles_select_var/qnc_22dec.nc'
var_name = ['nd','lwp', 'tau', 're']
var_name_model = ['qnc','lwp_dw', 'tau_dw', 're_dw']
titel_var = ['Nd', 'LWP', 'tau', 'reff']
titel_kind = [' (MODIS)', ' (ICON)']
cbar_label = ['$\mathrm{cm^{-3}}$', '$\mathrm{g\,m^{-2}}$', '','$\mathrm{{\mu}m}$']
fs_label = 20
fs_titel = 20
nd_bound = np.arange(0,202,2)
re_bound = np.arange(0,32,2)
lwp_bound = np.arange(0,201,2)
tau_bound = np.arange(0,110,10)
limit = np.zeros(4)
limit[0] = 6
limit[1] = 34
limit[2] = -178
limit[3] = -142
#gridSize = 0.02
#limit[0] = 15
#limit[1] = 25
#limit[2] = -165
#limit[3] = -145
gridSize = 0.01
titel = '22 DEC 2020'


fig, ((ax0,ax1)) = plt.subplots(2, 1, figsize=(17, 12))
fig.suptitle(titel, fontsize=fs_titel)
nd_modis =  postpro_func.read_nc(file_name, var_name[0])
postpro_func.visulize_sat_revised(ax0, nd_modis, postpro_func.read_nc(file_name, 'lat'),
                                  postpro_func.read_nc(file_name, 'lon'), nd_bound, cbar_label[0],
                                  titel_var[3] + titel_kind[0], limit)
#2end subplot
qnc = postpro_func.read_nc(file_name_model, var_name_model[0])[0, 60:73,:,:]* 1e-6
qnc_mean = np.ma.mean(qnc, axis=0)
postpro_func.visulize_model(ax1, qnc_mean
, nd_bound, cbar_label[0],titel_var[3]+ titel_kind[1], limit)
lat = postpro_func.read_nc(file_name_model,'lat')
#print(lat)
#3rd subplot
#postpro_func.visulize_sat(ax2, postpro_func.read_nc(file_name, var_name[0]),
#postpro_func.read_nc(file_name,'lat') , postpro_func.read_nc(file_name,'lon') , nd_bound, cbar_label[0],
#titel_var[0] + titel_kind[0], limit)
#4th subplot
#postpro_func.visulize_sat(ax3, postpro_func.read_nc(file_name, var_name[3]),
#postpro_func.read_nc(file_name,'lat') , postpro_func.read_nc(file_name,'lon') , re_bound, cbar_label[3],
#titel_figure[3], limit)
#postpro_func.visulize_model(ax3, postpro_func.read_nc(file_name_model, var_name_model[0])
#, nd_bound, cbar_label[0],titel_var[0] + titel_kind[1], limit)
#plt.savefig('nd_re_25_0.5km.png')
#plt.savefig('28Dec.pdf')
plt.show()

re_modis_1dim = nd_modis.flatten()
re_modis_1dim_n = re_modis_1dim[re_modis_1dim >5.]
re_modis_1dim_c = re_modis_1dim_n.compressed()
control_nd_mean_1dim =qnc_mean.flatten()
control_nd_mean_1dim_n = control_nd_mean_1dim[control_nd_mean_1dim > 5.]
control_nd_mean_1dim_c = control_nd_mean_1dim_n.compressed()

def weight(var):
    #weight = (1 + np.zeros(len(var))) / len(var)
    weight = np.zeros_like(var) + 1. / (var.size)
    return weight
def lable_hist(var):
    median = str(np.median(var))
    mean = str((np.mean(var)))

    std = str(np.std(var))
    lable = '('+'mean = ' + mean +')'
    return lable
fig, axs0 = plt.subplots(figsize= (20, 15))
numbin = np.arange(0, 200, 1)
#numbin = np.arange(0, 500, 10)
#numbin = np.arange(0,200,5)
font_tick = 30
font_legend = 30
font_lable = 30
line_width = 4
font_tick = 20
name = '$ \mathrm{N_d}$ ($\mathrm{cm^{-3}}$)'
#name =  '$ \mathrm{LWP}$ ($\mathrm{g m^{-2}}$)'
#name = '\u03BC'+'m'
#name = ''
#axs0.hist(,bins=numbin, weights=weight(var_per_test) , histtype='step',
#         linewidth=line_width, color='red', label='perturbed '+lable_hist(var_per_test))

axs0.hist(control_nd_mean_1dim_c,  bins=numbin, weights=weight(control_nd_mean_1dim_c),  histtype='step',
         linewidth=line_width, color='blue', label='per '+lable_hist(control_nd_mean_1dim_c))
axs0.hist(re_modis_1dim_c, bins=numbin, weights=weight(re_modis_1dim_c),  histtype='step',
         linewidth=line_width, color='black',  label='MODIS '+lable_hist(re_modis_1dim_c))
axs0.legend(loc='upper right', fontsize=font_legend, frameon=True)
#ticks = np.arange(0, 0.12, 0.02)
#ticks_x = np.arange(0,1200,200)
#axs0.set_yticks(ticks)
axs0.tick_params(axis='x', labelsize=font_tick)  # to Set Matplotlib Tick Labels Font Size
axs0.tick_params(axis='y', labelsize=font_tick)
axs0.set_xlabel(name, fontsize=font_lable)
axs0.set_ylabel('Relative Frequency', fontsize=font_lable)

axs0.set_title('Nd', fontsize= font_lable)
# axs0.annotate('(a)',xy=(-30,0.079),size=font_lable)
# plt.tight_layout()
axs0.grid(True)
plt.show()
plt.savefig('old_version_modis_geo.png')