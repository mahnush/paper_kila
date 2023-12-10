import matplotlib.pyplot as plt
import postpro_func
import numpy as np
import numpy.ma as ma
omps_path_22Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_22DEC_icon_domain_smooth.nc'
omps_path_23Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_23DEC_icon_domain_smooth.nc'
omps_path_24Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_24DEC.nc'
#omps_path_25Dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/OMPS_so2_25DEC.nc'

modis_22dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/modis_kila_2km_22Dec.nc'
modis_23dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/modis_kila_2km_23Dec.nc'
modis_24dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/modis_kila_2km_24Dec.nc'
#modis_25dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/modis_kila_25Dec.nc'

model_per_22dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/ncfiles_select_var/per_22dec.nc'
model_per_23dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/ncfiles_select_var/per_23dec.nc'
model_per_24dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/ncfiles_select_var/per_24dec.nc'
#model_per_25dec = '/home/mhaghigh/kilaue/model_kila/ncfiles_onetime/per_25dec.nc'

model_con_22dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/ncfiles_select_var/con_22dec.nc'
model_con_23dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/ncfiles_select_var/con_23dec.nc'
model_con_24dec = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/ncfiles_select_var/con_24dec.nc'
#model_con_25dec = '/home/mhaghigh/kilaue/model_kila/ncfiles_onetime/con_25dec.nc'

omps_data = [omps_path_22Dec, omps_path_23Dec, omps_path_24Dec]#, omps_path_25Dec]
model_per = [model_per_22dec, model_per_23dec, model_per_24dec]#, model_per_25dec]
model_con = [model_con_22dec, model_con_23dec, model_con_24dec]#, model_con_25dec]
modis_data = [modis_22dec, modis_23dec, modis_24dec, modis_24dec]#, modis_25dec]

var_name_model = ['nd_dw', 'lwp_dw', 'tau_dw', 're_dw']
var_name_modis = ['nd', 'lwp', 'tau', 're']
numdays = 3
#fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(30, 20))
#i =0
#for ax in axes.flatten():
def get_in_out(i):
    import scipy.interpolate as sci
    control_var = postpro_func.read_nc(model_con[i], var_name_model[1])
    perturbed_var = postpro_func.read_nc(model_per[i], var_name_model[1])
    control_var = control_var[0, :, :]*1e3
    perturbed_var = perturbed_var[0, :, :]*1e3
    modis_var = postpro_func.read_nc(modis_data[i], var_name_modis[1])
    lat_fine = postpro_func.read_nc(model_con[i], 'lat')
    #lat_fine = lat_fine[450 :950]
    lon_fine = postpro_func.read_nc(model_con[i],'lon')
    so2 = postpro_func.read_nc(omps_data[i], 'so2_PBL')
    so2 = np.transpose(so2)
    #so2 = so2[9 :20, :]
    modis_var = np.transpose(modis_var)
    #modis_var = modis_var[450:950, :]
    lat_so2 = postpro_func.read_nc(omps_data[i], 'lat')
    lon_so2 = postpro_func.read_nc(omps_data[i], 'lon')
    lon_coarse = lon_so2[:, 0]
    lat_coarse = lat_so2[0, :]
    so2_mask = np.ma.filled(so2, fill_value=0)
    f = sci.RectBivariateSpline(lat_coarse, lon_coarse, so2_mask)
    scale_interp = f(lat_fine, lon_fine)
    con_inside = control_var[scale_interp > 1.0]
    con_outside = control_var[scale_interp < 1.0]
    per_inside = perturbed_var[scale_interp > 1.0]
    per_outside = perturbed_var[scale_interp < 1.0]
    modis_inside = modis_var[scale_interp > 1.0]
    modis_outside = modis_var[scale_interp < 1.0]
    con_inside_mask = ma.masked_where(con_inside == 0, con_inside)
    #con_inside_mask = ma.masked_where(con_inside_mask_1 >= 1000, con_inside_mask_1)
    con_inside_mask = con_inside_mask.compressed()
    con_outside_mask = ma.masked_where(con_outside == 0, con_outside)
    con_outside_mask = con_outside_mask.compressed()
    per_inside_mask = ma.masked_where(per_inside == 0, per_inside)
    #per_inside_mask = ma.masked_where(per_inside_mask_1 >= 1000, per_inside_mask_1)
    per_inside_mask = per_inside_mask.compressed()
    per_outside_mask = ma.masked_where(per_outside == 0, per_outside)
    per_outside_mask = per_outside_mask.compressed()
    modis_inside_mask = ma.masked_where(modis_inside <= 0, modis_inside)
    modis_inside_mask = modis_inside_mask.compressed()
    modis_outside_mask = ma.masked_where(modis_outside <= 0., modis_outside)
    modis_outside_mask = modis_outside_mask.compressed()
    return con_inside_mask, con_outside_mask, per_inside_mask, per_outside_mask, modis_inside_mask, modis_outside_mask

allvar_in_per = []
allvar_out_per = []
allvar_in_con = []
allvar_out_con = []
allvar_in_modis = []
allvar_out_modis = []
for i in range (numdays):
    con_inside, con_outside, per_inside, per_outside, modis_inside, modis_outside = get_in_out(i)
    allvar_in_per= np.concatenate((allvar_in_per, per_inside ), axis=0)
    allvar_out_per = np.concatenate((allvar_out_per, per_outside), axis=0)
    allvar_in_con = np.concatenate((allvar_in_con, con_inside), axis=0)
    allvar_out_con = np.concatenate((allvar_out_con, con_outside), axis=0)
    allvar_in_modis = np.concatenate((allvar_in_modis, modis_inside), axis=0)
    allvar_out_modis = np.concatenate((allvar_out_modis, modis_outside), axis=0)


def weight(var):
    #weight = (1 + np.zeros(len(var))) / len(var)
    weight = np.zeros_like(var) + 1. / (var.size)
    return weight


def lable_hist(var):
    median = str(round(np.median(var)))
    mean = str(round(np.mean(var)))
    std = str(np.std(var))
    lable = '('+'median = ' + median + ' ,' + 'mean = ' + mean +')'
    return lable


fig, (axs0, axs1) = plt.subplots(1, 2, figsize = (12, 6))
numbin = np.arange(0., 1000, 10)
font_legend = 12
font_lable = 15
line_width = 2
font_tick = 13
#name = '$ \mathrm{N_d}$ ($\mathrm{cm^{-3}}$)'
name =  '$ \mathrm{LWP}$ ($\mathrm{g m^{-2}}$)'
axs0.hist(allvar_in_per,bins=numbin, weights=weight(allvar_in_per) , histtype='step',
          linewidth=line_width, color='red', label='Volcano '+lable_hist(allvar_in_per))

axs0.hist(allvar_in_con,  bins=numbin, weights=weight(allvar_in_con),  histtype='step',
          linewidth=line_width, color='blue', label='No-Volcano '+lable_hist(allvar_in_con))
axs0.hist(allvar_in_modis, bins=numbin, weights=weight(allvar_in_modis),  histtype='step',
         linewidth=line_width, color='black',  label='MODIS '+lable_hist(allvar_in_modis))
axs0.legend(loc='upper right', fontsize=font_legend, frameon=True)
ticks = np.arange(0, 0.18, 0.02)
#ticks_x = np.arange(0,1200,200)
axs0.set_yticks(ticks)
axs0.tick_params(axis='x', labelsize=font_tick)  # to Set Matplotlib Tick Labels Font Size
axs0.tick_params(axis='y', labelsize=font_tick)
axs0.set_xlabel(name, fontsize=font_lable)
axs0.set_ylabel('Relative Frequency', fontsize=font_lable)
axs1.hist(allvar_out_per, bins=numbin, weights=weight(allvar_out_per),  histtype='step',
          linewidth=line_width, color = 'red', label='Volcano '+lable_hist(allvar_out_per))
axs1.hist(allvar_out_con, bins=numbin,weights=weight(allvar_out_con),  histtype='step',
          linewidth=line_width, color = 'blue', label='No-Volcano '+lable_hist(allvar_out_con))
axs1.hist(allvar_out_modis, bins=numbin, weights=weight(allvar_out_modis), histtype='step',
          linewidth=line_width, color = 'black',label='MODIS ' + lable_hist(allvar_out_modis))
#temp = np.histogram(allre_in, bins=numbin, weights = weight(allre_in) )
#temp_1= temp[0]
#print(sum(temp_1))
axs1.legend(loc='upper right', fontsize=font_legend, frameon=True)
axs1.set_yticks(ticks)
#axs0.set_xticks(ticks_x)
#axs1.set_xticks(ticks_x)
axs1.tick_params(axis='x', labelsize=font_tick)  # to Set Matplotlib Tick Labels Font Size
axs1.tick_params(axis='y', labelsize=font_tick)

axs1.set_yticklabels([])
axs1.set_xlabel(name, fontsize=font_lable)
#axs1.set_ylabel('probability density function', fontsize=font_lable)

axs0.set_title('Inside Plume', fontsize= font_lable)
axs1.set_title('Outside Plume', fontsize=font_lable)
#axs0.annotate('(a)', xy=(-30, 0.079), size=font_lable)
#axs1.annotate('(b)', xy=(-30, 0.0785), size=font_lable)
plt.tight_layout()
axs0.grid(True)
axs1.grid(True)
plt.savefig('pdf_lwp_22_23_24Dec_fig5.png')
plt.savefig('pdf_lwp_22_23_24Dec_fig5.pdf')

plt.show()