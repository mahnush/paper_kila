import xarray as xr
import numpy as np
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LinearSegmentedColormap
import glob
import matplotlib
#matplotlib.use('TkAgg')
#re_name =  '\u03BC'+'m'
#lwp_name = '$ \mathrm{LWP}$ ($\mathrm{g m^{-2}}$)'
nd_name = '$ \mathrm{N_d}$ ($\mathrm{cm^{-3}}$)'
def read_data(var_path, var_name):
    ds = xr.open_dataset(var_path)
    var = ds[var_name][0, :, :]
    lon = ds['lon']
    lat = ds['lat']
    return var, lat, lon
def read_multiple_data(files, var_name):
    files_list =glob.glob(files)
    print(files_list)
    ds = xr.open_mfdataset(files_list)
    ds.time
    var = ds[var_name][:, :, :]
    lon = ds['lon']
    lat = ds['lat']
    mean_var = var.mean(dim={'time'})
    return mean_var, lat, lon

def plot_data_one_panel(var, lat, lon):
    fig = plt.figure(figsize=(12,6))
    ax = plt.axes(projection=ccrs.PlateCarree())
     #-- set projection
    #projection = ccrs.PlateCarree()
    #fig, ax = plt.subplots(figsize=(10,10), subplot_kw=dict(projection=projection))
    ax.coastlines(zorder=5)
    ax.set_extent([-178, -142, 6, 34], crs=ccrs.PlateCarree())
    ax.set_title( title, fontsize=12, fontweight='bold')
    ax.gridlines(draw_labels=True,
                     linewidth=0.5,
                     color='gray',
                     zorder=3,
                     xlocs=range(-180, 180, 5),
                     ylocs=range(-90, 90, 5))
    cmap_2 = LinearSegmentedColormap.from_list("",["white", "lightskyblue", "steelblue", "green", "yellowgreen", "yellow",
                                                "gold", "red", "firebrick", "darkred"])
    #-- set contour levels, labels
    varMin, varMax, varInt = -500, 500, 100
    levels = np.arange(varMin, varMax+varInt, varInt)
    nlevs  = levels.size
    labels = ['{:.2f}'.format(x) for x in levels]
    #-- create contour line plot
    cnplot = ax.contourf(lon, lat, var,
                             vmin=varMin,
                             vmax=varMax,
                             cmap='seismic',
                             levels=levels,
                             zorder=0,
                             transform=ccrs.PlateCarree())

    #-- add colorbar
    cbar = plt.colorbar(cnplot, orientation='horizontal', pad=0.1, shrink=0.6)
    cbar.set_label(nd_name)
    return
nrows = 1
ncols =2
fig, ax = plt.subplots(nrows=nrows, ncols=ncols,
                    subplot_kw = {'projection': ccrs.PlateCarree()},
                    figsize = (11, 8.5))
ax = ax.flatten()
def plot_data_multi_panel(var, lon, lat,i, title):

    ax[i].coastlines(zorder=5)
    ax[i].set_extent([-178, -142, 6, 34], crs=ccrs.PlateCarree())
    ax[i].set_title(title, fontsize=12, fontweight='bold')
    gl = ax[i].gridlines(draw_labels=True,
                     linewidth=0.5,
                     color='gray',
                     zorder=3,
                     xlocs=range(-180,180,5),
                     ylocs=range(-90,90,5))
    gl.ylabels_left = False
    gl.xlabels_top = False
    cmap_2 = LinearSegmentedColormap.from_list("",["white", "lightskyblue", "steelblue", "green", "yellowgreen", "yellow",
                                                "gold", "red", "firebrick", "darkred"])
    #-- set contour levels, labels
    varMin, varMax, varInt = 0, 200, 2
    levels = np.arange(varMin, varMax+varInt, varInt)
    nlevs  = levels.size
    labels = ['{:.2f}'.format(x) for x in levels]
    #-- create contour line plot
    cnplot = ax[i].contourf(lon, lat, var,
                             vmin=varMin,
                             vmax=varMax,
                             cmap=cmap_2,
                             levels=levels,
                             zorder=0,
                             transform=ccrs.PlateCarree())

    #-- add colorbar
    #cbar =plt.colorbar(cnplot, orientation='horizontal', pad=0.1, shrink=0.8)
    cbar = fig.colorbar(cnplot, ax= ax[i], orientation='horizontal', pad=0.1, shrink=0.8)
    cbar.set_label(nd_name)
    return cnplot

files_per =  '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/ncfiles_select_var/per_22dec.nc'
files_con = '/home/mhaghigh/kilaue/modis_kil/KV_ncfile/ncfiles_select_var/con_22dec.nc'
var_name = 'nd_dw'
data = [files_per, files_con]

title = ['volcano', 'no-volcano']

for i, exp in enumerate(data):
    var, lat, lon = read_data(exp, var_name)
    var = var
    var = np.ma.masked_where(var<1, var)
    plot_data_multi_panel(var, lon, lat, i ,title[i])

plt.savefig('nd_22dec.png')
plt.show()
var_con, lat, lon = read_data(files_con, var_name)
var_per, lat, lon = read_data(files_per, var_name)
var_con = var_con
var_per = var_per
var_diff = (var_per-var_con)
plot_data_one_panel(var_diff, lat, lon)
plt.savefig('nd_diff_22dec.png')
plt.show()