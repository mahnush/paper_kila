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
matplotlib.use('TkAgg')
def read_data(var_path, var_name):
    ds = xr.open_dataset(var_path)
    var = ds[var_name][0, :, :]*1e3
    lon = ds['lon']
    lat = ds['lat']
    return var, lat, lon
def read_multiple_data(files, var_name):
    files_list =glob.glob(files)
    print(files_list)
    ds = xr.open_mfdataset(files_list)
    ds.time
    var = ds[var_name][:,:,:]*1e3
    lon = ds['lon']
    lat = ds['lat']
    mean_var = var.mean(dim={'time'})
    return mean_var, lat, lon

def plot_data_one_panel(var, lon, lat):
    fig = plt.figure(figsize=(12,6))
    ax = plt.axes(projection=ccrs.PlateCarree())
     #-- set projection
    #projection = ccrs.PlateCarree()
    #fig, ax = plt.subplots(figsize=(10,10), subplot_kw=dict(projection=projection))
    ax.coastlines(zorder=5)
    ax.set_extent([-165, -145, 15, 25], crs=ccrs.PlateCarree())
    ax.set_title( title, fontsize=12, fontweight='bold')
    ax.gridlines(draw_labels=True,
                     linewidth=0.5,
                     color='gray',
                     zorder=3,
                     xlocs=range(-180,180,5),
                     ylocs=range(-90,90,5))
    cmap_2 = LinearSegmentedColormap.from_list("",["white", "lightskyblue", "steelblue", "green", "yellowgreen", "yellow",
                                                "gold", "red", "firebrick", "darkred"])
    #-- set contour levels, labels
    varMin, varMax, varInt = 0, 202, 2
    levels = np.arange(varMin, varMax+varInt, varInt)
    nlevs  = levels.size
    labels = ['{:.2f}'.format(x) for x in levels]
    #-- create contour line plot
    cnplot = ax.contourf(lon, lat, var,
                             vmin=varMin,
                             vmax=varMax,
                             cmap=cmap_2,
                             levels=levels,
                             zorder=0,
                             transform=ccrs.PlateCarree())

    #-- add colorbar
    cbar = plt.colorbar(cnplot, orientation='horizontal', pad=0.1, shrink=0.6)
    cbar.set_label('$\mathrm{gm^{-2}}$')
    plt.savefig('test.png')
    return
nrows = 1
ncols =2
fig, ax = plt.subplots(nrows=nrows,ncols=ncols,
                    subplot_kw={'projection': ccrs.PlateCarree()},
                    figsize=(11,8.5))
ax = ax.flatten()
def plot_data_multi_panel(var, lon, lat,i, title):

    ax[i].coastlines(zorder=5)
    ax[i].set_extent([-165, -145, 15, 25], crs=ccrs.PlateCarree())
    ax[i].set_title(title, fontsize=12, fontweight='bold')
    ax[i].gridlines(draw_labels=True,
                     linewidth=0.5,
                     color='gray',
                     zorder=3,
                     xlocs=range(-180,180,5),
                     ylocs=range(-90,90,5))
    cmap_2 = LinearSegmentedColormap.from_list("",["white", "lightskyblue", "steelblue", "green", "yellowgreen", "yellow",
                                                "gold", "red", "firebrick", "darkred"])
    #-- set contour levels, labels
    varMin, varMax, varInt = 0, 202, 2
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
    cbar = plt.colorbar(cnplot, orientation='horizontal', pad=0.1, shrink=0.6)
    cbar.set_label('$\mathrm{g m^{-2}}$')
    return cnplot

files_per =  '/work/bb1093/b380900/experiments/icon_lam_lim_kila_1dom/remap_NWP_LAM_DOM01_20201222T230000Z_0011.nc'
files_con = '/work/bb1093/b380900/experiments/icon_lam_lim_kila_1dom_con/remap_NWP_LAM_DOM01_20201222T230000Z_0011.nc'
var_name = 'modis_Liquid_Water_Path_Mean'
data = [files_per, files_con]

title = ['volcano', 'no-volcano']
for i,exp in enumerate(data):
    var, lat, lon = read_multiple_data(exp, var_name)
    plot_data_multi_panel(var, lon, lat, i ,title[i])

#fig.subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9,
 #                   wspace=0.02, hspace=0.02)

# Add a colorbar axis at the bottom of the graph
#cbar_ax = fig.add_axes([0.2, 0.2, 0.6, 0.02])

# Draw the colorbar
#cbar=fig.colorbar(cnplot, cax=cbar_ax, orientation='horizontal')

#plt.savefig('test_2.png')
plt.show()