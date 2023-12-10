def read_cloud_var(file_i, var_name):

    var_2D = file_i.select(var_name)
    #get scale factor and fill value for reff field
    attributes = var_2D.attributes()
    scale_factor = attributes['scale_factor']
    add_offset = attributes['add_offset']
    fillvalue = attributes['_FillValue']
    range = attributes['valid_range']
    #get variable
    var = var_2D.get()

    #get max and min range
    min_range = min(range)*scale_factor
    max_range = max(range)*scale_factor
    #get the valid var
    valid_var = (var-add_offset)*scale_factor
    #p.where(valid_var <min_range , valid_var,valid_var == fillvalue)
    #np.where(valid_var >max_range, valid_var,valid_var == fillvalue)
    valid_var = valid_var.flatten()
    return valid_var

def read_coordinate(file):
    lat_2D = file.select('Latitude')
    lon_2D=file.select('Longitude')
    lat=lat_2D.get()
    lon=lon_2D.get()
    latitude=lat.flatten()
    longitude=lon.flatten()
    return latitude, longitude

import numpy as np
def grid_coordinate(limit,gridSize ):
    minlat = float(limit[0])
    maxlat = float(limit[1])
    minlon = float(limit[2])
    maxlon = float(limit[3])
    dx = gridSize
    xdim=int(1+((maxlon-minlon)/dx))
    ydim=int(1+((maxlat-minlat)/dx))
    grdlat=np.full([xdim, ydim], -1.0)
    grdlon=np.full([xdim, ydim], -1.0)
    for i in range(xdim):
        for j in range(ydim):
            grdlon[i,j]=dx*i+minlon
            grdlat[i,j]=dx*j+minlat
    nlon = np.size(grdlon[:,0])
    nlat = np.size(grdlon[0,:])

    return grdlat,grdlon, nlon, nlat


def grid(limit, gsize, indata, inlat, inlon):
    dx = gsize
    dy = gsize
    minlat = float(limit[0])
    maxlat = float(limit[1])
    minlon = float(limit[2])
    maxlon = float(limit[3])
    xdim = int(1 + ((maxlon - minlon) / dx))
    ydim = int(1 + ((maxlat - minlat) / dy))
    print(xdim,ydim)
    sum_var = np.zeros((xdim, ydim))
    count = np.zeros((xdim, ydim))
    avg_var = np.full([xdim, ydim], -1.0)
    indata[indata < 0] = 0
    mask_re = np.where(indata != 0, 1, 0)
    for ii in range(len(indata)):

        if (inlat[ii] >= minlat and inlat[ii] <= maxlat and inlon[ii] >= minlon and inlon[ii] <= maxlon):
            i = round((inlon[ii] - minlon) / dx)
            i = int(i)
            j = round((inlat[ii] - minlat) / dy)
            j = int(j)
            sum_var[i, j] = sum_var[i, j] + indata[ii]
            count[i, j] += mask_re[ii]
            #count[i, j] += 1
    count = np.ma.masked_equal(count, 0)
    avg_var = sum_var / count

    avg_var = np.ma.masked_equal(avg_var, -1)

    return (avg_var)

def read_level1_var(file_i, var_name):

    var_2D = file_i.select(var_name)
    #get scale factor and fill value for reff field
    attributes = var_2D.attributes()
    print(attributes)
    scale_factor = attributes['radiance_scales'][9]
    add_offset  = attributes['radiance_offsets'][9]
    #get variable
    var = var_2D.get()
    var = var[9,:,:]
    print(np.shape(var))
    #get the valid var
    valid_var = (var-add_offset)*scale_factor
    valid_var = valid_var.flatten()
    return valid_var