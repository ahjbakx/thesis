"""
Created on Tue Feb 16 09:22:24 2021

@author: aaron
"""

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from cartopy import crs as ccrs
from palettable import scientific as scm


data_dir = '/Users/aaron/thesis/Data/polarisation/'
lats = fits.open(data_dir + 'latitude.fits')[0].data
lons = fits.open(data_dir + 'longitude.fits')[0].data
albedo = fits.open(data_dir + 'Av.fits')[0].data

mask = lats == -99
lats[mask] = np.nan
lons[mask] = np.nan
albedo[mask] = np.nan

levels = np.linspace(np.nanmin(albedo),np.nanmax(albedo),100)
proj = ccrs.Orthographic(central_longitude=0)

# x, y, _ = proj.transform_points(ccrs.PlateCarree(), lons, lats).T
# mask = np.invert(np.logical_or(np.isinf(x), np.isinf(y)))
# x = np.compress(mask, x)
# y = np.compress(mask, y)
# lon = lon[~np.array(mask)]
# lat = lat[~np.array(mask)]
# albedo = albedo[~np.array(mask)]


#%% Plot

ax = plt.axes(projection=proj)
cf = ax.contourf(lons, lats, albedo,
        cmap='Greys_r', levels=levels)
cb = plt.colorbar(cf)
plt.show()
