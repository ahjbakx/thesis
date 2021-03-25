#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 11:32:33 2021

@author: aaron
"""

import pandas as pd
import numpy as np
import pyshtools as pysh
from cartopy import crs as ccrs
from palettable import scientific as scm
import shapely
import geopandas as gpd
from shapely.geometry import Point

sf = gpd.read_file("/Users/aaron/thesis/Data/mare_shape/LROC_GLOBAL_MARE_180.shp")
# tolerance = 1.
# geoms = sf.simplify(tolerance, preserve_topology=True)

geoms = sf['geometry']
polies = []
for geo in geoms:
    if type(geo) == shapely.geometry.MultiPolygon:
        for poly in list(geo):
            # if poly.area>100:
            polies.append(poly)
    else:
        # if geo.area>100:
        polies.append(geo)


df = pd.read_excel(r'/Users/aaron/thesis/Data/mare_basalts/jgre20523-sup-0001-supinfo.xlsx')

lats = df['lat'].to_numpy()[1:]
lons = df['lon'].to_numpy()[1:]

Tb = df['tb'].to_numpy()[1:]

latrange = [90, -90]
lonrange = [-180, 180] 
gridres = 1
latgrid = np.arange(latrange[0], latrange[1]-gridres, -gridres)
longrid = np.arange(lonrange[0], lonrange[1]+gridres, gridres)


#%% Fill
from scipy import interpolate

basalt_thickness = np.empty((np.int((np.sum(np.abs(latrange))/gridres+1)), np.int((np.sum(np.abs(lonrange))/gridres+1))))
basalt_thickness[:] = np.nan

for px in range(len(Tb)):

    lat_index = np.where( latgrid==np.round(lats[px]) )[0][0]
    lon_index = np.where( longrid==np.round(lons[px]) )[0][0]
    
    local_thickness = np.float64(Tb[px]) 
    
    if local_thickness == 0.0:
        print(local_thickness)
        continue
    else:
        basalt_thickness[lat_index, lon_index] = local_thickness

array = np.ma.masked_invalid( basalt_thickness )
xx, yy = np.meshgrid(longrid, latgrid)
#get only the valid values
x1 = xx[~array.mask]
y1 = yy[~array.mask]
newarr = array[~array.mask]
basalt_thickness = interpolate.griddata((x1, y1), 
                                    newarr.ravel(),
                                    (xx, yy),
                                    method='nearest')

for lat in latgrid:
    print(lat)
    for lon in longrid:
        contains=False
        p = Point(lon, lat)
        for poly in polies:
            if poly.contains(p):
                contains=True
                continue
        if not contains:
            lat_index = np.where(latgrid==lat)[0][0]
            lon_index = np.where(longrid==lon)[0][0]
            basalt_thickness[lat_index, lon_index]= 0
    
#np.save("/Users/aaron/thesis/Data/mare_basalts/basalt-thickness", basalt_thickness)

#%% Plot map
pysh.utils.figstyle(rel_width=0.5)

basalt_thickness_grid=pysh.SHGrid.from_array(basalt_thickness/1000)

# fig1, ax1 = basalt_thickness_grid.plot(projection=ccrs.Mollweide(central_longitude=90.),
#                         cmap = scm.sequential.Oslo_20.mpl_colormap,
#                         colorbar='bottom',
#                         cb_triangles='both',
#                         cmap_limits = [0, 5000],
#                         # cb_tick_interval = 150,
#                         # cb_minor_tick_interval=50,
#                         cb_label = 'Basalt thickness, m',
#                         cmap_reverse=True,
#                         grid = True,
#                         show = False)

fig2, ax2 = basalt_thickness_grid.plot(projection=ccrs.Orthographic(central_longitude=180.),
                        cmap = scm.sequential.Oslo_20.mpl_colormap,
                        colorbar='bottom',
                        cb_triangles='max',
                        cmap_limits = [-.0001, 3],
                        cb_tick_interval = 1,
                        cb_minor_tick_interval=.5,
                        cb_label = 'Basalt thickness, km',
                        cmap_reverse=True,
                        grid = True,
                        show = False)
ax2.add_geometries(geoms=polies, crs=ccrs.PlateCarree(central_longitude=-180.),
                    linewidth=0.1, edgecolor='black', facecolor='none')

#fig2.savefig("/Users/aaron/thesis/Figures/WP4/basalt-thickness.png", format='png', dpi=300)
    