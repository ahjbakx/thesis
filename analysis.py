#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 10:28:39 2021

@author: aaron
"""

import pyshtools as pysh
import numpy as np
from cartopy import crs as ccrs
from palettable import scientific as scm
from scipy import interpolate
from matplotlib import pyplot as plt

# Width of image with respect to (journal) page
pysh.utils.figstyle(rel_width=0.75)

save_figures=True


#%% Import and prepare data

path = "/Users/aaron/thesis/Results/"
dirpath = path + "result_" + "04-03-21_19-19-16/"
lingrad = np.load(dirpath + "lingrad.npy")
linsurf = np.load(dirpath + "linsurf.npy")
dlingrad = np.load(dirpath + "dlingrad.npy")
dlinsurf = np.load(dirpath + "dlinsurf.npy")

with open(dirpath + "README.txt", "r") as f:
    lines = f.readlines()
    latrange=[ int(lines[4].split(' ')[-2][1:3]), int(lines[4].split(' ')[-1][0:3]) ]
    lonrange=[ int(lines[5].split(' ')[-2][1:5]), int(lines[5].split(' ')[-1][0:3]) ]
    fillres=int(lines[6].split(' ')[-1])
    gridres=int(lines[7].split(' ')[-1])


lats = np.arange(latrange[0], latrange[1]-gridres, -gridres)
lons = np.arange(lonrange[0], lonrange[1]+gridres, gridres)

def my_interpolate(array, lons, lats, method):
    
    array = np.ma.masked_invalid( array)
    xx, yy = np.meshgrid(lons, lats)
    #get only the valid values
    x1 = xx[~array.mask]
    y1 = yy[~array.mask]
    newarr = array[~array.mask]
    interpolated = interpolate.griddata((x1, y1), 
                                        newarr.ravel(),
                                        (xx, yy),
                                        method=method)
    return interpolated

lingrad_interpolated = my_interpolate(lingrad, lons, lats, 'cubic')
linsurf_interpolated = my_interpolate(linsurf, lons, lats, 'cubic')
dlingrad_interpolated = my_interpolate(dlingrad, lons, lats, 'cubic')
dlinsurf_interpolated = my_interpolate(dlinsurf, lons, lats, 'cubic')

# Grain density and porosity
longrid, latgrid = np.meshgrid(lons, lats)
# glm = pysh.SHCoeffs.from_file('/Users/aaron/thesis/Data/moon_gravity/grain_density_310.sh')
# grain = glm.expand(lat=latgrid, lon=longrid, lmax_calc=310, degrees=True)
grain = np.load("/Users/aaron/thesis/Data/moon_gravity/grain_density_310.npy")
porosity = 1 - linsurf_interpolated / grain

#%% Create maria and highlands mask

from cartopy.io import shapereader
from shapely.ops import cascaded_union
from shapely.geometry import Point

path = '/Users/aaron/thesis/Data/mare_shape/'
shape = shapereader.Reader(path + 'LROC_GLOBAL_MARE_180.shp')
polygon = cascaded_union(list(shape.geometries()))

def inpolygon(polygon, lons, lats):
    return np.array([Point(lon, lat).intersects(polygon) for lon, lat in zip(lons, lats)],
                    dtype=np.bool)

mask = inpolygon(polygon, longrid.ravel(), latgrid.ravel())

np.save(path + "mask.npy",  mask)

#%% Split analysis for maria and highlands
path = '/Users/aaron/thesis/Data/mare_shape/'
mask = np.load(path + "mask.npy")


#%% Maps
lingrad_grid = pysh.SHGrid.from_array(lingrad_interpolated)
fig1, ax1 = lingrad_grid.plot(ccrs.Orthographic(central_longitude=180.),
                               cmap=scm.diverging.Broc_20.mpl_colormap,
                               cmap_limits = [-80, 80],
                               colorbar='bottom',
                               cb_label='Linear density gradient, kg/m$^3$/km',
                               cb_tick_interval = 40,
                               cb_minor_tick_interval = 20,
                               cb_triangles='both',
                               grid=True,
                               show=False
                               )
    
linsurf_grid = pysh.SHGrid.from_array(linsurf_interpolated)
fig2, ax2 = linsurf_grid.plot(ccrs.Orthographic(central_longitude=180.),
                                cmap=scm.sequential.Bilbao_20.mpl_colormap,
                                cmap_limits = [1999, 2800],
                                colorbar='bottom',
                                cb_tick_interval = 200,
                                cb_minor_tick_interval = 100,
                                cb_label='Linear surface density, kg/m$^3$',
                                cb_triangles='both',
                                grid=True,
                                show=False
                                )
  
dlingrad_grid = pysh.SHGrid.from_array(dlingrad_interpolated)
fig3, ax3 = dlingrad_grid.plot(ccrs.Orthographic(central_longitude=180.),
                               cmap=scm.diverging.Broc_20.mpl_colormap,
                               cmap_limits = [-0.001, 5.],
                               colorbar='bottom',
                               cb_label='Linear density gradient, kg/m$^3$/km',
                               cb_tick_interval = 1,
                               cb_minor_tick_interval = 0.5,
                               cb_triangles='both',
                               grid=True,
                               show=False
                               )
   
dlinsurf_grid = pysh.SHGrid.from_array(dlinsurf_interpolated)
fig4, ax4 = dlinsurf_grid.plot(ccrs.Orthographic(central_longitude=180.),
                                cmap=scm.sequential.Bilbao_20.mpl_colormap,
                                cmap_limits = [-0.001, 20.],
                                colorbar='bottom',
                                cb_tick_interval = 5,
                                cb_minor_tick_interval = 2.5,
                                cb_label='Linear surface density, kg/m$^3$',
                                cb_triangles='both',
                                grid=True,
                                show=False
                                )

porosity_grid = pysh.SHGrid.from_array(porosity)
fig5, ax5 = porosity_grid.plot(ccrs.Orthographic(central_longitude=180.),
                               cmap=scm.sequential.Acton_20.mpl_colormap,
                               cmap_limits = [-0.001, 0.301],
                               colorbar='bottom',
                               cb_label='Porosity, -',
                               cb_tick_interval = 0.1,
                               cb_minor_tick_interval = 0.05,
                               cb_triangles='both',
                               grid=True,
                               cmap_reverse=True,
                               show=False
                               )

# grain_density_grid = pysh.SHGrid.from_array(grain)
# fig6, ax6 = grain_density_grid.plot(ccrs.Orthographic(central_longitude=180.),
#                                cmap=scm.sequential.Bilbao_20.mpl_colormap,
#                                cmap_limits = [2879, 2980],
#                                colorbar='bottom',
#                                cb_label='Grain density, kg/m$^3$',
#                                cb_tick_interval = 20,
#                                cb_minor_tick_interval = 10,
#                                cb_triangles='both',
#                                grid=True,
#                                show=False
#                                )

if save_figures:
    fig1.savefig(dirpath + "lindensgrad.pdf", format='pdf', dpi=150)
    fig2.savefig(dirpath + "linsurfdens.pdf", format='pdf', dpi=150)
    fig3.savefig(dirpath + "uncertainty_lindensgrad.pdf", format='pdf', dpi=150)
    fig4.savefig(dirpath + "uncertainty_linsurfdens.pdf", format='pdf', dpi=150)
    fig5.savefig(dirpath + "porosity.pdf", format='pdf', dpi=150)
    # fig6.savefig(dirpath + "grain_density.pdf", format='pdf', dpi=150)
    
#%% Plot 

# latmesh, lonmesh = np.meshgrid(latgrid, longrid)

# plt.figure(figsize=(3,3))       
# plt.contourf(lonmesh, latmesh, np.transpose(lingrad_interpolated),
#              cmap = scm.sequential.Devon_20.mpl_colormap) 
# plt.colorbar()
# plt.xlabel('Linear density gradient, kg/m^3/km')
# plt.ylabel('Linear surface density, kg/m^3')



#%% Plot local spectrum

# capwin.plot_windows(1, loss=True, show=False)
# ctud=[0,0.65,0.84]

# fig, ax = plt.subplots(1,1)
# ax.plot(degrees[lmin:lmax+1], global_eff_dens[lmin:lmax+1,0], '-k', label='global', linewidth=0.5, linestyle='dotted')
# ax.plot(degrees[lmin:lmax+1], local_eff_dens[lmin:lmax+1], '-k', label='local', linewidth=1)
# ax.plot(degrees[lmin:lmax+1], eff_dens_th, label='theoretical fit', color=ctud)
# ax.set(xlabel='Spherical harmonic degree', ylabel='Effective density, kg/m$^3$', 
#        xlim=(lmin,lmax), ylim=(2200, 2600))
# ax.legend()
# ax.grid()