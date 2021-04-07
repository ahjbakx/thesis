#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 09:36:46 2021

@author: aaron
"""

import pyshtools as pysh
import numpy as np
from cartopy import crs as ccrs
from palettable import scientific as scm
from scipy import interpolate
import matplotlib
from matplotlib import pyplot as plt

folder = "result_robust3_25-03-21_16-00-10"

plot_maps = False

#%% Import data

def make_empty_matrix(latrange, lonrange, gridres):
    
    empty_matrix = np.empty((np.int((np.sum(np.abs(latrange))/gridres+1)), np.int((np.sum(np.abs(lonrange))/gridres+1))))
    empty_matrix[:] = np.nan
    
    return empty_matrix

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


latrange = [90, -90]
lonrange = [-180, 180]
gridres = 1

lats = np.arange(latrange[0], latrange[1]-gridres, -gridres)
lons = np.arange(lonrange[0], lonrange[1]+gridres, gridres)

data_folder = "/Users/aaron/thesis/Data/"

albedo = make_empty_matrix(latrange, lonrange, gridres)
with open(data_folder + "albedo.txt", "r") as f:
        lines = f.readlines()
        l = 0
        for line in lines:
            albedo[l,:] = line.split(",")
            l += 1
            
grain_size = make_empty_matrix(latrange, lonrange, gridres)
with open(data_folder + "grain_size.txt", "r") as f:
        lines = f.readlines()
        l = 0
        for line in lines:
            grain_size[l,:] = line.split(",")
            l += 1
            
path = "/Users/aaron/thesis/Results/"
dirpath = path + folder + "/"
rho_surf = np.load(dirpath + "linsurf.npy")
rho_surf_interpolated = my_interpolate(rho_surf, lons, lats, 'cubic')
grain_density = np.load("/Users/aaron/thesis/Data/grain_density/grain_density.npy")
porosity = 1 - rho_surf_interpolated / grain_density

#%% Data preparation
longrid, latgrid = np.meshgrid(lons, lats)
mask_near = np.logical_and( abs(longrid)<70, abs(latgrid)<70 )

maria_mask = np.load("/Users/aaron/thesis/Data/mare_basalts/maria-mask.npy") 
mr_mask = maria_mask.astype(bool)
hl_mask = np.logical_not(mr_mask)

albedo_mr = [a for a in albedo[ np.logical_and(mr_mask, mask_near) ] if ( a != 0 and not np.isnan(a) ) ]
albedo_hl = [a for a in albedo[ np.logical_and(hl_mask, mask_near) ] if ( a != 0 and not np.isnan(a) ) ]

grain_size_mr = [g for g in grain_size[ np.logical_and(mr_mask, mask_near) ] if ( g != 0 and not np.isnan(g) ) ]
grain_size_hl = [g for g in grain_size[ np.logical_and(hl_mask, mask_near) ] if ( g != 0 and not np.isnan(g) ) ]

porosity_mr = [p for p in porosity[ np.logical_and(mr_mask, mask_near) ] if ( p != 0 and not np.isnan(p) ) ]
porosity_hl = [p for p in porosity[ np.logical_and(hl_mask, mask_near) ] if ( p != 0 and not np.isnan(p) ) ]

print(len(albedo_mr), len(grain_size_mr), len(porosity_mr))
print(len(albedo_hl), len(grain_size_hl), len(porosity_hl))



#%% 3D Histogram

matplotlib.rcParams.update({'font.size': 60})
plt.rc('axes', axisbelow=True)


fig, ax = plt.subplots()
ax.set_xlabel("Porosity, -")
ax.set_ylabel("Median grain size, $\mu$m")
ax.set_title('Highlands')
im = ax.scatter(porosity_hl, grain_size_hl, s=1, c=albedo_hl, marker = 'o', 
           cmap = scm.sequential.Imola_20.mpl_colormap,
           vmin=10, vmax=30)
cbar = plt.colorbar(im)
cbar.set_label("Albedo, %")
ax.grid()
plt.show()

fig, ax = plt.subplots()
ax.set_xlabel("Porosity, -")
ax.set_ylabel("Median grain size, $\mu$m")
ax.set_title('Maria')
im = ax.scatter(porosity_mr, grain_size_mr, s=1, c=albedo_mr, marker = 'o', 
           cmap = scm.sequential.Imola_20.mpl_colormap,
           vmin=7, vmax=18)
cbar = plt.colorbar(im)
cbar.set_label("Albedo, %")
ax.grid()
plt.show()

#%% Maps
if plot_maps:
    pysh.utils.figstyle(rel_width=0.75)
    
    albedo_grid = pysh.SHGrid.from_array(albedo)
    fig1, ax1 = albedo_grid.plot(ccrs.Orthographic(central_longitude=180.),
                                    cmap=scm.sequential.GrayC_20.mpl_colormap,
                                    cmap_limits = [6, 36],
                                    colorbar='bottom',
                                    cb_label='Albedo, %',
                                    # cb_tick_interval = 40,
                                    # cb_minor_tick_interval = 20,
                                    # cb_triangles='both',
                                    cmap_reverse=True,
                                    grid=True,
                                    show=False
                                    )
        
    grain_size_grid = pysh.SHGrid.from_array(grain_size)
    fig2, ax2 = grain_size_grid.plot(ccrs.Orthographic(central_longitude=180.),
                                    cmap=scm.sequential.Nuuk_20.mpl_colormap,
                                    cmap_limits = [50, 120],
                                    colorbar='bottom',
                                    cb_label='Median grain size, $\mu$m',
                                    # cb_tick_interval = 40,
                                    # cb_minor_tick_interval = 20,
                                    # cb_triangles='both',
                                    grid=True,
                                    show=False
                                    )
    
    grain_density_grid = pysh.SHGrid.from_array(grain_density)
    fig3, ax3 = grain_density_grid.plot(projection=ccrs.Orthographic(central_longitude=180.),
                            cmap = scm.sequential.Bilbao_20.mpl_colormap,
                            colorbar='bottom',
                            cb_triangles='both',
                            cmap_limits = [2799, 3600],
                            cb_tick_interval = 200,
                            cb_minor_tick_interval=50,
                            cb_label = 'Grain density, kg/m$^3$',
                            grid = True,
                            show = False)
    
    porosity_grid = pysh.SHGrid.from_array(porosity)
    fig4, ax4 = porosity_grid.plot(ccrs.Orthographic(central_longitude=180.),
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
    
        