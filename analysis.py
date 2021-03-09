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

# if not fill:
#     dirpath = path + "result_" + "04-03-21_19-19-16/"
#     lingrad = np.load(dirpath + "lingrad.npy")
#     linsurf = np.load(dirpath + "linsurf.npy")
#     dlingrad = np.load(dirpath + "dlingrad.npy")
#     dlinsurf = np.load(dirpath + "dlinsurf.npy")


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

lingrad_interpolated = my_interpolate(lingrad, longrid, latgrid, 'cubic')
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
    
linsurf_interpolated = my_interpolate(linsurf, longrid, latgrid, 'cubic')
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
  
dlingrad_interpolated = my_interpolate(dlingrad, longrid, latgrid, 'cubic')
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
   
dlinsurf_interpolated = my_interpolate(dlinsurf, longrid, latgrid, 'cubic')
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

if save_figures:
    fig1.savefig(dirpath + "/lindensgrad.pdf", format='pdf', dpi=150)
    fig2.savefig(dirpath + "linsurfdens.pdf", format='pdf', dpi=150)
    fig3.savefig(dirpath + "/uncertainty_lindensgrad.pdf", format='pdf', dpi=150)
    fig4.savefig(dirpath + "uncertainty_linsurfdens.pdf", format='pdf', dpi=150)
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