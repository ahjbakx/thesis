#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 09:36:46 2021

@author: aaron
"""


#%% Input

plot_maps = True
plot_scatters = False
plot_correlation = True
folder = "result_robust2_16-03-21_06-03-48"

#%% Configuration

import pyshtools as pysh
import numpy as np
from cartopy import crs as ccrs
from palettable import scientific as scm
from scipy import interpolate
import matplotlib
from matplotlib import pyplot as plt
import dcor


matplotlib.rcParams.update({'font.size': 20})
plt.rc('axes', axisbelow=True)

matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rc('font', family='serif', serif='CMU Serif')

def make_empty_grid(latrange, lonrange, gridres):
    
    """
    Constructs an empty DH2 grid filled with NaNs
    """
    
    empty_grid = np.empty((np.int((np.sum(np.abs(latrange))/gridres+1)), np.int((np.sum(np.abs(lonrange))/gridres+1))))
    empty_grid[:] = np.nan
    
    return empty_grid

def my_interpolate(array, lons, lats, method):
    
    """
    Interpolates data onto 1-deg res DH2 grid
    """
    
    array = np.ma.masked_invalid( array)
    xx, yy = np.meshgrid(lons, lats)

    x1 = xx[~array.mask]
    y1 = yy[~array.mask]
    newarr = array[~array.mask]
    interpolated = interpolate.griddata((x1, y1), newarr.ravel(),
                                        (xx, yy), method=method)
    return interpolated


#%% Import data
print("Importing data...")

path = "/Users/aaron/thesis/Results/"
dirpath = path + folder + "/"

with open(dirpath + "README.txt", "r") as f:
        lines = f.readlines()
        latrange=[ int(lines[4].split(' ')[-2][1:3]), int(lines[4].split(' ')[-1][0:3]) ]
        lonrange=[ int(lines[5].split(' ')[-2][1:5]), int(lines[5].split(' ')[-1][0:3]) ]
        fillres=int(lines[6].split(' ')[-1])
        gridres=int(lines[7].split(' ')[-1])

lats = np.arange(latrange[0], latrange[1]-gridres, -gridres)
lons = np.arange(lonrange[0], lonrange[1]+gridres, gridres)

""" Import polarisation data """
data_folder = "/Users/aaron/thesis/Data/"

albedo = make_empty_grid(latrange, lonrange, gridres)
with open(data_folder + "albedo.txt", "r") as f:
        lines = f.readlines()
        l = 0
        for line in lines:
            albedo[l,:] = line.split(",")
            l += 1
            
grain_size = make_empty_grid(latrange, lonrange, gridres)
with open(data_folder + "grain_size.txt", "r") as f:
        lines = f.readlines()
        l = 0
        for line in lines:
            grain_size[l,:] = line.split(",")
            l += 1
            
""" Import grain density and calculate porosity """
path = "/Users/aaron/thesis/Results/"
dirpath = path + folder + "/"
rho_surf = np.load(dirpath + "linsurf.npy")
rho_surf_interpolated = my_interpolate(rho_surf, lons, lats, 'cubic')
grain_density = np.load("/Users/aaron/thesis/Data/grain_density/grain_density.npy")
porosity = (1 - rho_surf_interpolated / grain_density ) * 100

#%% Split data into maria and highlands
print("Splitting...")

""" Set up latlon meshgrid """
longrid, latgrid = np.meshgrid(lons, lats)

""" Construct data masks """
mask_near = np.logical_and( abs(longrid)<70, abs(latgrid)<70 )
maria_mask = np.load("/Users/aaron/thesis/Data/mare_basalts/maria-mask.npy") 
mr_mask = maria_mask.astype(bool)
hl_mask = np.logical_not(mr_mask)

""" Split into maria and highlands """
#TODO: functionality to select individual maria
albedo_mr = np.array([a for a in albedo[ np.logical_and(mr_mask, mask_near) ] if ( a != 0 and not np.isnan(a) ) ])
albedo_hl = np.array([a for a in albedo[ np.logical_and(hl_mask, mask_near) ] if ( a != 0 and not np.isnan(a) ) ])

grain_size_mr = np.array([g for g in grain_size[ np.logical_and(mr_mask, mask_near) ] if ( g != 0 and not np.isnan(g) ) ])
grain_size_hl = np.array([g for g in grain_size[ np.logical_and(hl_mask, mask_near) ] if ( g != 0 and not np.isnan(g) ) ])

porosity_mr = np.array([p for p in porosity[ np.logical_and(mr_mask, mask_near) ] if ( p != 0 and not np.isnan(p) ) ])
porosity_hl = np.array([p for p in porosity[ np.logical_and(hl_mask, mask_near) ] if ( p != 0 and not np.isnan(p) ) ])

""" Check if arrays have same length """
if not len(albedo_mr) == len(grain_size_mr) == len(porosity_mr):
    print("Maria arrays have unequal length")
elif not len(albedo_hl) == len(grain_size_hl) ==  len(porosity_hl):
    print("Highlands arrays have unequal length")
else:
    print("All arrays have equal length")

#%% Construct fixed parameter arrays for correlation calculation
# plot_correlation= False
""" Fixed grain size """
g_cor = []
g_lengths = []
g_res = 2
g_min = round( min(grain_size_hl), 0 )
g_max = round( max(grain_size_hl), 0 )
g_vals = np.arange(g_min, g_max, g_res)
for g_val in g_vals:
    # g_val = 75
    fx_gs = [g for g in grain_size_hl if ( g > g_val-g_res/2 and g < g_val+g_res/2  ) ]
    
    sort_ind = np.argsort(grain_size_hl)
    g_sort = grain_size_hl[sort_ind]
    fx_gs_ind = np.searchsorted(g_sort, fx_gs)
    
    fx_gs_a = albedo_hl[sort_ind][fx_gs_ind]
    fx_gs_p = porosity_hl[sort_ind][fx_gs_ind]
    
    g_cor.append( dcor.distance_correlation(fx_gs_a, fx_gs_p) )
    g_lengths.append( len(fx_gs) )
    
    # fig, ax = plt.subplots()
    # im = ax.scatter(fx_gs_p, fx_gs_a, s=10, c=fx_gs, marker = 'o', 
    #             cmap = scm.sequential.Imola_20.mpl_colormap)
    # ax.set_xlabel("Porosity, %")
    # ax.set_ylabel("Albedo, %")
    # cbar = plt.colorbar(im)
    # cbar.set_label("Median grain size, $\mu$m")
    # ax.grid()
    # ax.set(xlim=(0, 40), ylim=(0, 30))
    # plt.show()
    # print(g_cor)
    # break


if plot_correlation:
    ctud=[0,0.65,0.84]
    fig, ax = plt.subplots()
    ax.plot(g_vals, g_cor, color='black')
    ax.set_xlabel("Median grain size, $\mu m$")
    ax.set_ylabel("Distance correlation, -")
    ax.grid()
    
    ax2 = ax.twinx()
    ax2.set_ylabel('Number of datapoints, -', color=ctud)
    ax2.plot(g_vals, g_lengths, color=ctud)
    ax2.tick_params(axis='y', labelcolor=ctud)
    
    fig.tight_layout()
    ax.set_title('Fixed median grain size')
    plt.show()

""" Fixed porosity """
p_cor = []
p_lengths = []
p_res = 0.35
p_min = round( min(porosity_hl), 0 )
p_max = round( max(porosity_hl), 0 )
p_vals = np.arange(p_min, p_max, p_res)
for p_val in p_vals:
    # p_val = 21.5
    fx_ps = [p for p in porosity_hl if ( p > p_val-p_res/2 and p < p_val+p_res/2  ) ]
    
    sort_ind = np.argsort(porosity_hl)
    p_sort = porosity_hl[sort_ind]
    fx_ps_ind = np.searchsorted(p_sort, fx_ps)
    
    fx_gs_a = albedo_hl[sort_ind][fx_ps_ind]
    fx_gs_g = grain_size_hl[sort_ind][fx_ps_ind]
    
    p_cor.append( dcor.distance_correlation(fx_gs_a, fx_gs_g) )
    p_lengths.append( len(fx_ps) )
    
    # fig, ax = plt.subplots()
    # im = ax.scatter(fx_gs_g, fx_gs_a, s=10, c=fx_ps, marker = 'o', 
    #             cmap = scm.sequential.Imola_20.mpl_colormap)
    # ax.set_xlabel("Median grain size, $\mu$m")
    # ax.set_ylabel("Albedo, %")
    # cbar = plt.colorbar(im)
    # cbar.set_label("Porosity, %")
    # ax.grid()
    # ax.set(xlim=(40, 160), ylim=(0, 30))
    # plt.show()

    # print(p_cor)
    # break


if plot_correlation:
    ctud=[0,0.65,0.84]
    fig, ax = plt.subplots()
    ax.plot(p_vals, p_cor, color='black')
    ax.set_xlabel("Porosity, %")
    ax.set_ylabel("Distance correlation, -")
    ax.grid()
    
    ax2 = ax.twinx()
    ax2.set_ylabel('Number of datapoints, -', color=ctud)
    ax2.plot(p_vals, p_lengths, color=ctud)
    ax2.tick_params(axis='y', labelcolor=ctud)
    
    fig.tight_layout()
    ax.set_title('Fixed porosity')
    plt.show()
#%% 3D Histogram

if plot_scatters:
    
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
    ax.set_xlabel("Porosity, %")
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
                                    cmap_limits = [-0.001, 30],
                                    colorbar='bottom',
                                    cb_label='Porosity, %',
                                    cb_tick_interval = 10,
                                    cb_minor_tick_interval = 5,
                                    cb_triangles='both',
                                    grid=True,
                                    cmap_reverse=True,
                                    show=False
                                    )
    
        