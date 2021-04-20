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
save_figures = True
only_map=True
duo=False

#%% Import and prepare data

folder = "result_robust4_15-04-21_09-42-50"
path = "/Users/aaron/thesis/Server/"
dirpath = path + folder + "/"

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

def load_files(dirpath):
    
    lingrad = np.load(dirpath + "lingrad.npy")
    linsurf = np.load(dirpath + "linsurf.npy")
    dlingrad = np.load(dirpath + "dlingrad.npy")
    dlinsurf = np.load(dirpath + "dlinsurf.npy")
    
    
    if duo:
        lincrust = np.load(dirpath + "lincrust.npy")
        dlincrust = np.load(dirpath + "dlincrust.npy")

    with open(dirpath + "README.txt", "r") as f:
        lines = f.readlines()
        latrange=[ int(lines[4].split(' ')[-2][1:3]), int(lines[4].split(' ')[-1][0:3]) ]
        lonrange=[ int(lines[5].split(' ')[-2][1:5]), int(lines[5].split(' ')[-1][0:3]) ]
        fillres=int(lines[6].split(' ')[-1])
        gridres=int(lines[7].split(' ')[-1])
    
    
    lats = np.arange(latrange[0], latrange[1]-gridres, -gridres)
    lons = np.arange(lonrange[0], lonrange[1]+gridres, gridres)
    
    if duo:
        lingrad_interpolated = my_interpolate(lingrad, lons, lats, 'cubic')
        linsurf_interpolated = my_interpolate(linsurf, lons, lats, 'cubic')
        lincrust_interpolated = my_interpolate(lincrust, lons, lats, 'cubic')
        dlingrad_interpolated = my_interpolate(dlingrad, lons, lats, 'cubic')
        dlinsurf_interpolated = my_interpolate(dlinsurf, lons, lats, 'cubic')
        dlincrust_interpolated = my_interpolate(dlincrust, lons, lats, 'cubic')
        
        return lats, lons, lingrad_interpolated, linsurf_interpolated, lincrust_interpolated, dlingrad_interpolated, dlinsurf_interpolated, dlincrust_interpolated


    else:

        lingrad_interpolated = my_interpolate(lingrad, lons, lats, 'cubic')
        linsurf_interpolated = my_interpolate(linsurf, lons, lats, 'cubic')
        dlingrad_interpolated = my_interpolate(dlingrad, lons, lats, 'cubic')
        dlinsurf_interpolated = my_interpolate(dlinsurf, lons, lats, 'cubic')
        
        return lats, lons, lingrad_interpolated, linsurf_interpolated, dlingrad_interpolated, dlinsurf_interpolated

if duo:
    lats, lons, lingrad_interpolated, linsurf_interpolated, lincrust_interpolated, dlingrad_interpolated, dlinsurf_interpolated, dlincrust_interpolated = load_files(dirpath)
else:
    lats, lons, lingrad_interpolated, linsurf_interpolated, dlingrad_interpolated, dlinsurf_interpolated = load_files(dirpath)

# Grain density and porosity
longrid, latgrid = np.meshgrid(lons, lats)
# glm = pysh.SHCoeffs.from_file('/Users/aaron/thesis/Data/moon_gravity/grain_density_310.sh')
# grain = glm.expand(lat=latgrid, lon=longrid, lmax_calc=310, degrees=True)
grain = np.load("/Users/aaron/thesis/Data/grain_density/grain_density.npy")
porosity = 1 - linsurf_interpolated / grain

#%%

# from netCDF4 import Dataset

# validation = True
# val_folder = "/Users/aaron/thesis/Data/localisation_validation/"

# ds = Dataset(val_folder + "LIN_L250-650_TC40_rho.grd", "r", format="NETCDF4")

# linsurf_G19 = np.flipud(ds['z'][:])
# linsurf_dummy = linsurf_G19
# temp = linsurf_dummy[:,0:180]
# linsurf_dummy[:,0:180] = linsurf_dummy[:,181:]
# linsurf_dummy[:,181:] = temp

# dif_linsurf = linsurf_dummy - linsurf_interpolated


# linsurf_grid = pysh.SHGrid.from_array(dif_linsurf)
# fig2, ax2 = linsurf_grid.plot(ccrs.Orthographic(central_longitude=180.),
#                                 cmap=scm.sequential.Bilbao_20.mpl_colormap,
#                                 cmap_limits = [-200, 200],
#                                 colorbar='bottom',
#                                 cb_label='Surface density, kg/m$^3$',
#                                 cb_triangles='both',
#                                 grid=True,
#                                 show=False
#                                 )

# linsurf_grid = pysh.SHGrid.from_array(linsurf_G19)
# fig2, ax2 = linsurf_grid.plot(ccrs.Mollweide(central_longitude=90.),
#                                 cmap=scm.sequential.Bilbao_20.mpl_colormap,
#                                 cmap_limits = [1999, 2800],
#                                 colorbar='bottom',
#                                 cb_label='Surface density, kg/m$^3$',
#                                 cb_triangles='both',
#                                 grid=True,
#                                 show=False
#                                 )

#%% Create maria and highlands mask

# from cartopy.io import shapereader
# from shapely.ops import cascaded_union
# from shapely.geometry import Point
# import fiona

# save_mask = False
# path = '/Users/aaron/thesis/Data/mare_shape/'
# shape = shapereader.Reader(path + 'LROC_GLOBAL_MARE_180.shp')
# polygon = cascaded_union(list(shape.geometries()))

# mask = np.zeros(shape=longrid.shape, dtype=bool)

# for i in range(361):
#     lon = lons[i]
#     for j in range(181):
#         lat = lats[j]
#         within = Point(lon, lat).within(polygon)
#         print('Longitude: ', lon, 'Latitude: ', lat, within)
#         mask[j,i] = within
  
# # mask = inpolygon(polygon, lons, lats)
# # mask = mask.reshape(longrid.shape)
# if save_mask:
#     np.save(path + "mask.npy",  mask)


#%% Split analysis for maria and highlands

if not only_map:
    pysh.utils.figstyle(rel_width=0.5)
    path = '/Users/aaron/thesis/Data/mare_shape/'
    # mask_m= np.load(path + "mask.npy")
    lines = np.loadtxt(path+'mask.txt', delimiter=',')
    mask_m = lines.astype(bool)
    mask_h = np.logical_not(mask_m)
    mask_near = abs(longrid)<90
    
    linsurf_m = linsurf_interpolated[np.logical_and(mask_m, mask_near)]
    linsurf_h = linsurf_interpolated[np.logical_and(mask_h, mask_near)]
    
    lingrad_m = lingrad_interpolated[np.logical_and(mask_m, mask_near)]
    lingrad_h = lingrad_interpolated[np.logical_and(mask_h, mask_near)]
    
    porosity_m = porosity[np.logical_and(mask_m, mask_near)]
    porosity_h = porosity[np.logical_and(mask_h, mask_near)]
    
    graindens_m = grain[np.logical_and(mask_m, mask_near)]
    graindens_h = grain[np.logical_and(mask_h, mask_near)]
    
    fig1, ax = plt.subplots()
    nbins=50
    ax.hist(linsurf_interpolated[mask_near].ravel(), nbins, label='total', color=[0.75,0.75,0.75])
    ax.hist(linsurf_h, nbins, label='highlands', color=[0, 0.4470, 0.7410], alpha=0.75)
    ax.hist(linsurf_m, nbins, label='maria', color=[0.8500, 0.3250, 0.0980], alpha=0.75)
    ax.set(xlim=(1000,4000), yticklabels=[], xlabel='Surface density, kg/m$^3$')
    ax.legend()
    
    fig2, ax = plt.subplots()
    ax.hist(lingrad_interpolated[mask_near].ravel(), nbins, label='total', color=[0.75,0.75,0.75])
    ax.hist(lingrad_h, nbins, label='highlands', color=[0, 0.4470, 0.7410], alpha=0.75)
    ax.hist(lingrad_m, nbins, label='maria', color=[0.8500, 0.3250, 0.0980], alpha=0.75)
    ax.set(xlim=(-200,200), yticklabels=[], xlabel='Density gradient, kg/m$^3$/km')
    
    fig3, ax = plt.subplots()
    ax.hist(porosity[mask_near].ravel(), nbins, label='total', color=[0.75,0.75,0.75])
    ax.hist(porosity_h, nbins, label='highlands', color=[0, 0.4470, 0.7410], alpha=0.75)
    ax.hist(porosity_m, nbins, label='maria', color=[0.8500, 0.3250, 0.0980], alpha=0.75)
    ax.set(yticklabels=[], xlabel='Porosity, -')
    ax.legend()
    
    fig4, ax = plt.subplots()
    ax.hist(grain[mask_near].ravel(), nbins, label='total', color=[0.75,0.75,0.75])
    ax.hist(graindens_h, nbins, label='highlands', color=[0, 0.4470, 0.7410], alpha=0.75)
    ax.hist(graindens_m, nbins, label='maria', color=[0.8500, 0.3250, 0.0980], alpha=0.75)
    ax.set(yticklabels=[], xlabel='Grain density, kg/m$^3$')
    ax.legend()
    
    if save_figures:
        # fig1.savefig(dirpath + "linsurfdens-hist.png", format='png', dpi=300)
        # fig2.savefig(dirpath + "lindensgrad-hist.png", format='png', dpi=300)
        # fig3.savefig(dirpath + "porosity-hist.png", format='png', dpi=300)
        fig4.savefig(dirpath + "graindens-hist.png", format='png', dpi=300)
    
#%% Latitudinal dependency
if not only_map:
    mask_0 = abs(latgrid)<15
    mask_15 = np.logical_and( abs(latgrid)>15, abs(latgrid)<30 )
    mask_30 = np.logical_and( abs(latgrid)>30, abs(latgrid)<50 )
    mask_50 = abs(latgrid)>50
    
    mask_0 = np.logical_and(mask_0, mask_near)
    mask_15 = np.logical_and(mask_15, mask_near)
    mask_30 = np.logical_and(mask_30, mask_near)
    mask_50 = np.logical_and(mask_50, mask_near)
    
    linsurf_0 = linsurf_interpolated[mask_0]
    linsurf_15 = linsurf_interpolated[mask_15]
    linsurf_30 = linsurf_interpolated[mask_30]
    linsurf_50 = linsurf_interpolated[mask_50]
    linsurf_m_0 = linsurf_interpolated[np.logical_and(mask_m, mask_0)]
    linsurf_m_15 = linsurf_interpolated[np.logical_and(mask_m, mask_15)]
    linsurf_m_30 = linsurf_interpolated[np.logical_and(mask_m, mask_30)]
    linsurf_m_50 = linsurf_interpolated[np.logical_and(mask_m, mask_50)]
    linsurf_h_0 = linsurf_interpolated[np.logical_and(mask_h, mask_0)]
    linsurf_h_15 = linsurf_interpolated[np.logical_and(mask_h, mask_15)]
    linsurf_h_30 = linsurf_interpolated[np.logical_and(mask_h, mask_30)]
    linsurf_h_50 = linsurf_interpolated[np.logical_and(mask_h, mask_50)]
    
    lingrad_0 = lingrad_interpolated[mask_0]
    lingrad_15 = lingrad_interpolated[mask_15]
    lingrad_30 = lingrad_interpolated[mask_30]
    lingrad_50 = lingrad_interpolated[mask_50]
    lingrad_m_0 = lingrad_interpolated[np.logical_and(mask_m, mask_0)]
    lingrad_m_15 = lingrad_interpolated[np.logical_and(mask_m, mask_15)]
    lingrad_m_30 = lingrad_interpolated[np.logical_and(mask_m, mask_30)]
    lingrad_m_50 = lingrad_interpolated[np.logical_and(mask_m, mask_50)]
    lingrad_h_0 = lingrad_interpolated[np.logical_and(mask_h, mask_0)]
    lingrad_h_15 = lingrad_interpolated[np.logical_and(mask_h, mask_15)]
    lingrad_h_30 = lingrad_interpolated[np.logical_and(mask_h, mask_30)]
    lingrad_h_50 = lingrad_interpolated[np.logical_and(mask_h, mask_50)]
    
    porosity_0 = porosity[mask_0]
    porosity_15 = porosity[mask_15]
    porosity_30 = porosity[mask_30]
    porosity_50 = porosity[mask_50]
    porosity_m_0 = porosity[np.logical_and(mask_m, mask_0)]
    porosity_m_15 = porosity[np.logical_and(mask_m, mask_15)]
    porosity_m_30 = porosity[np.logical_and(mask_m, mask_30)]
    porosity_m_50 = porosity[np.logical_and(mask_m, mask_50)]
    porosity_h_0 = porosity[np.logical_and(mask_h, mask_0)]
    porosity_h_15 = porosity[np.logical_and(mask_h, mask_15)]
    porosity_h_30 = porosity[np.logical_and(mask_h, mask_30)]
    porosity_h_50 = porosity[np.logical_and(mask_h, mask_50)]
    
    
    def boxplots(data, label, legend=False):
    
        positions = [1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15]
        colors = [[0.75,0.75,0.75], [0.75,0.75,0.75], [0.75,0.75,0.75], [0.75,0.75,0.75],
                  [0.8500, 0.3250, 0.0980], [0.8500, 0.3250, 0.0980], [0.8500, 0.3250, 0.0980], [0.8500, 0.3250, 0.0980],
                  [0, 0.4470, 0.7410], [0, 0.4470, 0.7410], [0, 0.4470, 0.7410], [0, 0.4470, 0.7410]]
    
        fig, ax = plt.subplots()
        box1=ax.boxplot(data, sym='', positions=positions)
        ax.set(xticklabels=['', '', '', '', r'$|\beta|<15\degree$', r'$15\degree<|\beta|<30\degree$', r'$30\degree<|\beta|<50\degree$', r'$|\beta|>50\degree$', '' ,'' ,'', ''],
               ylabel=label)
        
        for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
            line=0
            skip=False
            if item == 'whiskers' or item == 'caps':
                for box in box1[item]:
                    plt.setp(box, color=colors[line])
                    if skip:
                        line+=1
                        skip=False
                    else:
                        skip=True
            else:
                for box in box1[item]:
                    plt.setp(box, color=colors[line])
                    line+=1
                    
        if legend:
            ax.legend([box1['boxes'][0], box1['boxes'][4], box1['boxes'][8]], ['total', 'maria', 'highlands'], loc='upper right')
        
        return fig, ax
    
    box_linsurf = [linsurf_0, linsurf_15, linsurf_30, linsurf_50,
                    linsurf_m_0, linsurf_m_15, linsurf_m_30, linsurf_m_50,
                    linsurf_h_0, linsurf_h_15, linsurf_h_30, linsurf_h_50]
    box_lingrad = [lingrad_0, lingrad_15, lingrad_30, lingrad_50,
                    lingrad_m_0, lingrad_m_15, lingrad_m_30, lingrad_m_50,
                    lingrad_h_0, lingrad_h_15, lingrad_h_30, lingrad_h_50]
    box_porosity = [porosity_0, porosity_15, porosity_30, porosity_50,
                porosity_m_0, porosity_m_15, porosity_m_30, porosity_m_50,
                porosity_h_0, porosity_h_15, porosity_h_30, porosity_h_50]
       
    fig1, ax1 = boxplots(box_linsurf, 'Surface density, kg/m$^3$', legend=True)
    fig2, ax2 = boxplots(box_lingrad, 'Density gradient, kg/m$^3$/km')
    fig3, ax3 = boxplots(box_porosity, 'Porosity, -')
    
    if save_figures:
        fig1.savefig(dirpath + "linsurfdens-box.png", format='png', dpi=300)
        fig2.savefig(dirpath + "lindensgrad-box.png", format='png', dpi=300)
        fig3.savefig(dirpath + "porosity-box.png", format='png', dpi=300)
    
    
#%% Boxplots cap radius
if not only_map:
    pysh.utils.figstyle(rel_width=0.5)
    ctud=[0,0.65,0.84]
    path = "/Users/aaron/thesis/Results/"
    
    folder30='result_robust2_16-03-21_06-03-48'
    folder15='result_val-1-G19_10-03-21_02-36-51'
    folder75='result_val-2-G19_10-03-21_13-49-06'
    folder5='result_robust1_17-03-21_12-21-43'
    folder25='result_robust1_17-03-21_12-21-43'
    
    lats, lons, lingrad30, linsurf30, dlingrad_interpolated, dlinsurf_interpolated = load_files(path + folder30 + "/")
    lats, lons, lingrad15, linsurf15, dlingrad_interpolated, dlinsurf_interpolated = load_files(path + folder15 + "/")
    lats, lons, lingrad75, linsurf75, dlingrad_interpolated, dlinsurf_interpolated = load_files(path + folder75 + "/")
    lats, lons, lingrad5, linsurf5, dlingrad_interpolated, dlinsurf_interpolated = load_files(path + folder5 + "/")
    lats, lons, lingrad25, linsurf25, dlingrad_interpolated, dlinsurf_interpolated = load_files(path + folder25 + "/")
    
    linsurf30 = linsurf30[mask_near]
    linsurf15 = linsurf15[mask_near]
    linsurf75 = linsurf75[mask_near]
    linsurf5 = linsurf5[mask_near]
    linsurf25 = linsurf25[mask_near]
    
    lingrad30 = lingrad30[mask_near]
    lingrad15 = lingrad15[mask_near]
    lingrad75 = lingrad75[mask_near]
    lingrad5 = lingrad5[mask_near]
    lingrad25 = lingrad25[mask_near]
    
    data_linsurf=[linsurf25, linsurf5,linsurf75,linsurf15,linsurf30]
    data_lingrad=[lingrad25, lingrad5,lingrad75,lingrad15,lingrad30]
    flierprops = dict(marker='.', color=ctud, markersize=1,
                      linestyle='none')
    
    fig1, ax1 = plt.subplots()
    box1=ax1.boxplot(data_linsurf, flierprops=flierprops)
    ax1.set(xticklabels=[r'$\theta_{\mathrm{cap}}=2.5\degree$',r'$5.0\degree$', r'$7.5\degree$', r'$15.0\degree$', r'$30.0\degree$' ],
           ylabel='Surface density, kg/m$^3$')
    
    for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        line=0
        skip=False
        if item == 'whiskers' or item == 'caps':
            for box in box1[item]:
                plt.setp(box, color=ctud)
                if skip:
                    line+=1
                    skip=False
                else:
                    skip=True
        else:
            for box in box1[item]:
                plt.setp(box, color=ctud)
                line+=1
                
    ax1.grid()
    
    fig2, ax2 = plt.subplots()
    box1=ax2.boxplot(data_lingrad, flierprops=flierprops)
    ax2.set(xticklabels=[r'$\theta_{\mathrm{cap}}=2.5\degree$',r'$5.0\degree$', r'$7.5\degree$', r'$15.0\degree$', r'$30.0\degree$' ],
           ylabel='Density gradient, kg/m$^3$/km')
    
    for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        line=0
        skip=False
        if item == 'whiskers' or item == 'caps':
            for box in box1[item]:
                plt.setp(box, color=ctud)
                if skip:
                    line+=1
                    skip=False
                else:
                    skip=True
        else:
            for box in box1[item]:
                plt.setp(box, color=ctud)
                line+=1
                
    ax2.grid()
    
    if save_figures:
        fig1.savefig("/Users/aaron/thesis/Figures/WP4/linsurfdens-theta-box.png", format='png', dpi=300)
        fig2.savefig("/Users/aaron/thesis/Figures/WP4/lindensgrad-theta-box.png", format='png', dpi=300)
    
       
#%% Maps
import shapely
import geopandas as gpd

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
            
# sf = shapefile.Reader("")
pysh.utils.figstyle(rel_width=0.75)

# lingrad = xr.DataArray(lingrad_interpolated)

# fig = pygmt.Figure()
# fig.grdimage(lingrad,
#              projection="G0/0/12c")
# fig.show()
# fig.savefig('hoi.pdf')

# for shape in sf.shapeRecords():
#     shape = shape.shape
#     points = np.array(shape.points)
#     intervals = list(shape.parts) + [len(shape.points)]
#     for (i, j) in zip(intervals[:-1], intervals[1:]):
#         ax1.plot(*zip(*points[i:j]), color='black', linewidth=1)

lingrad_grid = pysh.SHGrid.from_array(lingrad_interpolated)
fig1, ax1 = lingrad_grid.plot(ccrs.Orthographic(central_longitude=180.),
                                cmap=scm.diverging.Broc_20.mpl_colormap,
                                cmap_limits = [-80, 80],
                                colorbar='bottom',
                                cb_label='Density gradient, kg/m$^3$/km',
                                cb_tick_interval = 40,
                                cb_minor_tick_interval = 20,
                                cb_triangles='both',
                                grid=True,
                                show=False
                                )

linsurf_grid = pysh.SHGrid.from_array(linsurf_interpolated)
fig2, ax2 = linsurf_grid.plot(ccrs.Orthographic(central_longitude=180.),
                                cmap=scm.sequential.Bilbao_20.mpl_colormap,
                                cmap_limits = [1799, 3000],
                                colorbar='bottom',
                                # cb_tick_interval = 300,
                                # cb_minor_tick_interval = 150,
                                cb_label='Surface density, kg/m$^3$',
                                cb_triangles='both',
                                grid=True,
                                show=False
                                )
  
# dlingrad_grid = pysh.SHGrid.from_array(dlingrad_interpolated)
# fig3, ax3 = dlingrad_grid.plot(ccrs.Orthographic(central_longitude=180.),
#                                cmap=scm.diverging.Broc_20.mpl_colormap,
#                                cmap_limits = [-0.001, 5.],
#                                colorbar='bottom',
#                                cb_label='Density gradient, kg/m$^3$/km',
#                                cb_tick_interval = 1,
#                                cb_minor_tick_interval = 0.5,
#                                cb_triangles='both',
#                                grid=True,
#                                show=False
#                                )
   
# dlinsurf_grid = pysh.SHGrid.from_array(dlinsurf_interpolated)
# fig4, ax4 = dlinsurf_grid.plot(ccrs.Orthographic(central_longitude=180.),
#                                 cmap=scm.sequential.Bilbao_20.mpl_colormap,
#                                 cmap_limits = [-0.001, 20.],
#                                 colorbar='bottom',
#                                 cb_tick_interval = 5,
#                                 cb_minor_tick_interval = 2.5,
#                                 cb_label='Surface density, kg/m$^3$',
#                                 cb_triangles='both',
#                                 grid=True,
#                                 show=False
#                                 )

porosity_grid = pysh.SHGrid.from_array(porosity*100)
fig5, ax5 = porosity_grid.plot(ccrs.Orthographic(central_longitude=180.),
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

# grain_density_grid = pysh.SHGrid.from_array(grain)
# fig6, ax6 = grain_density_grid.plot(ccrs.Mollweide(central_longitude=0.),
#                                 cmap=scm.sequential.Bilbao_20.mpl_colormap,
#                                 cmap_limits = [2849, 3500],
#                                 colorbar='bottom',
#                                 cb_label='Grain density, kg/m$^3$',
#                                 cb_tick_interval = 150,
#                                 cb_minor_tick_interval = 50,
#                                 cb_triangles='both',
#                                 grid=True,
#                                 show=False
#                                 )

# ax1.add_geometries(geoms=polies, crs=ccrs.PlateCarree(central_longitude=-180.),
#                     linewidth=0.2, edgecolor='black', facecolor='none')

# ax2.add_geometries(geoms=polies, crs=ccrs.PlateCarree(central_longitude=-180.),
#                     linewidth=0.2, edgecolor='white', facecolor='none')

# ax3.add_geometries(geoms=polies, crs=ccrs.PlateCarree(central_longitude=-180.),
#                    linewidth=0.2, edgecolor='white', facecolor='none')

# ax4.add_geometries(geoms=polies, crs=ccrs.PlateCarree(central_longitude=-180.),
#                    linewidth=0.2, edgecolor='white', facecolor='none')

ax5.add_geometries(geoms=polies, crs=ccrs.PlateCarree(central_longitude=-180.),
                    linewidth=0.2, edgecolor='white', facecolor='none')

# ax6.add_geometries(geoms=polies, crs=ccrs.PlateCarree(central_longitude=-180.),
#                     linewidth=0.2, edgecolor='white', facecolor='none')

if save_figures:
    fig1.savefig(dirpath + "lindensgrad.png", format='png', dpi=300)
    fig2.savefig(dirpath + "linsurfdens.png", format='png', dpi=300)
    # fig3.savefig(dirpath + "uncertainty_lindensgrad.png", format='png', dpi=300)
    # fig4.savefig(dirpath + "uncertainty_linsurfdens.png", format='png', dpi=300)
    fig5.savefig(dirpath + "porosity.png", format='png', dpi=300)
    # fig6.savefig(dirpath + "grain_density.png", format='png', dpi=300)
    
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