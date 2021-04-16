#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 09:36:46 2021

@author: aaron
"""


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

def import_data(folder):
    """ Import data """
    path = "/Users/aaron/thesis/Server/"
    data_folder = "/Users/aaron/thesis/Data/"
    
    dirpath = path + folder + "/"
    
    config = get_data_config(dirpath)
    
    lats = np.arange(config['latrange'][0], config['latrange'][1]-config['gridres'], -config['gridres'])
    lons = np.arange(config['lonrange'][0], config['lonrange'][1]+config['gridres'], config['gridres'])
    
    data_folder = "/Users/aaron/thesis/Data/"
    albedo = make_empty_grid(config['latrange'], config['lonrange'], config['gridres'])
    with open(data_folder + "albedo.txt", "r") as f:
            lines = f.readlines()
            l = 0
            for line in lines:
                albedo[l,:] = line.split(",")
                l += 1
                
    grainsize = make_empty_grid(config['latrange'], config['lonrange'], config['gridres'])
    with open(data_folder + "grain_size.txt", "r") as f:
            lines = f.readlines()
            l = 0
            for line in lines:
                grainsize[l,:] = line.split(",")
                l += 1
    
    """ Calculate porosity """
    porosity = calculate_porosity(dirpath, lons, lats)
    
    """ Set up latlon meshgrid """
    longrid, latgrid = np.meshgrid(lons, lats)
    
    """ Split into maria and highlands """
    albedo = split_maria_highlands(albedo, longrid, latgrid)
    grainsize = split_maria_highlands(grainsize, longrid, latgrid)
    porosity = split_maria_highlands(porosity, longrid, latgrid)
    
    """ Check if arrays have same length """
    if not len(albedo['mr']) == len(grainsize['mr']) == len(porosity['mr']):
        print("Maria arrays have unequal length")
    elif not len(albedo['hl']) == len(grainsize['hl']) ==  len(porosity['hl']):
        print("Highlands arrays have unequal length")
    else:
        return albedo, grainsize, porosity, config

def calculate_porosity(folder_path, lons, lats):
    grain_density = np.load("/Users/aaron/thesis/Data/grain_density/grain_density.npy")
    rho_surf = np.load(folder_path + "linsurf.npy")
    rho_surf_interpolated = my_interpolate(rho_surf, lons, lats, 'cubic')
    porosity = (1 - rho_surf_interpolated / grain_density ) * 100
    return porosity

def split_maria_highlands(data, longrid, latgrid):
    #TODO: functionality to select individual maria

    maria_mask = np.load("/Users/aaron/thesis/Data/mare_basalts/maria-mask.npy") 
    mask_near = np.logical_and( abs(longrid)<70, abs(latgrid)<70 )
    mr_mask = maria_mask.astype(bool)
    hl_mask = np.logical_not(mr_mask)
    
    data_split = dict()
    
    data_split['tot'] = np.array([d for d in data[ mask_near ] if ( d != 0 and not np.isnan(d) ) ])
    data_split['mr'] = np.array([d for d in data[ np.logical_and(mr_mask, mask_near) ] if ( d != 0 and not np.isnan(d) ) ])
    data_split['hl'] = np.array([d for d in data[ np.logical_and(hl_mask, mask_near) ] if ( d != 0 and not np.isnan(d) ) ])

    return data_split

def get_data_config(folder_path):
    config=dict()
    with open(folder_path + "README.txt", "r") as f:
        lines = f.readlines()
        config['latrange']=[ int(lines[4].split(' ')[-2][1:3]), int(lines[4].split(' ')[-1][0:3]) ]
        config['lonrange']=[ int(lines[5].split(' ')[-2][1:5]), int(lines[5].split(' ')[-1][0:3]) ]
        config['fillres']=int(lines[6].split(' ')[-1])
        config['gridres']=int(lines[7].split(' ')[-1])
        config['caprad'] = lines[3].split(' ')[-1][0:-1]
        
    return config

def fixed_correlation_analysis(fixdata, otherdata, res, roi, xlabel='', plot=False):

    """ Fixed grain size """
    dcors = []
    pcors = []
    lengths = []
    min_val = round( min(fixdata[roi]), 0 )
    max_val = round( max(fixdata[roi]), 0 )
    vals = np.arange(min_val, max_val, res)
    for val in vals:

        fx_vals = [i for i in fixdata[roi] if ( i > val-res/2 and i < val+res/2  ) ]
        
        sort_ind = np.argsort(fixdata[roi]) 
        fx_vals_ind = np.searchsorted(fixdata[roi][sort_ind], fx_vals)
        
        fx_rest_1 = otherdata[0][roi][sort_ind][fx_vals_ind]
        fx_rest_2 = otherdata[1][roi][sort_ind][fx_vals_ind]
        
        dcors.append( dcor.distance_correlation(fx_rest_1, fx_rest_2) )
        pcors.append( np.corrcoef(fx_rest_1, fx_rest_2)[0,1] )
        lengths.append( len(fx_vals) )
        
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
    
    if plot:
        if np.nanmean(pcors) < 0:
            c='red'
        else:
            c='green'
        
        ctud=[0,0.65,0.84]
        fig, ax = plt.subplots()
        ax.plot(vals, dcors, color='black', label='SRB')
        ax.plot(vals, np.abs(pcors), color=c, label='Pearson')
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Correlation coefficient, -")
        ax.grid()
        
        ax2 = ax.twinx()
        ax2.set_ylabel('Number of datapoints, -', color=ctud)
        ax2.plot(vals, lengths, color=ctud)
        ax2.tick_params(axis='y', labelcolor=ctud)
        
        fig.tight_layout()
        ax.legend()
        plt.show()
        
    coeffs = dict()
    coeffs['SRB'] = np.array(dcors)
    coeffs['Pearson'] = np.array(pcors)
    
    return coeffs, np.array(lengths)
#%% Construct fixed parameter arrays for correlation calculation """

folder = "result_robust2_16-03-21_06-03-48"
albedo, grainsize, porosity, config = import_data(folder)

""" Fixed grain size """
fixed_correlation_analysis(fixdata=grainsize, 
                           otherdata=[albedo, porosity], 
                           res=2, roi='hl', xlabel="Median grain size, $\mu m$",
                           plot=True)


""" Fixed porosity """
fixed_correlation_analysis(fixdata=porosity, 
                           otherdata=[albedo, grainsize], 
                           res=0.35, roi='hl', xlabel="Porosity, %",
                           plot=True)

#%% 3D scatter

fig, ax = plt.subplots()
ax.set_xlabel("Porosity, -")
ax.set_ylabel("Median grain size, $\mu$m")
ax.set_title('Highlands')
im = ax.scatter(porosity['hl'], grainsize['hl'], s=2, c=albedo['hl'], marker = 'o', 
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
im = ax.scatter(porosity['mr'], grainsize['mr'], s=2, c=albedo['mr'], marker = 'o', 
           cmap = scm.sequential.Imola_20.mpl_colormap,
           vmin=7, vmax=18)
cbar = plt.colorbar(im)
cbar.set_label("Albedo, %")
ax.grid()
plt.show()



#%% Maps

pysh.utils.figstyle(rel_width=0.75)

albedo_grid = pysh.SHGrid.from_array(albedo['tot'])
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
    
grain_size_grid = pysh.SHGrid.from_array(grainsize['tot'])
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

porosity_grid = pysh.SHGrid.from_array(porosity['tot'])
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

 
#%% Correlation Study Multiple Models       


folders = ["result_robust2_16-03-21_06-03-48",
           "result_val-1-G19_10-03-21_02-36-51",
           "result_val-2-G19_10-03-21_13-49-06",
           "result_robust1_17-03-21_12-21-43"]

threshold = 400
roi = 'hl'
correlations = dict()
coefficients = ['SRB', 'Pearson']
statistics = ['min', 'max', 'mean']


for coefficient in coefficients:
    correlations[coefficient] = dict()
    for statistic in statistics:
        correlations[coefficient][statistic] = []
correlations['caprad'] = []

for folder in folders:
    
    albedo, grainsize, porosity, config = import_data(folder)
    
    coeffs, lengths = fixed_correlation_analysis(fixdata=grainsize, 
                           otherdata=[albedo, porosity], 
                           res=2, roi=roi)
    
    correlations['caprad'].append( float(config['caprad']) )
    
    for coefficient in coefficients:
        valid = coeffs[coefficient][lengths>threshold]
        for statistic in statistics:
            if statistic == 'min':
                correlations[coefficient][statistic].append( np.min(valid) )
            elif statistic == 'max':
                correlations[coefficient][statistic].append( np.max(valid) )
            elif statistic == 'mean':
                correlations[coefficient][statistic].append( np.mean(valid) )
            else:
                print('Unknown statistic: ', statistic)
              
for coefficient in coefficients:
    for statistic in statistics:
        correlations[coefficient][statistic] = np.array(correlations[coefficient][statistic])
                
        
plt.figure()
for coeff in coefficients:
    if coeff == 'SRB':
        c = 'black'
        e = np.array([np.abs(correlations[coeff]['min']-correlations[coeff]['mean']), 
                      np.abs(correlations[coeff]['max']-correlations[coeff]['mean'])])
    elif coeff == 'Pearson':
        if np.mean( correlations[coeff]['mean'] ) < 0:
            c = 'red'
            e = np.array([np.abs(correlations[coeff]['max']-correlations[coeff]['mean']), 
                          np.abs(correlations[coeff]['min']-correlations[coeff]['mean'])])
        else:
            c = 'green'
            e = np.array([np.abs(correlations[coeff]['min']-correlations[coeff]['mean']), 
                          np.abs(correlations[coeff]['max']-correlations[coeff]['mean'])])
   
            
    eb = plt.errorbar(correlations['caprad'], np.abs(correlations[coeff]['mean']), 
                 yerr=e, fmt='o', capsize=10, capthick=2, elinewidth=2,
                 color=c, label=coeff)

plt.xlabel('Cap radius, $^\circ$')
plt.ylabel('Correlation coefficient, -')
plt.legend()
plt.grid()
plt.title('Fixed grain size  (' + roi + ')')
plt.show()









