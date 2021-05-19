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
                
    albedo[albedo==0] = np.nan
    grainsize[grainsize==0] = np.nan
    albedo = my_interpolate(albedo, lons, lats, 'cubic')
    grainsize = my_interpolate(grainsize, lons, lats, 'cubic')
    
    """ Calculate porosity """
    porosity = calculate_porosity(dirpath, lons, lats)
    
    """ Set up latlon meshgrid """
    longrid, latgrid = np.meshgrid(lons, lats)
    
    """ Split into maria and highlands """
    albedo = split(albedo, longrid, latgrid, porosity)
    grainsize = split(grainsize, longrid, latgrid, porosity)
    porosity = split(porosity, longrid, latgrid, porosity)
    
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

def split(data, longrid, latgrid, porosity):
    #TODO: functionality to select individual maria

    por_mask = np.logical_and(porosity>0, porosity<100)
    maria_mask = np.load("/Users/aaron/thesis/Data/mare_basalts/maria-mask.npy") 
    mask_near = np.logical_and( abs(longrid)<70, abs(latgrid)<70 )
    mask_near = np.logical_and( mask_near, por_mask)
    mr_mask = maria_mask.astype(bool)
    hl_mask = np.logical_not(mr_mask)
    
    mask_0 = abs(latgrid)<15
    mask_15 = np.logical_and( abs(latgrid)>15, abs(latgrid)<30 )
    mask_30 = np.logical_and( abs(latgrid)>30, abs(latgrid)<50 )
    mask_50 = abs(latgrid)>50
    
    mask_0 = np.logical_and(mask_0, mask_near)
    mask_15 = np.logical_and(mask_15, mask_near)
    mask_30 = np.logical_and(mask_30, mask_near)
    mask_50 = np.logical_and(mask_50, mask_near)
    
    def mask_data(data, mask):
        return np.array([d for d in data[ mask ] if ( d != 0 and not np.isnan(d) ) ])
    
    data_split = dict()
    
    data_split['tot'] = mask_data(data, mask_near)
    data_split['mr'] = mask_data(data, np.logical_and(mr_mask, mask_near))
    data_split['hl'] = mask_data(data, np.logical_and(hl_mask, mask_near))
    
    data_split['tot_0lat15'] = mask_data(data, mask_0)
    data_split['tot_15lat30'] = mask_data(data, mask_15)
    data_split['tot_30lat50'] = mask_data(data, mask_30)
    data_split['tot_50lat70'] = mask_data(data, mask_50)
    
    data_split['mr_0lat15'] = mask_data(data, np.logical_and(mr_mask, mask_0))
    data_split['mr_15lat30'] = mask_data(data, np.logical_and(mr_mask, mask_15))
    data_split['mr_30lat50'] = mask_data(data, np.logical_and(mr_mask, mask_30))
    data_split['mr_50lat70'] = mask_data(data, np.logical_and(mr_mask, mask_50))
    
    data_split['hl_0lat15'] = mask_data(data, np.logical_and(hl_mask, mask_0))
    data_split['hl_15lat30'] = mask_data(data, np.logical_and(hl_mask, mask_15))
    data_split['hl_30lat50'] = mask_data(data, np.logical_and(hl_mask, mask_30))
    data_split['hl_50lat70'] = mask_data(data, np.logical_and(hl_mask, mask_50))

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
"""
result_robust2_16-03-21_06-03-48
result_val-1-G19_10-03-21_02-36-51
result_val-2-G19_10-03-21_13-49-06
result_robust1_17-03-21_12-21-43
"""
folder = "result_duo20_15-04-21_14-45-44"
albedo, grainsize, porosity, config = import_data(folder)

roi = 'mr'

fig, ax = plt.subplots()
ax.set_xlabel("Porosity, %")
ax.set_ylabel("Median grain size, $\mu$m")
ax.set_title('Cap radius = ' + config['caprad'] + '$^\circ$ (' + roi + ')')
im = ax.scatter(porosity[roi], grainsize[roi], s=2, c=albedo[roi], marker = 'o', 
           cmap = scm.sequential.Imola_20.mpl_colormap,
           vmin=10, vmax=30)
cbar = plt.colorbar(im)
cbar.set_label("Albedo, %")
ax.grid()
plt.show()

# fig, ax = plt.subplots()
# ax.set_xlabel("Porosity, %")
# ax.set_ylabel("Median grain size, $\mu$m")
# ax.set_title('Maria')
# im = ax.scatter(porosity['mr'], grainsize['mr'], s=2, c=albedo['mr'], marker = 'o', 
#            cmap = scm.sequential.Imola_20.mpl_colormap,
#            vmin=7, vmax=18)
# cbar = plt.colorbar(im)
# cbar.set_label("Albedo, %")
# ax.grid()
# plt.show()



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

# calculate correlation between albedo and ...
between = 'porosity'

# in region of interest ... 
roi = 'tot'

# with model ...
duo = True

if duo:
    folders = ["result_duo18_15-04-21_01-34-38",
            "result_duo19_15-04-21_02-10-54",
            "result_duo20_15-04-21_14-45-44"]
    model = 'duo'
else:
    folders = ["result_robust2_16-03-21_06-03-48",
                "result_val-1-G19_10-03-21_02-36-51",
                "result_val-2-G19_10-03-21_13-49-06",
                "result_robust1_17-03-21_12-21-43"]
    model = 'lin'



threshold = 400

correlations = dict()
coefficients = ['SRB', 'Pearson']
statistics = ['Q1', 'Q3', 'median']


for coefficient in coefficients:
    correlations[coefficient] = dict()
    for statistic in statistics:
        correlations[coefficient][statistic] = []
correlations['caprad'] = []


for folder in folders:
    
    albedo, grainsize, porosity, config = import_data(folder)
    
    if between == 'porosity':
        fixdata = grainsize
        other = porosity
        res = 2
    elif between == 'grain size':
        fixdata = porosity
        other = grainsize
        res = 2
    else:
        print('Unknown fixdata ' + between )
        
    coeffs, lengths = fixed_correlation_analysis(fixdata=fixdata, 
                           otherdata=[albedo, other], 
                           res=res, roi=roi)
    
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
            elif statistic == 'Q1':
                correlations[coefficient][statistic].append( np.percentile(valid, 25) )
            elif statistic == 'median':
                correlations[coefficient][statistic].append( np.percentile(valid, 50) )
            elif statistic == 'Q3':
                correlations[coefficient][statistic].append( np.percentile(valid, 75) )
            else:
                print('Undefined statistic: ', statistic)
              
for coefficient in coefficients:
    for statistic in statistics:
        correlations[coefficient][statistic] = np.array(correlations[coefficient][statistic])
                
        
plt.figure()
for coeff in coefficients:
    if coeff == 'SRB':
        c = 'black'
        e = np.array([np.abs(correlations[coeff]['Q1']-correlations[coeff]['median']), 
                      np.abs(correlations[coeff]['Q3']-correlations[coeff]['median'])])
    elif coeff == 'Pearson':
        if np.mean( correlations[coeff]['median'] ) < 0:
            c = 'red'
            e = np.array([np.abs(correlations[coeff]['Q3']-correlations[coeff]['median']), 
                          np.abs(correlations[coeff]['Q1']-correlations[coeff]['median'])])
        else:
            c = 'green'
            e = np.array([np.abs(correlations[coeff]['Q1']-correlations[coeff]['median']), 
                          np.abs(correlations[coeff]['Q3']-correlations[coeff]['median'])])
   
            
    eb = plt.errorbar(correlations['caprad'], np.abs(correlations[coeff]['median']), 
                 yerr=e, fmt='o', capsize=10, capthick=2, elinewidth=2,
                 color=c, label=coeff)


plt.xlabel('Cap radius, $^\circ$')
plt.ylabel('Correlation coefficient, -')
plt.ylim(0.0, 0.8)
plt.legend(loc='lower right')
plt.grid()
plt.title('Between albedo and ' + between + '  (' + roi + ', ' + model + ')')
plt.show()

#%% Latitudinal dependency

folders = ["result_robust2_16-03-21_06-03-48",
           "result_val-1-G19_10-03-21_02-36-51",
           "result_val-2-G19_10-03-21_13-49-06",
           "result_robust1_17-03-21_12-21-43"]

latitudes = dict()
caprads = []
# Read cap radii of above models
for folder in folders:
    
    albedo, grainsize, porosity, config = import_data(folder)
    latitudes[config['caprad']] = dict()
    caprads.append(config['caprad'])

threshold = 200
statistics = ['Q1', 'Q3', 'median']
roi = 'hl'
lats = [roi + '_0lat15', roi + '_15lat30', roi + '_30lat50', roi + '_50lat70']

for caprad in caprads:
    for statistic in statistics:
        latitudes[caprad][statistic] = []
latitudes['lats'] = [7.5, 22.5, 40, 60]

for folder in folders:
    
    albedo, grainsize, porosity, config = import_data(folder)
    
    caprad = config['caprad']
    
    for lat in lats:
        coeffs, lengths = fixed_correlation_analysis(fixdata=porosity, 
                               otherdata=[albedo, grainsize], 
                               res=2, roi=lat)
        
        valid = coeffs['SRB'][lengths>threshold]
        for statistic in statistics:
            if statistic == 'min':
                latitudes[caprad][statistic].append( np.min(valid) )
            elif statistic == 'max':
                latitudes[caprad][statistic].append( np.max(valid) )
            elif statistic == 'mean':
                latitudes[caprad][statistic].append( np.mean(valid) )
            elif statistic == 'Q1':
                latitudes[caprad][statistic].append( np.percentile(valid, 25) )
            elif statistic == 'median':
                latitudes[caprad][statistic].append( np.percentile(valid, 50) )
            elif statistic == 'Q3':
                latitudes[caprad][statistic].append( np.percentile(valid, 75) )
            else:
                print('Undefined statistic: ', statistic)
              
for caprad in caprads:
    for statistic in statistics:
        latitudes[caprad][statistic] = np.array(latitudes[caprad][statistic])
                


        
plt.figure()
plt.xticks([0.11, 0.32, 0.57,0.86 ], 
           [r'$|\beta|$<15$^\circ$', r'15$^\circ$<$|\beta$|<30$^\circ$', r'30$^\circ$<$|\beta|$<50$^\circ$', r'$|\beta|$>50$^\circ$'], 
           rotation=20)
for caprad in caprads:
    e = np.array([np.abs(latitudes[caprad]['Q1']-latitudes[caprad]['median']), 
                  np.abs(latitudes[caprad]['Q3']-latitudes[caprad]['median'])])
            
    eb = plt.errorbar([0.11, 0.32, 0.57,0.86], latitudes[caprad]['median'],
                      yerr=e, fmt='o', capsize=10, capthick=2, elinewidth=2, label=r'$\theta_{cap}$=' + caprad + '$^\circ$')

plt.xlabel('Selenographic latitude')
plt.ylabel('Correlation coefficient, -')
# plt.ylim(-0.1, 0.7)
plt.legend( bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid()
plt.ylim(0.0, 0.6)
plt.title('Fixed porosity (' + roi + ')')
plt.show()

#%% Correlation total overview


folders = dict()

folders['duo'] = ["result_duo18_15-04-21_01-34-38",
        "result_duo19_15-04-21_02-10-54",
        "result_duo20_15-04-21_14-45-44",
        "result_duo21_22-04-21_03-05-16"]

folders['lin'] = ["result_robust2_16-03-21_06-03-48",
        "result_val-1-G19_10-03-21_02-36-51",
        "result_val-2-G19_10-03-21_13-49-06",
        "result_robust1_17-03-21_12-21-43"]


models = ['duo', 'lin']
rois = ['tot', 'mr', 'hl']
betweens = ['porosity', 'grain size']

threshold = 200
res = 2

correlations = dict()
statistics = ['Q1', 'Q3', 'median']

for model in models:
    correlations[model] = dict()
    for roi in rois:
        correlations[model][roi] = dict()
        for between in betweens:
            correlations[model][roi][between] = dict()
            for statistic in statistics:
                correlations[model][roi][between][statistic] = []
    correlations[model]['caprad'] = []

for model in models:
    for folder in folders[model]:
        
        albedo, grainsize, porosity, config = import_data(folder)
        
        correlations[model]['caprad'].append( float(config['caprad']) )
        
        for roi in rois:
            for between in betweens:
                
                if between == 'porosity':
                    fixdata = grainsize
                    other = porosity
                    
                elif between == 'grain size':
                    fixdata = porosity
                    other = grainsize
                else:
                    print('Unknown parameter: ' + between )
                    
                
            
                coeffs, lengths = fixed_correlation_analysis(fixdata=fixdata, 
                                       otherdata=[albedo, other], 
                                       res=res, roi=roi)
        
                valid = coeffs['SRB'][lengths>threshold]
                for statistic in statistics:
                    if statistic == 'min':
                        correlations[model][roi][between][statistic].append( np.min(valid) )
                    elif statistic == 'max':
                        correlations[model][roi][between][statistic].append( np.max(valid) )
                    elif statistic == 'mean':
                        correlations[model][roi][between][statistic].append( np.mean(valid) )
                    elif statistic == 'Q1':
                        correlations[model][roi][between][statistic].append( np.percentile(valid, 25) )
                    elif statistic == 'median':
                        correlations[model][roi][between][statistic].append( np.percentile(valid, 50) )
                    elif statistic == 'Q3':
                        correlations[model][roi][between][statistic].append( np.percentile(valid, 75) )
                    else:
                        print('Undefined statistic: ', statistic)
for model in models:
    for roi in rois:
        for between in betweens:
            for statistic in statistics:
                correlations[model][roi][between][statistic] = np.array(correlations[model][roi][between][statistic])
                    

#%% 
model = 'duo'
for between in betweens:
    plt.figure()
    
    for roi in rois:
        
        if roi == 'tot':
            c = 'black'
            label = 'total'
        elif roi == 'mr':
            c = (0.8500, 0.3250, 0.0980)
            label = 'maria'
        elif roi == 'hl':
            c = (0, 0.4470, 0.7410)
            label = 'highlands'
        else:
            print('Undefined region of interest: ', roi)

        e = np.array([np.abs(correlations[model][roi][between]['Q1']-correlations[model][roi][between]['median']), 
                              np.abs(correlations[model][roi][between]['Q3']-correlations[model][roi][between]['median'])])
                
        eb = plt.errorbar(correlations[model]['caprad'], np.abs(correlations[model][roi][between]['median']), 
                     yerr=e, fmt='-o', capsize=10, capthick=2, elinewidth=2,
                     color=c, label=label)
    
    
    plt.xlabel('Cap radius, $^\circ$')
    plt.ylabel('Correlation coefficient, -')
    plt.ylim(0.0, 0.8)
    plt.legend(loc='upper right', ncol=3)
    plt.title('Between albedo and ' + between)
    plt.grid()
    plt.show()


plt.figure()

for roi in rois:
    
    if roi == 'tot':
        c = 'black'
        label = 'total'
    elif roi == 'mr':
        c = (0.8500, 0.3250, 0.0980)
        label = 'maria'
    elif roi == 'hl':
        c = (0, 0.4470, 0.7410)
        label = 'highlands'
    else:
        print('Undefined region of interest: ', roi)
    
    
    mediandif = ( correlations['duo'][roi]['porosity']['median'] + correlations['duo'][roi]['grain size']['median'])/2 - \
        ( correlations['lin'][roi]['porosity']['median'] + correlations['lin'][roi]['grain size']['median'])/2
      
    plt.plot(correlations[model]['caprad'], np.abs(mediandif), color=c, marker='o',
                  markersize=10, label=label)

plt.xlabel('Cap radius, $^\circ$')
plt.ylabel('Difference in correlation, -')
plt.ylim(-0.02, 0.15)
plt.legend(loc='upper right', ncol=3)
plt.grid()
plt.show()
