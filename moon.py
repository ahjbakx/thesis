"""
Created on Thu Feb  4 09:41:53 2021

@author: aaron

Calculate & plot:
    * power spectrum
    * selenoid
    * gravity field
    * topography
    * Bouguer correction
    * Bouguer anomaly
    
Input:
    * maximum degree
    * type of projection (global/nearside)
    * option to auto-save figures
"""

import pyshtools as pysh
import numpy as np
from pyshtools import constants
from cartopy import crs as ccrs
from palettable import scientific as scm

# Width of image with respect to (journal) page
pysh.utils.figstyle(rel_width=0.75)

# Maximum degree
lmax=1200

# Type of projection 
projection = ccrs.Mollweide(central_longitude=270.) # Global
# projection = ccrs.Orthographic(central_longitude=0) # Nearside

save_figures = False
#%% Import SH gravity coefficients

# clm = pysh.datasets.Moon.GRGM1200B_RM1_1E0(lmax=lmax)
clm = pysh.SHGravCoeffs.from_file('/Users/aaron/thesis/Data/moon_gravity/sha.grgm1200b_rm1_1e1_sigma.txt',
                                  r0_index=1,
                                  gm_index=0,
                                  errors=True)

# print('------ clm info ------')
# clm.info()

# set angular rotation rate (for centripetal force)
clm.set_omega(constants.Moon.omega.value) 

#%% Power spectrum
# fig, ax = clm.plot_spectrum(function='geoid', show=False)

# 2D power spectrum
# fig, ax = clm.plot_spectrum2d(function='total', cmap_rlimits=(1.e-10, 0.01), errors = True, show=False)
 
#%% Selenoid
# print('Calculate selenoid')

geoid = clm.geoid(clm.gm/clm.r0, lmax=lmax)

fig, ax = geoid.plot(projection=projection,
                        cmap =scm.sequential.Davos_20.mpl_colormap,
                        cmap_limits = [-550, 550],
                        colorbar='bottom',
                        cb_triangles='both',
                        cb_label='Geoid, m',
                        cb_tick_interval=100,
                        cb_minor_tick_interval=50,
                        grid = True,
                        show = False)

if save_figures:
    fig.savefig('/Users/aaron/thesis/Figures/WP2/selenoid.pdf', 
                format='pdf', 
                dpi=150)

#%% Gravity field

grav = clm.expand(lmax=lmax, a=clm.r0, f=0.)

fig, ax = grav.plot_total(projection=projection,
                        cmap = scm.diverging.Vik_20.mpl_colormap,
                        cmap_limits = [-400, 400],
                        colorbar='bottom',
                        cb_triangles='both',
                        grid = True,
                        show = False)
if save_figures:
    fig.savefig('/Users/aaron/thesis/Figures/WP2/free-air-anomaly.pdf', 
                format='pdf', 
                dpi=150)
    
#%% Topography

shape = pysh.datasets.Moon.MoonTopo2600p(lmax=lmax)
shape_grid = shape.expand(grid='DH2')
topo_grid = (shape_grid - clm.r0 - geoid.geoid) / 1.e3

fig, ax = topo_grid.plot(projection= projection,
                        cmap=scm.sequential.Davos_20.mpl_colormap,
                        cmap_limits = [-6, 7],
                        cb_tick_interval=1,
                        cb_minor_tick_interval=0.5,
                        colorbar='bottom',
                        cb_label='Topography, km',
                        cb_triangles='both',
                        grid = True,
                        show = False)

if save_figures:
    fig.savefig('/Users/aaron/thesis/Figures/WP2/topography.pdf', 
                format='pdf', 
                dpi=150)

#%% Bouguer correction

bc = pysh.SHGravCoeffs.from_shape(shape, rho=2500., gm=clm.gm, lmax=lmax)
bc = bc.change_ref(r0=clm.r0)


# print('------ bc info ------')
# bc.info()

#%% Bouguer anomaly

ba = clm - bc

ba_grid = ba.expand(lmax=lmax, a=clm.r0, f=0.)

fig, ax = ba_grid.plot_total(projection=projection,
                        cmap=scm.diverging.Vik_20.mpl_colormap,
                        cmap_limits = [-650, 650],
                        colorbar='bottom',
                        cb_triangles='both',
                        grid = True,
                        show = False)

if save_figures:
    fig.savefig('/Users/aaron/thesis/Figures/WP2/bouguer-anomaly.pdf', 
                format='pdf', 
                dpi=150)


#%% Isostatic anomaly
ia = ba.copy()

co_deg = 20

values = np.zeros(shape=(1, co_deg * (co_deg + 1) + co_deg ), dtype=float)
ls = np.zeros(shape=(1, co_deg * (co_deg + 1) + co_deg ), dtype=int)
ms = np.zeros(shape=(1, co_deg * (co_deg + 1) + co_deg ), dtype=int)
for deg in range(co_deg + 1):
    ls[ :, (deg - 1) * deg + deg - 1 : deg * (deg + 1) + deg ] = deg * np.ones( 2 * deg + 1, dtype=int)
    ms[ :, (deg - 1) * deg + deg - 1 : deg * (deg + 1) + deg ] = np.arange(-deg, deg + 1, 1, dtype=int)


ia.set_coeffs(values=values, ls=ls, ms=ms)

ia_grid = ia.expand(lmax=lmax, a=clm.r0, f=0.)

fig, ax = ia_grid.plot_total(projection=projection,
                        cmap=scm.diverging.Vik_20.mpl_colormap,
                        cmap_limits = [-300, 300],
                        colorbar='bottom',
                        cb_triangles='both',
                        grid = True,
                        show = False)
if save_figures:
    fig.savefig('/Users/aaron/thesis/Figures/WP2/bouguer-anomaly_7-1200.pdf', 
                format='pdf', 
                dpi=150)

