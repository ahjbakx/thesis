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
print('Load GRGM1200B RM1 1.0 gravity coefficients')

clm = pysh.datasets.Moon.GRGM1200B_RM1_1E0(lmax=lmax)

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

# selenoid = clm.geoid(u0, lmax=lmax)

# fig, ax = selenoid.plot(projection=projection,
#                         cmap = 'RdBu_r',
#                         colorbar='bottom',
#                         grid = True,
#                         show = False)


#%% Gravity field
print('Expand & plot gravity field')

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
print('Expand & plot topography')

shape = pysh.datasets.Moon.MoonTopo2600p(lmax=lmax)
shape_grid = shape.expand(grid='DH2')
topo_grid = (shape_grid - clm.r0) / 1.e3

fig, ax = topo_grid.plot(projection= projection,
                        cmap=scm.sequential.Davos_20.mpl_colormap,
                        cmap_limits = [-7, 7],
                        colorbar='bottom',
                        cb_label='Topography, km',
                        cb_triangles='both',
                        grid = True,
                        show = False)

if save_figures:
    fig.savefig('/Users/aaron/thesis/Figures/WP2/topography.pdf', 
                format='pdf', 
                dpi=300)

#%% Bouguer correction
print('Calculate Bouguer correction')

bc = pysh.SHGravCoeffs.from_shape(shape,
                                  rho=2500.,
                                  gm=clm.gm,
                                  lmax=lmax)
bc = bc.change_ref(r0=clm.r0)


# print('------ bc info ------')
# bc.info()

#%% Bouguer anomaly
print('Calculate, expand & plot Bouguer anomaly')

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
                dpi=250)

ba_filtered = ba.copy()
ba_filtered.set_coeffs(values=[0., 0., 0., 
                                0., 0., 0., 0., 0., 
                                0., 0., 0., 0., 0., 0., 0.,
                                0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                        ls=[1, 1, 1, 
                            2, 2, 2, 2, 2, 
                            3, 3, 3, 3, 3, 3, 3,
                            4, 4, 4, 4, 4, 4, 4, 4, 4,
                            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                        ms=[-1, 0, 1, 
                            -2, -1, 0, 1, 2, 
                            -3, -2, -1, 0, 1, 2, 3,
                            -4, -3, -2, -1, 0, 1, 2, 3, 4,
                            -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5,
                            -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6])

ba_filtered_grid = ba_filtered.expand(lmax=lmax, a=clm.r0, f=0.)

fig, ax = ba_filtered_grid.plot_total(projection=projection,
                        cmap=scm.diverging.Vik_20.mpl_colormap,
                        cmap_limits = [-300, 300],
                        colorbar='bottom',
                        cb_triangles='both',
                        grid = True,
                        show = False)
if save_figures:
    fig.savefig('/Users/aaron/thesis/Figures/WP2/bouguer-anomaly_7-1200.pdf', 
                format='pdf', 
                dpi=250)

