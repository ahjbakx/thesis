"""
Created on Fri Feb  5 13:39:54 2021

@author: aaron

plot a local image of a pre-specified region of interest

"""

import pyshtools as pysh
from pyshtools import constants
from cartopy import crs as ccrs


# Width of image with respect to (journal) page
pysh.utils.figstyle(rel_width=0.5)

#%% Import data
lmax=1200

# clm = pysh.datasets.Moon.GRGM1200B_RM1_1E0(lmax=lmax)
clm = pysh.SHGravCoeffs.from_file('/Users/aaron/thesis/Data/moon_gravity/sha.grgm1200b_rm1_1e1_sigma.txt',
                                  r0_index=1,
                                  gm_index=0,
                                  errors=True)

clm.set_omega(constants.Moon.omega.value) 

#%% Plot region of interest

# Apollo 15 landing site
lat_min = 9.4308
lat_max = 42.7708
lon_min = -13.0174
lon_max = 20.3226

roi = (lon_min, lon_max, lat_min, lat_max)


grav = clm.expand(lmax=lmax, a=clm.r0, f=0.)

fig, ax = grav.plot_total(projection= ccrs.Mollweide(central_longitude=270.),
                        cmap = 'RdBu_r',
                        cmap_limits = [-400, 400],
                        colorbar='bottom',
                        cb_triangles='both',
                        show = False)
ax.set_extent(roi)