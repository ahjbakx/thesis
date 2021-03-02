"""
Created on Tue Mar  2 10:15:14 2021

@author: aaron
"""

import pyshtools as pysh
import numpy as np
from matplotlib import pyplot as plt
from pyshtools import constants
from cartopy import crs as ccrs
from palettable import scientific as scm

# Width of image with respect to (journal) page
pysh.utils.figstyle(rel_width=0.5)

#%% Load data

# Maximum degree
lmin=250
lmax=650

# clm = pysh.datasets.Moon.GRGM1200B_RM1_1E0(lmax=lmax)
clm = pysh.SHGravCoeffs.from_file('/Users/aaron/thesis/Data/moon_gravity/sha.grgm1200b_rm1_1e1_sigma.txt',
                                  r0_index=1,
                                  gm_index=0,
                                  errors=True)

hlm = pysh.datasets.Moon.MoonTopo2600p(lmax=1200)

# set angular rotation rate (for centripetal force)
clm.set_omega(constants.Moon.omega.value) 

spectrum = clm.spectrum(function='total', lmax=1200)
degrees = np.arange(lmin, lmax+1)

bc = pysh.SHGravCoeffs.from_shape(hlm, rho=1., gm=clm.gm, lmax=1200)
bc = bc.change_ref(r0=clm.r0)

ghat = pysh.SHCoeffs.from_array(bc.coeffs)
gobs = pysh.SHCoeffs.from_array(clm.coeffs)

global_eff_dens = gobs.admittance(ghat)

global_grid = pysh.SHCoeffs.from_array(clm.coeffs).expand()

#%% Multitaper approach

lwin = 58
caprad = 15.

capwin = pysh.SHWindow.from_cap(theta=caprad, lwin=lwin)

k = capwin.number_concentrated(0.99)

clat = 70.
clon = 270.
capwin.rotate(clat=clat, clon=clon, nwinrot=k)

mts_topo, mt_topo_sd = capwin.multitaper_spectrum(hlm, k,
                                                    clat=clat,
                                                    clon=clon)

mts_topograv, mt_topograv_sd = capwin.multitaper_cross_spectrum(bc, hlm, k,
                                                    clat=clat,
                                                    clon=clon)

# local_spectrum = capwin.biased_spectrum(spectrum[0], k)
# local_eff_dens = local_topograv_spectrum / local_topography_spectrum


#%% Plot local spectrum

capwin.plot_windows(1, loss=True, show=False)

fig, ax = plt.subplots(1,1)
ax.plot(np.arange(1143), mts_topo)
ax.plot(np.arange(1201), hlm.spectrum(lmax=1200))
ax.set(yscale='log', xlim=(lmin, lmax))

bcx = pysh.SHCoeffs.from_array(bc.coeffs)
fig, ax = plt.subplots(1,1)
ax.plot(np.arange(1143), mts_topograv)
ax.plot(np.arange(1201), hlm.cross_spectrum(bcx,lmax=1200))
ax.set(yscale='log', xlim=(lmin, lmax))



# fig, ax = plt.subplots(1,1)
# ax.plot(degrees, spectrum[0][lmin:lmax+1], '-k', label='global')
# ax.plot(degrees, local_spectrum[lmin:lmax+1], '-b', label='local')
# ax.set(yscale='log', xlabel='Spherical harmonic degree', ylabel='Power', xlim=(lmin,lmax))
# ax.legend()

# ctud=[0,0.65,0.84]
# eff_dens_th = xhat[1] + xhat[0]/k
# fig, ax = plt.subplots(1,1)
# ax.plot(degrees, eff_dens[lmin:lmax+1,0], '-k', label='global', linewidth=0.5, linestyle='dotted')
# ax.plot(degrees, local_eff_dens[lmin:lmax+1], '-k', label='local', linewidth=1)
# ax.plot(degrees, eff_dens_th, label='theoretical', color=ctud)
# ax.set(xlabel='Spherical harmonic degree', ylabel='Effective density, kg/m$^3$', xlim=(lmin,lmax))
# ax.legend()
# ax.grid()

#%% Determine lons&lats of tapers

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


grid = capwin.to_shgrid(itaper=4)

array = grid.to_array()

lons = np.linspace(0,360,119)
lats = np.linspace(90, -90, 119)
lon, lat = np.meshgrid(lons, lats)

test = find_nearest(lats, clat)
roilats = (lats > clat - caprad) & (lats < clat + caprad)
roilons  = (lons > clon - caprad) & (lons < clon + caprad)
roilon, roilat = np.meshgrid(roilons, roilats)

plt.figure(figsize=(3,3))       
plt.contourf(lon, lat, roilon&roilat,
              cmap = scm.sequential.Bilbao_20.mpl_colormap) 
plt.colorbar()





