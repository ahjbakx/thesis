"""
Created on Tue Mar  2 10:15:14 2021

@author: ahjbakx

Localised spectral analysis - multitaper approach

    * spherical caps of given radius and bandwith
    * local effective density spectrum
    * linear least-squares fit of compaction models

"""

import pyshtools as pysh
import numpy as np
from matplotlib import pyplot as plt
from pyshtools import constants
from cartopy import crs as ccrs
from palettable import scientific as scm
from scipy import interpolate

# Width of image with respect to (journal) page
pysh.utils.figstyle(rel_width=0.5)

#%% Load data

# Maximum degree
lmin = 250
lmax = 650
lwin = 58


# clm = pysh.datasets.Moon.GRGM1200B_RM1_1E0(lmax=lmax)
clm = pysh.SHGravCoeffs.from_file('/Users/aaron/thesis/Data/moon_gravity/sha.grgm1200b_rm1_1e1_sigma.txt',
                                  r0_index=1,
                                  gm_index=0,
                                  errors=True, lmax=lmax+lwin)

hlm = pysh.datasets.Moon.MoonTopo2600p(lmax=lmax+lwin)

# set angular rotation rate (for centripetal force)
clm.set_omega(constants.Moon.omega.value) 

spectrum = clm.spectrum(function='total', lmax=lmax+lwin)
degrees = np.arange(lmax+lwin+1)

bc = pysh.SHGravCoeffs.from_shape(hlm, rho=1., gm=clm.gm, lmax=lmax+lwin)
bc = bc.change_ref(r0=clm.r0)

ghat = pysh.SHCoeffs.from_array(bc.coeffs)
gobs = pysh.SHCoeffs.from_array(clm.coeffs)

global_eff_dens = gobs.admittance(ghat)


#%% Define multitaper spherical caps

clat = 0.
clon = 90.
caprad = 15.
concentration_threshold = 0.99

# Construct spherical caps with certain radius and bandwith
capwin = pysh.SHWindow.from_cap(theta=caprad, lwin=lwin)

# Use the best caps above a concentration threshold
k = capwin.number_concentrated(concentration_threshold)
print('Number of best spherical caps: ', k)

# Rotate best spherical caps to lonlat of interest
capwin.rotate(clat=clat, clon=clon, nwinrot=k)

# Determine multitaper (local) spectrum of Bouguer correction
mts_bc, mt_bc_sd = capwin.multitaper_spectrum(ghat, k,
                                                    clat=clat,
                                                    clon=clon)

# Determine multitaper (local) cross spectrum of gravity and Bouguer correction
mtxs_topograv, mt_topograv_sd = capwin.multitaper_cross_spectrum(ghat, gobs, k,
                                                    clat=clat,
                                                    clon=clon)

# Calculate local effective density
local_eff_dens = mtxs_topograv / mts_bc

#%% Least squares fit of local spectrum

k = np.sqrt(degrees[lmin:lmax+1]*(degrees[lmin:lmax+1]+1))/(clm.r0/1000)
x = 1/k
y = local_eff_dens[lmin:lmax+1]

H = np.ones((len(x), 2));
for i in range(len(x)):
    H[i,0] = x[i]

xhat = np.matmul(np.matmul(np.linalg.inv(np.matmul(np.transpose(H),H)), np.transpose(H)),y)
yhat = np.matmul(H, xhat)

Px = np.linalg.inv( 1/np.cov(y)*np.matmul(np.transpose(H), H) )

print('Initial guess:\n', 'a=', xhat[0], '+-', np.sqrt(Px[0,0]), '\n', 
      'rho=', xhat[1], '+-', np.sqrt(Px[1,1]), '\n',
      'effective density +- is', np.sqrt(Px[0,0] + Px[1,1]))

eff_dens_th = xhat[1] + xhat[0]/k

#%% Add to grid

def my_interpolate(array, lons, lats, method):
    
    array = np.ma.masked_invalid(lingrad)
    xx, yy = np.meshgrid(longrid, latgrid)
    
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
lonrange = [90, -270]
fillres = 10
gridres = 1

latgrid = np.arange(latrange[0], latrange[1]-gridres, -gridres)
longrid = np.arange(lonrange[0], lonrange[1]-gridres, -gridres)

lingrad = np.empty((np.int((np.sum(np.abs(latrange))/gridres+1)), np.int((np.sum(np.abs(lonrange))/gridres+1))))
lingrad[:] = np.nan
linsurf = np.empty((np.int((np.sum(np.abs(latrange))/gridres+1)), np.int((np.sum(np.abs(lonrange))/gridres+1))))
linsurf[:] = np.nan

# clat_index = np.where(latgrid == clat)[0][0]
# clon_index = np.where(longrid == clon)[0][0]

# fill dummy matrix
latfills = np.arange(0, np.int(180/gridres+1), np.int(fillres/gridres))
lonfills = np.arange(0, np.int(360/gridres+1), np.int(fillres/gridres))

np.random.seed(seed=0)
for lat in latfills:
    for lon in lonfills:
        lingrad[lat, lon] = np.random.rand()


# lingrad[clat_index, clon_index] = xhat[0]
# linsurf[clat_index, clon_index] = xhat[1]

lingrad_interpolated = my_interpolate(lingrad, longrid, latgrid, 'cubic')


lingrad_grid = pysh.SHGrid.from_array(lingrad_interpolated)


fig, ax = lingrad_grid.plot(ccrs.Mollweide(central_longitude=180.),
                        cmap=scm.diverging.Vik_20.mpl_colormap,
                        cmap_limits = [0, 1],
                        colorbar='bottom',
                        cb_triangles='both',
                        grid = False,
                        show = False)

#%% Plot 

latmesh, lonmesh = np.meshgrid(latgrid, longrid)

plt.figure(figsize=(3,3))       
plt.contourf(lonmesh, latmesh, np.transpose(lingrad_interpolated),
             cmap = scm.sequential.Devon_20.mpl_colormap) 
plt.colorbar()
plt.xlabel('Linear density gradient, kg/m^3/km')
plt.ylabel('Linear surface density, kg/m^3')



#%% Plot local spectrum

capwin.plot_windows(1, loss=True, show=False)
ctud=[0,0.65,0.84]

fig, ax = plt.subplots(1,1)
ax.plot(degrees[lmin:lmax+1], global_eff_dens[lmin:lmax+1,0], '-k', label='global', linewidth=0.5, linestyle='dotted')
ax.plot(degrees[lmin:lmax+1], local_eff_dens[lmin:lmax+1], '-k', label='local', linewidth=1)
ax.plot(degrees[lmin:lmax+1], eff_dens_th, label='theoretical fit', color=ctud)
ax.set(xlabel='Spherical harmonic degree', ylabel='Effective density, kg/m$^3$', 
       xlim=(lmin,lmax), ylim=(2200, 2600))
ax.legend()
ax.grid()




