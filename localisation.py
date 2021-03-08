"""
Created on Tue Mar  2 10:15:14 2021

@author: ahjbakx

Localised spectral analysis - multitaper approach

    * spherical caps of given radius and bandwith
    * local effective density spectrum
    * linear least-squares fit of compaction models
    * append fitting parameters to grid and interpolate to finer grid
    
"""

import pyshtools as pysh
import numpy as np
import os
import time
from datetime import datetime
from matplotlib import pyplot as plt
from pyshtools import constants
from cartopy import crs as ccrs
from palettable import scientific as scm
from scipy import interpolate

path = "/Users/aaron/thesis/Results/"

# Width of image with respect to (journal) page
pysh.utils.figstyle(rel_width=0.5)

#%% Load data

""" Degree inputs """
lmin = 250
lmax = 650
lwin = 58


# clm = pysh.datasets.Moon.GRGM1200B_RM1_1E0(lmax=lmax)
""" Import GRAIL GRGM1200B RM1 10 model """
clm = pysh.SHGravCoeffs.from_file('/Users/aaron/thesis/Data/moon_gravity/sha.grgm1200b_rm1_1e1_sigma.txt',
                                  r0_index=1,
                                  gm_index=0,
                                  errors=True, lmax=lmax+lwin)

hlm = pysh.datasets.Moon.MoonTopo2600p(lmax=lmax+lwin)

""" set angular rotation rate (for centripetal force) """
clm.set_omega(constants.Moon.omega.value) 

spectrum = clm.spectrum(function='total', lmax=lmax+lwin)
degrees = np.arange(lmax+lwin+1)

bc = pysh.SHGravCoeffs.from_shape(hlm, rho=1., gm=clm.gm, lmax=lmax+lwin)
bc = bc.change_ref(r0=clm.r0)

ghat = pysh.SHCoeffs.from_array(bc.coeffs)
gobs = pysh.SHCoeffs.from_array(clm.coeffs)

global_eff_dens = gobs.admittance(ghat)



#%% Localised Spectral Analysis: Multitaper Spherical Cap Approach

""" Define spherical cap radius and concentration threshold """
caprad = 15.
concentration_threshold = 0.99
fill = False

""" Additional comments to README file """
comments = "linear"

start_time = time.time()

def make_empty_matrix(latrange, lonrange, gridres):
    
    empty_matrix = np.empty((np.int((np.sum(np.abs(latrange))/gridres+1)), np.int((np.sum(np.abs(lonrange))/gridres+1))))
    empty_matrix[:] = np.nan
    
    return empty_matrix

latrange = [90, -90]
lonrange = [-180, 180] # Not 0-360 due to interpolation discontinuity on nearside
fillres = 10 # grid to be filled by LS parameters
gridres = 1 # finer grid to be interpolated

""" Set up lonlat values for fine grid """
latgrid = np.arange(latrange[0], latrange[1]-gridres, -gridres)
longrid = np.arange(lonrange[0], lonrange[1]+gridres, gridres)

""" Construct matrices to be filled """
lingrad = make_empty_matrix(latrange, lonrange, gridres) # linear density gradient
linsurf = make_empty_matrix(latrange, lonrange, gridres) # linear surface density
dlingrad = make_empty_matrix(latrange, lonrange, gridres) # uncertainty linear gradient
dlinsurf = make_empty_matrix(latrange, lonrange, gridres) # uncertainty surface density

""" Fill first with dummy data """
latdummies = np.arange(0, np.int(180/gridres+1), np.int(fillres/gridres))
londummies = np.arange(0, np.int(360/gridres+1), np.int(fillres/gridres))
np.random.seed(seed=1)
for lat in latdummies:
    for lon in londummies:
        lingrad[lat, lon] = np.random.rand()
        linsurf[lat, lon] = np.random.rand()
        dlingrad[lat, lon] = np.random.rand()
        dlinsurf[lat, lon] = np.random.rand()
        
if fill:

    """ Construct latlon **indices** of fill grid w.r.t. fine grid """
    latfills = np.arange(-90, 90+1, fillres)
    lonfills = np.arange(-90, 90+1, fillres)
    
    """ Fill nearside with actual data """
    for lat in latfills:
        for lon in lonfills:
            
            clat_index = np.where(np.isclose(latgrid, lat))[0][0]
            clon_index = np.where(np.isclose(longrid, lon))[0][0]
            
            clat = np.float( latgrid[clat_index] )
            clon = np.float( longrid[clon_index] )
            
            print('Latitude = ', clat, 'Longitude = ', clon)
            
            # Dummy data
            # lingrad[clat_index, clon_index] = 10*np.random.rand()
            # linsurf[clat_index, clon_index] = 10*np.random.rand()
            # dlingrad[clat_index, clon_index] = 10*np.random.rand()
            # dlinsurf[clat_index, clon_index] = 10*np.random.rand()
            # continue
        
            """ Construct spherical caps with certain radius and bandwith """
            capwin = pysh.SHWindow.from_cap(theta=caprad, lwin=lwin)
            
            """ Use the best caps above a concentration threshold """
            k = capwin.number_concentrated(concentration_threshold)
            #print('Number of best spherical caps: ', k)
        
            """ Rotate best spherical caps to lonlat of interest """
            capwin.rotate(clat=clat, clon=clon, nwinrot=k)
            
            """ Determine multitaper (local) spectrum of Bouguer correction """
            mts_bc, mt_bc_sd = capwin.multitaper_spectrum(ghat, k,
                                                                clat=clat,
                                                                clon=clon)
            
            """ Determine multitaper (local) cross spectrum of gravity and Bouguer correction """
            mtxs_topograv, mt_topograv_sd = capwin.multitaper_cross_spectrum(ghat, gobs, k,
                                                                clat=clat,
                                                                clon=clon)
            
            """ Calculate local effective density """
            local_eff_dens = mtxs_topograv / mts_bc
    
            """ Least squares fit of local spectrum """
            k = np.sqrt(degrees[lmin:lmax+1]*(degrees[lmin:lmax+1]+1))/(clm.r0/1000)
            x = 1/k
            y = local_eff_dens[lmin:lmax+1]
            
            H = np.ones((len(x), 2));
            for i in range(len(x)):
                H[i,0] = x[i]
            
            xhat = np.matmul(np.matmul(np.linalg.inv(np.matmul(np.transpose(H),H)), np.transpose(H)),y)
            yhat = np.matmul(H, xhat)
            
            Px = np.linalg.inv( 1/np.cov(y)*np.matmul(np.transpose(H), H) )
            
            """ Add data to grid """
            lingrad[clat_index, clon_index] = xhat[0]
            linsurf[clat_index, clon_index] = xhat[1]
            dlingrad[clat_index, clon_index] = np.sqrt(Px[0,0])
            dlinsurf[clat_index, clon_index] = np.sqrt(Px[1,1])
            
            # print('Initial guess:\n', 'a=', xhat[0], '+-', np.sqrt(Px[0,0]), '\n', 
            #       'rho=', xhat[1], '+-', np.sqrt(Px[1,1]), '\n',
            #       'effective density +- is', np.sqrt(Px[0,0] + Px[1,1]))
            
            # eff_dens_th = xhat[1] + xhat[0]/k
        
print("Finished! Runtime: ", round((time.time() - start_time)/60,2), "minutes" )
      
""" Save arrays to files """
try:
    """ Make directory for results """
    dirpath = path + "result_" + datetime.now().strftime('%d-%m-%y_%H-%M-%S') + "/"
    os.mkdir(dirpath)
    
    """ Write input parameters to file """
    f = open(dirpath + "README.txt", "x")
    f.write("lmin = " + str(lmin) + "\n" + "lmax = " + str(lmax) + "\n" +
            "lwin = " + str(lwin) + "\n" + "caprad = " + str(caprad) + "\n" +
            "latrange = " + str(latrange) + "\n" + "lonrange = " + str(lonrange) + "\n" +
            "fillres = " + str(fillres) + "\n" + "gridres = " + str(gridres) + "\n" +
            "comments = " + comments)
    f.close()
    
    """ Save arrays """
    np.save(dirpath + "lingrad", lingrad)
    np.save(dirpath + "linsurf", linsurf)
    np.save(dirpath + "dlingrad", dlingrad)
    np.save(dirpath + "dlinsurf", dlinsurf)
    
except OSError:
    print ("Creation of the directory %s failed" % dirpath)   
else:
    print ("Successfully written results to %s" % dirpath)



# Interpolate data to finer grid and plot
import_data=False
save_figures=True

if import_data:
    dirpath = path + "result_" + "04-03-21_19-19-16/"
    lingrad = np.load(dirpath + "lingrad.npy")
    linsurf = np.load(dirpath + "linsurf.npy")
    dlingrad = np.load(dirpath + "dlingrad.npy")
    dlinsurf = np.load(dirpath + "dlinsurf.npy")


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
                               cb_triangles='both',
                               grid=True,
                               show=False
                               )
    
linsurf_interpolated = my_interpolate(linsurf, longrid, latgrid, 'cubic')
linsurf_grid = pysh.SHGrid.from_array(linsurf_interpolated)
fig2, ax2 = linsurf_grid.plot(ccrs.Orthographic(central_longitude=180.),
                                cmap=scm.sequential.Davos_20.mpl_colormap,
                                cmap_limits = [1999, 2800],
                                colorbar='bottom',
                                cb_tick_interval = 200,
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
                               cb_triangles='both',
                               grid=True,
                               show=False
                               )
   
dlinsurf_interpolated = my_interpolate(dlinsurf, longrid, latgrid, 'cubic')
dlinsurf_grid = pysh.SHGrid.from_array(dlinsurf_interpolated)
fig4, ax4 = dlinsurf_grid.plot(ccrs.Orthographic(central_longitude=180.),
                                cmap=scm.sequential.Davos_20.mpl_colormap,
                                cmap_limits = [-0.001, 20.],
                                colorbar='bottom',
                                cb_tick_interval = 5,
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




