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
import argparse
from datetime import datetime
from pyshtools import constants

path = os.getcwd() + "/"

parser = argparse.ArgumentParser()
parser.add_argument("--lwin", type=int, required=True)
parser.add_argument("--caprad", type=float, required=True)
parser.add_argument("--res", type=int, required=True)
parser.add_argument("--savename", type=str, required=True)
args = parser.parse_args()

# status generator
def range_with_status(total):
    """ iterate from 0 to total and show progress in console """
    n=0
    while n<total:
        # done = '#'*(n+1)
        todo = '-'*(total-n-1)
        # s = '<{0}>'.format(done+todo)
        s = 'Progress: ' + str(np.round((n+1)/total*100,1)) + '%'
        if not todo:
            s+='\n'        
        if n>0:
            s = '\r'+s
        print(s , end='')
        yield n
        n+=1

#%% Load data

print("Loading data...")

""" Degree inputs """
lmin = 250
lmax = 650
lwin = args.lwin

# clm = pysh.datasets.Moon.GRGM1200B_RM1_1E0(lmax=lmax)
""" Import GRAIL GRGM1200B RM1 10 model """
clm = pysh.SHGravCoeffs.from_file(path + 'sha.grgm1200b_rm1_1e1_sigma.txt',
                                  r0_index=1,
                                  gm_index=0,
                                  errors=True, lmax=lmax+lwin)

hlm = pysh.datasets.Moon.MoonTopo2600p(lmax=lmax+lwin)

""" set angular rotation rate (for centripetal force) """
clm.set_omega(constants.Moon.omega.value) 

spectrum = clm.spectrum(function='total', lmax=lmax+lwin)
degrees = np.arange(lmax+lwin+1)

R = clm.r0/1000

bc = pysh.SHGravCoeffs.from_shape(hlm, rho=1., gm=clm.gm, lmax=lmax+lwin)
bc = bc.change_ref(r0=clm.r0)

ghat = pysh.SHCoeffs.from_array(bc.coeffs)
gobs = pysh.SHCoeffs.from_array(clm.coeffs)

global_eff_dens = gobs.admittance(ghat)

basalt_thickness = np.load(path + 'basalt-thickness.npy')

#%% Localised Spectral Analysis: Multitaper Spherical Cap Approach

""" Define spherical cap radius and concentration threshold """
caprad = args.caprad
concentration_threshold = 0.99

""" Additional comments to README file """
comments = "linear + mare models"

start_time = time.time()

def make_empty_matrix(latrange, lonrange, gridres):
    
    empty_matrix = np.empty((np.int((np.sum(np.abs(latrange))/gridres+1)), np.int((np.sum(np.abs(lonrange))/gridres+1))))
    empty_matrix[:] = np.nan
    
    return empty_matrix

latrange = [90, -90]
lonrange = [-180, 180] # Not 0-360 due to interpolation discontinuity on nearside
fillres = args.res # grid to be filled by LS parameters
gridres = 1 # finer grid to be interpolated

""" Set up lonlat values for fine grid """
latgrid = np.arange(latrange[0], latrange[1]-gridres, -gridres)
longrid = np.arange(lonrange[0], lonrange[1]+gridres, gridres)

""" Construct matrices to be filled """
lingrad = make_empty_matrix(latrange, lonrange, gridres) # linear density gradient
linsurf = make_empty_matrix(latrange, lonrange, gridres) # linear surface density
lincrust = make_empty_matrix(latrange, lonrange, gridres) # upper crustal density
dlingrad = make_empty_matrix(latrange, lonrange, gridres) # uncertainty linear gradient
dlinsurf = make_empty_matrix(latrange, lonrange, gridres) # uncertainty surface density
dlincrust = make_empty_matrix(latrange, lonrange, gridres) # uncertainty upper crustal density

print("Localisation...")
""" Fill first with dummy data """
latdummies = np.arange(0, np.int(180/gridres+1), np.int(fillres/gridres))
londummies = np.arange(0, np.int(360/gridres+1), np.int(fillres/gridres))
np.random.seed(seed=1)
for lat in latdummies:
    for lon in londummies:
        lingrad[lat, lon] = 160*np.random.rand()-80
        linsurf[lat, lon] = 800*np.random.rand()+2000
        lincrust[lat, lon] = 800*np.random.rand()+2000
        dlingrad[lat, lon] = 5*np.random.rand()
        dlinsurf[lat, lon] = 20*np.random.rand()
        dlincrust[lat, lon] = 20*np.random.rand()
        
""" Construct latlon **indices** of fill grid w.r.t. fine grid """
latfills = np.arange(-90, 90+1, fillres)
lonfills = np.arange(-90, 90+1, fillres)

""" Fill nearside with actual data """
for l in range_with_status(len(latfills)):
    lat = latfills[l]
    for lon in lonfills:
        clat_index = np.where(np.isclose(latgrid, lat))[0][0]
        clon_index = np.where(np.isclose(longrid, lon))[0][0]
        
        clat = np.float( latgrid[clat_index] )
        clon = np.float( longrid[clon_index] )
    
        """ Construct spherical caps with certain radius and bandwith """
        capwin = pysh.SHWindow.from_cap(theta=caprad, lwin=lwin)
        
        """ Use the best caps above a concentration threshold """
        num = capwin.number_concentrated(concentration_threshold)
        #print('Number of best spherical caps: ', num)
    
        """ Rotate best spherical caps to lonlat of interest """
        capwin.rotate(clat=clat, clon=clon, nwinrot=num)
        
        """ Determine multitaper (local) spectrum of Bouguer correction """
        mts_bc, mt_bc_sd = capwin.multitaper_spectrum(ghat, num,
                                                            clat=clat,
                                                            clon=clon)
        
        """ Determine multitaper (local) cross spectrum of gravity and Bouguer correction """
        mtxs_topograv, mt_topograv_sd = capwin.multitaper_cross_spectrum(ghat, gobs, num,
                                                            clat=clat,
                                                            clon=clon)
        
        """ Calculate local effective density """
        local_eff_dens = mtxs_topograv / mts_bc

        
        """ Least squares fit of local spectrum """
        
        """ Linear model """
        y = local_eff_dens[lmin:lmax+1]
        k = np.sqrt(degrees[lmin:lmax+1]*(degrees[lmin:lmax+1]+1)) / R
        x = 1/k

        H = np.ones((len(x), 2));
        for i in range(len(x)):
            H[i,0] = x[i]
  
        xhat = np.matmul(np.matmul(np.linalg.inv(np.matmul(np.transpose(H),H)), np.transpose(H)),y)
        yhat = np.matmul(H, xhat)
        
        Px = np.linalg.inv( 1/np.cov(y)*np.matmul(np.transpose(H), H) )
        
        """ Add data to grid """
        lingrad[clat_index, clon_index] = xhat[0]
        linsurf[clat_index, clon_index] = xhat[1]
        lincrust[clat_index, clon_index] = xhat[1]
        dlingrad[clat_index, clon_index] = np.sqrt(Px[0,0])
        dlinsurf[clat_index, clon_index] = np.sqrt(Px[1,1])
        dlincrust[clat_index, clon_index] = 0
     
        if xhat[0] < 0: # if density gradient is negative
            """ Two-layered Mare model """
            
            """ Read basalt thickness """
            if clat_index - caprad < 0 or clat_index + caprad > basalt_thickness.shape[0] or clon_index - caprad < 0  or clon_index + caprad > basalt_thickness.shape[1]:
                Tb = basalt_thickness[clat_index, clon_index] / 1000
            else:                
                Tb = np.average( basalt_thickness[clat_index-caprad:clat_index+caprad, 
                                  clon_index-caprad:clon_index+caprad] ) / 1000
            
            a = 21
            rho_0 = 2390
            
            N = 1000 # amount of infinitesimal layers
            rN = 20 # depth of analysis

            r_0 = R - Tb # spherical radius of upper crustal layer
            r_n = np.linspace(r_0, r_0 - rN, N) # spherical radii of crustal layers
            
            rhosum = 0
            for n in range(N):
                if n == 0:
                    term = 0
                elif n == 1:
                    term = a * ( r_0 - r_n[n] ) * ( r_n[n] / R ) ** ( degrees[lmin:lmax+1] + 2 )
                else:
                    term = a * ( r_n[ n-1 ] - r_n[ n ] )  * ( r_n[n] / R ) ** ( degrees[lmin:lmax+1] + 2 )
                    
                rhosum += term
                    
            y = local_eff_dens[lmin:lmax+1] - rhosum - rho_0 * ( r_0 / R ) ** ( degrees[lmin:lmax+1]+ 2 )
            
            x = (r_0 / R) ** ( degrees[lmin:lmax+1] + 2 )
            H = np.ones((len(x), 1));
            for i in range(len(x)):
                H[i,0] = 1 - x[i]
      
            xhat = np.matmul(np.matmul(np.linalg.inv(np.matmul(np.transpose(H),H)), np.transpose(H)),y)
            yhat = np.matmul(H, xhat)
                            
            """ Add data to grid """
            lingrad[clat_index, clon_index] = a
            linsurf[clat_index, clon_index] = xhat[0]
            lincrust[clat_index, clon_index] = rho_0
            dlingrad[clat_index, clon_index] = 0
            dlinsurf[clat_index, clon_index] = np.sqrt( Px[0,0] )
            dlincrust[clat_index, clon_index] = 0
    
print("Finished: " + args.savename + "! Runtime: ", round((time.time() - start_time)/60,2), "minutes" )

""" Save arrays to files """
try:
    """ Make directory for results """
    dirpath = path + "result_" + args.savename + "_" + datetime.now().strftime('%d-%m-%y_%H-%M-%S') + "/"
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
    np.save(dirpath + "lincrust", lincrust)
    np.save(dirpath + "dlingrad", dlingrad)
    np.save(dirpath + "dlinsurf", dlinsurf)
    np.save(dirpath + "dlincrust", dlincrust)
    
except OSError:
    print ("Creation of the directory %s failed" % dirpath)   
else:
    print ("Successfully written results to %s" % dirpath)



