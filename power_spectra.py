"""
Created on Tue Feb  9 12:40:48 2021

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

# Maximum degree
lmax=1200


#%% Import SH gravity coefficients

# clm = pysh.datasets.Moon.GRGM1200B_RM1_1E0(lmax=lmax)
clm = pysh.SHGravCoeffs.from_file('/Users/aaron/thesis/Data/moon_gravity/sha.grgm1200b_rm1_1e1_sigma.txt',
                                  r0_index=1,
                                  gm_index=0,
                                  errors=True)

hlm = pysh.datasets.Moon.MoonTopo2600p(lmax=lmax)

# set angular rotation rate (for centripetal force)
clm.set_omega(constants.Moon.omega.value) 

#%% Calculate spectra

# bc = pysh.SHGravCoeffs.from_shape(hlm, rho=2500., gm=clm.gm, lmax=lmax)
# bc = bc.change_ref(r0=clm.r0)
# ba = clm - bc

degrees = np.arange(0, lmax+1)


clm_spectrum = clm.spectrum(lmax=lmax)
hlm_spectrum = hlm.spectrum(lmax=lmax)
# bc_spectrum = bc.spectrum(lmax=lmax)
# ba_spectrum = ba.spectrum(lmax=lmax)

clm_spectrum[0][1]=1.0e4


#%% Plot power spectra

plt.figure(figsize=(6,2))
plt.plot(degrees[1:len(degrees)], clm_spectrum[0][1: len(clm_spectrum[0])], label='GRAIL1200B RM1 10')
plt.plot(degrees[500:len(degrees)], clm_spectrum[1][500: len(clm_spectrum[1])], label='Error')
# plt.plot(degrees, hlm_spectrum, label='topography')
plt.yscale('log')
plt.grid()
plt.ylim(1.0e-6, 1.0e6)
plt.ylabel('Power, m$^2$')
plt.xlabel('Spherical harmonic degree')
plt.legend()

