"""
Created on Tue Feb  9 12:40:48 2021

@author: ahjbakx
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

#%% Calculate effective density

bc = pysh.SHGravCoeffs.from_shape(hlm, rho=1., gm=clm.gm, lmax=lmax)
bc = bc.change_ref(r0=clm.r0)

ghat = pysh.SHCoeffs.from_array(bc.coeffs)
gobs = pysh.SHCoeffs.from_array(clm.coeffs)

degs = np.arange(lmin, lmax+1)
eff_dens = gobs.admittance(ghat)
corr = gobs.correlation(ghat)

#%% Least squares
k = np.sqrt(degs*(degs+1))/(clm.r0/1000)
x = 1/k
y = eff_dens[lmin:lmax+1,0]

H = np.ones((len(x), 2));
for i in range(len(x)):
    H[i,0] = x[i]


xhat = np.matmul(np.matmul(np.linalg.inv(np.matmul(np.transpose(H),H)), np.transpose(H)),y)
yhat = np.matmul(H, xhat)

Px = np.linalg.inv( 1/np.cov(y)*np.matmul(np.transpose(H), H) )

print('Initial guess:\n', 'a=', xhat[0], '+-', np.sqrt(Px[0,0]), '\n', 
      'rho=', xhat[1], '+-', np.sqrt(Px[1,1]), '\n',
      'effective density +- is', np.sqrt(Px[0,0] + Px[1,1]))

#%% Grid search

# Linear model
num = 1000


# Define linear model search space
a_vec = np.linspace(0.8*xhat[0], 1.2*xhat[0], num)
rho_vec = np.linspace(0.8*xhat[1], 1.2*xhat[1], num)

a_mesh, rho_mesh = np.meshgrid(a_vec, rho_vec)
misfits = np.zeros(a_mesh.shape)

i,j=0,0
for a, rho in np.nditer([a_mesh, rho_mesh]):
    current_eff_dens = rho + a/k
    misfit = np.sum((eff_dens[lmin:lmax+1,0] - current_eff_dens)**2)
    misfits[i,j]=misfit
    
    if i < num - 1:
        i += 1
    else:
        i = 0
        j += 1



plt.figure(figsize=(3,3))       
plt.contourf(a_mesh, rho_mesh, np.log10(misfits),
             cmap = scm.sequential.Devon_20.mpl_colormap) 
plt.colorbar()
plt.xlabel('Linear density gradient, kg/m^3/km')
plt.ylabel('Linear surface density, kg/m^3')





#%% Plots
ctud=[0,0.65,0.84]
fig, ax = plt.subplots()
eff_dens_th = xhat[1] + xhat[0]/k
d_eff_dens_th = np.sqrt(Px[0,0] + Px[1,1])
ax.plot(degs, eff_dens[lmin:lmax+1,0], label='observed', color='black', linewidth=0.5, linestyle='dotted')
ax.plot(degs, eff_dens_th, label='theoretical', color=ctud)
ax.plot(degs, eff_dens_th + d_eff_dens_th, linewidth=0.75, color=ctud, linestyle='--', label='uncertainty fit')
ax.plot(degs, eff_dens_th - d_eff_dens_th, linewidth=0.75, color=ctud, linestyle='--')
ax.set_ylabel('Effective density, kg/m$^3$')
ax.set_xlabel('Spherical harmonic degree')
ax.grid()
ax.set_ylim(2350, 2600)
ax.set_xlim(lmin, lmax)
ax.legend()

# plt.figure(figsize=(6,2))
# plt.plot(degs, corr[lmin:lmax+1])
# plt.ylabel('Correlation')
# plt.xlabel('Spherical harmonic degree')
# plt.grid()
# plt.ylim(0.8, 1)
# plt.xlim(lmin, lmax)



# #%% Calculate spectra

# # bc = pysh.SHGravCoeffs.from_shape(hlm, rho=2500., gm=clm.gm, lmax=lmax)
# # bc = bc.change_ref(r0=clm.r0)
# # ba = clm - bc

# degrees = np.arange(0, lmax+1)


# clm_spectrum = clm.spectrum(lmax=lmax)
# hlm_spectrum = hlm.spectrum(lmax=lmax)
# # bc_spectrum = bc.spectrum(lmax=lmax)
# # ba_spectrum = ba.spectrum(lmax=lmax)

# clm_spectrum[0][1]=1.0e4


# #%% Plot power spectra

# plt.figure(figsize=(6,2))
# plt.plot(degrees[1:len(degrees)], clm_spectrum[0][1: len(clm_spectrum[0])], label='GRAIL1200B RM1 10')
# plt.plot(degrees[500:len(degrees)], clm_spectrum[1][500: len(clm_spectrum[1])], label='Error')
# # plt.plot(degrees, hlm_spectrum, label='topography')
# plt.yscale('log')
# plt.grid()
# plt.ylim(1.0e-6, 1.0e6)
# plt.ylabel('Power, m$^2$')
# plt.xlabel('Spherical harmonic degree')
# plt.legend()

