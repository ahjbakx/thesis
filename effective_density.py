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
lmin=0
lmax=1200

# clm = pysh.datasets.Moon.GRGM1200B_RM1_1E0(lmax=lmax)
clm = pysh.SHGravCoeffs.from_file('/Users/aaron/thesis/sha.grgm1200b_rm1_1e1_sigma.txt',
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
admittance = clm.admittance(hlm)
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



#%% Theoretical effective density spectra

lmin = 250
lmax = 650
lwin = 58
degrees = np.arange(lmax+lwin+1)

R = clm.r0/1000
rho_b = 2963
Tb = 0.3
rho_0 = 2390
a = 21

N = 1000
rN = 20

# latrange = [90, -90]
# lonrange = [-180, 180] # Not 0-360 due to interpolation discontinuity on nearside
# # fillres = args.res # grid to be filled by LS parameters #TODO
# fillres = 10.
# gridres = 1 # finer grid to be interpolated

# latgrid = np.arange(latrange[0], latrange[1]-gridres, -gridres)
# longrid = np.arange(lonrange[0], lonrange[1]+gridres, gridres)

# longrid, latgrid = np.meshgrid(longrid, latgrid)
# topo = hlm.expand(lat=latgrid, lon=longrid, lmax_calc=lmax+lwin, degrees=True)


def rho_eff(degrees, R, N, rN, rho_0, a, rho_b, Tb ):
    
    r_0 = R - Tb
    r_n = np.linspace(r_0, r_0 - rN, N)
    rho_n = rho_0 + a * ( R-r_n )
    
    rhosum = 0
    for n in range(N):
        if n == 0:
            term = 0
        elif n == 1:
            term = ( rho_n[n] - rho_0 ) * ( r_n[n] / R ) ** ( degrees + 2 )
        else:
            term = ( rho_n[n] - rho_n[ n-1 ] )  * ( r_n[n] / R ) ** ( degrees + 2 )
            
        rhosum += term
            
    rho_eff = rho_b + (rho_0 - rho_b ) * ( r_0 / R) ** (degrees + 2) + rhosum
    
    return rho_eff


# rho_eff_test = rho_eff(degrees, R, N, rN, rho_0, a, rho_b, Tb)


#%% Validation Gong et al. (2016)
fig1, ax1 = plt.subplots()

lw = 1.5

rho_eff_Tb0 = rho_eff(degrees, R, N, rN, rho_0, a, rho_b, 0)
ax1.plot(degrees[lmin:lmax+1], rho_eff_Tb0[lmin:lmax+1], label='$T_b$ = 0 km', linewidth=lw)

rho_eff_Tb03 = rho_eff(degrees, R, N, rN, rho_0, a, rho_b, 0.3)
ax1.plot(degrees[lmin:lmax+1], rho_eff_Tb03[lmin:lmax+1], label='$T_b$ = 0.3 km', linewidth=lw)

rho_eff_Tb08 = rho_eff(degrees, R, N, rN, rho_0, a, rho_b, 0.8)
ax1.plot(degrees[lmin:lmax+1], rho_eff_Tb08[lmin:lmax+1], label='$T_b$ = 0.8 km', linewidth=lw)

rho_eff_Tb16 = rho_eff(degrees, R, N, rN, rho_0, a, rho_b, 1.6)
ax1.plot(degrees[lmin:lmax+1], rho_eff_Tb16[lmin:lmax+1], label='$T_b$ = 1.6 km', linewidth=lw)

ax1.set_ylabel('Effective density, kg/m$^3$')
ax1.set_xlabel('Spherical harmonic degree')
ax1.grid()
ax1.set_ylim(2250, 3000)
ax1.set_xlim(lmin, lmax)
ax1.legend(loc='upper left')


fig2, ax2 = plt.subplots()

rho_eff_rhob2800 = rho_eff(degrees, R, N, rN, rho_0, a, 2800, Tb)
ax2.plot(degrees[lmin:lmax+1], rho_eff_rhob2800[lmin:lmax+1], label=r'$\rho_b$ = 2800 kg/m$^3$', linewidth=lw)

rho_eff_rhob3000 = rho_eff(degrees, R, N, rN, rho_0, a, 3000, Tb)
ax2.plot(degrees[lmin:lmax+1], rho_eff_rhob3000[lmin:lmax+1], label=r'$\rho_b$ = 3000 kg/m$^3$', linewidth=lw)

rho_eff_rhob3200 = rho_eff(degrees, R, N, rN, rho_0, a, 3200, Tb)
ax2.plot(degrees[lmin:lmax+1], rho_eff_rhob3200[lmin:lmax+1], label=r'$\rho_b$ = 3200 kg/m$^3$', linewidth=lw)

rho_eff_rhob3400 = rho_eff(degrees, R, N, rN, rho_0, a, 3400, Tb)
ax2.plot(degrees[lmin:lmax+1], rho_eff_rhob3400[lmin:lmax+1], label=r'$\rho_b$ = 3400 kg/m$^3$', linewidth=lw)

ax2.set_ylabel('Effective density, kg/m$^3$')
ax2.set_xlabel('Spherical harmonic degree')
ax2.grid()
ax2.set_ylim(2250, 3000)
ax2.set_xlim(lmin, lmax)
ax2.legend(loc='upper left')


fig3, ax3 = plt.subplots()

rho_eff_rho02200 = rho_eff(degrees, R, N, rN, 2200, a, rho_b, Tb)
ax3.plot(degrees[lmin:lmax+1], rho_eff_rho02200[lmin:lmax+1], label=r'$\rho_0$ = 2200 kg/m$^3$', linewidth=lw)

rho_eff_rho02400 = rho_eff(degrees, R, N, rN, 2400, a, rho_b, Tb)
ax3.plot(degrees[lmin:lmax+1], rho_eff_rho02400[lmin:lmax+1], label=r'$\rho_0$ = 2400 kg/m$^3$', linewidth=lw)

rho_eff_rho02600 = rho_eff(degrees, R, N, rN, 2600, a, rho_b, Tb)
ax3.plot(degrees[lmin:lmax+1], rho_eff_rho02600[lmin:lmax+1], label=r'$\rho_0$ = 2600 kg/m$^3$', linewidth=lw)

rho_eff_rho02800 = rho_eff(degrees, R, N, rN, 2800, a, rho_b, Tb)
ax3.plot(degrees[lmin:lmax+1], rho_eff_rho02800[lmin:lmax+1], label=r'$\rho_0$ = 2800 kg/m$^3$', linewidth=lw)


ax3.set_ylabel('Effective density, kg/m$^3$')
ax3.set_xlabel('Spherical harmonic degree')
ax3.grid()
ax3.set_ylim(2250, 3000)
ax3.set_xlim(lmin, lmax)
ax3.legend(loc='upper left')


fig4, ax4 = plt.subplots()

rho_eff_amin20 = rho_eff(degrees, R, N, rN, rho_0, -20, rho_b, Tb)
ax4.plot(degrees[lmin:lmax+1], rho_eff_amin20[lmin:lmax+1], label='$a$ = -20 kg/m$^3$/km', linewidth=lw)

rho_eff_a0 = rho_eff(degrees, R, N, rN, rho_0, 0, rho_b, Tb)
ax4.plot(degrees[lmin:lmax+1], rho_eff_a0[lmin:lmax+1], label='$a$ = 0 kg/m$^3$/km', linewidth=lw)

rho_eff_a20 = rho_eff(degrees, R, N, rN, rho_0, 20, rho_b, Tb)
ax4.plot(degrees[lmin:lmax+1], rho_eff_a20[lmin:lmax+1], label='$a$ = 20 kg/m$^3$/km', linewidth=lw)

rho_eff_a40 = rho_eff(degrees, R, N, rN, rho_0, 40, rho_b, Tb)
ax4.plot(degrees[lmin:lmax+1], rho_eff_a40[lmin:lmax+1], label='$a$ = 40 kg/m$^3$/km', linewidth=lw)


ax4.set_ylabel('Effective density, kg/m$^3$')
ax4.set_xlabel('Spherical harmonic degree')
ax4.grid()
ax4.set_ylim(2250, 3000)
ax4.set_xlim(lmin, lmax)
ax4.legend(loc='upper right')

fig1.savefig("/Users/aaron/thesis/Figures/WP4/maria-model-Tb.png", format='png', dpi=300)
fig2.savefig("/Users/aaron/thesis/Figures/WP4/maria-model-rhob.png", format='png', dpi=300)
fig3.savefig("/Users/aaron/thesis/Figures/WP4/maria-model-rho0.png", format='png', dpi=300)
fig4.savefig("/Users/aaron/thesis/Figures/WP4/maria-model-a.png", format='png', dpi=300)

#%% Plots
ctud=[0,0.65,0.84]

# fig, ax = plt.subplots()
# ax.plot(degs[11:1201], admittance[11:1201,0], 'k', linewidth=0.5)
# ax.grid()
# ax.set(xlim=(lmin, lmax), ylim=(0, 500),
#         ylabel="Admittance, mGal/km", xlabel="Spherical harmonic degree")
# ax2=ax.twinx()
# ax2.plot(degs[41:1200], corr[41:1200],color=ctud, linewidth=0.5)
# ax2.set_ylabel("Correlation [-]",color=ctud)
# ax2.set(ylim=(0.8, 1), yticks=[0.8, 0.85, 0.9, 0.95, 1.0])


# fig, ax = plt.subplots()
# ax.plot(degs[40:1193], eff_dens[40:1193,0], color=ctud, linewidth=0.5)
# ax.grid()
# ax.set(xlim=(lmin, lmax), ylim=(2350, 3000),
#        ylabel="Effective density, kg/m$^3$", xlabel="Spherical harmonic degree")


# ctud=[0,0.65,0.84]
# fig, ax = plt.subplots()
# eff_dens_th = xhat[1] + xhat[0]/k
# d_eff_dens_th = np.sqrt(Px[0,0] + Px[1,1])
# ax.plot(degs, eff_dens[lmin:lmax+1,0], label='observed', color='black', linewidth=0.5, linestyle='dotted')
# ax.plot(degs, eff_dens_th, label='theoretical', color=ctud)
# ax.plot(degs, eff_dens_th + d_eff_dens_th, linewidth=0.75, color=ctud, linestyle='--', label='uncertainty fit')
# ax.plot(degs, eff_dens_th - d_eff_dens_th, linewidth=0.75, color=ctud, linestyle='--')
# ax.set_ylabel('Effective density, kg/m$^3$')
# ax.set_xlabel('Spherical harmonic degree')
# ax.grid()
# ax.set_ylim(2350, 2600)
# ax.set_xlim(lmin, lmax)
# ax.legend()

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

