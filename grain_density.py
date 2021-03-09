"""
Created on Mon Mar  8 13:31:11 2021

@author: aaron
"""

import pyshtools as pysh
from cartopy import crs as ccrs
from palettable import scientific as scm

save_figures = False

glm = pysh.SHCoeffs.from_file('/Users/aaron/thesis/Data/moon_gravity/grain_density_310.sh')

grain = glm.expand()

fig1, ax1 = grain.plot(projection=ccrs.Mollweide(central_longitude=270.),
                        cmap = scm.sequential.Bilbao_20.mpl_colormap,
                        colorbar='bottom',
                        cb_triangles='both',
                        cmap_limits = [2879, 3020],
                        cb_tick_interval = 20,
                        cb_label = 'Grain density, kg/m$^3$',
                        grid = True,
                        show = False)

fig2, ax2 = grain.plot(projection=ccrs.Orthographic(central_longitude=0.),
                        cmap = scm.sequential.Bilbao_20.mpl_colormap,
                        colorbar='bottom',
                        cb_triangles='both',
                        cmap_limits = [2879, 2980],
                        cb_tick_interval = 20,
                        cb_label = 'Grain density, kg/m$^3$',
                        grid = True,
                        show = False)


if save_figures:
    fig1.savefig("/Users/aaron/thesis/Figures/WP4/global-grain-density.pdf", format='pdf', dpi=150)
    fig2.savefig("/Users/aaron/thesis/Figures/WP4/nearside-grain-density.pdf", format='pdf', dpi=150)
    