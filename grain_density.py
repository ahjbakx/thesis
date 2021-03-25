"""
Created on Mon Mar  8 13:31:11 2021

@author: aaron
"""

import pyshtools as pysh
from cartopy import crs as ccrs
from palettable import scientific as scm
import numpy as np
from scipy import interpolate

pysh.utils.figstyle(rel_width=0.5)
save_figures = False

# glm = pysh.SHCoeffs.from_file('/Users/aaron/thesis/Data/moon_gravity/grain_density_310.sh')
# data = np.genfromtxt('/Users/aaron/thesis/Data/moon_gravity/grain_density_310.dat',
#                      dtype=None,
#                      delimiter=',')
# data = data.reshape((721,1441))
#%% Import LP elemental abundances data

from astropy.io import ascii
elem = ascii.read('/Users/aaron/thesis/Data/grain_density/lpgrs_high1_elem_abundance_2deg.tab.txt')

min_lats = elem['col2']
max_lats = elem['col3']
min_lons = elem['col4']
max_lons = elem['col5']
wts_TiO2 = elem['col12']
wts_FeO = elem['col13']

latrange = [90, -90]
lonrange = [-180, 180] 
gridres = 1
latgrid = np.arange(latrange[0], latrange[1]-gridres, -gridres)
longrid = np.arange(lonrange[0], lonrange[1]+gridres, gridres)
grain_density = np.empty((np.int((np.sum(np.abs(latrange))/gridres+1)), np.int((np.sum(np.abs(lonrange))/gridres+1))))
grain_density[:] = np.nan

#%% Fill grid

for px in range(len(min_lats)):

    gdens = 27.3*wts_FeO[px]*100. + 11.0*wts_TiO2[px]*100. + 2773.
    
    min_lat_index = np.where( latgrid==min_lats[px] )[0][0]
    max_lat_index = np.where( latgrid==max_lats[px] )[0][0]
    min_lon_index = np.where( longrid==np.round(min_lons[px]) )[0][0]
    max_lon_index = np.where( longrid==np.round(max_lons[px]) )[0][0]
    

    grain_density[np.int((max_lat_index+min_lat_index+1)/2), np.int((min_lon_index+max_lon_index+1)/2)] = gdens
    

#%% plot map

array = np.ma.masked_invalid( grain_density )
xx, yy = np.meshgrid(longrid, latgrid)
#get only the valid values
x1 = xx[~array.mask]
y1 = yy[~array.mask]
newarr = array[~array.mask]
grain_density = interpolate.griddata((x1, y1), 
                                    newarr.ravel(),
                                    (xx, yy),
                                    method='cubic')

import shapely
import geopandas as gpd

sf = gpd.read_file("/Users/aaron/thesis/Data/mare_shape/LROC_GLOBAL_MARE_180.shp")
# tolerance = 1.
# geoms = sf.simplify(tolerance, preserve_topology=True)

geoms = sf['geometry']
polies = []
for geo in geoms:
    if type(geo) == shapely.geometry.MultiPolygon:
        for poly in list(geo):
            if poly.area>100:
                polies.append(poly)
    else:
        if geo.area>100:
            polies.append(geo)

grain_density_grid=pysh.SHGrid.from_array(grain_density)

fig1, ax1 = grain_density_grid.plot(projection=ccrs.Mollweide(central_longitude=270.),
                        cmap = scm.sequential.Bilbao_20.mpl_colormap,
                        colorbar='bottom',
                        cb_triangles='both',
                        cmap_limits = [2799, 3600],
                        cb_tick_interval = 150,
                        cb_minor_tick_interval=50,
                        cb_label = 'Grain density, kg/m$^3$',
                        grid = True,
                        show = False)

fig2, ax2 = grain_density_grid.plot(projection=ccrs.Orthographic(central_longitude=180.),
                        cmap = scm.sequential.Bilbao_20.mpl_colormap,
                        colorbar='bottom',
                        cb_triangles='both',
                        cmap_limits = [2799, 3600],
                        cb_tick_interval = 200,
                        cb_minor_tick_interval=50,
                        cb_label = 'Grain density, kg/m$^3$',
                        grid = True,
                        show = False)
ax2.add_geometries(geoms=polies, crs=ccrs.PlateCarree(central_longitude=-180.),
                    linewidth=0.2, edgecolor='white', facecolor='none')


# if save_figures:
#     fig1.savefig("/Users/aaron/thesis/Figures/WP4/global-grain-density.pdf", format='pdf', dpi=150)
#     fig2.savefig("/Users/aaron/thesis/Figures/WP4/nearside-grain-density.png", format='png', dpi=300)
    