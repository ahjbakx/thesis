clear all; close all; clc;

fid  = fopen('/Users/aaron/thesis/Data/moon_topography/LDEM_16.IMG');
data = fread(fid, [5760 2880], 'int16');

% height in km, 0.5 is offset
height = 0.5*data/1000;

latitude_vec = linspace(90, -90, 180*16);
longitude_vec = linspace(0, 360, 360*16);

[latitude, longitude] = meshgrid(latitude_vec, longitude_vec);
clear latitude_vec longitude_vec


figure('Position', [500 500 900 900])
% axesm('mollweid', 'Grid','on','Frame','on', 'PLabelLocation',30, 'GLineWidth', 1)
axesm('vperspec', 'Grid','on','Frame','on', 'PLabelLocation',30, 'GLineWidth', 1, 'Origin',[0 0])
geoshow(latitude, longitude, height, 'DisplayType','texturemap')
% geoshow('/Users/aaron/thesis/Data/mare_shape/LROC_GLOBAL_MARE_180.shp', 'DisplayType','polygon','FaceColor','none','EdgeColor','k', 'LineWidth',0.5)

colormap(load("topocmap.mat").cmap);

% Position colorbar below Moon
cbar = colorbar('southoutside');

mid = (cbar.Limits(1) + cbar.Limits(2))/2;
set(cbar, 'Ticks', [-10 -8 -6 -4 -2 0 2 4 6 8 10])

caxis([-10 10])

xlabel(cbar, 'height (km)','Interpreter','latex' )

% Adjust width of colorbar
x1=get(gca,'position');
x=get(cbar,'Position');
x(1) = x(1) + 0.1250;
x(3) = 0.5;
set(cbar,'Position',x)
set(gca,'position',x1)

axis off
set(gca, 'FontSize', 20)

exportgraphics(gca,'FiguresWP2/topography_nearside.pdf','Resolution', 300)