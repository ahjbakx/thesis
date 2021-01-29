% plot results
clear;
close all;
clc;

%% load data
load('GSH/Results/data_moon_full_1_420_28-Jan-2021 20:27:48')


%%
lon = data.grd.lon(1,:);
lats = data.grd.lat(:,1);

[longitude, latitude] = meshgrid(lon, lats);

figure('Position', [500 500 900 900])
axesm('vperspec', 'Grid','on','Frame','on', 'PLabelLocation',30, 'GLineWidth', 1, 'Origin',[0 0])
geoshow(latitude, longitude, ((data.vec.R)).*1e5, 'DisplayType','texturemap')
geoshow('/Users/aaron/thesis/Data/mare_shape/LROC_GLOBAL_MARE_180.shp', 'DisplayType','polygon','FaceColor','none','EdgeColor','k', 'LineWidth',0.5)

% % Gray colormap
% colormap gray

colormap(load("gravicmap.mat").cmap);

% Position colorbar below Moon
cbar = colorbar('southoutside');

mid = (cbar.Limits(1) + cbar.Limits(2))/2;
% set(cbar, 'Ticks', [-10 -8 -6 -4 -2 0 2 4 6 8 10])
% 
caxis([-600 1000])

xlabel(cbar, 'Free-air gravity (mGal)','Interpreter','latex' )

% Adjust width of colorbar
x1=get(gca,'position');
x=get(cbar,'Position');
x(1) = x(1) + 0.1250;
x(3) = 0.5;
set(cbar,'Position',x)
set(gca,'position',x1)

axis off
set(gca, 'FontSize', 20)

exportgraphics(gca,'FiguresWP2/gravity_nearside_with_mare.png','Resolution', 300)