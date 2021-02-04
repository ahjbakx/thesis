% plot results
clear;
close all;
clc;

%% load data
% topography = load('/Users/aaron/thesis/GSH/Results/data_moon_topography_layer_1_100_03-Feb-2021 13:56:16');
% gravity = load('/Users/aaron/thesis/GSH/Results/(moon_full) 1-420 Moon Gravity field');
bouguer = load('/Users/aaron/thesis/GSH/Results/data_moon_topography_layer_7_100_03-Feb-2021 14:06:44');

%% global map
lon = bouguer.data.grd.lon(1,:);
lats = bouguer.data.grd.lat(:,1);

[longitude, latitude] = meshgrid(lon, lats);

figure('Position', [500 500 900 900])
axesm('mollweid', 'Grid','on','Frame','on', 'PLabelLocation',30, 'GLineWidth', 1, 'Origin',[0 270])
geoshow(latitude, longitude, ((bouguer.data.vec.R)).*1e5, 'DisplayType','texturemap')
% geoshow('/Users/aaron/thesis/Data/mare_shape/LROC_GLOBAL_MARE_180.shp', 'DisplayType','polygon','FaceColor','none','EdgeColor','k', 'LineWidth',0.5)

% % Gray colormap
% colormap gray

colormap(load("/Users/aaron/thesis/Colormaps/longwavelength.mat").cmap);

% Position colorbar below Moon
cbar = colorbar('southoutside');

mid = (cbar.Limits(1) + cbar.Limits(2))/2;
% set(cbar, 'Ticks', [-10 -8 -6 -4 -2 0 2 4 6 8 10])
% 
caxis([-400 700])

xlabel(cbar, 'Bouguer anomaly (mGal)','Interpreter','latex' )

% Adjust width of colorbar
x1=get(gca,'position');
x=get(cbar,'Position');
x(1) = x(1) + 0.1250;
x(3) = 0.5;
set(cbar,'Position',x)
set(gca,'position',x1)

axis off
set(gca, 'FontSize', 20)

% exportgraphics(gca,'FiguresWP2/gravity_4_5.png','Resolution', 100)
% exportgraphics(gca,'FiguresWP2/gravity_nearside_with_mare.png','Resolution', 300)

%% local map
lon = gravity.data.grd.lon(1,:);
lats = gravity.data.grd.lat(:,1);

[longitude, latitude] = meshgrid(lon, lats);

figure('Position', [500 500 900 900])
axesm mollweid
worldmap([22.1 30], [-0.35 7.61])
geoshow(latitude, longitude, ((gravity.data.vec.R)).*1e5, 'DisplayType','texturemap')
% geoshow('/Users/aaron/thesis/Data/mare_shape/LROC_GLOBAL_MARE_180.shp', 'DisplayType','polygon','FaceColor','none','EdgeColor','k', 'LineWidth',0.5)
setm(gca, 'FontSize', 20);
setm(gca, 'GColor', 'k', 'GLineWidth', 1.5);
plotm(26.1008, 3.6527, 'ok', 'MarkerSize', 50)
plotm(26.1008, 3.6527, '.k', 'MarkerSize', 10)

% % Gray colormap
% colormap gray

colormap(load("/users/aaron/thesis/Colormaps/longwavelength.mat").cmap);

% Position colorbar below Moon
cbar = colorbar('southoutside');

mid = (cbar.Limits(1) + cbar.Limits(2))/2;
% set(cbar, 'Ticks', [-10 -8 -6 -4 -2 0 2 4 6 8 10])
% 
caxis([-300 400])

xlabel(cbar, 'Free-air gravity ($$1\leq l \leq 1200$$) (mGal)','Interpreter','latex' )

% Adjust width of colorbar
x1=get(gca,'position');
x=get(cbar,'Position');
x(1) = x(1) + 0.1250;
x(3) = 0.5;
x(2) = 0.06;
set(cbar,'Position',x)
set(gca,'position',x1)

axis off
set(gca, 'FontSize', 50)

%exportgraphics(gca,'/users/aaron/thesis/FiguresWP2/gravity_A15_1200.png','Resolution', 100)