clear all; close all; clc;


roundn = @(x,n) round(x*10^n)./10^n;

fid  = fopen('/Users/aaron/thesis/Data/moon_topography/LDEM_16.IMG');
data = fread(fid, [5760 2880], 'int16');

% height in km, 0.5 is offset
height_km = roundn(0.5*data/1000,4);

latitude_vec = roundn(linspace(-89.5, 89.5, 2880),4);
longitude_vec = roundn(linspace(0.5, 359.5, 5760),4);


[latitude, longitude] = meshgrid(latitude_vec, longitude_vec);
clear latitude_vec longitude_vec

%% Dummy

latitude_dummy = linspace(-89.5, 89.5, 2880/16);
longitude_dummy = linspace(0.5, 359.5, 5760/16);


%% Rewrite to boundary and density file
constant_density = 2.5; % g/cm^3

% Layer 1
bd1 = zeros(180*360,3);
rho1 = zeros(180*360,3);

% Layer 2
bd2 = zeros(180*360,3);
rho2 = zeros(180*360,3);

height_dummy = zeros(360, 180);

index = 1;
col_count = 1;
row_count = 1;
for col = 1:16:2880
    for row = 1:16:5760
        bd1(index, 1) = longitude_dummy(row_count); 
        bd1(index, 2) = latitude_dummy(col_count); 
        bd1(index, 3) = height_km(row, col);
        
        rho1(index, 1) = longitude_dummy(row_count); 
        rho1(index, 2) = latitude_dummy(col_count); 
        rho1(index, 3) = constant_density;
        
        bd2(index, 1) = longitude_dummy(row_count); 
        bd2(index, 2) = latitude_dummy(col_count); 
        
        rho2(index, 1) = longitude_dummy(row_count); 
        rho2(index, 2) = latitude_dummy(col_count); 
        rho2(index, 3) = constant_density;
        
        height_dummy(row_count, col_count) = height_km(row, col);
        
        index = index + 1;
        row_count = row_count + 1;
    end
    row_count = 1;
    col_count = col_count + 1;
end

%% Plot dummy
[latitude_dummy, longitude_dummy] = meshgrid(latitude_dummy, longitude_dummy);

figure('Position', [500 500 900 900])
% axesm('mollweid', 'Grid','on','Frame','on', 'PLabelLocation',30, 'GLineWidth', 1)
axesm('mollweid', 'Grid','on','Frame','on', 'PLabelLocation',30, 'GLineWidth', 1, 'Origin',[0 270])
geoshow(-latitude_dummy, longitude_dummy, height_dummy, 'DisplayType','texturemap')
% geoshow('/Users/aaron/thesis/Data/mare_shape/LROC_GLOBAL_MARE_180.shp', 'DisplayType','polygon','FaceColor','none','EdgeColor','k', 'LineWidth',0.5)

colormap(load("/Users/aaron/thesis/Colormaps/topography.mat").cmap);

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


%% Save

writematrix(bd1, 'Data/moon_topography/topography_boundaries_layer1.dat','Delimiter','tab');
writematrix(rho1, 'Data/moon_topography/topography_densities_layer1.dat','Delimiter','tab');


writematrix(bd2, 'Data/moon_topography/topography_boundaries_layer2.dat','Delimiter','tab');
writematrix(rho2, 'Data/moon_topography/topography_densities_layer2.dat','Delimiter','tab');

%% Plot

figure('Position', [500 500 900 900])
% axesm('mollweid', 'Grid','on','Frame','on', 'PLabelLocation',30, 'GLineWidth', 1)
axesm('vperspec', 'Grid','on','Frame','on', 'PLabelLocation',30, 'GLineWidth', 1, 'Origin',[0 0])
geoshow(latitude, longitude, height_km, 'DisplayType','texturemap')
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

% exportgraphics(gca,'FiguresWP2/topography_nearside.pdf','Resolution', 300)