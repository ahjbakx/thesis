clear all; close all; clc

latitude = load_fits("latitude");
longitude = load_fits("longitude");

fits_albedo_v = load_fits("Av");
fits_albedo_r = load_fits("Ar");

fits_pmax_v = load_fits("Pv");
fits_pmax_r = load_fits("Pr");

% interpolate albedo to 650 nm with units log(percentage)
median_albedo_v = median(median(fits_albedo_v, 'omitnan'), 'omitnan');
median_albedo_r = median(median(fits_albedo_r, 'omitnan'), 'omitnan');
albedo_650 = log10( 100 * (fits_albedo_v + (median_albedo_r - median_albedo_v) / (676.3 - 558.6) * (650 - 558.6)) );
% size_albedo = size(fits_albedo_v);
% albedo_650 = zeros(size_albedo);
% for i = 1:size_albedo(1)
%     for j = 1:size_albedo(2)
%         slope = (fits_albedo_r(i,j) - fits_albedo_v(i,j) ) / (676.3 - 558.6);
%         albedo_650(i,j) = log10( 100 * (fits_albedo_v(i,j) + slope * (650 - 558.6)) );
%     end
% end


% interpolate pmax to 650 nm with units log(permilli)
size_pmax = size(fits_pmax_v);
pmax_650 = zeros(size_pmax);
for i = 1:size_pmax(1)
    for j = 1:size_pmax(2)
        slope = (fits_pmax_r(i,j) - fits_pmax_v(i,j) ) / (676.3 - 558.6);
        pmax_650(i,j) = log10( 1000 * (fits_pmax_v(i,j) + slope * (650 - 558.6)) );
    end
end

a = 0.845;

%% Calculate

size_pixels = size(latitude);
grain_size_SO92 = zeros(size_pixels);
for i = 1:size_pixels(1)
    for j = 1:size_pixels(2)
        grain_size_SO92(i,j) = log10(0.03 * exp( 2.9 * ( albedo_650(i,j) + a * pmax_650(i,j) ) ) );
    end
end

figure('Position', [500 500 900 900])
axesm vperspec
geoshow(latitude, longitude, grain_size_SO92, 'DisplayType','texturemap')

colormap(load("cmap.mat").cmap);
cbar = colorbar('southoutside');

mid = (cbar.Limits(1) + cbar.Limits(2))/2;
%set(cbar, 'Ticks', round([cbar.Limits(1), mid, cbar.Limits(2)], 2))
%set(cbar, 'Ticks', round([1.80, 1.93, 2.06], 2))
%caxis(round(cbar.Limits, 2))
caxis([1.80 2.16])
xlabel(cbar, strcat('log d ($\mu$m)'),'Interpreter','latex' )

% geoshow('/Users/aaron/thesis/Data/mare_shape/LROC_GLOBAL_MARE_180.shp', 'DisplayType','polygon','FaceColor','none','EdgeColor','w');

axis off
set(gca, 'FontSize', 20)
