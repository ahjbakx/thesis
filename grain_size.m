clear all; close all; clc

%% Load data

latitude = load_fits("latitude");
longitude = load_fits("longitude");

fits_albedo_v = load_fits("Av");
fits_albedo_r = load_fits("Ar");

fits_pmax_v = load_fits("Pv");
fits_pmax_r = load_fits("Pr");


%% Interpolate to 650 nm

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

%% Plot 650 nm maps with maria shape

% plot_map(albedo_650, ' $A$', latitude, longitude, true, 'Albedo650', true)
%plot_map(pmax_650, ' $A$', latitude, longitude, true, 'Albedo650', true)

%% Calculate grain size

a = 0.845;

size_pixels = size(latitude);
grain_size_SO92 = zeros(size_pixels);
for i = 1:size_pixels(1)
    for j = 1:size_pixels(2)
        grain_size_SO92(i,j) = log10(0.03 * exp( 2.9 * ( albedo_650(i,j) + a * pmax_650(i,j) ) ) );
    end
end

%% Plot grain size

figure('Position', [500 500 900 900])
axesm vperspec
geoshow(latitude, longitude, grain_size_SO92, 'DisplayType','texturemap')

crameri lajolla
cbar = colorbar('southoutside');

mid = (cbar.Limits(1) + cbar.Limits(2))/2;
%set(cbar, 'Ticks', round([cbar.Limits(1), mid, cbar.Limits(2)], 2))
%set(cbar, 'Ticks', round([1.80, 1.93, 2.06], 2))
%caxis(round(cbar.Limits, 2))
caxis([1.80 2.10])
xlabel(cbar, strcat('log d ($\mu$m)'),'Interpreter','latex' )

% Adjust width of colorbar
x1=get(gca,'position');
x=get(cbar,'Position');
x(1) = x(1) + 0.1250;
x(3) = 0.5;
set(cbar,'Position',x)
set(gca,'position',x1)

% plotm(0.6875, 23.4333, 'wd', 'MarkerSize', 20, 'LineWidth', 5)
% geoshow('/Users/aaron/thesis/Data/mare_shape/LROC_GLOBAL_MARE_180.shp', 'DisplayType','polygon','FaceColor','none','EdgeColor','k', 'LineWidth', 0.5);

axis off
set(gca, 'FontSize', 20)

exportgraphics(gca,'Figures/WP3/grain_size_crameri.png','Resolution', 300)

%% Visualise Umov's law in 650 nm images

[maria_mask, highlands_mask] = get_maria_and_highlands_mask();

albedo_maria = albedo_650(maria_mask);
% albedo_maria = albedo_maria(1:1000,:);
albedo_highlands = albedo_650(highlands_mask);
% albedo_highlands = albedo_highlands(1:1000,:);

pmax_maria = pmax_650(maria_mask);
% pmax_maria = pmax_maria(1:1000,:);
pmax_highlands = pmax_650(highlands_mask);
% pmax_highlands = pmax_highlands(1:1000,:);

logA = linspace(0.7, 1.6, 10);
logPmax = log10( (10.^((1.81 - logA)/0.724))*10); % convert to percent

figure('Position', [500 500 500 500], 'Renderer', 'painters')

plot(albedo_highlands, pmax_highlands, '.', 'Color', [0, 0.4470, 0.7410]); hold on;
plot(albedo_maria, pmax_maria, '.', 'Color',  [0.8500, 0.3250, 0.0980]);
plot(logA, logPmax, 'k', 'LineWidth', 2)
% Create a legend with 3 entries
[h,icons] = legend({'highlands', 'maria', 'Umov law'}, 'FontSize', 20);
icons = findobj(icons, 'type', 'line');
set(icons,'MarkerSize',20);

xlim([0.7 1.6])

xlabel('$\log{A}$ (\%)', 'Interpreter', 'latex')
ylabel('$\log{P_{\textrm{max}}}$ (\%)', 'Interpreter', 'latex')

grid on
set(gca, 'FontSize', 20)



figure('Position', [500 500 500 500], 'Renderer', 'painters')

plot(albedo_highlands, pmax_highlands - log10( (10.^((1.81 - albedo_highlands)/0.724))*10), '.', 'Color', [0, 0.4470, 0.7410]); hold on;
plot(albedo_maria, pmax_maria - log10( (10.^((1.81 - albedo_maria)/0.724))*10), '.', 'Color',  [0.8500, 0.3250, 0.0980]);
% Create a legend with 3 entries
[h,icons] = legend({'highlands', 'maria'}, 'FontSize', 20);
icons = findobj(icons, 'type', 'line');
set(icons,'MarkerSize',20);

xlim([0.7 1.6])

xlabel('$\log{A}$ (\%)', 'Interpreter', 'latex')
ylabel('$\Delta\log{P_{\textrm{max}}}$ (\%)', 'Interpreter', 'latex')

grid on
set(gca, 'FontSize', 20)



%% Compute grain size for different locations

A11_latitude = 0.6875;
A11_longitude = 23.4333;

% latmask = latitude == min(latitude(latitude >= A11_latitude));
% lonmask = longitude == min(longitude(longitude >= A11_longitude));

latmask = latitude >= A11_latitude;
lonmask = longitude >= A11_longitude;
pos = latmask & lonmask;
max(max(pos));


