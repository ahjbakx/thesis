clear all; close all; clc

%% Load data

latitude = load_fits("latitude");
longitude = load_fits("longitude");

fits_albedo_v = load_fits("Av");
fits_albedo_r = load_fits("Ar");

fits_pmax_v = load_fits("Pv");
fits_pmax_r = load_fits("Pr");

albedo_b = log10(100*load_fits("Ab"));
pmax_b = log10(1000*load_fits("Pb"));

%% Interpolate to 630 nm

% interpolate albedo to 630 nm with units log(percentage)
median_albedo_v = median(median(fits_albedo_v, 'omitnan'), 'omitnan');
median_albedo_r = median(median(fits_albedo_r, 'omitnan'), 'omitnan');
albedo_630 = log10( 100 * (fits_albedo_v + (median_albedo_r - median_albedo_v) / (676.3 - 558.6) * (630 - 558.6)) );

% interpolate pmax to 630 nm with units log(permilli)
size_pmax = size(fits_pmax_v);
pmax_630 = zeros(size_pmax);
for i = 1:size_pmax(1)
    for j = 1:size_pmax(2)
        slope = (fits_pmax_r(i,j) - fits_pmax_v(i,j) ) / (676.3 - 558.6);
        pmax_630(i,j) = log10( 1000 * (fits_pmax_v(i,j) + slope * (630 - 558.6)) );
    end
end

%% Plot 630 nm maps with maria shape

% plot_map(albedo_630, ' $A$', latitude, longitude, true, 'Albedo650', true)
%plot_map(pmax_630, ' $A$', latitude, longitude, true, 'Albedo650', true)

%% Calculate grain size for 630 nm

a = 0.845;

size_pixels = size(latitude);
grain_size_SO92 = zeros(size_pixels);
polarimetric_anomaly = zeros(size_pixels);
for i = 1:size_pixels(1)
    for j = 1:size_pixels(2)
        polarimetric_anomaly(i,j) = 10^(pmax_b(i,j)^0.795)*10^albedo_b(i,j); % for B band
        grain_size_SO92(i,j) = 0.03 * exp( 2.9 * ( albedo_630(i,j) + a * pmax_630(i,j) ) );
    end
end

%% Plot polarimetric anomaly
close all;
figure('Position', [500 500 900 900])
axesm vperspec
geoshow(latitude, longitude, polarimetric_anomaly, 'DisplayType','texturemap')

colormap gray
% cbar = colorbar('southoutside');

% mid = (cbar.Limits(1) + cbar.Limits(2))/2;
%set(cbar, 'Ticks', round([cbar.Limits(1), mid, cbar.Limits(2)], 2))
%set(cbar, 'Ticks', round([1.80, 1.93, 2.06], 2))
%caxis(round(cbar.Limits, 2))
% xlabel(cbar, '$b=(P_{\textrm{max}})^aA$','Interpreter','latex' )
caxis([300 800])
% Adjust width of colorbar
% x1=get(gca,'position');
% x=get(cbar,'Position');
% x(1) = x(1) + 0.1250;
% x(3) = 0.5;
% set(cbar,'Position',x)
% set(gca,'position',x1)

% plotm(0.6875, 23.4333, 'wd', 'MarkerSize', 20, 'LineWidth', 5)
% geoshow('/Users/aaron/thesis/Data/mare_shape/LROC_GLOBAL_MARE_180.shp', 'DisplayType','polygon','FaceColor','none','EdgeColor','k', 'LineWidth', 0.5);

axis off
set(gca, 'FontSize', 20)

% exportgraphics(gca,'Figures/WP3/polarimetric_anomaly.png','Resolution', 300)

%% Visualise Umov's law in 630 nm images

[maria_mask, highlands_mask] = get_maria_and_highlands_mask();

albedo_maria = albedo_630(maria_mask);
% albedo_maria = albedo_maria(1:1000,:);
albedo_highlands = albedo_630(highlands_mask);
% albedo_highlands = albedo_highlands(1:1000,:);

pmax_maria = pmax_630(maria_mask);
% pmax_maria = pmax_maria(1:1000,:);
pmax_highlands = pmax_630(highlands_mask);
% pmax_highlands = pmax_highlands(1:1000,:);

logA = linspace(0.7, 1.6, 10);
logPmax = log10( (10.^((1.871 - logA)/0.7950))*10); % convert to percent

% figure('Position', [500 500 500 500], 'Renderer', 'painters')
% 
% plot(albedo_highlands, pmax_highlands, '.', 'Color', [0, 0.4470, 0.7410]); hold on;
% plot(albedo_maria, pmax_maria, '.', 'Color',  [0.8500, 0.3250, 0.0980]);
% plot(logA, logPmax, 'k', 'LineWidth', 2)
% % Create a legend with 3 entries
% [h,icons] = legend({'highlands', 'maria', 'Umov law'}, 'FontSize', 20);
% icons = findobj(icons, 'type', 'line');
% set(icons,'MarkerSize',20);
% 
% xlim([0.7 1.6])
% 
% xlabel('$\log{A}$ (\%)', 'Interpreter', 'latex')
% ylabel('$\log{P_{\textrm{max}}}$ (\%)', 'Interpreter', 'latex')
% 
% grid on
% set(gca, 'FontSize', 20)


%% Deviation from Umov law
% figure('Position', [500 500 500 500], 'Renderer', 'painters')
% 
% plot(albedo_highlands, pmax_highlands - log10( (10.^((1.81 - albedo_highlands)/0.845))*10), '.', 'Color', [0, 0.4470, 0.7410]); hold on;
% plot(albedo_maria, pmax_maria - log10( (10.^((1.81 - albedo_maria)/0.845))*10), '.', 'Color',  [0.8500, 0.3250, 0.0980]);
% % Create a legend with 3 entries
% % [h,icons] = legend({'highlands', 'maria'}, 'FontSize', 20, 'Location', 'northwest');
% % icons = findobj(icons, 'type', 'line');
% % set(icons,'MarkerSize',20);
% 
% xlim([0.7 1.6])
% ylim([-0.4 0.7])
% 
% xlabel('$\log{A}$ (\%)', 'Interpreter', 'latex')
% ylabel('$\Delta\log{P_{\textrm{max}}}$ (\%)', 'Interpreter', 'latex')
% 
% grid on
% set(gca, 'FontSize', 20)

%% 3D histogram
close all;

maria = [albedo_maria, pmax_maria];
highlands = [albedo_highlands, pmax_highlands];

% Specify number of bins
nbins=[200 200];

figure('Position', [500 500 500 500], 'Renderer', 'painters')
[N_highlands, C_highlands] = hist3(highlands, nbins);
hist3(highlands, nbins, 'CDataMode','auto','FaceColor','interp', 'EdgeColor', 'none');
view(0,90)
crameri bilbao
xlabel('$\log{A}$ (\%)', 'Interpreter', 'latex')
ylabel('$\log{P_{\textrm{max}}}$ (\%)', 'Interpreter', 'latex')
set(gca, 'FontSize', 20)
xlim([0.8 1.6])
ylim([1.5 2.3])
xl = xlim;
yl = ylim;
xticks([0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7])
grid off; box on;
hold on;
plot3(logA, logPmax, linspace(max(max(N_highlands)),max(max(N_highlands)),length(logA)), '--k')

figure('Position', [500 500 500 500], 'Renderer', 'painters')
[N_maria, C_maria] = hist3(maria, nbins);
hist3(maria, nbins, 'CDataMode','auto','FaceColor','interp', 'EdgeColor', 'none');
view(0,90)
crameri bilbao
xlabel('$\log{A}$ (\%)', 'Interpreter', 'latex')
ylabel('$\log{P_{\textrm{max}}}$ (\%)', 'Interpreter', 'latex')
set(gca, 'FontSize', 20)
xlim(xl)
ylim(yl)
xticks([0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7])
grid off; box on;
hold on;
plot3(logA, logPmax, linspace(max(max(N_maria)),max(max(N_maria)),length(logA)), '--k')



%% Compute grain size for different locations

% Apollo 11
target_longitude = 23.4333;
target_latitude = 0.6875;

% % Apollo 12
% target_longitude = -23.3856;
% target_latitude = -3.1975;
% 
% % Apollo 14
% target_longitude = -17.4653;
% target_latitude = -3.6733;

% % Apollo 15
% target_longitude = 3.6527;
% target_latitude = 26.1008;

% Point of interest
P = [target_longitude, target_latitude];

% Radius of Moon
RM = 1737;

% Convert the array of lat/lon coordinates to Cartesian vectors
% ph2cart expects radians
% use radius 1, so normalisation of vectors not needed
[X,Y,Z] = sph2cart( longitude*pi/180,  latitude*pi/180, 1);

% Same for point of interest    
[xP,yP,zP] = sph2cart(P(1)*pi/180, P(2)*pi/180, 1);

% The minimum distance, and the linear index where that distance was found
% force the dot product into the interval [-1 +1]. This prevents 
% slight overshoots due to numerical artifacts
dotProd = xP*X(:) + yP*Y(:) + zP*Z(:);
[minDist, index] = min( RM*acos( min(max(-1,dotProd),1) ) );

% Convert that linear index to 2D subscripts
[lon,lat] = ind2sub(size(longitude), index);

grain_size_SO92(lon,lat);

%% Plot grain size

figure('Position', [500 500 900 900])
axesm vperspec
geoshow(latitude, longitude, grain_size_SO92, 'DisplayType','texturemap')


% Apollo landing sites
plotm(0.6875, 23.4333, 'xw', 'MarkerSize', 20);
textm(0.6875, 23.4333,'  A11', 'FontSize', 20, 'Color', 'w')

plotm(-3.1975, -23.3856, 'xw', 'MarkerSize', 20);
textm(-3.1975, -23.3856,'A12  ', 'FontSize', 20, 'HorizontalAlignment', 'right', 'Color', 'w')

plotm(-3.6733, -17.4653, 'xw', 'MarkerSize', 20);
textm(-3.6733, -17.4653,'  A14', 'FontSize', 20, 'Color', 'w')

plotm(26.1008, 3.6527, 'xw', 'MarkerSize', 20);
textm(26.1008, 3.6527,'  A15', 'FontSize', 20, 'Color', 'w')

crameri lajolla
cbar = colorbar('southoutside');

mid = (cbar.Limits(1) + cbar.Limits(2))/2;
%set(cbar, 'Ticks', round([cbar.Limits(1), mid, cbar.Limits(2)], 2))
%set(cbar, 'Ticks', round([1.80, 1.93, 2.06], 2))
%caxis(round(cbar.Limits, 2))
% caxis([1.80 2.10])
caxis([60 130])
xlabel(cbar, strcat('Median grain size, $\mu$m'),'Interpreter','latex' )

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

% exportgraphics(gca,'Figures/WP3/grain_size.png','Resolution', 300)