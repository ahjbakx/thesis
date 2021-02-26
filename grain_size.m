clear all; close all; clc

%% Load data

latitude = load_fits("latitude");
longitude = load_fits("longitude");

fits_albedo_v = load_fits("Av");
fits_albedo_r = load_fits("Ar");

fits_pmax_v = load_fits("Pv");
fits_pmax_r = load_fits("Pr");

albedo_b = log10(100*load_fits("Ab")); % units of log% percent
pmax_b = log10(1000*load_fits("Pb")); % units of log%% permilli

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

%% Estimate regression line with least-squares fit

x = -albedo_630(longitude>-100);
y = pmax_630(longitude>-100);

H = ones(length(x), 2);
for i = 1:length(x)
    H(i,1) = x(i);
end

xhat = inv((transpose(H)*H))*transpose(H)*y;
yhat = H*xhat;

%plot(-x, yhat)

a = 1/xhat(1);
b = a * xhat(2);

py = cov(y);
px = inv(transpose(H)*inv(py)*H);

dx1 = sqrt(px(1,1));
da = xhat(1)^2*dx1;

dx2 = sqrt(px(2,2));
db = a*xhat(2)*sqrt( (da/a)^2 + (dx2/xhat(2))^2 );



%% Split into total, maria and highlands

[maria_mask, highlands_mask] = get_maria_and_highlands_mask();

latitude_total = latitude(longitude>-99);
latitude_maria = latitude(maria_mask);
latitude_highlands = latitude(highlands_mask);

longitude_total = longitude(longitude>-99);
longitude_maria = longitude(maria_mask);
longitude_highlands = longitude(highlands_mask);

albedo_total = albedo_630(longitude>-99);
albedo_maria = albedo_630(maria_mask);
albedo_highlands = albedo_630(highlands_mask);

pmax_total = pmax_630(longitude>-99);
pmax_maria = pmax_630(maria_mask);
pmax_highlands = pmax_630(highlands_mask);

maria = [albedo_maria, pmax_maria];
highlands = [albedo_highlands, pmax_highlands];

%% Visualise 630 nm data with regression line

logA = linspace(0.7, 1.6, 10);
logPmax = log10( (10.^((1.871 - logA)/0.7950))*10); % convert to percent


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
plot3(logA, logPmax, linspace(max(max(N_highlands)),max(max(N_highlands)),length(logA)),'--k')
plot3(-x, yhat, linspace(max(max(N_highlands)),max(max(N_highlands)),length(x)), 'k')


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
p1 = plot3(logA, logPmax, linspace(max(max(N_maria)),max(max(N_maria)),length(logA)), '--k');
p2 = plot3(-x, yhat, linspace(max(max(N_maria)),max(max(N_maria)),length(x)), 'k');
legend([p1 p2],'Shkuratov & Opanasenko (1992)', 'Best LLS fit')



%% Plot 630 nm maps with maria shape

% plot_map(albedo_630, ' $A$', latitude, longitude, true, 'Albedo650', true)
%plot_map(pmax_630, ' $A$', latitude, longitude, true, 'Albedo650', true)

%% Calculate grain size for 630 nm

a_SO92 = 0.845;
grain_size_SO92 = get_grain_size(albedo_630, pmax_630, a_SO92);

grain_size_self = get_grain_size(albedo_630, pmax_630, a);

grain_size_total = grain_size_self(longitude>-99);
grain_size_maria = grain_size_self(maria_mask);
grain_size_highlands = grain_size_self(highlands_mask);

%% Plot polarimetric anomaly

size_pixels = size(latitude);
polarimetric_anomaly = zeros(size_pixels);
for i = 1:size_pixels(1)
    for j = 1:size_pixels(2)
        polarimetric_anomaly(i,j) = 10^(pmax_b(i,j)^0.795)*10^albedo_b(i,j); % for B band
    end
end

close all;
figure('Position', [500 500 900 900])
axesm vperspec
geoshow(latitude, longitude, polarimetric_anomaly, 'DisplayType','texturemap')

colormap gray
caxis([300 710])

axis off
set(gca, 'FontSize', 20)

% exportgraphics(gca,'Figures/WP3/polarimetric_anomaly.png','Resolution', 300)

%% Visualise Umov's law in 630 nm images



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

%% Visualise Grain size as function of selenographic latitude
close all;

grain_size_latitude_total = [grain_size_total(abs(latitude_total) < 15); grain_size_total(abs(latitude_total) > 15 & abs(latitude_total) < 30 ); grain_size_total(abs(latitude_total) > 30 & abs(latitude_total) < 50); grain_size_total(abs(latitude_total) > 50);
    grain_size_maria(abs(latitude_maria) < 15); grain_size_maria( abs(latitude_maria) > 15 & abs(latitude_maria) < 30); grain_size_maria(abs(latitude_maria) > 30 & abs(latitude_maria) < 50); grain_size_maria(abs(latitude_maria) > 50);
    grain_size_highlands(abs(latitude_highlands) < 15); grain_size_highlands(abs(latitude_highlands) > 15 & abs(latitude_highlands) < 30); grain_size_highlands(abs(latitude_highlands) > 30 & abs(latitude_highlands) < 50); grain_size_highlands(abs(latitude_highlands) > 50)];

latitudes_total = [repmat([1], length(grain_size_total(abs(latitude_total) < 15)), 1);
                   repmat([5], length(grain_size_total(abs(latitude_total) > 15 & abs(latitude_total) < 30  )), 1);
                   repmat([9], length(grain_size_total(abs(latitude_total) > 30 & abs(latitude_total) < 50 )), 1);
                   repmat([13], length(grain_size_total(abs(latitude_total) > 50)), 1);
                   
                   repmat([2], length(grain_size_maria(abs(latitude_maria) < 15)), 1);
                   repmat([6], length(grain_size_maria(abs(latitude_maria) > 15 & abs(latitude_maria) < 30 )), 1);
                   repmat([10], length(grain_size_maria(abs(latitude_maria) > 30 & abs(latitude_maria) < 50 )), 1);
                   repmat([14], length(grain_size_maria(abs(latitude_maria) > 50)), 1);
                   
                   repmat([3], length(grain_size_highlands(abs(latitude_highlands) < 15)), 1);
                   repmat([7], length(grain_size_highlands(abs(latitude_highlands) > 15 & abs(latitude_highlands) < 30 )), 1);
                   repmat([11], length(grain_size_highlands(abs(latitude_highlands) > 30 & abs(latitude_highlands) < 50 )), 1);
                   repmat([15], length(grain_size_highlands(abs(latitude_highlands) > 50)), 1)];
            
positions = [1 1.25 1.5 2 2.25 2.5 3 3.25 3.5 4 4.25 4.5];

figure('Position', [500 500 900 450])

b = boxplot(grain_size_latitude_total,latitudes_total, 'positions', positions, 'Notch','on', 'Symbol','', 'Whisker', 0);

set(gca,'XTickLabel',{'',[char(946) ' < 15', char(176)], '', '', ['15',char(176),' < ', char(946) ' < 30',char(176)], '', '', ['30',char(176),' < ', char(946) ' < 50',char(176)], '', '', [char(946) ' > 50',char(176)], ''} )

color = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0 0 0]; [0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0 0 0]; [0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0 0 0]; [0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0 0 0]];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   set(h(j), 'Color', color(j,:))
end
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
for j = 1:length(lines)
    set(lines(j), 'Color', color(j,:));
end

% c = get(gca, 'Children');
% 
[~, hobj, ~, ~] = legend([b(4,4), b(5,5), b(5,6)], {'total', 'maria', 'highlands'}, 'Position',[0.15 0.7 0.2 0.2]  );
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
ht = findobj(hobj,'type','text');
set(ht,'FontSize',20);

ylabel(['Median grain size, ', char(181), 'm'])
set(gca, 'FontSize', 20)
ylim([65 100])

% exportgraphics(gca,'Figures/WP3/latitude-dependency.png','Resolution', 300)

%% Compute grain size for different locations

% Apollo 11
A11_longitude = 23.4333;
A11_latitude = 0.6875;
[g11, dg11, A11, Pmax11, dA11, dPmax11] = get_grain_size_at_lat_lon(grain_size_self, a, da, A11_latitude, A11_longitude, albedo_630, pmax_630);

% Apollo 12
A12_longitude = -23.3856;
A12_latitude = -3.1975;
[g12, dg12, A12, Pmax12, dA12, dPmax12] = get_grain_size_at_lat_lon(grain_size_self, a, da, A12_latitude, A12_longitude, albedo_630, pmax_630);

% Apollo 14
A14_longitude = -17.4653;
A14_latitude = -3.6733;
[g14, dg14, A14, Pmax14, dA14, dPmax14] = get_grain_size_at_lat_lon(grain_size_self, a, da, A14_latitude, A14_longitude, albedo_630, pmax_630);

% Apollo 15
A15_longitude = 3.6527;
A15_latitude = 26.1008;
[g15, dg15, A15, Pmax15, dA15, dPmax15] = get_grain_size_at_lat_lon(grain_size_self, a, da, A15_latitude, A15_longitude, albedo_630, pmax_630);

% Apollo 16
A16_longitude = 15.5144;
A16_latitude = -8.9913;
[g16, dg16, A16, Pmax16, dA16, dPmax16] = get_grain_size_at_lat_lon(grain_size_self, a, da, A16_latitude, A16_longitude, albedo_630, pmax_630);

% Apollo 17
A17_longitude = 30.7658;
A17_latitude = 20.1653;
[g17, dg17, A17, Pmax17, dA17, dPmax17] = get_grain_size_at_lat_lon(grain_size_self, a, da, A17_latitude, A17_longitude, albedo_630, pmax_630);


%% Plot grain size

figure('Position', [500 500 900 900])
axesm vperspec
geoshow(latitude, longitude, grain_size_self, 'DisplayType','texturemap')


% Apollo landing sites
plotm(A11_latitude, A11_longitude, '.w', 'MarkerSize', 20);
textm(A11_latitude, A11_longitude,'  A11', 'FontSize', 20, 'Color', 'w')

plotm(A12_latitude, A12_longitude, '.w', 'MarkerSize', 20);
textm(A12_latitude, A12_longitude,'A12  ', 'FontSize', 20, 'HorizontalAlignment', 'right', 'Color', 'w')

plotm(A14_latitude, A14_longitude, '.w', 'MarkerSize', 20);
textm(A14_latitude, A14_longitude,'  A14', 'FontSize', 20, 'Color', 'w')

plotm(A15_latitude, A15_longitude, '.w', 'MarkerSize', 20);
textm(A15_latitude, A15_longitude,'  A15', 'FontSize', 20, 'Color', 'w')


plotm(A16_latitude, A16_longitude, '.w', 'MarkerSize', 20);
textm(A16_latitude, A16_longitude,'  A16', 'FontSize', 20, 'Color', 'w')

plotm(A17_latitude, A17_longitude, '.w', 'MarkerSize', 20);
textm(A17_latitude, A17_longitude,'  A17', 'FontSize', 20, 'Color', 'w')

plotm(23.7, -47.4, 'ow', 'MarkerSize', 20);
textm(23.7, -47.4,'   Aristarchus', 'FontSize', 20, 'Color', 'w')

plotm(9.7, -20.1, 'ow', 'MarkerSize', 30);
textm(9.7, -20.1,'   Copernicus', 'FontSize', 20, 'Color', 'w')

plotm(-43.4, -11.1, 'ow', 'MarkerSize', 30);
textm(-43.4, -11.1,'   Tycho', 'FontSize', 20, 'Color', 'w')

plotm(8.1,-38, 'ow', 'MarkerSize', 30);
textm(8.1,-38,'Kepler   ', 'FontSize', 20, 'HorizontalAlignment', 'right', 'Color', 'w')

plotm(-24.7, -65.3, 'ow', 'MarkerSize', 30);
textm(-24.7, -65.3, '   Byrgius', 'FontSize', 20, 'Color', 'w')


plotm(-8.9, 61.1, 'ow', 'MarkerSize', 30);
textm(-8.9, 61.1, 'Langrenus   ', 'FontSize', 20, 'HorizontalAlignment', 'right', 'Color', 'w')


crameri nuuk
cbar = colorbar('southoutside');

caxis([50 110])
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

%% Local plots
close all;

% Apollo 15 landing site
lat_min = 9.4308;
lat_max = 42.7708;
lon_min = -13.0174;
lon_max = 20.3226;

local_mask = latitude > lat_min & latitude < lat_max & longitude > lon_min & longitude < lon_max;

figure('Position', [500 500 900 900])
axesm('mollweid', 'MapLatLimit',[lat_min lat_max]', 'MapLonLimit', [lon_min lon_max])
geoshow(latitude, longitude, grain_size_self, 'DisplayType','texturemap')

crameri nuuk
cbar = colorbar('southoutside');


caxis([50 110])
xlabel(cbar, strcat('Median grain size, $\mu$m'),'Interpreter','latex' )

% Adjust width of colorbar
x1=get(gca,'position');
x=get(cbar,'Position');
x(1) = x(1)+0.08;
x(3) = 0.5;
set(cbar,'Position',x)
set(gca,'position',x1)

axis off
set(gca, 'FontSize', 20)

%% Functions

function grain_size = get_grain_size(albedo_630, pmax_630, a)
      
    % albedo and pmax in natural units
    
    latitude = load_fits("latitude");
    longitude = load_fits("longitude");
    
    size_pixels = size(latitude);
    grain_size = zeros(size_pixels);
    for i = 1:size_pixels(1)
        for j = 1:size_pixels(2)
            grain_size(i,j) = 0.03 * exp( 2.9 * ( albedo_630(i,j) + a * pmax_630(i,j) ) );
        end
    end

end

function [g, dg, A, Pmax, dA, dPmax] = get_grain_size_at_lat_lon(grain_size, a, da, target_latitude, target_longitude, albedo, pmax)
    
    latitude = load_fits("latitude");
    longitude = load_fits("longitude");

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
    
    A = 10^albedo(lon, lat);
    dA = 0.8/100*A;
    Pmax = 10^pmax(lon, lat);
    dPmax = 0.8/100*Pmax;
    g = grain_size(lon,lat);
    dg = get_uncertainty(10^(latitude(lon, lat)), 10^(longitude(lon,lat)), g, a, da);
end

function dd = get_uncertainty(A, Pmax, d, a, da)
    % A, Pmax in %
    % d in m
    dA = 0.8/100*A;
    dLogA = dA/(log(10)*A);

    dPmax = 0.8/100*Pmax;
    dLogPmax = dPmax/(log(10)*Pmax);
    daLogPmax = a * log10(Pmax)*sqrt( (da/a)^2 + (dLogPmax/log10(Pmax))^2);

    dd = d * 2.9 * sqrt( (dLogA)^2 + (daLogPmax)^2);

end