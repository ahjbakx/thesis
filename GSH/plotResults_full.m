% plot results
clear;
close all;
clc;

%% load data
load('Results/data_moon_test_1_1200_27-Jan-2021 17:13:11')


%%
lon = data.grd.lon(1,:);
lats = data.grd.lat(:,1);

load coast;

%%

[longitude, latitude] = meshgrid(lon, lats);

figure('Position', [500 500 900 900])
axesm mollweid
geoshow(latitude, longitude, data.pot, 'DisplayType','texturemap')
%geoshow('/Users/aaron/thesis/Data/mare_shape/LROC_GLOBAL_MARE_180.shp', 'DisplayType','polygon','FaceColor','none','EdgeColor','k', 'LineWidth',1.5)

% % Gray colormap
% colormap gray

% Position colorbar below Moon
cbar = colorbar('southoutside');

axis off
set(gca, 'FontSize', 20)

%%

figure;
subplot(2,2,1)
hold on
imagesc(lon,lats,((data.pot)));c=colorbar; 
%plot(long,lat,'k','LineWidth',1.5);
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Potential gravity field'])
ylabel(c,'m*m/s/s') 

subplot(2,2,2)
hold on
imagesc(lon,lats,((data.vec.R)).*1e5);c=colorbar; 
%plot(long,lat,'k','LineWidth',1.5);
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['R-component of gravity vector'])
ylabel(c,'mGal') 

subplot(2,2,3)
hold on
imagesc(lon,lats,((data.vec.T)).*1e5);c=colorbar; 
%plot(long,lat,'k','LineWidth',1.5);
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['T-component of gravity vector (North-South)'])
ylabel(c,'mGal') 

subplot(2,2,4)
hold on
imagesc(lon,lats,((data.vec.L)).*1e5);c=colorbar; 
%plot(long,lat,'k','LineWidth',1.5);
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['L-component of gravity vector (East-West)'])
ylabel(c,'mGal') 

%% Tensor

figure;
subplot(3,3,1)
hold on
imagesc(lon,lats,((data.ten.Trr).*1e9));c=colorbar; 
%plot(long,lat,'k','LineWidth',1.5);
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Trr-component of gravity gradient tensor'])
ylabel(c,'Eotvos') 

subplot(3,3,2)
hold on
imagesc(lon,lats,((data.ten.Trt).*1e9));c=colorbar; 
%plot(long,lat,'k','LineWidth',1.5);
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Trt-component of gravity gradient tensor'])
ylabel(c,'Eotvos') 

subplot(3,3,3)
hold on
imagesc(lon,lats,((data.ten.Trl).*1e9));c=colorbar; 
%plot(long,lat,'k','LineWidth',1.5);
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Trl-component of gravity gradient tensor'])
ylabel(c,'Eotvos') 

subplot(3,3,5)
hold on
imagesc(lon,lats,((data.ten.Ttt).*1e9));c=colorbar; 
%plot(long,lat,'k','LineWidth',1.5);
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Ttt-component of gravity gradient tensor'])
ylabel(c,'Eotvos') 

subplot(3,3,6)
hold on
imagesc(lon,lats,((data.ten.Ttl).*1e9));c=colorbar; 
%plot(long,lat,'k','LineWidth',1.5);
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Ttl-component of gravity gradient tensor'])
ylabel(c,'Eotvos') 

subplot(3,3,9)
hold on
imagesc(lon,lats,((data.ten.Tll).*1e9));c=colorbar; 
%plot(long,lat,'k','LineWidth',1.5);
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Tll-component of gravity gradient tensor'])
ylabel(c,'Eotvos') 

figure
hold on
imagesc(lon,lats,((data.vec.R).*1e5));c=colorbar; 
%plot(long,lat,'k','LineWidth',1.5);
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Radial vector component of gravity field'])
ylabel(c,'mGal') 

%%

figure
hold on
imagesc(lon,lats,((data.ten.Txx+data.ten.Tyy+data.ten.Tzz).*1e9));c=colorbar; 
%plot(long,lat,'k','LineWidth',1.5);
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['TraceTensor'])
ylabel(c,'Eotvos') 

