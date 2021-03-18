clear all; close all; clc;

roi = shaperead('mare_shape/LROC_GLOBAL_MARE_180.shp');

% latitude = fitsread("polarisation_data/latitude.fits");
% longitude = fitsread("polarisation_data/longitude.fits");

latitude = linspace(90, -90, 181);
longitude = linspace(-180, 180, 361);

[longitude, latitude] = meshgrid(longitude, latitude);

maria_mask = zeros(size(latitude),'logical');
for i = 1:length(roi)
    i=i
    rx = roi(i).X;
    ry = roi(i).Y;
    mask_temp = inpolygon(longitude,latitude,rx,ry);
    maria_mask = maria_mask | mask_temp;
end

% % CHECK
% save('maria_mask2.txt', 'maria_mask')
csvwrite('mask.txt',maria_mask)