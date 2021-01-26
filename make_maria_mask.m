clear all; close all; clc;

roi = shaperead('mare_shape/LROC_GLOBAL_MARE_180.shp');

latitude = fitsread("polarisation_data/latitude.fits");
longitude = fitsread("polarisation_data/longitude.fits");

maria_mask = zeros(size(latitude),'logical');
for i = 1:length(roi)
    i=i
    rx = roi(i).X;
    ry = roi(i).Y;
    mask_temp = inpolygon(longitude,latitude,rx,ry);
    maria_mask = maria_mask | mask_temp;
end

% CHECK
save('bin/maria_mask.mat', 'maria_mask')