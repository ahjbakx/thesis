clear all; close all; clc;

% load moonalb20c
% 
% axesm ortho
% geoshow(moonalb20c,moonalb20cR,'DisplayType','texturemap')
% colormap gray
% axis off

% t = Tiff("Lunar_Clementine_UVVIS_FeO_ClrBinned_70S70N_1km.tif");
% data = read(t);



% Get lat/lon of a single pixel
% [x,y] = pix2map(info.RefMatrix, 1000, 1000);
% [lat,lon] = projinv(info, x,y);






% % get lat/lon of geotiff
% [AY, AX]= size(data);
% lon = zeros(size(data));
% lat = zeros(size(data));
% x = zeros(1, AX);
% y = zeros(1, AY);
% height = info.Height; % Integer indicating the height of the image in pixels
% width = info.Width; % Integer indicating the width of the image in pixels
% [rows,cols] = meshgrid(1:height,1:width);
% % Getting the latitude and longitude of the points
% [x,y] = pix2map(info.RefMatrix, rows, cols);
% [lat,lon] = projinv(info, x,y);
% lat = lat';
% lon = lon';
% 
% [test1, test2] = mfwdtran(lat, lon)


% figure
% geoshow(latitude, longitude, data);
% %axesm eqdcylin
% 
% axesm('vperspec','FLatLimit', [],'MapLonLimit',[-70 70])
% axis off
% set(gca, 'FontSize', 20)


% Vq = zeros(size(latitude));
% [Xq,Yq]=meshgrid(-70.0:0.5:70.0,-70.0:0.5:70.0); 
% 
% for i=1:3
%     Vq(:,:,i) = interp2(lat,lon,double(data(:,:,i)), Xq,Yq);
% end

% tif = 'Lunar_Clementine_UVVIS_FeO_ClrBinned_70S70N_1km.tif';
tif = 'Lunar_Clementine_UVVIS_750nm_Global_Mosaic_118m_v2.1';

info = geotiffinfo( tif );
[data, R] = geotiffread( tif );

data_reduced = imresize(data, size(data)/10, 'bicubic');

R_reduced = R;
R_reduced.RasterSize = size(data_reduced);

figure
mapshow(data_reduced, R_reduced)

% latitude = fitsread("data/latitude.fits");
% longitude = fitsread("data/longitude.fits");
% 
% % Get pixel centers
% [x,y] = pixcenters(R,size(data));
% [X,Y] = meshgrid(x,y);
% 
% % Get geocoordinates
% mstruct = geotiff2mstruct(info);
% [lat,lon] = minvtran(mstruct,X,Y);
% 
% % P = 'vperspec';
% % axesm('MapProjection',P);
% % h = pcolorm(lat,lon,data);
% % axis off
% % colorbar

% data = inpaintn(data); 
% 
% FeO = linspace(7, 25); 
% [i,j,k] = meshgrid(1:10918,1:4246, 1:3);
% FeOmap = interp3(double(data(:,:,1)),double(data(:,:,2)),double(data(:,:,3)), FeO, i(:,:,1),j(:,:,2),k(:,:,3));

% [IND,map] = rgb2ind(data,32);
% figure
% imagesc(IND)
% colormap(map)
% axis image
% colorbar
% 

