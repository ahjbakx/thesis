clear all; close all; clc

latitude = fitsread("polarisation_data/latitude.fits");
longitude = fitsread("polarisation_data/longitude.fits");


% plot_map(albedo, '$A$', latitude, longitude, false, 'AlbedoR')
% plot_map(pmax, '$P_{\textrm{max}}$', latitude, longitude, false, 'PmaxR')

boxplot_albedo()
boxplot_pmax()
% contour_plots()

%% Plot maria & highlands
% [maria_mask, highlands_mask] = get_maria_and_highlands_mask();
% 
% albedo = load_fits_log("Av");
% size_albedo = size(albedo);
% albedo_maria = zeros(size_albedo);
% albedo_highlands = zeros(size_albedo);
% for i = 1:size_albedo(1)
%     for j = 1:size_albedo(2)
%         albedo_maria(i,j) = albedo(i,j)*maria_mask(i,j);
%         albedo_highlands(i,j) = albedo(i,j)*highlands_mask(i,j);
%     end
% end
% 
% plot_map(albedo_maria, '$A$', latitude, longitude, false, '')
% plot_map(albedo_highlands, '$A$', latitude, longitude, false, '')



%% Functions
 
function fitslog = load_fits_log(name)

    fits = fitsread( strcat("polarisation_data/", name, ".fits") );

    fits_size = size(fits);
    for i = 1:fits_size(1)
        for j = 1:fits_size(2)
           if fits(i,j) == -99 || fits(i,j) == 0
               fits(i,j) = 1.0e-101;
           end
        end
    end
    
   fitslog = log10(100*fits);

end

function fits = load_fits(name)

    fits = fitsread( strcat("polarisation_data/", name, ".fits") );

    fits_size = size(fits);
    for i = 1:fits_size(1)
        for j = 1:fits_size(2)
           if fits(i,j) == -99
               fits(i,j) = 1.0e-101;
           end
        end
    end
    
   fits = 100*fits;

end

function plot_map(quantity, label, latitude, longitude, save, savename)
    
    
    figure('Position', [500 500 900 900])
    axesm vperspec
    geoshow(latitude, longitude, quantity, 'DisplayType','texturemap')

    % Gray colormap
    colormap gray

    % Position colorbar below Moon
    cbar = colorbar('southoutside');

    % Colorbar labels: low, mid, high value, with log(%)
    mid = (cbar.Limits(1) + cbar.Limits(2))/2;
    set(cbar, 'Ticks', round([cbar.Limits(1), mid, cbar.Limits(2)], 2))
    caxis(round(cbar.Limits, 2))
    xlabel(cbar, strcat('log', label,' (\%)'),'Interpreter','latex' )

    % Adjust width of colorbar
    x1=get(gca,'position');
    x=get(cbar,'Position');
    x(1) = x(1) + 0.1250;
    x(3) = 0.5;
    set(cbar,'Position',x)
    set(gca,'position',x1)
    
    axis off
    set(gca, 'FontSize', 20)
    %geoshow('mare_shape/LROC_GLOBAL_MARE_180.shp', 'DisplayType','polygon','FaceColor','none','EdgeColor','w');
    
    if save == true
        tightmap
        set(gca,'LooseInset',get(gca,'TightInset'));
        saveas(gcf, strcat('Figures/', savename, '.eps'))
    end

end

function [maria_mask, highlands_mask] = get_maria_and_highlands_mask()
    latitude = fitsread("polarisation_data/latitude.fits");
    longitude = fitsread("polarisation_data/longitude.fits");
    albedo = load_fits_log("Av");
    maria_mask = load("maria_mask.mat").maria_mask & (longitude > 15 | longitude<-15) ;
    highlands_mask = (latitude>-99 & (longitude > 15 | longitude<-15) & albedo>-99) & not(maria_mask); % filter NaN and infs
end

function contour_plots(save, savename)
    % Get masks for maria and highlands
    [maria_mask, highlands_mask] = get_maria_and_highlands_mask();

    % Load albedo and pmax data
    albedo_u = load_fits_log("Au");
    pmax_u = load_fits_log("Pu");

    albedo_b = load_fits_log("Ab");
    pmax_b = load_fits_log("Pb");

    albedo_v = load_fits_log("Av");
    pmax_v = load_fits_log("Pv");

    albedo_r = load_fits_log("Ar");
    pmax_r = load_fits_log("Pr");

    % Mask data for maria and highlands
    albedo_u_maria = albedo_u(maria_mask);
    pmax_u_maria = pmax_u(maria_mask);
    albedo_u_highlands = albedo_u(highlands_mask & albedo_u >= 0);
    pmax_u_highlands = pmax_u(highlands_mask & albedo_u >= 0);

    albedo_b_maria = albedo_b(maria_mask);
    pmax_b_maria = pmax_b(maria_mask);
    albedo_b_highlands = albedo_b(highlands_mask & albedo_u >= 0);
    pmax_b_highlands = pmax_b(highlands_mask & albedo_u >= 0);

    albedo_v_maria = albedo_v(maria_mask);
    pmax_v_maria = pmax_v(maria_mask);
    albedo_v_highlands = albedo_v(highlands_mask & albedo_u >= 0);
    pmax_v_highlands = pmax_v(highlands_mask & albedo_u >= 0);

    albedo_r_maria = albedo_r(maria_mask);
    pmax_r_maria = pmax_r(maria_mask);
    albedo_r_highlands = albedo_r(highlands_mask & albedo_u >= 0);
    pmax_r_highlands = pmax_r(highlands_mask & albedo_u >= 0);

    % Combine horizontal and vertical data for histogram
    u_maria = [albedo_u_maria, pmax_u_maria];
    u_highlands = [albedo_u_highlands, pmax_u_highlands];

    b_maria = [albedo_b_maria, pmax_b_maria];
    b_highlands = [albedo_b_highlands, pmax_b_highlands];

    v_maria = [albedo_v_maria, pmax_v_maria];
    v_highlands = [albedo_v_highlands, pmax_v_highlands];

    r_maria = [albedo_r_maria, pmax_r_maria];
    r_highlands = [albedo_r_highlands, pmax_r_highlands];

    % Specify number of bins
    nbins=[150 150];

    % Create histograms
    [N_u_highlands,C_u_highlands] = hist3(u_highlands, nbins);
    [N_u_maria,C_u_maria] = hist3(u_maria, nbins);

    [N_b_highlands,C_b_highlands] = hist3(b_highlands, nbins);
    [N_b_maria,C_b_maria] = hist3(b_maria, nbins);

    [N_v_highlands,C_v_highlands] = hist3(v_highlands, nbins);
    [N_v_maria,C_v_maria] = hist3(v_maria, nbins);

    [N_r_highlands,C_r_highlands] = hist3(r_highlands, nbins);
    [N_r_maria,C_r_maria] = hist3(r_maria, nbins);

    % Calculate percentiles for contours
    percentiles_u_highlands = [prctile(nonzeros(N_u_highlands), 25, 'all'), prctile(nonzeros(N_u_highlands), 50, 'all'), prctile(nonzeros(N_u_highlands), 75, 'all')];
    percentiles_u_maria = [prctile(nonzeros(N_u_maria), 25, 'all'), prctile(nonzeros(N_u_maria), 50, 'all'), prctile(nonzeros(N_u_maria), 75, 'all')];

    percentiles_b_highlands = [prctile(nonzeros(N_b_highlands), 25, 'all'), prctile(nonzeros(N_b_highlands), 50, 'all'), prctile(nonzeros(N_b_highlands), 75, 'all')];
    percentiles_b_maria = [prctile(nonzeros(N_b_maria), 25, 'all'), prctile(nonzeros(N_b_maria), 50, 'all'), prctile(nonzeros(N_b_maria), 75, 'all')];

    percentiles_v_highlands = [prctile(nonzeros(N_v_highlands), 25, 'all'), prctile(nonzeros(N_v_highlands), 50, 'all'), prctile(nonzeros(N_v_highlands), 75, 'all')];
    percentiles_v_maria = [prctile(nonzeros(N_v_maria), 25, 'all'), prctile(nonzeros(N_v_maria), 50, 'all'), prctile(nonzeros(N_v_maria), 75, 'all')];

    percentiles_r_highlands = [prctile(nonzeros(N_r_highlands), 25, 'all'), prctile(nonzeros(N_r_highlands), 50, 'all'), prctile(nonzeros(N_r_highlands), 75, 'all')];
    percentiles_r_maria = [prctile(nonzeros(N_r_maria), 25, 'all'), prctile(nonzeros(N_r_maria), 50, 'all'), prctile(nonzeros(N_r_maria), 75, 'all')];

    % Smooth contours with convolution
    K = (1/10)*ones(3);

    N_u_highlands = conv2(N_u_highlands,K,'same');
    N_u_maria = conv2(N_u_maria,K,'same');

    N_b_highlands = conv2(N_b_highlands,K,'same');
    N_b_maria = conv2(N_b_maria,K,'same');

    N_v_highlands = conv2(N_v_highlands,K,'same');
    N_v_maria = conv2(N_v_maria,K,'same');

    N_r_highlands = conv2(N_r_highlands,K,'same');
    N_r_maria = conv2(N_r_maria,K,'same');

    % Make contour plots
    figure('Position', [500 500 1200 900])

    subplot(2,3,1)
    contour(C_u_maria{1},C_u_maria{2},N_u_maria, percentiles_u_maria, 'LineColor', 'b'); hold on;
    contour(C_u_highlands{1},C_u_highlands{2},N_u_highlands, percentiles_u_highlands, 'LineColor', 'r');

    xlabel('log A (\%)', 'Interpreter', 'latex')
    ylabel('log $P_{max}$ (\%)', 'Interpreter', 'latex')

    xlim([0.4 1.5]); ylim([0.5 1.6])
    xticks([0.4 0.6 0.8 1.0 1.2 1.4])
    grid on
    set(gca, 'FontSize', 15)

    legend('maria', 'highlands')

    subplot(2,3,2)
    contour(C_b_maria{1},C_b_maria{2},N_b_maria, percentiles_b_maria, 'LineColor', 'b'); hold on;
    contour(C_b_highlands{1},C_b_highlands{2},N_b_highlands, percentiles_b_highlands, 'LineColor', 'r');

    xlabel('log A (\%)', 'Interpreter', 'latex')
    ylabel('log $P_{max}$ (\%)', 'Interpreter', 'latex')

    xlim([0.4 1.5]); ylim([0.5 1.6])
    xticks([0.4 0.6 0.8 1.0 1.2 1.4])
    grid on
    set(gca, 'FontSize', 15)

    subplot(2,3,4)
    contour(C_v_maria{1},C_v_maria{2},N_v_maria, percentiles_v_maria, 'LineColor', 'b'); hold on;
    contour(C_v_highlands{1},C_v_highlands{2},N_v_highlands, percentiles_v_highlands, 'LineColor', 'r');

    xlabel('log A (\%)', 'Interpreter', 'latex')
    ylabel('log $P_{max}$ (\%)', 'Interpreter', 'latex')

    xlim([0.4 1.5]); ylim([0.5 1.6])
    xticks([0.4 0.6 0.8 1.0 1.2 1.4])
    grid on
    set(gca, 'FontSize', 15)

    subplot(2,3,5)
    contour(C_r_maria{1},C_r_maria{2},N_r_maria, percentiles_r_maria, 'LineColor', 'b'); hold on;
    contour(C_r_highlands{1},C_r_highlands{2},N_r_highlands, percentiles_r_highlands, 'LineColor', 'r');

    xlabel('log A (\%)', 'Interpreter', 'latex')
    ylabel('log $P_{max}$ (\%)', 'Interpreter', 'latex')

    xlim([0.4 1.5]); ylim([0.5 1.6])
    xticks([0.4 0.6 0.8 1.0 1.2 1.4])
    grid on
    set(gca, 'FontSize', 15)

    subplot(2,3,6)
    contour(C_u_maria{1},C_u_maria{2},N_u_maria, [percentiles_u_maria(2) percentiles_u_maria(2)], 'LineColor', 'b'); hold on;
    contour(C_b_maria{1},C_b_maria{2},N_b_maria, [percentiles_b_maria(2) percentiles_b_maria(2)], 'LineColor', 'b');
    contour(C_v_maria{1},C_v_maria{2},N_v_maria, [percentiles_v_maria(2) percentiles_v_maria(2)], 'LineColor', 'b');
    contour(C_r_maria{1},C_r_maria{2},N_r_maria, [percentiles_r_maria(2) percentiles_r_maria(2)], 'LineColor', 'b');

    text(0.8,1.0, .1, 'U', 'Color', 'b', 'FontSize', 15)
    text(0.9,0.9, .1, 'B', 'Color', 'b', 'FontSize', 15)
    text(1.0,0.8, .1, 'V', 'Color', 'b', 'FontSize', 15)
    text(1.1,0.7, .1, 'R', 'Color', 'b', 'FontSize', 15)
    
    contour(C_u_highlands{1},C_u_highlands{2},N_u_highlands, [percentiles_u_highlands(2) percentiles_u_highlands(2)], 'LineColor', 'r');
    contour(C_b_highlands{1},C_b_highlands{2},N_b_highlands, [percentiles_b_highlands(2) percentiles_b_highlands(2)], 'LineColor', 'r');
    contour(C_v_highlands{1},C_v_highlands{2},N_v_highlands, [percentiles_v_highlands(2) percentiles_v_highlands(2)], 'LineColor', 'r');
    contour(C_r_highlands{1},C_r_highlands{2},N_r_highlands, [percentiles_r_highlands(2) percentiles_r_highlands(2)], 'LineColor', 'r');


    xlabel('log A (\%)', 'Interpreter', 'latex')
    ylabel('log $P_{max}$ (\%)', 'Interpreter', 'latex')

    xlim([0.4 1.5]); ylim([0.5 1.6])
    xticks([0.4 0.6 0.8 1.0 1.2 1.4])
    grid on
    set(gca, 'FontSize', 15)

end

function boxplot_albedo()

    latitude = fitsread("polarisation_data/latitude.fits");
    longitude = fitsread("polarisation_data/longitude.fits");

    [maria_mask, highlands_mask] = get_maria_and_highlands_mask();
    
    % Load pmax data
    albedo_u = load_fits("Au");
    albedo_b = load_fits("Ab");
    albedo_v = load_fits("Av");
    albedo_r = load_fits("Ar");

    % Mask maria pmax
    albedo_u_maria = albedo_u(maria_mask);
    albedo_b_maria = albedo_b(maria_mask);
    albedo_v_maria = albedo_v(maria_mask);
    albedo_r_maria = albedo_r(maria_mask);

    % Mask highlands pmax
    albedo_u_highlands = albedo_u(highlands_mask);
    albedo_b_highlands = albedo_b(highlands_mask);
    albedo_v_highlands = albedo_v(highlands_mask);
    albedo_r_highlands = albedo_r(highlands_mask);

    % Construct pmax matrices
    albedo_maria = [albedo_u_maria; albedo_b_maria; albedo_v_maria; albedo_r_maria];
    albedo_highlands = [albedo_u_highlands; albedo_b_highlands; albedo_v_highlands; albedo_r_highlands];

    % Calculate mean of albedo
    mean_albedo_maria = [mean(albedo_u_maria), mean(albedo_b_maria), mean(albedo_v_maria), mean(albedo_r_maria)]; 
    mean_albedo_highlands = [mean(albedo_u_highlands), mean(albedo_b_highlands), mean(albedo_v_highlands), mean(albedo_r_highlands)]; 

    % Construct wavelength matrices
    size_albedo_maria = size(albedo_u_maria);
    wavelengths_maria = [repmat([373.8], size_albedo_maria(1), 1);
                   repmat([443.5], size_albedo_maria(1), 1);
                   repmat([558.6], size_albedo_maria(1), 1);
                   repmat([676.6], size_albedo_maria(1), 1)];

    size_albedo_highlands = size(albedo_u_highlands);
    wavelengths_highlands = [repmat([373.8], size_albedo_highlands(1), 1);
                   repmat([443.5], size_albedo_highlands(1), 1);
                   repmat([558.6], size_albedo_highlands(1), 1);
                   repmat([676.6], size_albedo_highlands(1), 1)];

    figure('Position', [500 500 900 900], 'Renderer', 'painters')
    boxplot(albedo_maria, wavelengths_maria, 'Notch','on', 'Symbol','', 'Colors', 'b'); hold on
    plot([1 2 3 4], mean_albedo_maria, '-ob')

    boxplot(albedo_highlands, wavelengths_highlands, 'Notch','on', 'Symbol','', 'Colors', 'r'); 
    plot([1 2 3 4], mean_albedo_highlands, '-or')

    xlabel('Wavelength (nm)'); ylabel('$A$ (\%)', 'Interpreter', 'latex')
    grid on
    set(gca, 'FontSize', 20)

    ylim([0 35])

    boxes = findobj(gca, 'Tag', 'Box');
    legend(boxes([end 1]), 'maria', 'highlands')
end

function boxplot_pmax()

    latitude = fitsread("polarisation_data/latitude.fits");
    longitude = fitsread("polarisation_data/longitude.fits");

    [maria_mask, highlands_mask] = get_maria_and_highlands_mask();
    
    % Load pmax data
    pmax_u = load_fits("Pu");
    pmax_b = load_fits("Pb");
    pmax_v = load_fits("Pv");
    pmax_r = load_fits("Pr");

    % Mask maria pmax
    pmax_u_maria = pmax_u(maria_mask);
    pmax_b_maria = pmax_b(maria_mask);
    pmax_v_maria = pmax_v(maria_mask);
    pmax_r_maria = pmax_r(maria_mask);

    % Mask highlands pmax
    pmax_u_highlands = pmax_u(highlands_mask);
    pmax_b_highlands = pmax_b(highlands_mask);
    pmax_v_highlands = pmax_v(highlands_mask);
    pmax_r_highlands = pmax_r(highlands_mask);

    % Construct pmax matrices
    pmax_maria = [pmax_u_maria; pmax_b_maria; pmax_v_maria; pmax_r_maria];
    pmax_highlands = [pmax_u_highlands; pmax_b_highlands; pmax_v_highlands; pmax_r_highlands];

    % Calculate mean of albedo
    mean_pmax_maria = [mean(pmax_u_maria), mean(pmax_b_maria), mean(pmax_v_maria), mean(pmax_r_maria)]; 
    mean_pmax_highlands = [mean(pmax_u_highlands), mean(pmax_b_highlands), mean(pmax_v_highlands), mean(pmax_r_highlands)]; 

    % Construct wavelength matrices
    size_pmax_maria = size(pmax_u_maria);
    wavelengths_maria = [repmat([373.8], size_pmax_maria(1), 1);
                   repmat([443.5], size_pmax_maria(1), 1);
                   repmat([558.6], size_pmax_maria(1), 1);
                   repmat([676.6], size_pmax_maria(1), 1)];

    size_pmax_highlands = size(pmax_u_highlands);
    wavelengths_highlands = [repmat([373.8], size_pmax_highlands(1), 1);
                   repmat([443.5], size_pmax_highlands(1), 1);
                   repmat([558.6], size_pmax_highlands(1), 1);
                   repmat([676.6], size_pmax_highlands(1), 1)];

    figure('Position', [500 500 900 900], 'Renderer', 'painters')
    boxplot(pmax_maria, wavelengths_maria, 'Notch','on', 'Symbol','', 'Colors', 'b'); hold on
    plot([1 2 3 4], mean_pmax_maria, '-ob')

    boxplot(pmax_highlands, wavelengths_highlands, 'Notch','on', 'Symbol','', 'Colors', 'r'); 
    plot([1 2 3 4], mean_pmax_highlands, '-or')

    xlabel('Wavelength (nm)'); ylabel('$P_{max}$ (\%)', 'Interpreter', 'latex')
    grid on
    set(gca, 'FontSize', 20)

    ylim([0 35])

    boxes = findobj(gca, 'Tag', 'Box');
    legend(boxes([end 1]), 'maria', 'highlands')
end

function area = get_area_weights()

    latitude = fitsread("polarisation_data/latitude.fits");
    longitude = fitsread("polarisation_data/longitude.fits");

    size_area = size(latitude);
    area = zeros(size_area);

    for i = 1:size_area(1)
        for j = 1:size_area(2)
             if latitude(i,j) > -99
                 if latitude(i,j-1) > -99 & longitude(i-1,j) > -99
                     beta = deg2rad(latitude(i,j));
                     lambda = deg2rad(longitude(i,j)); 
                     area(i,j) = abs(cos(beta)*cos(lambda));                
                 else
                     area(i,j) = 0;       
                 end 
             else
                 area(i,j) = 0; 
             end
        end
    end

    % figure
    % axesm vperspec
    % geoshow(latitude, longitude, area, 'DisplayType','texturemap')
    % title('Area')
    % colorbar
end


function boxplot_albedo_weighted()

    latitude = fitsread("polarisation_data/latitude.fits");
    longitude = fitsread("polarisation_data/longitude.fits");

    [maria_mask, highlands_mask] = get_maria_and_highlands_mask();
    
    area = get_area_weights();
    
    albedo_u = load_fits("Au");
    albedo_b = load_fits("Ab");
    albedo_v = load_fits("Av");
    albedo_r = load_fits("Ar");

    % pmax_u = load_fits("Pu");
    % pmax_b = load_fits("Pb");
    % pmax_v = load_fits("Pv");
    % pmax_r = load_fits("Pr");

    % Mask maria albedo
    albedo_u_maria = albedo_u(maria_mask);
    albedo_b_maria = albedo_b(maria_mask);
    albedo_v_maria = albedo_v(maria_mask);
    albedo_r_maria = albedo_r(maria_mask);

    % Mask highlands albedo
    albedo_u_highlands = albedo_u(highlands_mask);
    albedo_b_highlands = albedo_b(highlands_mask);
    albedo_v_highlands = albedo_v(highlands_mask);
    albedo_r_highlands = albedo_r(highlands_mask);

    % Construct albedo matrices
    albedo_maria = [albedo_u_maria albedo_b_maria albedo_v_maria albedo_r_maria];
    albedo_highlands = [albedo_u_highlands albedo_b_highlands albedo_v_highlands albedo_r_highlands];

    % Calculate mean of albedo
    mean_albedo_maria = [mean(albedo_u_maria), mean(albedo_b_maria), mean(albedo_v_maria), mean(albedo_r_maria)]; 
    mean_albedo_highlands = [mean(albedo_u_highlands), mean(albedo_b_highlands), mean(albedo_v_highlands), mean(albedo_r_highlands)]; 

%     % Construct wavelength matrices
%     size_albedo_maria = size(albedo_u_maria);
%     wavelengths_maria = [repmat([373.8], size_albedo_maria(1), 1);
%                    repmat([443.5], size_albedo_maria(1), 1);
%                    repmat([558.6], size_albedo_maria(1), 1);
%                    repmat([676.6], size_albedo_maria(1), 1)];
% 
%     size_albedo_highlands = size(albedo_u_highlands);
%     wavelengths_highlands = [repmat([373.8], size_albedo_highlands(1), 1);
%                    repmat([443.5], size_albedo_highlands(1), 1);
%                    repmat([558.6], size_albedo_highlands(1), 1);
%                    repmat([676.6], size_albedo_highlands(1), 1)];
    
    maria_weights = [area(maria_mask) area(maria_mask) area(maria_mask) area(maria_mask)];
    highlands_weights = [area(highlands_mask) area(highlands_mask) area(highlands_mask) area(highlands_mask)];
    weights = zeros(size(maria_weights));
    
    
    wavelengths_single = [373.8 443.5 558.6 676.6];
    wavelengths = zeros(size(wavelengths_single), 2);
    wavelengths(:,:,1) = wavelengths_single;
    wavelengths(:,:,2) = wavelengths_single;
    
    albedos = zeros(size(albedo_maria), 2);
    albedos(:,:,1) = albedo_maria;
    albedos(:,:,2) = albedo_highlands;

    figure('Position', [500 500 900 900])
    iosr.statistics.boxPlot(wavelengths, albedos, 'symbolMarker', 'none', 'lineColor', {'b', 'r'}, 'medianColor','k', 'weights', maria_weights); hold on;
    
    %boxplot(albedo_maria, wavelengths_maria, 'Symbol','', 'Colors', 'b'); hold on
    plot([1 2 3 4], mean_albedo_maria, '-ob')
    
    
    iosr.statistics.boxPlot([373.8 443.5 558.6 676.6], albedo_highlands, 'symbolMarker', 'none', 'lineColor', 'r','medianColor','k', 'weights', highlands_weights)
    %boxplot(albedo_highlands, wavelengths_highlands, 'Notch','on', 'Symbol','', 'Colors', 'r'); 
    plot([1 2 3 4], mean_albedo_highlands, '-or')

    xlabel('Wavelength (nm)'); ylabel('$A$ (\%)', 'Interpreter', 'latex')
    grid on
    set(gca, 'FontSize', 20)
% 
%     ylim([0 30])
% 
%     boxes = findobj(gca, 'Tag', 'Box');
%     legend(boxes([end 1]), 'maria', 'highlands', 'Location', 'northwest')
end
