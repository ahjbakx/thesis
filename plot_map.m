function plot_map(quantity, label, latitude, longitude, save, savename, shape)
    
    
    figure('Position', [500 500 900 900])
    axesm('vperspec', 'Grid','on','Frame','on', 'PLabelLocation',30, 'GLineWidth', 1, 'Origin',[0 0])
    geoshow(latitude, longitude, quantity, 'DisplayType','texturemap')

    % Gray colormap
    colormap gray

    % Position colorbar below Moon
    cbar = colorbar('southoutside');

    % Colorbar labels: low, mid, high value, with log(%)
    mid = (cbar.Limits(1) + cbar.Limits(2))/2;
    %set(cbar, 'Ticks', round([cbar.Limits(1), mid, cbar.Limits(2)], 2))
    set(cbar, 'Ticks', [0.86 1.15 1.44])
    %caxis(round(cbar.Limits, 2))
    caxis([0.85 1.44])
    xlabel(cbar, strcat('log ', label,' (\%)'),'Interpreter','latex' )

    % Adjust width of colorbar
    x1=get(gca,'position');
    x=get(cbar,'Position');
    x(1) = x(1) + 0.1250;
    x(3) = 0.5;
    set(cbar,'Position',x)
    set(gca,'position',x1)
    
    axis off
    set(gca, 'FontSize', 20)
    
    if shape == true
        geoshow('/Users/aaron/thesis/Data/mare_shape/LROC_GLOBAL_MARE_180.shp', 'DisplayType','polygon','FaceColor','none','EdgeColor','w', 'LineWidth',0.5)
    end
    
    if save == true
        exportgraphics(gca,strcat('FiguresWP1/', savename ,'.png'),'Resolution', 300)
    end

end