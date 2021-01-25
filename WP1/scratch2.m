clear all; close all; clc;
% (requires Statistics or Machine Learning Toolbox)
load carbig

% arrange data
[y,x,g] = iosr.statistics.tab2box(Cylinders,MPG,when);

% sort
IX = [1 3 2]; % order
g = g{1}(IX);
y = y(:,:,IX);

% plot
figure
h = iosr.statistics.boxPlot(x,y, 'symbolColor','k','medianColor','k','symbolMarker','+', 'boxcolor',{[1 1 1],[.75 .75 .75],[.5 .5 .5]}, 'scalewidth',true,'xseparator',true, 'groupLabels',g,'showLegend',true);
box on
title('MPG by number of cylinders and period')
xlabel('Number of cylinders')
ylabel('MPG')