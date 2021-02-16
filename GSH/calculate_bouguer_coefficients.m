clear all; close all; clc;

file = '/Users/aaron/thesis/GSH/Results/moon_topography_layer.mat';
SHcoefficients = load(file);

Vtopo = SHcoefficients.V;

% Moon SH Coefficients
file = '/Users/aaron/thesis/Data/moon_gravity/SH_coefficients_GSHformat_1200B-RM1.mat';
SHcoefficients = load(file);
Vgrav = SHcoefficients.coefficients_new;

% Bouguer
Vgrav(Vgrav(:,1)>100|Vgrav(:,1)<1,:) = [];
Vtopo(Vtopo(:,1)>100|Vtopo(:,1)<1,:) = [];
V = zeros(size(Vtopo));

V(:,1:2) = Vtopo(:,1:2);
V(:,3:4) = Vgrav(:,3:4) - Vtopo(:,3:4);

% save(['GSH/Results/Bouguer_1ppd.mat'],'V')