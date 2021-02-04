% makefile for the complete GSH circle for a particular model
clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

% Model
% Load previous saved model

model_name = 'moon_topography_layer';
load(model_name);

% % Construct new model
% inputModel

%%%%%%%%%%%%%%%%%%% Computation area %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%

latLim =    [-89.5 89.5 0.50];  % [deg] min latitude, max latitude, resolution latitude (preferable similar to latitude)
lonLim =    [-180 180 0.50];% [deg] min longitude, max longitude, resolution longitude (preferable similar to latitude)
height =    0.0; % height of computation above spheroid
SHbounds =  [7 100]; % Truncation settings: lower limit, upper limit SH-coefficients used

%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%

%% Global Spherical Harmonic Analysis 

% tic;
% [V] = model_SH_analysis(Model);
% toc

% % Moon SH Coefficients
% file = '/Users/aaron/thesis/Data/moon_gravity/SH_coefficients_GSHformat_1200B-RM1.mat';
% SHcoefficients = load(file);
% [V] = SHcoefficients.coefficients_new;
% clear data

% Bouguer Coefficients
file = '/Users/aaron/thesis/GSH/Results/Bouguer_1ppd.mat';
SHcoefficients = load(file);
[V] = SHcoefficients.V;
clear data


% save(['Results/' Model.name '.mat'],'V')

%% Global Spherical Harmonic Synthesis

tic;
[data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V,Model);
toc

%% Save data

DATE = datestr(now);
save(['Results/data_' Model.name '_' num2str(SHbounds(1)) '_' num2str(SHbounds(2)) '_' DATE '.mat'],'data')