% ---------------------------------------------------- %
% generate the data set
%
% .mat: sigma  = std deviation of the error of y - G
%       data   = observations (values of y) vector
%       z_data = corresponding observation points
% ----------------------------------------------------- %
close all; clear all;

randn('seed',20)
rand('seed',20)

% generate data
z_data = linspace(0.1,1,10)';
sigma  = 0.2;
data   = -0.5*rand(1)*(z_data.^2-z_data) + sigma*randn(length(z_data),1);

[filepart,~,~] = fileparts(pwd);
savepath = fullfile(filepart,'1D_Toy','Results','observations.mat');
save(savepath,'sigma','data','z_data')
