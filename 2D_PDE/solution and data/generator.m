% ---------------------------------------------------- %
% generate the data
%
% observations.mat:  
%        sigma  = std deviation of the error of y - G
%        data   = observations (values of y) vector
%        z_data = corresponding observation points
%        x      = value of parameter
% ----------------------------------------------------- %
close all; clear all;

randn('seed',20)
rand('seed',20)

[filepart,~,~] = fileparts(pwd);
addpath(fullfile(filepart,'solver_pde_2D'));

sigma  = 0.5;
% observation points
z_data = [0.25, 0.25; 0.25, 0.75; 0.75, 0.25; 0.75, 0.75];
% generate data
x = 2*rand(2,1) - 1;
tic;
data =  data_generator(z_data,x,sigma);
toc;
savepath = fullfile(filepart, 'Results','observations.mat');
save(savepath,'sigma','data','z_data','x')

%% data generator
% ------------------------------------------------------ %
% function [data] = data_generator(z_data, x, sigma)
% generate data based on the model y = G(u) + v
% input:  z_data = observation points
%         x      = value of parameter
%         sigma  = std deviation of the error of y - G
%
% output: data   = observations (values of y) vector
% ------------------------------------------------------ %
function [data] = data_generator(z_data, x, sigma)
lx = 10;
ly = 10;

Gu = pde_solver_2D(lx,ly,z_data,x);

data = Gu + sigma*randn(length(z_data), 1);
end
