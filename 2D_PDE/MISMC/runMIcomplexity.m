% -------------------------------------------------------------------------- %
% run complexity test of MISMC with the self-normalised increments estimator 
% and ratio estimator
% -------------------------------------------------------------------------- %
close all; clear all;

randn('seed',20)
rand('seed',20)

[filepart,~,~] = fileparts(pwd);
addpath(fullfile(filepart,'solver_pde_2D'));
addpath(fullfile(filepart,'utils'));
 
Eps = [ 0.001 0.0025 0.005 0.01 0.025 ];
M = 200;
K = 2;

% load data
loadpath = fullfile(filepart, 'Results','observations.mat');
load(loadpath,'sigma','data','z_data','x')

% MISMC with TP index set
tic;
mi_complexity(@opre_mismc_l,Eps,M,sigma,data,z_data,K,1);
toc;

% MISMC with TD index set
tic;
mi_complexity(@opre_mismc_l,Eps,M,sigma,data,z_data,K,2);
toc;

%rmpath(fullfile(filepart,'pde_solver_2D'));