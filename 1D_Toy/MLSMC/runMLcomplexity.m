% -------------------------------------------------------------------------- %
% run complexity test of MLSMC with the self-normalised increments estimator 
% and ratio estimator
% -------------------------------------------------------------------------- %
close all; clear all;

[filepart,~,~] = fileparts(pwd); 
addpath(fullfile(filepart, 'solver_pde_1D'))
addpath(fullfile(filepart, 'utils'))

randn('seed',20)
rand('seed',20)

Eps = [ 0.0025 0.005 0.01 0.02 0.04 ];
M   = 100;

% load data
loadpath = fullfile(filepart,'Results','observations.mat');
load(loadpath,'sigma','data','z_data')

tic
ml_complexity(@opre_mlsmc_l,Eps,M,sigma,data,z_data);
toc


