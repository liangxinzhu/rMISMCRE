% ----------------------------------------------------------------- %
% run empirical test of MLSMC for normalised/unnormalized estimator 
% ----------------------------------------------------------------- %
close all; clear all;

[filepart,~,~] = fileparts(pwd); 
addpath(fullfile(filepart, 'solver_pde_1D'))
addpath(fullfile(filepart, 'utils'))

randn('seed',20)
rand('seed',20)

N    = 1000;
Lmax = 6;

% load data
loadpath = fullfile(filepart,'Results','observations.mat');
load(loadpath,'sigma','data','z_data')

tic
ml_empirical(@opre_mlsmc_l, N, Lmax, sigma, data, z_data);
toc
