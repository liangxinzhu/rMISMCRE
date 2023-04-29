% -------------------------------------------------------------------------- %
% run complexity test of rMISMC with the self-normalised increments estimator 
% and ratio estimator
% -------------------------------------------------------------------------- %
close all; clear all;

randn('seed',20)
rand('seed',20)

[filepart,~,~] = fileparts(pwd);
addpath(fullfile(filepart,'solver_pde_2D'));
addpath(fullfile(filepart,'utils'));

Nt = [10^2 10^3 10^4 10^5];
M = 200;
K = 2;

% load data
loadpath = fullfile(filepart, 'Results','observations.mat');
load(loadpath,'sigma','data','z_data','x')

tic;
rmi_complexity(@opre_rmismc_l,Nt,M,K,sigma,data,z_data);
toc;

%rmpath(fullfile(filepart,'solver2Dnlin'));