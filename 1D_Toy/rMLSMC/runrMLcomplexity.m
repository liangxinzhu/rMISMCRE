% --------------------------------------------------------------------------- %
% run complexity test of rMLSMC with the self-normalised increments estimator 
% and ratio estimator
% --------------------------------------------------------------------------- %
close all; clear all;

[filepart,~,~] = fileparts(pwd); 
addpath(fullfile(filepart, 'solver_pde_1D'))
addpath(fullfile(filepart, 'utils'))

randn('seed',20)
rand('seed',20)

Nt  = [10^2 4*10^2 16*10^2 64*10^2 256*10^2];
M  = 100; 
N0 = 2^3;

% load data
loadpath = fullfile(filepart,'Results','observations.mat');
load(loadpath,'sigma','data','z_data')

tic
rml_complexity(@opre_rmlsmc_l,M,Nt,N0,sigma,data,z_data);
toc


