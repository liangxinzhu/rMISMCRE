% ----------------------------------------------------------------- %
% run empirical test of MISMC for normalised/unnormalized estimator 
% ----------------------------------------------------------------- %
close all; clear all;

randn('seed',20)
rand('seed',20)

[filepart,~,~] = fileparts(pwd);
addpath(fullfile(filepart,'solver_pde_2D'));
addpath(fullfile(filepart,'utils'));

N      = 1000; 
Lx     = 6;
Ly     = 6;
K = 2;

% load data
loadpath = fullfile(filepart, 'Results','observations.mat');
load(loadpath,'sigma','data','z_data','x')

tic;
mi_empirical(@opre_mismc_l,N, Lx, Ly,sigma,data,z_data,K);
toc;
