% ----------------------------------------------------------------- %
% run empirical test of MLSMC for normalised/unnormalized estimator 
% ----------------------------------------------------------------- %
close all; clear all;

randn('seed',20)
rand('seed',20)

[filepart,~,~] = fileparts(pwd);
addpath(fullfile(filepart,'solver_pde_2D'));
addpath(fullfile(filepart,'utils'));

N = 1000;
M = 20;
K = 2;
L = 4;

% load data
loadpath = fullfile(filepart, 'Results','observations.mat');
load(loadpath,'sigma','data','z_data','x')

tic;
ml_empirical(@opre_mlsmc_l,N,M,L,sigma,data,z_data,K);
toc;