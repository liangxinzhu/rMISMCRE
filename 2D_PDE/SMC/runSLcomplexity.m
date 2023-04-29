% --------------------------- %
% run complexity test of SMC 
% --------------------------- %
close all; clear all;

randn('seed',20)
rand('seed',20)

[filepart,~,~] = fileparts(pwd);
addpath(fullfile(filepart,'solver_pde_2D'));
addpath(fullfile(filepart,'utils'));
 
Eps = [ 0.001 0.0025 0.005 0.01 0.025 ];
M = 200;

% load data
loadpath = fullfile(filepart, 'Results','observations.mat');
load(loadpath,'sigma','data','z_data','x')

tic;
sl_complexity(@opre_sl,Eps,M,sigma,data,z_data);
toc;