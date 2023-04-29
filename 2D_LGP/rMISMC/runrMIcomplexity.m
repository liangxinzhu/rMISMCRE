% -------------------------------------------------------------------------- %
% run complexity test of rMISMC with the self-normalised increments estimator 
% and ratio estimator
% -------------------------------------------------------------------------- %
close all; clear all;

randn('seed',20)
rand('seed',20)
randset = RandStream('mlfg6331_64'); 

[filepart,~,~] = fileparts(pwd);
addpath(fullfile(filepart,'utils'));
addpath(fullfile(filepart,'solution and data'));

s2=1.91; %\sigma^2
b=(33/pi/2)^2; %\theta_3
a = 1; %\theta_2
mu = 0; %\theta_1
rate = 2.6; %r+1

M    = 200; % number of realisations
Nt  = 2*[80 2^8 2^10 2^12]; % required total number of samples
k = 5; % starting level

% observations
load FinPine
xq=Like.x;yq=Like.y;
data = [xq yq]; 

tic;
rmi_complexity(@opre_rmismc_l,Nt,M,data,k,rate,a,b,mu,randset);
toc;
