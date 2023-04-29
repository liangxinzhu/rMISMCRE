% -------------------------------------------------------------------------- %
% run complexity test of MISMC with the self-normalised increments estimator 
% and ratio estimator
% -------------------------------------------------------------------------- %
close all; clear all;

randn('seed',20)
rand('seed',20)

[filepart,~,~] = fileparts(pwd);
addpath(fullfile(filepart,'utils'));
addpath(fullfile(filepart,'solution and data'));

s2=1.91; %\sigma^2
b=(33/pi/2)^2; %\theta_3
a = 1; %\theta_2
mu = 0; %\theta_1
rate = 2.6; %r+1

M    = 200; % number of realisations
Eps  = 2.^(-12:1:-8); % required accuracy
K = 5; % starting level

% observations
load FinPine
xq=Like.x;yq=Like.y;
data = [xq yq]; 

%% MISMC TD
tic;
mi_complexity(@opre_mismc_l,Eps,M,data,K,rate,a,b,mu,1);
toc;
%% MISMC TP
tic;
mi_complexity(@opre_mismc_l,Eps,M,data,K,rate,a,b,mu,2);
toc;