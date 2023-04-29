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
b=(33/pi)^2; %\theta_3
a = 1; %\theta_2
mu = 0; %\theta_1
rate = 2.6; %r+1

M = 200; %number of realisations
K = 5; %starting level

% observations
load FinPine
xq=Like.x;yq=Like.y;
data = [xq yq];
%% MISMC TD
Eps  = 2.^(-10:1:-7); %required accuracy 
tic;
mi_complexity(@opre_mismc_l,Eps,M,data,K,rate,a,b,mu,1);
toc;
%% MISMC TP
Eps  = 2.^(-6:1:-3); %required accuracy 
tic;
mi_complexity(@opre_mismc_l,Eps,M,data,K,rate,a,b,mu,2);
toc;