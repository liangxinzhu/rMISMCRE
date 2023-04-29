% ----------------------------------------------------------------- %
% run empirical test of MISMC for normalised/unnormalized estimator 
% ----------------------------------------------------------------- %
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

N = 1000; %number of samples
M = 20; %number of realisations
Lmax = 8; %maximum level for empirical tests
K = 5; %starting level

% observations
load FinPine
xq=Like.x;yq=Like.y;
data = [xq yq]; 

tic;
mi_empirical(@opre_mismc_l,N,M,Lmax,data,K,rate,a,b,mu);
toc;
