% --------------------------------------------------------------- %
% approximate reference solution by MLSMC with high acccuracy
%
% solutions.mat:  
%        solutions(1) = solution with self-normalized estimator
%        solutions(2) = solution with ratio estimator
% --------------------------------------------------------------- %
close all; clear all;

randn('seed',20)
rand('seed',20)

[filepart,~,~] = fileparts(pwd);
addpath(fullfile(filepart,'utils'));

s2=1.91;
b=(33/pi/2)^2;
a = 1;
mu = 0;
rate = 2.6;

eps = 2^-9;
K = 4;

% observations
load FinPine
xq=Like.x;yq=Like.y;
data = [xq yq]; 

tic;
[solutions] = mismc_solution(@opre_mlsmc_l,eps,data,K,rate,a,b,mu);
toc;

savepath = fullfile(filepart,'Results','solutions.mat');
save(savepath,'solutions')

%% solution generator
% --------------------------------------------------------------------- %
% function [solutions] = mismc_solution(mlsmc_l,eps,data,K,rate,a,b,mu)
%
% inputs: mlsmc_l = mlsmc single level l estimator
%         eps     = required accuracy
%         data    = observations (values of y) vector
%         K       = starting refinement level
%         rate    = r+1
%         a       = \theta_2
%         b       = \theta_3
%         mu      = \theta_1
%
% ouputs: solutions(1) = reference solution with sn
%         solutions(2) = reference solution with re 
% --------------------------------------------------------------------- %
function [solutions] = mismc_solution(mlsmc_l,eps,data,K,rate,a,b,mu)
s = 0.8;
beta  = 1.6;
L = max(ceil(log2(1/eps/(2^s-1))/s),K+1);
L_set = (K:1:L);
var   = 2.^(-beta*L_set);
cost  = 2.^(2*L_set);
K_l = sum(sum( sqrt(var.*cost) ));
M = ceil(  eps^(-2)*K_l.*sqrt(var./cost) );

sumlsn  = 0;
sumlret = 0;
sumlreb = 0;
costl   = 0;

for l = K:L
    N = M(l-K+1);
    [sums, cst] = mlsmc_l(l,l,N,data,rate,a,b,mu,K);
    sumlsn  = sumlsn  + sums(1);
    sumlret = sumlret + sums(2);
    sumlreb = sumlreb + sums(3);
    costl   = costl + N*cst;
end

solutions(1) = sumlsn;
solutions(2) = sumlret/sumlreb;

end