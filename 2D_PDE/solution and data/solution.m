% --------------------------------------------------------------- %
% approximate reference solution by MLSMC with high acccuracy
%
% solutions.mat:  
%        solutions(1) = solution with self-normalized estimator
%        solutions(2) = solution with ratio estimator
% --------------------------------------------------------------- %
[filepart,~,~] = fileparts(pwd);
addpath(fullfile(filepart,'solver_pde_2D'));
addpath(fullfile(filepart,'utils'));

eps = 1e-4;
K = 2;
tic;
[solutions] = mismc_solution(@opre_mlsmc_l,eps,sigma,data,z_data,K);
toc;
savepath = fullfile(filepart, 'Results','solutions.mat');
save(savepath,'solutions');
%rmpath(fullfile(filepart,'MLSMC'))
%rmpath(fullfile(filepart,'pde_solver_2D'));

%% solution generator
% --------------------------------------------------------------------- %
% function [solutions] = mismc_solution(mlsmc_l,eps,sigma,data,x_data,k)
%
% inputs: mlsmc_l = function for level l estimator
%         eps     = required accuracy
%         sigma   = std deviation of the error of y - G
%         data    = observations (values of y) vector
%         z_data  = corresponding observation points
%         K       = starting refinement level
%
% ouputs: solutions(1) = reference solution with sn
%         solutions(2) = reference solution with re 
% --------------------------------------------------------------------- %
function [solutions] = mismc_solution(mlsmc_l,eps,sigma,data,z_data,K)
% loading empirical results first
[filepart,~,~] = fileparts(pwd);
loadpath = fullfile(filepart, 'Results','empirical results','mlsmc_empirical.mat');
load(loadpath,'se','cost')
se1 = se{1};

solutions = zeros(1,2);

vc1 = 2*sum(sum(sqrt(se1.*cost))).*sqrt(se1./cost);

L = ceil(log2(sqrt(2)/eps/3)/2)-K;
M = ceil(vc1./eps/eps);

sumlsn  = 0;
sumlret = 0;
sumlreb = 0;

for l = 0:L
    N = M(l+1);
    [sums,~] = mlsmc_l(l,l,N,sigma,data,z_data,K);
    sumlsn  = sumlsn  + sums(1);
    sumlret = sumlret + sums(2);
    sumlreb = sumlreb + sums(3);
end
solutions(1) = sumlsn;
solutions(2) = sumlret/sumlreb;

end