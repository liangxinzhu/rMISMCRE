% ------------------------------------------------------------------ %
% function sl_complexity(opre_sl,Eps,M,sigma,data,z_data)
%
% single level Monte Carlo complexity test routine
%
% input: opre_sl = function for level l estimator
%        Eps     = vector of required accuracy
%        M       = number of realizations
%        sigma   = std deviation of the error of y - G
%        data    = observations (values of y) vector
%        z_data  = corresponding observation points
%
% .mat: sl           = results of single level estimator
%       slcosts      = total costs of sl
%       mse_sl_solsn = mse of sl estimator with solution with sn
%       mse_sl_solre = mse of sl estimator with solution with re
% ------------------------------------------------------------------ %
function sl_complexity(opre_sl,Eps,M,sigma,data,z_data)
% results
sl = zeros(M,length(Eps));
slcosts = zeros(1,length(Eps));

% loading emperical results first
[filepart,~,~] = fileparts(pwd);
loadpath = fullfile(filepart, 'Results','empirical results','mismc_empirical.mat');
load(loadpath,'se')
se1 = se{1};
loadpath = fullfile(filepart, 'Results','solutions.mat');
load(loadpath,'solutions')
solution_sn = solutions(1);
solution_re = solutions(2);

for i = 1:length(Eps)
    eps = Eps(i);
    L = ceil(log2(1/eps/3)/2);
    N = ceil(se1(1)/eps/eps);
    costl = 0;
    parfor j = 1:M    
        [sums, cst] = opre_sl(L,N,sigma,data,z_data);
        sl(j,i) = sums;
        costl = costl + cst;
    end
    slcosts(i) = costl/M;
end

mse_sl_solsn = sum((sl - solution_sn).^2)./M;
mse_sl_solre = sum((sl - solution_re).^2)./M;

savepath = fullfile(filepart, 'Results','complexity results','sl_complexity.mat');
save(savepath,'sl','slcosts','mse_sl_solsn','mse_sl_solre')  
end