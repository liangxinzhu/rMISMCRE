% -------------------------------------------------------------- %
% function sl_complexity(opre_sl,Eps,M,sigma,data,z_data)
%
% single level Monte Carlo complexity test routine
%
% Input: opre_sl  = function for level l estimator
%        Eps      = desired accuracy array for MLMC
%        M        = number of realizations for convergence tests
%        sigma    = std deviation of the error of y - G
%        data     = observations (values of y) vector
%        z_data   = corresponding observation points
%
% .mat: sl      = results of single level estimator
%       slcosts = total costs of sl
%       mse_sl  = mse of sl estimator
% --------------------------------------------------------------- %
function sl_complexity(opre_sl,Eps,M,sigma,data,z_data)

sl = zeros(M,length(Eps));
slcosts = zeros(1,length(Eps));

% weak convergence rate
s = 2;

[filepart,~,~] = fileparts(pwd);
loadpath = fullfile(filepart, 'Results','empirical results','mlsmc_empirical.mat');
load(loadpath,'se')
se1 = se{1};

% analytic solution of the quantity of interest (x^2)
x2 = analytic_x2(sigma,data,z_data);

for i = 1:length(Eps)
    eps = Eps(i);
    
    costl = 0; 
    
    % calculate MSE
    deno = (2^s - 1)*eps;
    L = max(1,ceil( log2(sqrt(2)/deno)/s )); 
    N = ceil(2*se1(1)/eps/eps);
    for j = 1:M
        [sums, costs] = opre_sl(L,N,sigma,data,z_data);
        sl(j,i) = sums;
        costl   = costl + costs*N;
    end
    slcosts(i) = costl/M;
end
mse_sl = sum((sl-x2).^2)./M;

savepath = fullfile(filepart,'Results','complexity results','sl_complexity.mat');
save(savepath,'mse_sl','slcosts')

end