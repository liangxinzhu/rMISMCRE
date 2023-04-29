% -------------------------------------------------------------- %
% function ml_complexity(mlsmc_l,Eps,M,sigma,data,z_data)
%
% multilevel Monte Carlo complexity test routine
%
% Input: mlmc_l = function for level l estimator
%        Eps    = desired accuracy array for MLMC
%        M      = number of realizations for convergence tests
%        sigma  = std deviation of the error of y - G
%        data   = observations (values of y) vector
%        z_data = corresponding observation points
%
% .mat: mlsn     = results of mlmc with self-noramlized estimator
%       mlre     = results of mlmc with ratio estimator
%       mlcosts  = total costs of mlmc
%       mse_mlsn = mse of self-normalised increments estimator
%       mse_mlre = mse of ratio estimator 
% --------------------------------------------------------------- %
function ml_complexity(mlsmc_l,Eps,M,sigma,data,z_data)

mlsn = zeros(M,length(Eps));
mlre = zeros(M,length(Eps));
mlcosts = zeros(1,length(Eps));

% weak convergence rate
s = 2;

% load empirical results
[filepart,~,~] = fileparts(pwd);
loadpath = fullfile(filepart, 'Results','empirical results','mlsmc_empirical.mat');
load(loadpath,'se','cost')
se1 = se{1};

% analytic solution of the quantity of interest (x^2)
x2 = analytic_x2(sigma,data,z_data);

for i = 1:length(Eps)
    eps = Eps(i);
    deno = (2^s - 1)*eps;
    L = max(1,ceil( log2(sqrt(2)/deno)/s ));    
    N_l = ceil(2*eps^(-2)*sum(sqrt(se1(1:L).*cost(1:L))).* ...
        sqrt(se1(1:L)./cost(1:L)) );
    costl = 0;
    parfor j = 1:M
        sumlsn = 0;
        sumlret = 0;
        sumlreb = 0;
        for l=1:L
             Nl = N_l(l);
            [sums, cost] = mlsmc_l(l,Nl,sigma,data,z_data);
            sumlsn  = sumlsn + sums(1);
            sumlret = sumlret + sums(2);
            sumlreb = sumlreb + sums(3);
            costl   = costl + cost*Nl;
        end
        mlsn(j,i) = sumlsn;
        mlre(j,i) = sumlret/sumlreb;
    end
    mlcosts(i) = costl/M;
end
mse_mlsn = sum((mlsn-x2).^2)./M;
mse_mlre = sum((mlre-x2).^2)./M;

savepath = fullfile(filepart,'Results','complexity results','mlsmc_complexity.mat');
save(savepath,'mlsn','mlre','mse_mlsn','mse_mlre','mlcosts')

end