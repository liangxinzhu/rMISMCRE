% -------------------------------------------------------------- %
% function rml_complexity(rmlsmc_l,M,N,N0,sigma,data,z_data)
%
% randomized multilevel Monte Carlo complexity test routine
%
% Input: rmlmc_l   = function for level l estimator
%        M         = number of realizations for convergence tests
%        N         = total number of sample array
%        N0        = minimum number of sample
%        sigma     = std deviation of the error of y - G
%        data      = observations (values of y) vector
%        z_data    = corresponding observation points
%
% .mat: rmlsn     = results of mlmc with self-noramlized estimator
%       rmlre     = results of mlmc with ratio estimator
%       rmlcosts  = total costs of mlmc
%       mse_rmlsn = mse of self-normalised increments estimator
%       mse_rmlre = mse of ratio estimator 
% --------------------------------------------------------------- %
function rml_complexity(rmlsmc_l,M,Nt,N0,sigma,data,z_data)

rmlsn = zeros(M,length(Nt));
rmlre = zeros(M,length(Nt));
rmlcosts = zeros(1,length(Nt));

beta  = 4;
gamma = 1;

% analytic solution of the quantity of interest (u^2)
x2 = analytic_x2(sigma,data,z_data);

for i = 1:length(Nt)
    N = Nt(i);
    
    % random samples from level distribution
    p = 1 - 2^(-(beta+gamma)/2);
    l_s = geornd(p,N/N0,1);

    % compute the number of samples for each level
    L = max(l_s);
    N_l = zeros(1,L+1);
    for l = 0:L
        N_l(l+1) = sum(l_s==l)*N0;
    end

    % compute the distribution of levels
    p_l(1,1:L+1) = (1-2^(-(beta+gamma)/2))*2.^((-(beta+gamma)/2).*(0:L));
    
    costl = 0;
    
    % calculate MSE
    parfor j = 1:M
        %[P,cl] = mlsmc(mlsmc_l,eps,L,N_l,sigma,data,x_data);
        sumlsn = 0;
        sumlret = 0;
        sumlreb = 0;
        for l=1:L+1
            if N_l(l) > 0
                Nl = N_l(l);
                pl = p_l(l);
                [sums, cost] = rmlsmc_l(l,Nl,sigma,data,z_data);
                sumlsn  = sumlsn + sums(1)*Nl/pl/N;
                sumlret = sumlret + sums(2)*Nl/pl/N;
                sumlreb = sumlreb + sums(3)*Nl/pl/N;
                costl   = costl + cost*Nl;
            end
        end
        rmlsn(j,i) = sumlsn;
        rmlre(j,i) = sumlret/sumlreb;
    end
    rmlcosts(i) = costl/M;
end
mse_rmlsn = sum((rmlsn-x2).^2)./M;
mse_rmlre = sum((rmlre-x2).^2)./M;

[filepart,~,~] = fileparts(pwd);
savepath = fullfile(filepart,'Results','complexity results','rmlsmc_complexity.mat');
save(savepath,'rmlre','rmlsn','rmlcosts','mse_rmlre','mse_rmlsn')

end