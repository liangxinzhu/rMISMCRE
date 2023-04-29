% ---------------------------------------------------------------------------------------- %
% function rmismc_complexity(rmismc_l,Nt,M,K,sigma,data,z_data)
%
% random multi-index Monte Carlo complexity test routine
%
% input:  rmismc_l = function for level l estimator
%         Nt       = total number of samples
%         M        = number of realizations
%         K        = starting refinement level
%         sigma    = std deviation of the error of y - G
%         data     = observations (values of y) vector
%         z_data   = corresponding observation points
%
% .mat: rmisn           = results of mimc with sn
%       rmire           = results of mimc with re
%       rmicosts        = total costs of mimc
%       mse_rmisn_solsn = mse of self-normalised increments estimator with solution with sn
%       mse_rmire_solsn = mse of ratio estimator with solution with sn
%       mse_rmisn_solre = mse of self-normalised increments estimator with solution with re
%       mse_rmire_solre = mse of ratio estimator with solution with re
%       rL             = actural finest level of refinement
%       K       = starting refinement level
% ------------------------------------------------------------------------------------------- %
function rmi_complexity(rmismc_l,Nt,M,K,sigma,data,x_data)
% results
rmisn = zeros(M,length(Nt));
rmire = zeros(M,length(Nt));
rmicosts = zeros(1,length(Nt));

% loading solution
[filepart,~,~] = fileparts(pwd);
loadpath = fullfile(filepart, 'Results','solutions.mat');
load(loadpath,'solutions')
solution_sn = solutions(1);
solution_re = solutions(2);

% random samples from level distribution
betax = 4;
betay = 4;
gammax = 1;
gammay = 1;
N0 = 1;
px = 1-2^(-(betax+gammax)/2);
py = 1-2^(-(betay+gammay)/2);

for i = 1:length(Nt)
    N = Nt(i);
    l_s = [geornd(px,N/N0,1), geornd(py,N/N0,1)];
    Lx = max(l_s(:,1));
    Ly = max(l_s(:,2));
    N_l = zeros(Lx+1,Ly+1); 
    p_l = zeros(Lx+1,Ly+1);
    for l1 = 0:Lx
        index_x = find(l_s(:,1) == l1);
        for l2 = 0:Ly
            index_y = find(l_s(:,2) == l2);
            N_l(l1+1,l2+1) = length(intersect(index_x, index_y)) * N0;
            p_l(l1+1, l2+1) = px*py*2^(-l1*(betax+gammax)/2 - l2*(betay+gammay)/2);
        end
    end
    costl = 0;
    parfor j = 1:M    
        sumlsn  = 0;
        sumlret = 0;
        sumlreb = 0;
        for lx = 0:Lx
            for ly = 0:Ly
                Nl = N_l(lx+1,ly+1);
                pl = p_l(lx+1,ly+1);
                if Nl ~= 0
                    [sums, cost] = rmismc_l(lx+K,ly+K,Nl,sigma,data,x_data,K);
                    sumlsn  = sumlsn  + sums(1);
                    sumlret = sumlret + sums(2)/(N*pl);
                    sumlreb = sumlreb + sums(3)/(N*pl);
                    costl   = costl + Nl*cost;
                end
            end
        end
        rmisn(j,i) = sumlsn;
        rmire(j,i) = sumlret/sumlreb;   
    end
    rmicosts(i) = costl/M;
end

mse_rmisn_solsn = sum((rmisn - solution_sn).^2)./M;
mse_rmire_solsn = sum((rmire - solution_sn).^2)./M;
mse_rmisn_solre = sum((rmisn - solution_re).^2)./M;
mse_rmire_solre = sum((rmire - solution_re).^2)./M;

rL = [Lx; Ly];

savepath = fullfile(filepart, 'Results','complexity results','rmismc_complexity.mat');
save(savepath,'rmisn','rmire','mse_rmisn_solsn','mse_rmire_solsn','mse_rmisn_solre','mse_rmire_solre','rmicosts','rL','K')  
end

