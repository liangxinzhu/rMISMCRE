% ---------------------------------------------------------------- %
% function rmi_complexity(rmismc_l,Nt,rel,data,K,rate,a,b,mu,randset)
%
% randomized multi-index Monte Carlo complexity test routine
%
% inputs: rmismc_l = mismc single level routine
%         Eps      = required accuracy
%         M        = number of realisations
%         K        = starting level
%         rate     = \beta^{\prime}
%         a        = \theta_2
%         b        = \theta_3
%         mu       = \theta_1
%         randset  = a random set
%
% .mat: rmisn           = results of mimc with sn
%       rmire           = results of mimc with re
%       rmicosts        = total costs of mimc
%       K               = starting refinement level
%       mse_rmisn_solsn = mse of self-normalised increments estimator with solution with sn
%       mse_rmire_solsn = mse of ratio estimator with solution with sn
%       mse_rmisn_solre = mse of self-normalised increments estimator with solution with re
%       mse_rmire_solre = mse of ratio estimator with solution with re
% ---------------------------------------------------------------- %

function rmi_complexity(rmismc_l,Nt,M,data,K,rate,a,b,mu,randset)
% results
rmisn = zeros(M,length(Nt));
rmire = zeros(M,length(Nt));
rmicosts = zeros(1,length(Nt));

[filepart,~,~] = fileparts(pwd);
loadpath = fullfile(filepart,'Results','solutions.mat');
load(loadpath,'solutions')
solution_sn = solutions(1);
solution_re = solutions(2);

beta  = 1.6;
gamma = 1;

% random samples from level distribution
N0 = 2^4;%2^2,2^6
L_max = 30;
pp = 2^(-(beta+gamma)/2);
px = [10 pp.^(0:L_max-1)];
py = px;
C_p = sum(px);


for i = 1:length(Nt)
    N = Nt(i);
    l_s1 = randsample(randset,0:L_max,N/N0,true,px)';
    l_s2 = randsample(randset,0:L_max,N/N0,true,py)';
    l_s = [l_s1, l_s2];
    Lx = max(l_s(:,1));
    Ly = max(l_s(:,2));
    Nlt = zeros(Lx+1,Ly+1); 
    p = zeros(Lx+1,Ly+1);
    for l1 = 0:Lx
        index_x = find(l_s(:,1) == l1);
        for l2 = 0:Ly
            index_y = find(l_s(:,2) == l2);
            Nlt(l1+1,l2+1) = length(intersect(index_x, index_y)) * N0;
            p(l1+1, l2+1) = 1/C_p/C_p*px(l1+1)*py(l2+1);
        end
    end
    costl = 0;
    parfor j = 1:M
        sumlsn  = 0;
        sumlret = 0;
        sumlreb = 0;
        for lx = 0:Lx
            for ly = 0:Ly
                Nl = Nlt(lx+1,ly+1);
                pl = p(lx+1,ly+1)
                if Nl ~= 0
                [sums, cst] = rmismc_l(lx+K,ly+K,Nl,data,rate,a,b,mu,K);
                sumlsn  = sumlsn  + sums(1);
                sumlret = sumlret + sums(2)/(N*pl);
                sumlreb = sumlreb + sums(3)/(N*pl);
                costl   = costl + Nl*cst;
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

savepath = fullfile(filepart,'Results','complexity results','rmismc_complexity.mat');
save(savepath,'rmisn','rmire','rmicosts','K','mse_rmisn_solsn','mse_rmire_solsn','mse_rmisn_solre','mse_rmire_solre')
        
end

