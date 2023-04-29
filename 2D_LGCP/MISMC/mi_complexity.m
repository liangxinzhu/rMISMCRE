% -------------------------------------------------------------------------------------------- %
% function mi_complexity(mismc_l,Eps,M,data,K,rate,a,b,mu,setType)
%
% multi-index Monte Carlo complexity test routine
%
% inputs: mismc_l = function for level l estimator
%         Eps     = vector of required accuracy
%         M       = number of realisations
%         data    = observation points
%         K       = starting level
%         rate    = r+1
%         a       = \theta_2
%         b       = \theta_3
%         mu      = \theta_1
%         setType = 1: TD;2: TP
%
% .mat: misn           = results of mimc with sn
%       mire           = results of mimc with re
%       micosts        = total costs of mimc
%       K              = starting refinement level
%       mse_misn_solsn = mse of self-normalised increments estimator with solution with sn
%       mse_mire_solsn = mse of ratio estimator with solution with sn
%       mse_misn_solre = mse of self-normalised increments estimator with solution with re
%       mse_mire_solre = mse of ratio estimator with solution with re
% ------------------------------------------------------------------------------------------ %

function mi_complexity(mismc_l,Eps,M,data,K,rate,a,b,mu,setType)
% results
misn = zeros(M,length(Eps));
mire = zeros(M,length(Eps));
micosts = zeros(1,length(Eps));

[filepart,~,~] = fileparts(pwd);
loadpath = fullfile(filepart,'Results','solutions.mat');
load(loadpath,'solutions')
solution_sn = solutions(1);
solution_re = solutions(2);

% weak and strong convergence rate
s = 0.8;
beta  = 1.6;

% total degree index set
switch setType
    case 1
        for i = 1:length(Eps)
             eps = Eps(i);
             temp = 0.5/log(2)/s;
             logEps = -log(eps);
             L = ceil(temp*(logEps + log(temp*logEps) + ...
                4*exp(log(2)*s*sqrt(2))*(temp^2+temp) )) -6 ;
%            L = max(ceil(log((log(1/eps))^2/eps)/log(2)/2/alpha)+4, 2*k+1);
            for al = K:L
                for bl = K:L-al
                    L_ma(al-K+1,bl-K+1) = al + bl;
                end
            end
            
            var   = 2.^(-beta*L_ma);
            var(1,1) = 1e-3;
            cost  = L_ma.*2.^L_ma;
            %cost  = 2.^L_ma;
            K_l = sum(sum( sqrt(var.*cost) ));
            N_l = max(2*ceil( eps^(-2)*K_l.*sqrt(var./cost) ),2);
            N_l = fliplr(triu(fliplr(N_l)));
            costl = 0;
            parfor j = 1:M
                sumlsn  = 0;
                sumlret = 0;
                sumlreb = 0;
                for lx = K:L
                    for ly = K:L-lx
                        Nl = N_l(lx-K+1,ly-K+1);
                        [sums, cst] = mismc_l(lx,ly,Nl,data,rate,a,b,mu,K);
                        sumlsn  = sumlsn  + sums(1);
                        sumlret = sumlret + sums(2);
                        sumlreb = sumlreb + sums(3);
                        costl   = costl + Nl*cst;
                    end
                end
                misn(j,i) = sumlsn;
                mire(j,i) = sumlret/sumlreb;
            end
            micosts(i) = costl/M;
        end

        mse_misn_solsn = sum((misn - solution_sn).^2)./M;
        mse_mire_solsn = sum((mire - solution_sn).^2)./M;
        mse_misn_solre = sum((misn - solution_re).^2)./M;
        mse_mire_solre = sum((mire - solution_re).^2)./M;

        savepath = fullfile(filepart,'Results','complexity results','mismc_complexity_td.mat');
        save(savepath,'misn','mire','micosts','K','mse_misn_solsn','mse_mire_solsn','mse_misn_solre','mse_mire_solre')

% tensor product index set
    case 2
        for i = 1:length(Eps)
            eps = Eps(i);

            L = ceil( log2(2/eps)/s ) + 2;
            for al = K:L
                for bl = K:L
                    L_ma(al-K+1,bl-K+1) = al + bl;
                end
            end
            
            var   = 2.^(-beta*L_ma);
            var(1,1) = 1e-3;
            cost  = 2.^L_ma;
            
            K_l = sum(sum( sqrt(var.*cost) ));
            N_l = max(ceil(  eps^(-2)*K_l.*sqrt(var./cost) ),2);
            costl = 0;
            parfor j = 1:M
                sumlsn  = 0;
                sumlret = 0;
                sumlreb = 0;
                % total degree index set
                for lx = K:L
                    for ly = K:L
                        Nl = N_l(lx-K+1,ly-K+1);
                        [sums, cst] = mismc_l(lx,ly,Nl,data,rate,a,b,mu,K);
                        sumlsn  = sumlsn  + sums(1);
                        sumlret = sumlret + sums(2);
                        sumlreb = sumlreb + sums(3);
                        costl   = costl + Nl*cst;
                    end
                end
                misn(j,i) = sumlsn;
                mire(j,i) = sumlret/sumlreb;
            end
            micosts(i) = costl/M;
        end

        mse_misn_solsn = sum((misn - solution_sn).^2)./M;
        mse_mire_solsn = sum((mire - solution_sn).^2)./M;
        mse_misn_solre = sum((misn - solution_re).^2)./M;
        mse_mire_solre = sum((mire - solution_re).^2)./M;

        savepath = fullfile(filepart,'Results','complexity results','mismc_complexity_tp.mat');
        save(savepath,'misn','mire','micosts','K','mse_misn_solsn','mse_mire_solsn','mse_misn_solre','mse_mire_solre')
end

end

