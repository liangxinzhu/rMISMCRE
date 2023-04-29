% ----------------------------------------------------------------------------------------- %
% function mismc_complexity(mismc_l,Eps,M,sigma,data,x_data,K,setType)
%
% multi-index Monte Carlo complexity test routine
%
% input: mimc_l  = function for level l estimator
%        Eps     = vector of required accuracy
%        M       = number of realizations
%        sigma   = std deviation of the error of y - G
%        data    = observations (values of y) vector
%        z_data  = corresponding observation points
%        K       = starting refinement level
%        setType = type of index set
%
% .mat: misn           = results of mimc with sn
%       mire           = results of mimc with re
%       micosts        = total costs of mimc
%       mse_misn_solsn = mse of self-normalised increments estimator with solution with sn
%       mse_mire_solsn = mse of ratio estimator with solution with sn
%       mse_misn_solre = mse of self-normalised increments estimator with solution with re
%       mse_mire_solre = mse of ratio estimator with solution with re
%       rL             = actural finest level of refinement
%       K       = starting refinement level
% -------------------------------------------------------------------------------------------- %
function mi_complexity(mismc_l,Eps,M,sigma,data,z_data,K,setType)
% results
misn = zeros(M,length(Eps));
mire = zeros(M,length(Eps));
micosts = zeros(1,length(Eps));

% loading emperical results and solutions
[filepart,~,~] = fileparts(pwd);
loadpath = fullfile(filepart, 'Results','empirical results','mismc_empirical.mat');
load(loadpath,'se','cost')
se1 = se{1};
loadpath = fullfile(filepart, 'Results','solutions.mat');
load(loadpath,'solutions')

solution_sn = solutions(1);
solution_re = solutions(2);

% tensor product index set
if setType == 1
for i = 1:length(Eps)
    eps = Eps(i);
    L = ceil(log2(2/eps)/2)-K;
    vc1 = sum(sum(sqrt(se1(1:L+1,1:L+1).*cost(1:L+1,1:L+1)))).*...
        sqrt(se1(1:L+1,1:L+1)./cost(1:L+1,1:L+1));
    N_l = ceil(vc1./eps/eps);
    costl = 0;
    parfor j = 1:M    
        sumlsn  = 0;
        sumlret = 0;
        sumlreb = 0;
        for lx = 0:L
            for ly = 0:L
                Nl = N_l(lx+1,ly+1);
                [sums, cst] = mismc_l(lx,ly,Nl,sigma,data,z_data,K);
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

rL = L + K;
savepath = fullfile(filepart, 'Results','complexity results','mismc_complexity_tp.mat');
save(savepath,'misn','mire','micosts','mse_misn_solsn','mse_mire_solsn','mse_misn_solre','mse_mire_solre','rL','K')  

% total degree index set
elseif setType == 2
for i = 1:length(Eps)
    eps = Eps(i);
    L = ceil( 0.25/log(2)*(log((log(1/eps))^2/eps)) )-K;
    se = fliplr(triu(fliplr(se1(1:L+1,1:L+1)))); 
    vc1 = sum(sum(sqrt(se.*cost(1:L+1,1:L+1)))).*...
        sqrt(se(1:L+1,1:L+1)./cost(1:L+1,1:L+1));
    N_l = ceil( vc1./eps/eps );
    costl = 0;
    parfor j = 1:M    
        sumlsn  = 0;
        sumlret = 0;
        sumlreb = 0;
        for lx = 0:L
            for ly = 0:L-lx
                Nl = N_l(lx+1,ly+1);
                [sums, cst] = mismc_l(lx,ly,Nl,sigma,data,z_data,K);
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

rL = L + K;
savepath = fullfile(filepart, 'Results','complexity results','mismc_complexity_td.mat');
save(savepath,'misn','mire','micosts','mse_misn_solsn','mse_mire_solsn','mse_misn_solre','mse_mire_solre','rL','K')    

else
    fprintf(1,'Wrong setType')
end
end

