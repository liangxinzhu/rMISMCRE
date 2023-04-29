% ------------------------------------------------------------------ %
% function mi_empirical(mismc_l,N,M,Lmax,data,k,rate,a,b,mu)
%
% mult-index Monte Carlo empirical test routine
%
% input: mismc_l = function for level l estimator
%        N       = number of samples for convergence tests
%        M       = number of realisations
%        Lmax    = max number of levels for convergence tests in x and y
%        data    = observation points
%        k       = starting level
%        rate    = r+1
%        a       = \theta_2
%        b       = \theta_3
%        mu      = \theta_1
%
% .mat: Lmax = maximum level we set
%       cost = cost at each level
%       K       = starting refinement level
%       we1  = difference with self-normalised increments estimator
%       we2  = difference unnormalised integral
%       we3  = difference normalising constant
%       se1  = variance of sn increments estimtor
%       se2  = unnormalised integral of squared target
%       se3  = unnormalised integral of 1
%       se4  = square of difference with self-normalised increments estimator
%       se5  = square of difference unnormalised integral
%       se6  = square of difference normalising constant
% ------------------------------------------------------------------ %

function mi_empirical(mismc_l,N,M,Lmax,data,K,rate,a,b,mu)
we1  = zeros(Lmax-K+1,Lmax-K+1);
we2  = zeros(Lmax-K+1,Lmax-K+1);
we3  = zeros(Lmax-K+1,Lmax-K+1);
se1  = zeros(Lmax-K+1,Lmax-K+1);
se2  = zeros(Lmax-K+1,Lmax-K+1);
se3  = zeros(Lmax-K+1,Lmax-K+1);
se4  = zeros(Lmax-K+1,Lmax-K+1);
se5  = zeros(Lmax-K+1,Lmax-K+1);
se6  = zeros(Lmax-K+1,Lmax-K+1);
cost = zeros(Lmax-K+1,Lmax-K+1);

for lx = K:Lmax
    for ly = K:Lmax

    sums = [];
    cst  = 0;
    
    fprintf(1,'lx = %.1d, ly = %.1d \n',lx,ly)
    
    parfor j = 1:M
        fprintf(1,'j = %.1d \n',j)
        [sums_j, cst_j] = mismc_l(lx,ly,N,data,rate,a,b,mu,K);
        sums(j,:) = sums_j;
        cst  = cst  + cst_j/M;      
    end
    
    idx = lx - K + 1;
    idy = ly - K + 1;
    
    cost(idx,idy)  = cst;
    we1(idx,idy)  = abs( sum(sums(:,1))/M );
    we2(idx,idy)  = abs( sum(sums(:,2))/M );
    we3(idx,idy)  = abs( sum(sums(:,3))/M );
    
    se1(idx,idy)  = N*sum( (sums(:,1) - sum(sums(:,1))/M).^2 )/M;
    se2(idx,idy)  = N*sum( (sums(:,2) - sum(sums(:,2))/M).^2 )/M;
    se3(idx,idy)  = N*sum( (sums(:,3) - sum(sums(:,3))/M).^2 )/M;
    
    se4(idx,idy)  = sum(sums(:,1).^2)/M;
    se5(idx,idy)  = sum(sums(:,2).^2)/M;
    se6(idx,idy)  = sum(sums(:,3).^2)/M;
    end 
end
we = {we1,we2,we3};
se = {se1,se2,se3,se4,se5,se6};

[filepart,~,~] = fileparts(pwd);
savepath = fullfile(filepart,'Results','empirical results','mismc_empirical.mat');
save(savepath,'Lmax','cost','K','we','se')  

end

