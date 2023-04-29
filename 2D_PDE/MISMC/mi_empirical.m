% ------------------------------------------------------------------ %
% function mismc_emperical(mimc_l,N,Lx,Ly,sigma,data,z_data,K)
%
% mult-index Monte Carlo empirical test routine
%
% input:  mimc_l = function for level l estimator
%         N      = number of samples for convergence tests
%         Lx     = number of levels for convergence tests in x
%         Ly     = number of levels for convergence tests in y
%         sigma  = std deviation of the error of y - G
%         data   = observations (values of y)
%         z_data = corresponding observation points
%         K       = starting refinement level
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
function mi_empirical(mismc_l,N, Lx, Ly,sigma,data,z_data,K)

we1  = zeros(Lx+1,Ly+1);
we2  = zeros(Lx+1,Ly+1);
we3  = zeros(Lx+1,Ly+1);
se1  = zeros(Lx+1,Ly+1);
se2  = zeros(Lx+1,Ly+1);
se3  = zeros(Lx+1,Ly+1);
se4  = zeros(Lx+1,Ly+1);
se5  = zeros(Lx+1,Ly+1);
se6  = zeros(Lx+1,Ly+1);
cost = zeros(Lx+1,Ly+1);

for lx = 0:Lx
    for ly = 0:Ly
    
    sums = [];
    cst  = 0;
    
    M = 20;
    
    parfor j = 1:M  
        [sums_j, cst_j] = mismc_l(lx,ly,N,sigma,data,z_data,K);
        sums(j,:) = sums_j;
        cst  = cst  + cst_j/M;   
    end
    
    idx = lx + 1;
    idy = ly + 1;
    
    cost(idx,idy) = cst;
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
savepath = fullfile(filepart, 'Results','empirical results','mismc_empirical.mat');
save(savepath,'Lx','Ly','K','cost','we','se')
end

