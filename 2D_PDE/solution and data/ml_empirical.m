% ------------------------------------------------------------------ %
% function ml_empirical(mlsmc_l,N,M,L,sigma,data,z_data,K)
%
% multlevel Monte Carlo empirical test routine
%
% input: mlsmc_l = function for level l estimator
%        N       = number of samples for convergence tests
%        M       = number of realizations
%        L       = number of levels for convergence tests
%        sigma   = std deviation of the error of y - G
%        data    = observations (values of y)
%        z_data  = corresponding observation points
%        K       = starting refinement level
%
% .mat: Lmax = maximum level we set
%       cost = cost at each level
%       K    = starting refinement level
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
function ml_empirical(mlsmc_l,N,M,L,sigma,data,z_data,K)
we1  = zeros(1,L+1);
we2  = zeros(1,L+1);
we3  = zeros(1,L+1);
se1  = zeros(1,L+1);
se2  = zeros(1,L+1);
se3  = zeros(1,L+1);
se4  = zeros(1,L+1);
se5  = zeros(1,L+1);
se6  = zeros(1,L+1);
cost = zeros(1,L+1);

tic;
for l = 0:L
    
    sums = [];
    cst  = 0;
    parfor j = 1:M
        [sums_j, cst_j] = mlsmc_l(l,l,N,sigma,data,z_data,K);
        sums(j,:) = sums_j;
        cst  = cst  + cst_j/M;  
    end
    id = l + 1;
    cost(1,id) = cst;
    
    we1(1,id)  = abs( sum(sums(:,1))/M );
    we2(1,id)  = abs( sum(sums(:,2))/M );
    we3(1,id)  = abs( sum(sums(:,3))/M );
    
    se1(1,id)  = N*sum( (sums(:,1) - sum(sums(:,1))/M).^2 )/M;
    se2(1,id)  = N*sum( (sums(:,2) - sum(sums(:,2))/M).^2 )/M;
    se3(1,id)  = N*sum( (sums(:,3) - sum(sums(:,3))/M).^2 )/M;

    se4(1,id)  = sum(sums(:,1).^2)/M;
    se5(1,id)  = sum(sums(:,2).^2)/M;
    se6(1,id)  = sum(sums(:,3).^2)/M;
end

we = {we1,we2,we3};
se = {se1,se2,se3,se4,se5,se6};

[filepart,~,~] = fileparts(pwd);
savepath = fullfile(filepart, 'Results','empirical results','mlsmc_empirical.mat');
save(savepath,'L','K','cost','we','se')
